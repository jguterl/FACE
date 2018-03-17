module modFACE_run
use modFACE_header
use modFACE_interface
 use modFACE_input
use modFACE_init
implicit none
contains
subroutine run_FACE(face_input,face_output)
        ! core routine of FACE:
        ! - read input file
        ! - initialize variable
        ! - restore state or restart
        ! - run the solver
        ! - dump final state files

    type(FACE_inputs) ,intent(in)  :: face_input
    type(FACE_outputs),intent(out) :: face_output
        ! get time at beginning of run
        call cpu_time(tcpustart)
        !read input
        call input_run(face_input)
        ! initialize variables
        call initialize
        ! restore from state or restart file
        call restore
        ! get initial inventory of species in material and surface
        call compute_init_inventory
        ! get info on initial wall temperature
        call compute_wall_temp_info(init_wall_temp)
        ! print info run
        call print_info_run
        ! execute time loop
        call time_loop
        ! dump final state
        call store_state(trim(final_state_file))
        ! get final inventory of species in material and surface
        call compute_final_inventory
        ! get info on final wall temperature
        call compute_wall_temp_info(final_wall_temp)
        ! compute the particle balance
        call compute_particle_balance
        ! compute outgassing flux info
        call compute_outgassing_flux
        ! get time at end of run
        call cpu_time(tcpufinish)
        ! generate output for FACE
        call output_run(face_output)
        ! print summary of the run
        call print_summary


    end subroutine run_FACE

    subroutine compute_init_inventory
        integer :: j,k
        real(DP):: s
        do k=1,nspc
            s=0d0
            do j=0,ngrd
                s=s+dens(ndt,j,k)*dx(j)
            enddo
            init_inventory(k)%Ntotbulk=s
            init_inventory(k)%Ntotsrf=dsrfl(ndt,k)+dsrfr(ndt,k)
            init_inventory(k)%Nnetbulk=0.d0
            init_inventory(k)%Nnetsrf=0.d0
        enddo

    end subroutine compute_init_inventory

        subroutine compute_final_inventory
        integer :: j,k
        real(DP):: s
        do k=1,nspc
            s=0d0
            do j=0,ngrd
                s=s+dens(ndt,j,k)*dx(j)
            enddo
            final_inventory(k)%Ntotbulk=s
            final_inventory(k)%Ntotsrf=dsrfl(ndt,k)+dsrfr(ndt,k)
            final_inventory(k)%Nnetbulk=final_inventory(k)%Ntotbulk-init_inventory(k)%Ntotbulk
            final_inventory(k)%Nnetsrf =final_inventory(k)%Ntotsrf-init_inventory(k)%Ntotsrf

        enddo



    end subroutine compute_final_inventory

    subroutine compute_particle_balance
    ! WARNING: this particle balance is performed for hydrogen assuming that free hydrogen are species k=1 and are trapped in species k=3,5,7,...
      integer :: k
      particle_balance%Ninflux=trace_flux(1)%sum_inflx
      particle_balance%Nnet=0d0
      particle_balance%Noutflux=0d0
      do k=1,nspc,2
      particle_balance%Nnet=particle_balance%Nnet+final_inventory(k)%Nnetbulk+final_inventory(k)%Nnetsrf
      particle_balance%Noutflux=particle_balance%Noutflux+trace_flux(k)%sum_Gdes_l+trace_flux(k)%sum_Gdes_r
      enddo



    end subroutine compute_particle_balance

    subroutine output_run(face_output)
    type(FACE_outputs),intent(out) :: face_output
    face_output%cpu_runtime=tcpufinish-tcpustart
    face_output%error_status=error_status

    call FACE2fluidcode(face_output%fluidcode_output)

    end subroutine output_run

    subroutine compute_outgassing_flux
    integer k
    real(DP)::sigma
    type(outgassing_fluxes) :: outgassing_flux
    outgassing_flux%min_Gdes=trace_flux(1)%min_Gdes_l
    outgassing_flux%max_Gdes=trace_flux(1)%max_Gdes_l
    outgassing_flux%Gpermeation=trace_flux(1)%max_Gdes_r
    if (iteration.gt.0) then
    outgassing_flux%ave_Gdes=trace_flux(1)%sum_Gdes_l/iteration
    sigma=trace_flux(1)%sig_Gdes_l/iteration-outgassing_flux%ave_Gdes**2
    else
     outgassing_flux%ave_Gdes=0d0
     sigma=0d0
    endif

    if (sigma.gt.0) then
    outgassing_flux%sig_Gdes=sqrt(sigma)
    else
    outgassing_flux%sig_Gdes=-sqrt(-sigma)
    endif

    ! we add a verifiction to throw an error if desorption of H occurs as desorption of filled traps.
    ! Filled traps usually do not diffuse in material but this check might be eased if necessary...
    do k=2,nspc
    if (trace_flux(k)%max_Gdes_l.gt.0d0.or.trace_flux(k)%max_Gdes_r.gt.0d0) then
    call face_error('Desorption of filled of emptry traps occurs. Please check consistancy of results, input and code')
    endif
    enddo

    end subroutine compute_outgassing_flux

    subroutine compute_wall_temp_info(wall_temp)
        type(wall_temperatures),intent(out):: wall_temp
        real(DP)::max_temp,min_temp
        integer::j
        wall_temp%sfr_temp_l=temp(ndt,0)
        wall_temp%sfr_temp_r=temp(ndt,ngrd)
        max_temp=0d0
        do j=0,ngrd
        max_temp=max(max_temp,temp(ndt,j))
        enddo
        wall_temp%max_temp=max_temp
        min_temp=1d99
        do j=0,ngrd
        min_temp=min(min_temp,temp(ndt,j))
        enddo
        wall_temp%min_temp=min_temp

        end subroutine compute_wall_temp_info

     subroutine input_run(face_input)

        type(FACE_inputs),intent(in):: face_input

        call init_input()
        if (face_input%read_input_file) then
            call read_inputfile(face_input%input_filename)
        else
            call set_input_parameters()
        endif



        call write_input_log

        if (couple_fluidcode) then
        call fluidcode2FACE(face_input%fluidcode_input)
        endif

    end subroutine input_run


subroutine time_loop
        ! iterative loop where R-D and heat flux equations are iteratively solved

        call print_milestone('starting time iteration')

        iteration=0
        do while (time .le. end_time)
            call print_timestep_info()                   ! print info on current time step
            call save()                                  ! dump data
            call step()                                  ! numerical solver
            call store_restart(trim(restart_filename))   ! store restart if requested
            iteration=iteration+1
        enddo

        call print_milestone('time iteration completed')

    end subroutine time_loop

    subroutine print_summary()
integer ::k
write(iout,*) 'Successful run of FACE !'
write(iout,*) '***************Summary***************'
write(iout,*) '# iteration = ', iteration
write(iout,'("cpu time of execution= ",es12.3," seconds.")') tcpufinish-tcpustart
write(iout,*) '*************************************'

do k=1,nspc
write(iout,'(a,i2,a,e9.2,a,e9.2)') "Final inventory k=",k ,"; bulk/surface:",final_inventory(k)%Nnetbulk,' / ',&
 final_inventory(k)%Nnetsrf
enddo
write(iout,*) " "
      write(iout,*) " *** Particle balance ***"
      write(iout,'(a,es12.3,a,es12.3,a,es12.3)') " * Ninflux=",particle_balance%Ninflux,"; Nnet_material=",particle_balance%Nnet,&
      "; Noutflux=",particle_balance%Noutflux
      write(iout,'(a,es12.3)') " * Ninflux-Noutflux+Nnet=",particle_balance%Ninflux-particle_balance%Noutflux+particle_balance%Nnet
      write(iout,*) " *** "
end subroutine print_summary

    subroutine print_info_run
call print_milestone('info run')
write(iout,'(a,es12.3,a,es12.3,a,es12.3)') "start_time=",start_time,"; end_time=",end_time," ;dt_face=",dt_face
write(iout,'(a,i5)') 'estimated # iterations =', int((end_time-start_time)/dt_face)
write(iout,'(a,i4,a,i4,a,a4)') 'nspc = ', nspc,"; ngrid = ",ngrd,"; solve_heat_eq : ",solve_heat_eq
write(iout,*) ''
end subroutine print_info_run
end module modFACE_run
