module modFACE_run
!$     use omp_lib 
use modFACE_header
 use modFACE_input
use modFACE_init
    use modFACE_coupling
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
        walltime_start=tcpustart
        !$ walltime_start = omp_get_wtime ( )
        
        call print_headline('Start FACE run')
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
        !walltime_end=tcpufinish
        !$ walltime_end = omp_get_wtime ( )
       
! generate output for FACE
        call output_run(face_output)

        call print_headline('FACE run completed')
        ! print summary of the run
        call print_summary


    end subroutine run_FACE

    subroutine compute_init_inventory
        integer :: j,k
        real(DP):: s
        do k=1,nspc

            init_inventory(k)%Ntotbulk=integrale_dens(k)
            init_inventory(k)%Ntotsrf=0d0
            if (left_surface_model(k).eq."S") then
            init_inventory(k)%Ntotsrf=init_inventory(k)%Ntotsrf+dsrfl(ndt,k)
            endif
            if (right_surface_model(k).eq."S") then
            init_inventory(k)%Ntotsrf=init_inventory(k)%Ntotsrf+dsrfr(ndt,k)
            endif

            init_inventory(k)%Nnetbulk=0.d0
            init_inventory(k)%Nnetsrf=0.d0
        enddo

    end subroutine compute_init_inventory

        subroutine compute_final_inventory
        integer :: j,k
        real(DP):: s
        do k=1,nspc

            final_inventory(k)%Ntotbulk=integrale_dens(k)
            final_inventory(k)%Ntotsrf=0d0
            if (left_surface_model(k).eq."S") then
            final_inventory(k)%Ntotsrf=final_inventory(k)%Ntotsrf+dsrfl(ndt,k)
            endif
             if (right_surface_model(k).eq."S") then
            final_inventory(k)%Ntotsrf=final_inventory(k)%Ntotsrf+dsrfr(ndt,k)
            endif
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

      particle_balance%Nnet=particle_balance%Nnet+final_inventory(1)%Nnetbulk+final_inventory(1)%Nnetsrf
      particle_balance%Noutflux=particle_balance%Noutflux+trace_flux(1)%sum_Gdes_l+trace_flux(1)%sum_Gdes_r
      particle_balance%Noutflux_l=trace_flux(1)%sum_Gdes_l
      particle_balance%Noutflux_r=trace_flux(1)%sum_Gdes_r
      if (nspc.gt.2) then
      do k=3,nspc,2
      particle_balance%Nnet=particle_balance%Nnet+final_inventory(k)%Nnetbulk+final_inventory(k)%Nnetsrf
      particle_balance%Noutflux=particle_balance%Noutflux+trace_flux(k)%sum_Gdes_l+trace_flux(k)%sum_Gdes_r
      enddo
      endif

      particle_balance%p_net=particle_balance%Ninflux-particle_balance%Noutflux-particle_balance%Nnet
      particle_balance%p_max=max(max(abs(particle_balance%Ninflux),abs(particle_balance%Noutflux)),abs(particle_balance%Nnet))
      if (particle_balance%p_max.ne.0d0) then
      particle_balance%f_lost=abs(particle_balance%p_net)/particle_balance%p_max
      else
      particle_balance%f_lost=1d99
      endif

    end subroutine compute_particle_balance

    subroutine output_run(face_output)
    type(FACE_outputs),intent(out) :: face_output
    face_output%cpu_runtime=tcpufinish-tcpustart
    face_output%error_status=error_status
    face_output%nspc=nspc

    allocate(FACE_output%init_inventory(face_output%nspc))
    allocate(FACE_output%final_inventory(face_output%nspc))
    FACE_output%final_inventory=final_inventory
    FACE_output%init_inventory=init_inventory
    FACE_output%particle_balance=particle_balance
    FACE_output%outgassing_flux=outgassing_flux
    FACE_output%init_wall_temp=init_wall_temp
    FACE_output%final_wall_temp=final_wall_temp

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
    if (trace_flux(k)%max_Gdes_l.gt.1d-4*trace_flux(1)%max_Gdes_l) then
    call face_warning('Desorption of traps occurs and is larger than 0.1% of H desorption.Check consistancy(k)='&
    ,trace_flux(k)%max_Gdes_l,"Gdes_l(1)=",trace_flux(1)%max_Gdes_l)
    endif
    enddo

    end subroutine compute_outgassing_flux

    subroutine compute_wall_temp_info(wall_temp)
        type(wall_temperatures),intent(out):: wall_temp
        real(DP)::max_temp,min_temp
        integer::j
        wall_temp%srf_temp_l=temp(ndt,0)
        wall_temp%srf_temp_r=temp(ndt,ngrd)
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





        if (face_input%couple_fluidcode) then
        call input_fluidcode(face_input%fluidcode_input)
        endif

    end subroutine input_run


subroutine time_loop
        ! iterative loop where R-D and heat flux equations are iteratively solved

        call print_milestone('starting time iteration')

        iteration=0

        do while ((time .lt. end_time) .AND. (iteration.lt.(int(max_iter))))

            call save                                  ! dump data
            call step                                  ! numerical solver

            call store_restart(trim(restart_filename))   ! store restart if requested
            iteration=iteration+1
        enddo

        call print_milestone('time iteration completed')

    end subroutine time_loop

    subroutine print_summary()
integer ::k,nthreads
character(string_length) :: str

call print_section('Summary')
write(str,*) '# iteration = ', iteration
call print_line(str)
write(str,'("cpu time of execution= ",es12.3," seconds.")') tcpufinish-tcpustart
call print_line(str)
write(str,'("wall time of execution= ",es12.3," seconds.")') walltime_end-walltime_start
call print_line(str)
nthreads=1
!$nthreads = OMP_GET_NUM_THREADS();
write(str,*) "nthreads  = ",nthreads 
call print_line(str)

call print_end_section('Summary')
call print_section('Final inventory')
do k=1,nspc
 write(str,'(a,i2,a,e12.3)') "Final net inventory k= ",k    ," - bulk    : ",final_inventory(k)%Nnetbulk
 call print_line(str)
 write(str,'(a,a,a,e12.3)') "                       ","  " ," - surface : ", final_inventory(k)%Nnetsrf
call print_line(str)
write(str,*) " "
call print_line(str)
enddo
call print_end_section('Final inventory')
call print_section('Particle balance')
 write(str,'(a,es12.3,a,es12.3,a,es12.3,a,es12.3,a,es12.3,a)') "Ninflux=",particle_balance%Ninflux,&
 "; Nnet_material=",particle_balance%Nnet,"; Noutflux=",particle_balance%Noutflux,&
      "(" , particle_balance%Noutflux_l,"/", particle_balance%Noutflux_r,")"
     call print_line(str)
      write(str,'(a,es12.3)') "Ninflux-Noutflux-Nnet=",particle_balance%p_net
      call print_line(str)
      write(str,'(a,es12.3)') "Fraction of numerically lost particles =",particle_balance%f_lost
call print_line(str)
call print_end_section('Particle balance')
end subroutine print_summary

    subroutine print_info_run
    character(string_length)::str
call print_section('Run features')
write(str,'(a,es12.3,a,es12.3,a,es12.3)') "start_time=",start_time,"; end_time=",end_time," ;dt_face=",dt_face
call print_line(str)
write(str,'(a,i5)') 'estimated # iterations =', int((end_time-start_time)/dt_face)
call print_line(str)
write(str,'(a,i4,a,i4,a,a4)') 'nspc = ', nspc,"; ngrid = ",ngrd,"; solve_heat_eq : ",solve_heat_eq_string
call print_line(str)
call print_end_section('Run features')
end subroutine print_info_run
end module modFACE_run
