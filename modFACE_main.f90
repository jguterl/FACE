module modFACE_main
    use modFACE_header
    use modFACE_parser
    use modFACE_output
    use modFACE_help
    use modFACE_input
    use modFACE_init
    use modFACE_IO
    use modFACE_interface
    implicit none

contains

    subroutine FACE_main(face_input)
        type(FACE_inputs) :: face_input

        ! Initialization of path, casename, logfile,...
        call initialization_main(face_input)

        !  Execute requested run mode
        face_mode:select case(face_input%run_mode)

            case("help")
                call display_help

            case ("version")
                call write_version

            case ("print_default_input")
                call write_version()
                call write_default_inputfile(default_inputfile,"H+Tr")

            case ("default_H") ! default modes-> print out default input file then read file
                call write_version()
                call write_default_inputfile(default_inputfile,"H")
                call write_header_log()
                call run_FACE(face_input)

            case ("default_H+Tr") ! default modes-> print out default input file then read file
                call write_version()
                call write_default_inputfile(default_inputfile,"H+Tr")
                call write_header_log()
                call run_FACE(face_input)

            case default ! all other modes which request reading of input file
                call write_short_version()
                call write_header_log()
                call run_FACE(face_input)

        end select face_mode

        ! Final message if execution was successful

        call cpu_time(tcpufinish)

        call print_summary

        call finalize
        stop
    end    subroutine FACE_main

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

    subroutine run_FACE(face_input)
        ! core routine of FACE:
        ! - read input file
        ! - initialize variable
        ! - restore state or restart
        ! - run the solver
        ! - dump final state files

    type(FACE_inputs) :: face_input

        !read input
        call get_input(face_input)

        ! initialize variables
        call initialize()
        ! restore from state or restart file
        call restore()
        ! get initial inventory of species in material and surface
        call get_init_inventory
        ! execute time loop
        call time_loop()
        ! dump final state
        call store_state(final_state_file)
        ! get final inventory of species in material and surface
        call get_final_inventory

        call particle_balance

        if (face_input%couple_fluidcode) then
        call FACE2fluidcode()
        endif

    end subroutine run_FACE

    subroutine initialization_main(face_input)

        type(FACE_inputs) :: face_input

        call cpu_time(tcpustart)

        call init_casename(face_input%casename)

        call init_path(face_input%path)
        ! initialize help: required to check keywords in input file
        call init_help()
        ! Initialize log output (file or standard output (unit=6))
        call init_log(trim(face_input%logfile))

    end subroutine initialization_main

        subroutine finilization_main()

        if (iout.gt.10) then
        close(iout)
        endif

    end subroutine finilization_main

    subroutine get_init_inventory
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

    end subroutine get_init_inventory

        subroutine get_final_inventory
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
        write(iout,'(a,e9.2)') "Final inveotry net surface:",final_inventory(k)%Nnetsrf
        enddo

    end subroutine get_final_inventory

    subroutine particle_balance
    ! WARNING: this particle balance is performed for hydrogen assuming that free hydrogen are species k=1 and are trapped in species k=3,5,7,...
      real(DP):: Nnet,Ninflux,Noutflux
      integer :: k
      Ninflux=trace_flux(1)%sum_inflx
      Nnet=0d0
      do k=1,nspc,2
      Nnet=Nnet+final_inventory(k)%Nnetbulk+final_inventory(k)%Nnetsrf
      Noutflux=trace_flux(k)%sum_Gdes_l+trace_flux(k)%sum_Gdes_r
      enddo
      write(iout,*) " *** Particle balance ***"
      write(iout,*) " * Ninflux=",Ninflux,"; Nnet_material=",Nnet,"; Noutflux=",Noutflux
      write(iout,*) " * Ninflux+Noutflux-Nnet",Ninflux-Noutflux+Nnet
      write(iout,*) " *** "


    end subroutine

end module modFACE_main
