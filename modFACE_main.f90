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
        ! execute time loop
        call time_loop()
        ! dump final state
        call output_final_state()

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

end module modFACE_main
