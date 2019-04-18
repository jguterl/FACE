module modFACE_main
    use modFACE_header
    use modFACE_help
    use modFACE_IO

    use modFACE_run
    implicit none

contains

 subroutine FACE_main(face_input,face_output)
        type(FACE_inputs),intent(in) :: face_input
        type(FACE_outputs),intent(out) :: face_output

        ! Initialization of path, casename, logfile,...
        call initialization_main(face_input)

        !  Execute requested run mode
        face_mode:select case(face_input%run_mode)

            case ("help")
                call print_help

            case ("version")
                call print_version

            case ("list-keyword")
                call print_version
                call print_list_keyword
            case ("print_default_input")
                call print_version()
                call write_default_inputfile(default_inputfile,"H+Tr")
                case ("print_default_input_nohelp")
                call print_version()
                call write_default_inputfile_nohelp(default_inputfile,"H+Tr")

            case ("default_H") ! default modes-> print out default input file then read file
                call print_version()
                call write_default_inputfile(default_inputfile,"H")
                call write_default_inputfile_nohelp(default_inputfile,"H")
                call write_header_log()
                call run_FACE(face_input,face_output)

            case ("default_H+Tr") ! default modes-> print out default input file then read file
                call print_version()
                call write_default_inputfile(default_inputfile,"H+Tr")
                call write_default_inputfile_nohelp(default_inputfile,"H+Tr")
                call write_header_log()
                call run_FACE(face_input,face_output)

            case default ! all other modes which request reading of input file
                call write_header_log()
                call run_FACE(face_input,face_output)

        end select face_mode

        call finalization_main
    end    subroutine FACE_main





    subroutine initialization_main(face_input)

        type(FACE_inputs) :: face_input



        call init_casename(face_input%casename)

        call init_path(face_input%path)
        ! initialize help: required to check keywords in input file
        call init_help()
        ! Initialize log output (file or standard output (unit=6))
        call init_log(trim(face_input%logfile))

    end subroutine initialization_main

 subroutine finalization_main()
        call close_timedata_files
        call deallocate_variables()
        call close_log
        write(iout,*) 'Quitting FACE after complete normal execution...'
    end subroutine finalization_main
        subroutine finilization_main

        if (iout.gt.10) then
        close(iout)
        endif

    end subroutine finilization_main



end module modFACE_main
