module modFACE_interface
    use modFACE_precision
    use modFACE_error
    use modFACE_coupling
    use modFACE_header
    use modFACE_main
    implicit none





contains
    subroutine read_arguments(face_input)
        type(FACE_inputs)::face_input
        CHARACTER(string_length) :: arg
        integer narg,k
        !      running -pg --print-grid : only print the grid correspondign to the input file. No run
        !      running -ps --print-species:   dump vol and surface after reading restart file or history file if requested
        !      running -pdi --print-default-input:   print default_input.face
        !      running -h  --help:
        !      running -dh --default-H: run with default input file with only H
        !      running -dhtr --default-H+Tr: run with default input file with only H+Tr
        !      running  -c --couple-fluidcode
        !      running -in --inputfile inputfile
        !              -l --log logfile
        !default:
        face_input%logfile="no"
        face_input%input_filename=default_inputfile
        face_input%run_mode='default'
        face_input%path='run_FACE'
        face_input%casename='casename'
        narg=iargc()
        !write(*,*) '# arguments=',narg

        k=1
        do while (k.le.narg)
            CALL getarg(k, arg)
            if (verbose_interface) write(*,*)arg
            ! check inputfile flag -in or --inputfile
            if (arg.eq."-in" .or. arg.eq."--inputfile") then
                if (k.le.narg) then
                    k=k+1
                    call getarg(k,arg)
                    if (verbose_interface) write(*,*)arg
                    if (arg(1:1).eq.'-') then
                        call face_error('Invalid input filename :',arg)
                    endif
                else
                    call facE_error('Missing input filename')
                endif
                face_input%input_filename=arg
                k=k+1

            ! check inputfile flag -in or --inputfile
            elseif (arg.eq."-dh" .or. arg.eq."--default-H") then
                if (face_input%run_mode.ne."default") then
                    write(*,*) 'ERROR: Current flag incompatible with previous flags (current flag: ',arg,')'
                    stop 'Exiting FACE...'
                else
                    face_input%run_mode="default_H"
                    face_input%input_filename=default_inputfile
                    k=k+1
                endif

            elseif (arg.eq."-dhtr" .or. arg.eq."--default-H+Tr") then
                if (face_input%run_mode.ne."default") then
                    write(*,*) 'ERROR: Current flag incompatible with previous flags (current flag: ',arg,')'
                    stop 'Exiting FACE...'
                else
                    face_input%run_mode="default_H+Tr"
                    face_input%input_filename=default_inputfile
                    k=k+1
                endif

            elseif (arg.eq."-pd" .or. arg.eq."--print-default-input") then
                if (face_input%run_mode.ne."default") then
                    write(*,*) 'ERROR: Current flag incompatible with previous flags (current flag: ',arg,')'
                    stop 'Exiting FACE...'
                else
                    face_input%run_mode="print_default_input"
                    k=k+1
                endif
                    elseif (arg.eq."-vi" .or. arg.eq."--verbose-input") then
                    verbose_input=.true.
                    k=k+1
                    elseif (arg.eq."-vs" .or. arg.eq."--verbose-surface") then
                    verbose_surface=.true.
                    k=k+1
                    elseif (arg.eq."-vp" .or. arg.eq."--verbose-parser") then
                    verbose_parser=.true.
                    k=k+1
                    elseif (arg.eq."-vd" .or. arg.eq."--verbose-debug") then
                    verbose_debug=.true.
                    k=k+1
                    elseif (arg.eq."-vit" .or. arg.eq."--verbose-init") then
                    verbose_init=.true.
                    k=k+1

                    elseif (arg.eq."-vc" .or. arg.eq."--verbose-couple") then
                    verbose_couple=.true.
                    k=k+1
                    elseif (arg.eq."-vv" .or. arg.eq."--verbose-version") then
                    verbose_version=.true.
                    k=k+1
                    elseif (arg.eq."-vh" .or. arg.eq."--verbose-header") then
                    verbose_header=.true.
                    k=k+1


            ! help
            elseif (arg.eq."-h" .or. arg.eq."--help") then
                if (face_input%run_mode.ne."default") then
                    call face_error('Current flag incompatible with previous flags (current flag: ',arg,')')
                else
                    face_input%run_mode="help"
                    k=k+1
                endif

            ! version
            elseif (arg.eq."-v" .or. arg.eq."--version") then
                if (face_input%run_mode.ne."default") then
                    call face_error('Current flag incompatible with previous flags (current flag: ',arg,')')
                else
                    face_input%run_mode="version"
                    k=k+1
                endif
                ! keyword_list
                elseif (arg.eq."-lk" .or. arg.eq."--list-keyword") then
                if (face_input%run_mode.ne."default") then
                    call face_error('Current flag incompatible with previous flags (current flag: ',arg,')')
                else
                    face_input%run_mode="list-keyword"
                    k=k+1
                endif
                 ! output folder
            elseif (arg.eq."-p" .or. arg.eq."--path") then
                if (k.le.narg) then
                    k=k+1
                    call getarg(k,arg)
                    if (verbose_interface) write(*,*)arg
                    if (arg(1:1).eq.'-') then
                        write(*,*) 'ERROR: Invalid folder path :',arg
                        stop 'Exiting FACE...'
                    endif
                else
                    write(*,*) 'ERROR: Missing folder path'
                    stop 'Exiting FACE...'
                endif
                face_input%path=arg
                k=k+1

            !
            elseif (arg.eq."-n" .or. arg.eq."--name-case") then
                if (k.le.narg) then
                    k=k+1
                    call getarg(k,arg)
                    if (verbose_interface) write(*,*)arg
                    if (arg(1:1).eq.'-') then
                        write(*,*) 'ERROR: Invalid case name :',arg
                        stop 'Exiting FACE...'
                    endif
                else
                    write(*,*) 'ERROR: Missing name case'
                    stop 'Exiting FACE...'
                endif
                face_input%casename=arg
                k=k+1

            elseif (arg.eq."-c" .or. arg.eq."--couple-fluidcode") then
                face_input%run_mode="coupled"
                k=k+1
            elseif (arg.eq."-l" .or. arg.eq."--log") then
                if (k.le.narg) then
                    k=k+1
                    call getarg(k,arg)
                    if (verbose_interface) write(*,*)arg
                    if (arg(1:1).eq.'-') then
                        write(*,*) 'ERROR: Invalid log filename :',arg
                        stop 'Exiting FACE...'
                    endif
                else
                    write(*,*) 'ERROR: Missing log filename'
                    stop 'Exiting FACE...'
                endif
                face_input%logfile=arg
                k=k+1
            else
                write(*,*) 'Unknow flag for FACE execution: ',arg
                stop 'Exiting FACE...'
            endif
        END DO
        if (face_input%input_filename.ne.default_inputfile.and.face_input%run_mode.eq.'default') then
            face_input%run_mode="normal"
        endif
        if (face_input%input_filename.eq.default_inputfile.and.face_input%run_mode.eq.'default') then
            face_input%run_mode="default_H"
        endif
        if (verbose_interface) then
            WRITE (*,*) "filename=",face_input%input_filename
            WRITE (*,*) "run_mode=",face_input%run_mode
            WRITE (*,*) "logfile=",face_input%logfile
        endif
    end subroutine read_arguments



    subroutine FACE_from_fluidcode

character(string_length)::str
integer:: iter
        type(fluidcode_inputs)::fluidcode_input
        type(fluidcode_outputs)::fluidcode_output
        ! folder
        fluidcode_input%path='solps_test'
        ! casename
        fluidcode_input%casename_base='solps1'
        ! first restore state file
        fluidcode_input%read_state_file="no"
        ! species information and allocate
        fluidcode_input%nspc_fluid=1
        call alloc_fluidcode_input(fluidcode_input)
        fluidcode_input%namespc(1:fluidcode_input%nspc_fluid)="D"
        fluidcode_input%indexspc(1:fluidcode_input%nspc_fluid)=1
        ! index of wall element
        fluidcode_input%wall_idx=1
        fluidcode_input%Ndump_space=100
        fluidcode_input%Ndump_time=100
        ! solps time step

        ! time step FACE
        fluidcode_input%dt0_face=1e-4

        ! input file FACE
        fluidcode_input%input_file="test_solpsiter.face"
        !  iteration
        do iter=1,10
        write(iout,*) "SOLPS: iteration =", iter
        fluidcode_input%iter=iter

        write(str,'(a,a1,i0,a1,i0)') trim(fluidcode_input%casename_base),'_',(fluidcode_input%wall_idx),'_',fluidcode_input%iter

        fluidcode_input%casename=trim(str)

        ! save the final state in
        fluidcode_input%final_state_file=trim(fluidcode_input%path)//'/'//trim(fluidcode_input%casename)//"_final.state"

        !Gamma in
        if (mod(iter,2).eq.0) then
        fluidcode_input%inflx_in(1:fluidcode_input%nspc_fluid)=1e21
        fluidcode_input%tempwall=800        ! temperature of the wall from fluid code
        fluidcode_input%dt=1e-3
        else
        fluidcode_input%inflx_in(1:fluidcode_input%nspc_fluid)=1e20
        fluidcode_input%tempwall=500        ! temperature of the wall from fluid code
        fluidcode_input%dt=1e-2

        endif
        ! solps time
        fluidcode_input%time=(fluidcode_input%iter-1)*fluidcode_input%dt


        fluidcode_input%qflx_in=0           ! Heat flux from fluid code

        fluidcode_input%solve_heat_eq_string="no"   ! if solve_heat_eq then use Qin otherwise T=tempwall for the entire bulk
        fluidcode_input%log_file=trim(fluidcode_input%path)//'/'//trim(fluidcode_input%casename)//".log"
fluidcode_input%log_file="no"
        call wrapper_FACE(fluidcode_input,fluidcode_output)
        fluidcode_input%read_state_file=fluidcode_input%final_state_file
        enddo

        call dealloc_fluidcode_input(fluidcode_input)
        call dealloc_fluidcode_output(fluidcode_output)

    end subroutine FACE_from_fluidcode

    subroutine FACE_standalone
     type(FACE_inputs)::face_input
      type(FACE_outputs)::face_output

    call read_arguments(face_input)
    call FACE_main(face_input,face_output)
    end subroutine FACE_standalone


    subroutine wrapper_FACE(fluidcode_input,fluidcode_output)
        type(FACE_inputs)::face_input
        type(FACE_outputs)::face_output
        type(fluidcode_inputs), intent(in) :: fluidcode_input
        type(fluidcode_outputs), intent(out) :: fluidcode_output
        call read_arguments(face_input)
        call fluidcode2FACE(face_input,fluidcode_input)
        call FACE_main(face_input,face_output)
        call FACE2fluidcode(face_output,fluidcode_output)
        !call trace_FACE_flx(fluidcode_input,fluidcode_output)
    end subroutine wrapper_FACE

subroutine trace_FACE_flx(fluidcode_input,fluidcode_output)
        type(fluidcode_inputs)::fluidcode_input
        type(fluidcode_outputs)::fluidcode_output
        integer ::unit_trace_flx,proc=0,ios
        character(string_length):: filename

        call set_unit(unit_trace_flx)
        filename=trim(fluidcode_input%path)//"/FACE.flux"
        open (unit_trace_flx, file=trim(filename),status='old', position='APPEND',iostat=ios)
        if (ios.ne.0) then
        open (unit_trace_flx, file=trim(filename),status='new',iostat=ios)
        if (ios.ne.0) then
         call face_error('Cannot open trace_FACE file ', trim(filename))
         else
         write(unit_trace_flx,'(3a4,5a12)') "proc", &
         "idx_wall", &
         "iter_fc", &
         "G_in", &
         "G_des",&
         "aveG_des",&
         "maxG_des",&
          "G_perm"
        endif
        endif

           write(unit_trace_flx,'(3i4,5es12.3)') proc,&
           fluidcode_input%wall_idx,&
           fluidcode_input%iter,&
           fluidcode_input%inflx_in(1),&
           fluidcode_output%outgassing_flux%Gdes,&
           fluidcode_output%outgassing_flux%ave_Gdes,&
           fluidcode_output%outgassing_flux%max_Gdes,&
           fluidcode_output%outgassing_flux%Gpermeation

           close(unit_trace_flx)

        end subroutine trace_FACE_flx

end module modFACE_interface
