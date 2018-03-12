module modFACE_interface
    use modFACE_precision
    implicit none

    save
    logical :: verbose_interface=.true.
    character(Lfn) :: default_inputfile="default_inputfile.face"
    type fluidcode_inputs
    integer                  :: wall_idx            ! Index of the wall stratum
        integer                  :: iter                ! Fluid code iteration
        real(DP)                 :: time                ! Fluid code time
        real(DP)                 :: dt                  ! Time step of the fluid
        real(DP)                 :: dt_face             ! Time step of FACE
        integer                 :: nspc                 ! Number of incoming species from fluid code
        integer,allocatable          :: indexspc(:)             ! Index of species in FACE (usually "k" in FACE)
        character(Lname),allocatable :: namespc(:)      ! Name of the incoming species from fluid code
        real(DP),allocatable         :: Gammain(:)      ! Particle flux'
        real(DP),allocatable         :: Emean(:)        ! Average energy of incoming particle enrg'
        real(DP)                 :: Qin                 ! Heat flux from fluid code
        character(Lfn)           :: history_file        ! restart from file history
        real(DP)                 :: tempwall            ! temperature of the wall from fluid code
        character(15)            :: solve_heat_eq       ! if solve_heat_eq then use Qin otherwise T=tempwall for the entire bulk
        character(Lfn)            ::casename            ! casename

    end type fluidcode_inputs

    type fluidcode_out
        real(DP),allocatable   :: out_flux(:)            ! outgassing flux
        real(DP),allocatable   :: Ntot(:)                ! total species in material
        real(DP),allocatable   :: Nnet_dt_fc(:)          ! total species added during dt_fx
        real(DP),allocatable   :: Nnet_free_dt_fc(:)     ! total free species added during dt_fx
        real(DP),allocatable   :: Nnet_trapped_dt_fc(:)  ! total trapped species added during dt_fx
        real(DP),allocatable   :: Nnet_srf_dt_fc(:)      ! total species added on surface during dt_fx
        real(DP),allocatable   :: Nnet_vol_dt_fc(:)      ! total species added in bulk during dt_fx
    end type fluidcode_out

    type FACE_inputs
        character(Lfn):: run_mode
        character(Lfn):: input_filename
        logical   :: read_input_file=.true.
        character(Lfn)::logfile
        logical :: couple_fluidcode
        type(fluidcode_inputs) :: fluidcode_input
        character(Lfn)::path
        character(Lfn)::casename
    end type FACE_inputs

contains
    subroutine read_arguments(face_input)
        type(FACE_inputs)::face_input
        CHARACTER(lfn) :: arg
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
        write(*,*) '# arguments=',narg

        k=1
        DO while (k.le.narg)
            CALL getarg(k, arg)
            if (verbose_interface) write(*,*)arg
            ! check inputfile flag -in or --inputfile
            if (arg.eq."-in" .or. arg.eq."--inputfile") then
                if (k.le.narg) then
                    k=k+1
                    call getarg(k,arg)
                    if (verbose_interface) write(*,*)arg
                    if (arg(1:1).eq.'-') then
                        write(*,*) 'ERROR: Invalid input filename :',arg
                        stop 'Exiting FACE...'
                    endif
                else
                    write(*,*) 'ERROR: Missing input filename'
                    stop 'Exiting FACE...'
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

            elseif (arg.eq."-pg" .or. arg.eq."--print-grid") then
                if (face_input%run_mode.ne."default") then
                    write(*,*) 'ERROR: Current flag incompatible with previous flags (current flag: ',arg,')'
                    stop 'Exiting FACE...'
                else
                    face_input%run_mode="print_grid"
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

            ! help
            elseif (arg.eq."-h" .or. arg.eq."--help") then
                if (face_input%run_mode.ne."default") then
                    write(*,*) 'ERROR: Current flag incompatible with previous flags (current flag: ',arg,')'
                    stop 'Exiting FACE...'
                else
                    face_input%run_mode="help"
                    k=k+1
                endif

            ! version
            elseif (arg.eq."-v" .or. arg.eq."--version") then
                if (face_input%run_mode.ne."default") then
                    write(*,*) 'ERROR: Current flag incompatible with previous flags (current flag: ',arg,')'
                    stop 'Exiting FACE...'
                else
                    face_input%run_mode="version"
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


    subroutine set_default_mode(default_mode)
        character(*):: default_mode

        if (default_mode.eq."H") then
        elseif (default_mode.eq."H+5traps") then
        endif

    end subroutine



    subroutine wrapper_FACE(face_input,fluidcode_input)
        type(FACE_inputs),intent(out)::face_input
        type(fluidcode_inputs), intent(in) :: fluidcode_input
        face_input%logfile="no"
        face_input%input_filename=default_inputfile
        face_input%run_mode='default'
        face_input%path='run_FACE'
        face_input%casename='casename'
        face_input%fluidcode_input=fluidcode_input
    end subroutine wrapper_FACE


    subroutine set_fluidcode_input(fluidcode_input)
        type(fluidcode_inputs)::fluidcode_input
        character(200)::str
        !    type fluidcode_input
        !        integer                  :: wall_idx        ! Index of the wall stratum
        !        integer                  :: iter            ! Fluid code iteration
        !        real(DP)                 :: time            ! Fluid code time
        !        real(DP)                 :: dt              ! Time step of the fluid
        !        real(DP)                 :: nspc            ! Number of incoming species from fluid code
        !        character(Lname),allocatable :: namespc(:)  ! Name of the incoming species from fluid code
        !        real(DP),allocatable     :: Gammain(:)      ! Particle flux'
        !        real(DP),allocatable     :: Emean(:)        ! Average energy of incoming particle enrg'
        !        real(DP)                 :: Qin_            ! Heat flux from fluid code
        !        character(Lfn)           :: history_file    ! restart from file history
        !        real(DP)                 :: tempwall        ! temperature of the wall from fluid code
        !        character(15)            :: solve_heat_eq   ! if solve_heat_eq then use Qin otherwise T=tempwall for the entire bulk
        !    end type fluidcode_input

        fluidcode_input%wall_idx=1
        fluidcode_input%iter=1
        fluidcode_input%time=1e-5
        fluidcode_input%dt=1e-2
        fluidcode_input%nspc=1
        allocate(fluidcode_input%namespc(fluidcode_input%nspc))
        fluidcode_input%namespc(1:fluidcode_input%nspc)="D"
        allocate(fluidcode_input%indexspc(fluidcode_input%nspc))
        fluidcode_input%indexspc(1:fluidcode_input%nspc)=1
        allocate(fluidcode_input%Gammain(fluidcode_input%nspc))
        fluidcode_input%Gammain(1:fluidcode_input%nspc)=1e20
        allocate(fluidcode_input%Emean(fluidcode_input%nspc))
        fluidcode_input%Emean(1:fluidcode_input%nspc)=100
        fluidcode_input%Qin=0           ! Heat flux from fluid code
        fluidcode_input%history_file="no"    ! restart from file history
        fluidcode_input%tempwall=700        ! temperature of the wall from fluid code
        fluidcode_input%solve_heat_eq="no"   ! if solve_heat_eq then use Qin otherwise T=tempwall for the entire bulk
        fluidcode_input%casename='solps'
        write(str,'(a,a,i4.4,a,i6.6)') fluidcode_input%casename,'_',fluidcode_input%wall_idx,'_',fluidcode_input%iter
        fluidcode_input%casename=trim(str)
    !    end type fluidcode_input
    end subroutine set_fluidcode_input


end module modFACE_interface
