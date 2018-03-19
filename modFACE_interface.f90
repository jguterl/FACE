module modFACE_interface
    use modFACE_precision
    use modFACE_error
    use modFACE_header
    implicit none

    save
    logical:: couple_fluidcode=.false.
    logical :: verbose_interface=.true.
    character(string_length) :: default_inputfile="default_inputfile.face"
    type fluidcode_inputs
         integer                     :: wall_idx            ! Index of the wall stratum
        integer                      :: iter                ! Fluid code iteration
        integer                      :: Ndump_space         ! # space data files to be dumped
        integer                      :: Ndump_time          !  # of times time data are dumped
        logical                      :: append              ! append mode for dumping data
        real(DP)                     :: time                ! Fluid code time
        real(DP)                     :: dt                  ! Time step of the fluid
        real(DP)                     :: dt_face             ! Time step of FACE
        integer                      :: nspc_fluid                 ! Number of incoming species from fluid code
        integer,allocatable          :: indexspc(:)             ! Index of species in FACE (usually "k" in FACE)
        character(Lname),allocatable :: namespc(:)      ! Name of the incoming species from fluid code
        real(DP),allocatable         :: inflx_in(:)      ! Particle flux'
        real(DP),allocatable         :: Emean(:)        ! Average energy of incoming particle enrg'
        real(DP)                     :: qflx_in                 ! Heat flux from fluid code
        character(string_length)     :: restore_state_file  ! restart from this staste file
        character(string_length)     :: final_state_file    ! store final state in this state file
        real(DP)                     :: tempwall            ! temperature of the wall from fluid code
        character(15)                :: solve_heat_eq       ! if solve_heat_eq then use Qin otherwise T=tempwall for the entire bulk
        character(string_length)            ::casename            ! casename

    end type fluidcode_inputs

        type fluidcode_outputs
        integer              :: nspc_fluid                 ! Number of incoming species from fluid code
        integer              :: nspc_face                 ! Number of incoming species from fluid code
        integer,allocatable          :: indexspc(:)             ! Index of species in FACE (usually "k" in FACE)
             type(outgassing_fluxes):: outgassing_flux
             type(wall_temperatures) :: init_wall_temp,final_wall_temp
             type(particle_balances)  :: particle_balance
             type(inventories),allocatable       :: init_inventory(:),final_inventory(:)

           end type fluidcode_outputs






    type FACE_inputs
        character(string_length):: run_mode
        character(string_length):: input_filename
        logical   :: read_input_file=.true.
        character(string_length)::logfile
        logical :: couple_fluidcode
        type(fluidcode_inputs) :: fluidcode_input
        character(string_length)::path
        character(string_length)::casename
    end type FACE_inputs

    type FACE_outputs
        real(DP) :: cpu_runtime
        integer :: error_status

        type(fluidcode_outputs)   :: fluidcode_output

    end type FACE_outputs
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
        write(*,*) '# arguments=',narg

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
            elseif (arg.eq."-vi" .or. arg.eq."--verbose-input") then
                    verbose_input=.true.
                    k=k+1

                    elseif (arg.eq."-vp" .or. arg.eq."--verbose_parser") then
                    verbose_parser=.true.
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


    subroutine wrapper_FACE(face_input,fluidcode_input)
        type(FACE_inputs),intent(out)::face_input
        type(fluidcode_inputs), intent(in) :: fluidcode_input
        face_input%casename=fluidcode_input%casename
        face_input%logfile="no"!face_input%casename//".log"
        face_input%input_filename=default_inputfile
        face_input%run_mode='default'
        face_input%path='solps_iter'

        face_input%fluidcode_input=fluidcode_input
    end subroutine wrapper_FACE


    subroutine set_fluidcode_input(fluidcode_input)

        type(fluidcode_inputs),intent(out)::fluidcode_input
        character(string_length)::str

        fluidcode_input%wall_idx=1
        fluidcode_input%iter=2
        fluidcode_input%time=0d0
        fluidcode_input%Ndump_space=10
        fluidcode_input%Ndump_time=10
        fluidcode_input%dt=1e-3
        fluidcode_input%dt_face=1e-4
        fluidcode_input%nspc_fluid=1

        allocate(fluidcode_input%namespc(fluidcode_input%nspc_fluid))
        fluidcode_input%namespc(1:fluidcode_input%nspc_fluid)="D"
        allocate(fluidcode_input%indexspc(fluidcode_input%nspc_fluid))
        fluidcode_input%indexspc(1:fluidcode_input%nspc_fluid)=1
        allocate(fluidcode_input%inflx_in(fluidcode_input%nspc_fluid))
        fluidcode_input%inflx_in(1:fluidcode_input%nspc_fluid)=1e20
        allocate(fluidcode_input%Emean(fluidcode_input%nspc_fluid))
        fluidcode_input%qflx_in=0           ! Heat flux from fluid code
        fluidcode_input%restore_state_file="no"    ! restart from this state file
        fluidcode_input%final_state_file="no"    ! save from state file history
        fluidcode_input%tempwall=700        ! temperature of the wall from fluid code
        fluidcode_input%solve_heat_eq="no"   ! if solve_heat_eq then use Qin otherwise T=tempwall for the entire bulk
        fluidcode_input%casename='solps'
        write(str,'(a,a1,i0,a1,i0)') trim(fluidcode_input%casename),'_',(fluidcode_input%wall_idx),'_',fluidcode_input%iter
        fluidcode_input%casename=trim(str)
    !    output
    end subroutine set_fluidcode_input

!    subroutine allocate_fluidcode_output
!    type trace_outgassing_flx
!        real(DP),allocatable         :: Gdes(:)
!        real(DP),allocatable         :: min_Gdes(:)
!        real(DP),allocatable         :: max_Gdes(:)
!        real(DP),allocatable         :: ave_Gdes(:)    ! ave deviation of Gdes over FACE run
!        real(DP),allocatable         :: sig_Gdes(:)    ! sdt deviation of Gdes over FACE run
!        real(DP),allocatable         :: Gpermeation(:) ! =Gdes_r(:)
!     end type trace_outgassing_flx
!
!    type fluidcode_output
!        real(DP),allocatable   :: out_flux(:)            ! outgassing flux
!        real(DP),allocatable   :: Ntot(:)                ! total species in material
!        real(DP),allocatable   :: Nnet_dt_fc(:)          ! total species added during dt_fx
!        real(DP),allocatable   :: Nnet_free_dt_fc(:)     ! total free species added during dt_fx
!        real(DP),allocatable   :: Nnet_trapped_dt_fc(:)  ! total trapped species added during dt_fx
!        real(DP),allocatable   :: Nnet_srf_dt_fc(:)      ! total species added on surface during dt_fx
!        real(DP),allocatable   :: Nnet_vol_dt_fc(:)      ! total species added in bulk during dt_fx
!    end type fluidcode_output
!    end subroutine allocate_fluidcode_output


end module modFACE_interface
