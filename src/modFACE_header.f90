!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE: Module Name
!
!> @author
!> Module Author Name and Affiliation
!
! DESCRIPTION:
!> Brief description of module.
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
      module modFACE_header
      use modFACE_precision
      implicit none
!     ******************************************************************
!     * This file contains all common variables specifications         *
!     * for 1-dimensional First Wall simulation Code "FACE"            *
!     *                                                                *
!     * Author: Roman D. Smirnov                                       *
!     * E-mail: rsmirnov@ucsd.edu; rosmirnov@yahoo.com                 *
!     *                                                                *
!     ******************************************************************
!     ------------------------------------------------------------------
!      Array dimensions and grid parameters
!     ------------------------------------------------------------------
       logical::verbose_parser=.false.
           logical::verbose_input=.false.
    logical::verbose_init=.false.
    logical :: verbose_step=.false.
    logical :: verbose_debug=.false.
    logical :: verbose_couple=.false.
    logical :: verbose_restore=.false.
     logical :: verbose_maths=.false.
     logical :: verbose_interface=.false.
     logical :: verbose_help=.false.
     logical :: verbose_version=.false.
          logical :: verbose_header=.false.
          logical :: verbose_surface=.false.
    logical:: enforce_error=.true.


      integer::ngrd
      integer::nspc
      integer::neq
      integer:: nramp
    integer:: iout=6
      real(DP) :: tcpustart, tcpufinish,walltime_start,walltime_end
     integer :: Nprint_run_info=2 ! print info on current run every Nprint_run_info steps
     real(DP),parameter :: min_rate_surface=1.0d-20
      ! **  Some physical and mathematical constants
      real(DP),parameter ::ee=1.602176462d-19
      real(DP),parameter ::eps0=8.854187817d-12
      real(DP),parameter ::amass=1.66053886d-27
      real(DP),parameter ::pi=3.14159265358979d0
      real(DP),parameter ::twopi=6.28318530717959d0
      real(DP),parameter ::sqrt2=sqrt(2.d0)
      real(DP),parameter ::kb=1.3806504d-23
      real(DP),parameter ::eekb=1.160450595d+04
      real(DP),parameter ::sigma_sb=5.670400d-08

      ! ** coefficients for BDF
      !     --- 1st order BDF ---
      real(DP),parameter ::a11=1.d0
      real(DP),parameter ::a12=1.d0
      !     --- 2nd order BDF ---
      real(DP),parameter ::a21= 4.d0/3.d0
      real(DP),parameter ::a22=-1.d0/3.d0
      real(DP),parameter ::a23= 2.d0/3.d0
      !     --- 5th order BDF ---
      real(DP),parameter ::a51= 48.d0/25.d0
      real(DP),parameter ::a52=-36.d0/25.d0
      real(DP),parameter ::a53= 16.d0/25.d0
      real(DP),parameter ::a54=- 3.d0/25.d0
      real(DP),parameter ::a55= 12.d0/25.d0

      real(DP) ::solver_eps=3.d-3, solver_udspl=9.d-1, solver_fdspl=9.d0, solver_gdspl=1.d-3,jac_eps=1d-8
      real(DP) ::solver_fstp=1.d-1
      integer :: iter_solver_max=150
      logical :: finalcheck=.true.

      character(string_length) :: default_inputfile="default_inputfile.face"
!


!     ------------------------------------------------------------------
!      Spatial and temporal parameters
!     ------------------------------------------------------------------
      !numeric
      integer :: ndt ! >@var size of storage for time dependent variables
      ! grid
      real(DP):: length !length spatial domain
      real(DP),allocatable:: x(:)
      real(DP),allocatable:: dx(:)
      real(DP) ::alpha
      real(DP) ::grid_dx0
      character(string_length) ::grid_type
      character(string_length) ::grid_gen_mode
      ! time
      real(DP):: dt_face    !>@var current solver time step
      real(DP):: min_dt_face   !>@var current solver time step
      real(DP):: max_dt_face   !>@var current solver time step
      real(DP):: dt0_face ! nominal time step
      real(DP) :: reduction_factor_dt=1d0
      logical :: adjust_reduction_factor=.false.
      character(string_length):: adjust_reduction_factor_string
      integer:: Nstep_increase_dt=10
      integer :: solver_step_count=0
      real(DP) :: max_iter=1e9
      real(DP):: end_time  ! end time of simulations
      real(DP):: time  ! current time of simulations
      real(DP):: time_savevol  !
      real(DP):: time_savetime  !
      real(DP):: start_time ! start time  of simulation
!      real(DP):: cdt   ! factor for solver time step (dt)
      ! temp
      real(DP):: cero     ! erosion velocity
      real(DP):: cero_min
      real(DP):: cero_max
      real(DP):: gamero


      real(DP):: tramp0=0
      real(DP):: tramp1=0

      ! parameters controlling data dumping
      real(DP):: dump_space_dt=0
      real(DP):: dump_time_dt=0
      real(DP):: dump_restart_dt=0
      logical:: dump_space=.true.
      logical:: dump_time=.true.
      logical:: dump_restart=.true.

      real(DP):: nucut=1d99
      real(DP):: delta=0

      integer:: iter_solver=0 ! #of solver iterations at each time step
      integer::iteration=0
      integer:: order_solver

!       Save file numerations
      integer :: max_ifile
      parameter(max_ifile=10000)
      integer:: sfln_voldata=0
      integer:: sfln_srfdata=0
      integer:: sfln_heatdata=0
      real(DP):: normf=0.d0
      character(string_length)::restart_filename


!       Flags
      logical:: read_input_file=.true.
      logical:: restore_state_temp=.true.
      logical::dump_vol_append=.false.
      logical::dump_srf_append=.false.
      logical::dump_time_append=.false.
      logical:: solve_heat_eq=.false.
      logical:: variable_timestep
      integer::avr
      character(string_length):: input_filename
      character(string_length):: logfile
      character(string_length):: read_restart_file
      character(string_length):: read_state_file
      character(string_length):: steady_state
      character(string_length):: framp
      character(string_length):: solve_heat_eq_string
      character(string_length):: variable_timestep_string
      character(string_length):: final_state_file
      character(string_length):: casename

      type inventories
       real(DP)         ::  Nnetbulk
        real(DP)         :: Nnetsrf
        real(DP)         :: Ntotbulk
        real(DP)         :: Ntotsrf
     end type inventories

     type particle_balances
     real(DP):: Nnet=0d0,Ninflux=0d0,Noutflux=0d0,p_net=0d0,p_max=0d0,f_lost=0d0,Noutflux_l=0d0,Noutflux_r=0d0
     end type particle_balances

    type outgassing_fluxes
        real(DP)         :: Gdes
        real(DP)         :: min_Gdes
        real(DP)         :: max_Gdes
        real(DP)         :: ave_Gdes    ! ave deviation of Gdes over FACE run
        real(DP)         :: sig_Gdes   ! sdt deviation of Gdes over FACE run
        real(DP)         :: Gpermeation ! =Gdes_r
    end type outgassing_fluxes

     type wall_temperatures
        real(DP)         :: srf_temp_l,srf_temp_r
        real(DP)         :: mean_temp
        real(DP)         :: max_temp
        real(DP)         :: min_temp
    end type wall_temperatures

     type(inventories),allocatable :: init_inventory(:)
     type(inventories),allocatable :: final_inventory(:)
     type(particle_balances):: particle_balance
     type(outgassing_fluxes):: outgassing_flux
     type(wall_temperatures):: init_wall_temp
      type(wall_temperatures)::final_wall_temp

!
!     ------------------------------------------------------------------
!       Material parameters
!     ------------------------------------------------------------------

      real(DP):: temp0=300d0
      real(DP):: temp1=300d0
      real(DP):: dtemp=0d0
      real(DP):: lambda=0d0
      !
      real(DP):: cvlm=1d0
      real(DP):: csrf=1d0
      real(DP):: clng=1d0

!     implantation parameters
      character(lname),allocatable::implantation_model(:)
      real(DP),allocatable::implantation_depth(:)
      real(DP),allocatable::diagnostic_depth(:)
      integer,allocatable::j_implantation_depth(:)
      integer,allocatable::j_diagnostic_depth(:)
      real(DP),allocatable::implantation_width(:)
      real(DP),allocatable::enrg(:)
      real(DP),allocatable::inflx(:) ! influx of particles (may differ from nominal particle flux in pulsed_plasma mode)
      real(DP),allocatable::inflx_in_max(:) ! max influx of particles
      real(DP),allocatable::inflx_in(:) ! nominal influx of particles
      real(DP),allocatable::inflx_in_pulse_period(:) ! max influx of particles
      real(DP),allocatable::inflx_in_pulse_starttime(:) ! max influx of particles
      character(lname),allocatable::inflx_in_pulse(:) ! nominal influx of particles
      real(DP),allocatable::gas_pressure(:)
      real(DP),allocatable:: gas_temp(:)
      real(DP),allocatable:: mass(:)
      real(DP),allocatable:: temp(:,:)
      real(DP),allocatable:: rtime(:)
      real(DP),allocatable:: rtemp(:)
!
!     ------------------------------------------------------------------
!       Species parameters
!     ------------------------------------------------------------------
!
      character(lname),allocatable::namespc(:)
      character(lname),allocatable::left_surface_model(:)
      character(lname),allocatable::right_surface_model(:)

!      Volumetric species terms
      real(DP),allocatable::dens0(:)
      real(DP),allocatable::gxmax(:)
      real(DP),allocatable::gsigm(:)

      real(DP),allocatable::cdif0(:)
      real(DP),allocatable::edif (:)
      real(DP),allocatable::flx (:,:,:)
      real(DP),allocatable::rate_d (:,:,:)
      character(lname),allocatable:: gprof(:)
      real(DP),allocatable::jout(:,:)



      real(DP),allocatable::dsrfl0(:)
      real(DP),allocatable::dsrfr0(:)
      real(DP),allocatable::dsrfm(:)
      real(DP),allocatable::densm(:)
      real(DP),allocatable::etr   (:)
      real(DP),allocatable::edtr  (:)


!     Boundary species parameters

      real(DP),allocatable::order_desorption_left(:)
      real(DP),allocatable::order_desorption_right(:)
      real(DP),allocatable::Eabs_l(:)
      real(DP),allocatable::Edes_l(:)
      real(DP),allocatable::Eb_l(:)
      real(DP),allocatable::Eads_l(:)

      real(DP),allocatable::Eabs_r(:)
      real(DP),allocatable::Edes_r(:)
      real(DP),allocatable::Eb_r(:)
      real(DP),allocatable::Eads_r(:)

      real(DP),allocatable::Kabs_l(:)
      real(DP),allocatable::Kdes_l(:)
      real(DP),allocatable::Kb_l(:)
      real(DP),allocatable::Kads_l(:)
      real(DP),allocatable::Kabs_r(:)
      real(DP),allocatable::Kdes_r(:)
      real(DP),allocatable::Kb_r(:)
      real(DP),allocatable::Kads_r(:)


      real(DP),allocatable::K0abs_l(:)
      real(DP),allocatable::K0des_l(:)
      real(DP),allocatable::K0b_l(:)
      real(DP),allocatable::K0ads_l(:)

      real(DP),allocatable::K0abs_r(:)
      real(DP),allocatable::K0des_r(:)
      real(DP),allocatable::K0b_r(:)
      real(DP),allocatable::K0ads_r(:)

      real(DP),allocatable::nu (:)
      real(DP),allocatable::j0 (:)
      real(DP),allocatable::dens(:,:,:)
!
!     ------------------------------------------------------------------
!       Reaction parameters
!     ------------------------------------------------------------------
      real(DP),allocatable:: nuth (:,:)
      real(DP),allocatable:: kbin (:,:,:)
      real(DP),allocatable::nuth0(:,:)
      real(DP),allocatable:: kbin0(:,:,:)
      real(DP),allocatable::eth  (:,:)
      real(DP),allocatable:: ebin (:,:,:)
!
!     ------------------------------------------------------------------
!       Species variables
!     ------------------------------------------------------------------
!     Bulk
      real(DP),allocatable::srs(:,:,:)
      real(DP),allocatable::src_profile(:,:)
      real(DP),allocatable::srb (:,:,:,:)
      real(DP),allocatable::src (:,:,:)
      real(DP),allocatable::cdif(:,:,:)
      real(DP),allocatable::rct (:,:,:)
      real(DP),allocatable::ero_flx (:,:,:)  ! erosion flux
      real(DP),allocatable::dif_flx (:,:,:)  ! erosion flux
!     Surface
      integer :: error_status=0


      ! left surface at x(j=0)
      real(DP),allocatable:: dsrfl(:,:)   ! density on left surface
      real(DP),allocatable:: Gsrf_l(:,:) ! net flux of species onto the left surface
      real(DP),allocatable:: Gabs_l(:,:)
      real(DP),allocatable:: Gdes_l(:,:)
      real(DP),allocatable:: Gb_l(:,:)
      real(DP),allocatable:: Gads_l(:,:)

      !right surface at x(j=n+1)
      real(DP),allocatable:: dsrfr(:,:)   ! density on right surface
      real(DP),allocatable:: Gsrf_r(:,:) ! net flux of species onto the right surface
      real(DP),allocatable:: Gabs_r(:,:)
      real(DP),allocatable:: Gdes_r(:,:)
      real(DP),allocatable:: Gb_r(:,:)
      real(DP),allocatable:: Gads_r(:,:)

      type trace_fluxes
      real(DP):: sum_inflx
      real(DP):: sum_Gdes_l
      real(DP):: sig_Gdes_l
      real(DP):: min_Gdes_l
      real(DP):: max_Gdes_l
      real(DP):: sum_Gdes_r
      real(DP):: sig_Gdes_r
      real(DP):: min_Gdes_r
      real(DP):: max_Gdes_r
      end type trace_fluxes
      type(trace_fluxes),allocatable :: trace_flux(:)

         type onthefly_inventories
      real(DP):: int_dens
      real(DP):: int_dsrf
      real(DP):: net_int_dens
      real(DP):: net_int_dsrf
      real(DP):: int_des
      real(DP):: int_src
      end type onthefly_inventories
      type(onthefly_inventories),allocatable :: onthefly_inventory(:)
      logical :: active_cap=.false.
      logical :: print_onthefly_inventory=.false.
      character(string_length)::print_onthefly_inventory_string
      character(string_length)::active_cap_string
      character(string_length)::dump_vol_append_string
      character(string_length)::dump_srf_append_string
!     ------------------------------------------------------------------
!       Thermal variables
!     ------------------------------------------------------------------
      real(DP):: thcond  =0
      real(DP):: rho     =0
      real(DP):: cp      =0
      real(DP)::rhocp    =0
      real(DP):: emiss   =0
      real(DP)::qform    =0
      real(DP):: qflx_in =0 ! incoming heat flux from plasma
      real(DP):: rad     =0
      real(DP)::rad_min  =0
      real(DP):: rad_max =0
      real(DP)::t1       =0
      real(DP):: t2      =0
      real(DP)::t3       =0
!      real(DP)::tpulse   =0 ! period of plasma pulse
      real(DP),allocatable::rate_t(:,:)
      real(DP),allocatable::qflx(:,:)
      real(DP),allocatable::ero_qflx(:,:)

!
!     ------------------------------------------------------------------
!       Environment variables
!     ------------------------------------------------------------------
      character(string_length) path_folder ! top folder where simulations files anf folders are written in
      character(string_length) dat_folder ! top folder where vol,srf and heat data files are written in

type fluidcode_inputs
         integer                     :: wall_idx            ! Index of the wall stratum
        integer                      :: iter                ! Fluid code iteration
        integer                      :: Ndump_space         ! # space data files to be dumped
        integer                      :: Ndump_time          !  # of times time data are dumped
        logical                      :: append              ! append mode for dumping data
        real(DP)                     :: time                ! Fluid code time
        real(DP)                     :: dt                  ! Time step of the fluid
        real(DP)                     :: dt0_face             ! Time step of FACE
        integer                      :: nspc_fluid                 ! Number of incoming species from fluid code
        integer,allocatable          :: indexspc(:)             ! Index of species in FACE (usually "k" in FACE)
        character(Lname),allocatable :: namespc(:)      ! Name of the incoming species from fluid code
        real(DP),allocatable         :: inflx_in(:)      ! Particle flux'
        real(DP),allocatable         :: Emean(:)        ! Average energy of incoming particle enrg'
        real(DP)                     :: qflx_in                 ! Heat flux from fluid code
        character(string_length)     :: read_state_file  ! restart from this staste file
        character(string_length)     :: final_state_file    ! store final state in this state file
        real(DP)                     :: tempwall            ! temperature of the wall from fluid code
        character(15)                :: solve_heat_eq_string       ! if solve_heat_eq then use Qin otherwise T=tempwall for the entire bulk
        character(string_length)      ::casename            ! casename (see below)
        character(string_length)      ::casename_base       ! base to form casename=casename_base_iteration_idx_wall
        character(string_length)      ::input_file            ! casename
        character(string_length)      ::log_file            ! casename
        character(string_length)      ::path            ! casename
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
         integer :: nspc
        type(outgassing_fluxes):: outgassing_flux
             type(wall_temperatures) :: init_wall_temp,final_wall_temp
             type(particle_balances)  :: particle_balance
             type(inventories),allocatable       :: init_inventory(:),final_inventory(:)

    end type FACE_outputs



end module modFACE_header
