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
    logical :: verbose_step=.true.
    logical :: verbose_debug=.false.
    logical :: verbose_couple=.false.
    logical :: verbose_restore=.true.
    logical:: enforce_error=.true.
!      integer ngrdm, nspcm, ndt, nrampm
      integer::ngrd
      integer::nspc
      integer::neq
      integer:: nramp
    integer:: iout=6
      real :: tcpustart, tcpufinish

      ! coefficients for BDF
       real(DP):: a11, a12, a21,a22,a23,a51,a52, a53,a54, a55
        !     --- 1st order BDF ---
        parameter (a11=1.d0)
        parameter (a12=1.d0)
        !     --- 2nd order BDF ---
        parameter (a21= 4.d0/3.d0)
        parameter (a22=-1.d0/3.d0)
        parameter (a23= 2.d0/3.d0)
        !     --- 5th order BDF ---
              parameter (a51= 48.d0/25.d0)
              parameter (a52=-36.d0/25.d0)
              parameter (a53= 16.d0/25.d0)
              parameter (a54=- 3.d0/25.d0)
              parameter (a55= 12.d0/25.d0)
!
!      parameter ( ngrdm=1000, nspcm=40, ndt=5, nrampm=10000 )
!
!     ------------------------------------------------------------------
!      Some physical and mathematical constants
!     ------------------------------------------------------------------
      real(DP)::ee, eps0, pi, twopi, sqrt2, amass, kb, eekb, sigma_sb


      parameter (       ee=1.602176462d-19    )
      parameter (     eps0=8.854187817d-12    )
      parameter (    amass=1.66053886d-27     )
      parameter (       pi=3.14159265358979d0 )
      parameter (    twopi=6.28318530717959d0 )
      parameter (    sqrt2=sqrt(2.d0)         )
      parameter (       kb=1.3806504d-23      )
      parameter (     eekb=1.160450595d+04    )
      parameter ( sigma_sb=5.670400d-08       )
      logical :: finalcheck=.true.
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
      ! time
      real(DP):: dt_face    !>@var current solver time step
      real(DP):: dtmin ! minimun dt
      real(DP):: end_time  ! end time of simulations
      real(DP):: time  ! current time of simulations
      real(DP):: start_time ! start time  of simulation
      real(DP):: cdt   ! factor for solver time step (dt)
      ! temp
      real(DP):: cero     ! erosion velocity
      real(DP):: cero_min
      real(DP):: cero_max
      real(DP):: gamero


      real(DP):: tramp0=0
      real(DP):: tramp1=0

      real(DP):: tspc=0
      real(DP):: tstr=0
      real(DP):: ttm=0

      real(DP):: nucut=0
      real(DP):: delta=0

      integer:: iter_solver=0 ! #of solver iterations at each time step
      integer::iteration
      integer:: order_solver
!     ------------------------------------------------------------------
!       Save file numerations
!     ------------------------------------------------------------------
      integer :: max_ifile
      parameter(max_ifile=10000)
      integer:: sfln_voldata=0
      integer:: sfln_srfdata=0
      integer:: sfln_heatdata=0
      real(DP):: normf=0.d0
      character(string_length)::restart_filename
!
!     ------------------------------------------------------------------
!       Flags
!     ------------------------------------------------------------------
      logical:: read_input_file=.true.
      logical:: restore_state_temp=.true.
      integer::avr
      character(string_length):: input_filename
      character(string_length):: logfile
      character(string_length):: read_restart_file
      character(string_length):: read_state_file
      character(string_length):: steady_state
      character(string_length):: framp
      character(string_length):: solve_heat_eq
      character(string_length):: final_state_file
      character(string_length):: restore_state_file
      character(string_length):: casename
      character(string_length):: pulsed_flux

      type inventories
       real(DP)         ::  Nnetbulk
        real(DP)         :: Nnetsrf
        real(DP)         :: Ntotbulk
        real(DP)         :: Ntotsrf
     end type inventories

     type(inventories),allocatable :: init_inventory(:)
     type(inventories),allocatable :: final_inventory(:)

!
!     ------------------------------------------------------------------
!       Material parameters
!     ------------------------------------------------------------------

      real(DP):: temp0=0
      real(DP):: temp1=0
      real(DP):: dtemp=0
      real(DP):: lambda=0
      real(DP):: cvlm=0
      real(DP):: csrf=0
      real(DP):: clng=0
      real(DP):: lambda1c=0
      real(DP):: lambda2c=0
      real(DP):: lambda3c=0

!     species parameters
      real(DP),allocatable::enrg(:)
      real(DP),allocatable::inflx(:) ! influx of particles (may differ from nominal particle flux in pulsed_plasma mode)
      real(DP),allocatable::inflx_in_max(:) ! max influx of particles
      real(DP),allocatable::inflx_in(:) ! nominal influx of particles
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
      real(DP),allocatable::echl(:)
      real(DP),allocatable::qchl(:)
      real(DP),allocatable::ebl(:)
      real(DP),allocatable::esl(:)
      real(DP),allocatable::echr(:)
      real(DP),allocatable::qchr(:)
      real(DP),allocatable::ebr(:)
      real(DP),allocatable::esr(:)
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
      real(DP),allocatable::qchtl(:)
      real(DP),allocatable::qchtr(:)
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
      real(DP),allocatable::srb (:,:,:,:)
      real(DP),allocatable::src (:,:,:)
      real(DP),allocatable::cdif(:,:,:)
      real(DP),allocatable::rct (:,:,:)
      real(DP),allocatable::ero_flx (:,:,:)  ! erosion flux
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
!
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
      real(DP)::tpulse   =0 ! period of plasma pulse
      real(DP),allocatable::rate_t(:,:)
      real(DP),allocatable::qflx(:,:)
      real(DP),allocatable::ero_qflx(:,:)
!
!     ------------------------------------------------------------------
!       Environment variables
!     ------------------------------------------------------------------
      character(string_length) path_folder ! top folder where simulations files anf folders are written in
      character(string_length) dat_folder ! top folder where vol,srf and heat data files are written in



end module modFACE_header
