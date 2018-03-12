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
!      integer ngrdm, nspcm, ndt, nrampm
      integer::ngrd
      integer::nspc
      integer::neq
      integer:: nramp
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
      real(DP):: cero
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

      integer:: cnt=0
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
      real(DP):: normf=0d0
!
!     ------------------------------------------------------------------
!       Flags
!     ------------------------------------------------------------------
      logical:: read_input_file=.true.
      integer::avr
      character(Linput):: input_filename
      character(Linput):: logfile
      character(Linput):: restart_mode
      character(Linput):: stdst
      character(Linput):: framp
      character(Linput):: solve_heat_eq
      character(Linput):: store_history_file
      character(Linput):: restore_history_file
      character(Linput):: casename
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
      real(DP),allocatable::inflx(:)
      real(DP),allocatable::inflx_max(:)
      real(DP),allocatable::inflx_min(:)
      real(DP),allocatable::prg(:)
      real(DP),allocatable:: tg(:)
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
      real(DP),allocatable::rtd (:,:,:)
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
      real(DP),allocatable::k1l(:)
      real(DP),allocatable::k2l(:)
      real(DP),allocatable::k3l(:)
      real(DP),allocatable::k4l(:)
      real(DP),allocatable::k1r(:)
      real(DP),allocatable::k2r(:)
      real(DP),allocatable::k3r(:)
      real(DP),allocatable::k4r(:)
      real(DP),allocatable::r1l(:)
      real(DP),allocatable::r2l(:)
      real(DP),allocatable::r3l(:)
      real(DP),allocatable::r4l(:)
      real(DP),allocatable::r1r(:)
      real(DP),allocatable::r2r(:)
      real(DP),allocatable::r3r(:)
      real(DP),allocatable::r4r(:)
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
      real(DP),allocatable::ero (:,:,:)
!     Surface
      real(DP),allocatable::dsrfl(:,:)
      real(DP),allocatable:: rtsl(:,:)
      real(DP),allocatable::dsrfr(:,:)
      real(DP),allocatable:: rtsr(:,:)
      real(DP),allocatable::j1l(:,:)
      real(DP),allocatable:: j2l(:,:)
      real(DP),allocatable:: j3l(:,:)
      real(DP),allocatable:: j4l(:,:)
      real(DP),allocatable::j1r(:,:)
      real(DP),allocatable:: j2r(:,:)
      real(DP),allocatable::j3r(:,:)
      real(DP),allocatable::j4r(:,:)

!
!     ------------------------------------------------------------------
!       Thermal variables
!     ------------------------------------------------------------------
      real(DP):: thcond=0
      real(DP):: rho=0
      real(DP):: cp=0
      real(DP)::rhocp=0
      real(DP):: emiss=0
      real(DP)::qform=0
      real(DP):: qflx=0
      real(DP):: rad=0
      real(DP)::rad_min=0
      real(DP):: rad_max=0
      real(DP)::t1=0
      real(DP):: t2 =0
      real(DP)::t3=0
       real(DP)::tp=0
      real(DP),allocatable::rtt(:,:)
      real(DP),allocatable::flxt(:,:)
      real(DP),allocatable::erot(:,:)
!
!     ------------------------------------------------------------------
!       Environment variables
!     ------------------------------------------------------------------
      character *128 path_folder
end module modFACE_header
