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
call cpu_time(tcpustart)
call init_casename(face_input%casename)
call init_path(face_input%path)
! initialize help: required to check keywords in input file
call init_help()
! Initialize log output (file or standard output (unit=6))
call init_log(face_input%logfile)
! Apply requested run mode
if (face_input%run_mode.eq."help") then
call display_help
elseif (face_input%run_mode.eq."version") then
call display_version
elseif (face_input%run_mode.eq."print_default_input") then
call display_version()
call write_default_inputfile(default_inputfile,"H+Tr")
else ! all other mode request reading of input

      call display_version()
            ! default modes-> print out default input file then read file
      if (face_input%run_mode.eq."default_H") then
      call write_default_inputfile(default_inputfile,"H")
      elseif (face_input%run_mode.eq."default_H+Tr") then
      call write_default_inputfile(default_inputfile,"H+Tr")
      endif
      !read input
      call get_input(face_input)

      ! initialize variables
      call initialize()
      ! restore from previous state (restart or history)
      call restore()
      ! execute time loop
      call time_loop()
      ! dump final state
     call output_final_state()
     !
endif

! Final message if execution was successful

call cpu_time(tcpufinish)

call print_summary
call finalize

!      integer k, unt
!
!cc      open (unit=6, form='formatted', carriagecontrol='fortran')
!      open (unit=6, form='formatted')
!mnt      call getenv('FACE_PATH', path)
!      path='./'
!      call input()
!      call init()
!      call restore()
!      write (6, 1005) time, temp(ndt,0), temp(ndt,ngrd), dt
!1005  format ('  time=', 1pe19.9e4, ' s;'
!     +        '   T_l=', 1pe19.9e4, ' K;'
!     +        '   T_r=', 1pe19.9e4, ' K;'
!     +        '    dt=', 1pe19.9e4, ' s')
!1010  if (time .le. end_time) then
!       call save()
!       call step()
!       goto 1010
!      endif
!      call store()
!      do 1020, k=1,nspc
!       unt=100+k
!       close (unt)
!1020  continue
!      close (6)
      stop
  end    subroutine FACE_main

        subroutine time_loop
      call print_milestone('starting time iteration')
      iteration=0
      do while (time .le. end_time)
      call print_timestep_info()
      call save()
      call step()
      call store()
      iteration=iteration+1
      enddo
      call print_milestone('time iteration completed')
      end subroutine time_loop
end module modFACE_main
