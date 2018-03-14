    module modFACE_output
    use modFACE_header
     implicit none
     save
    integer:: iout=6
    integer:: ifile_heatdata
    integer :: ifile_voldata
    integer:: ifile_surfdata
    integer::ifile_log
    integer:: ifile_inputlog
    integer::current_ifile=10
    integer:: ifile_testpath

    public::iout,timestring
    contains

        subroutine init_log(logfile_)
        character(*)::logfile_
        integer ios
            if (logfile_.eq."no") then
                iout=6 !default output unit
            elseif (logfile_.eq."yes") then
                call set_ifile(iout)
                open (iout, file=trim(path_folder)//'log_FACE', status='replace', iostat=ios)
                if (ios.ne.0) then
                    write (*, '(2a)') 'ERROR: Cannot open log file:', 'log_FACE'
                    stop
                endif
            else
                call set_ifile(iout)
                open (iout, file=trim(logfile_), status='replace', iostat=ios)
                if (ios.ne.0) then
                    write (*, '(2a)') 'ERROR: Cannot open log file:', logfile_
                    stop
                endif
            endif
            if(verbose_init) write(iout,*) "iout=",iout
        end subroutine init_log

        ! set unit numbers for various files opened by FACE
        subroutine init_ifile()

            call set_ifile(ifile_heatdata)
            if(verbose_init) write(iout,*) "ifile_heatdata=",ifile_heatdata

            call set_ifile(ifile_voldata)
            if(verbose_init) write(iout,*) "ifile_voldata=",ifile_voldata

            call set_ifile(ifile_inputlog)
            if(verbose_init) write(iout,*) "ifile_inputlog=",ifile_inputlog

            call set_ifile(ifile_testpath)
            if(verbose_init) write(iout,*) "ifile_testpath=",ifile_testpath

        end subroutine init_ifile

subroutine set_ifile(ifile)
      integer ifile
      logical unit_open
      unit_open = .true.
      current_ifile=10
      Do While (unit_open.and.current_ifile.le.max_ifile)
         current_ifile=current_ifile + 1
         Inquire (Unit = current_ifile, Opened = Unit_open)
      End Do
      if (current_ifile.lt.max_ifile) then
      ifile=current_ifile
      else
      write(iout,*) "ERROR: cannot find available unit number <", max_ifile
      stop
      endif

      return


   end subroutine set_ifile



    subroutine timestamp ( )

  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine timestring ( string )

!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = * ) string
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( string, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

subroutine open_timedata_file
      ! must add single file for all species
!       write (name, '(a, i2.2, a)') 'dat/time.dat'
!       open (unit=unt, file=path(1:lnblnk(path))//name,
!     +       status='replace', err=2510)
!
!  write (unit=unt,'14(a)') '         01:time',&
!       '         02:Temp',&
!       '           03:densL',&
!       '           04:densR',&
!       '           05:NsrfL',&
!       '           06:NsrfR',&
!       '           07:GdesL',&
!       '           08:GdesR',&
!       '            09:Qnty',&
!       '             10:Src',&
!       '             11:Rct',&
!       '             12:Rtd',&
!       '           13:inflx',&
!       '            14:qrad'&
!      enddo
end subroutine open_timedata_file

subroutine print_timestep_info
character(200)::myfmt
write(myfmt,*) "('iter=', 1I6,' time=', 1pe12.4e2, ' s;   T_l=', 1pe12.4e2, ' K;   T_r=', 1pe12.4e2, ' K;"&
 ,"dt=', 1pe12.4e2, ' s, |f|=',1pe12.4e2)"
 write (iout, myfmt) iteration,time, temp(ndt,0), temp(ndt,ngrd), dt_face,normf

end subroutine print_timestep_info

subroutine print_milestone(str)
character(*)::str
write(iout,"('--- ',a,' ---')") str
end subroutine print_milestone

    end module modFACE_output
