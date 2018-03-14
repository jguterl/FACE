module modFACE_IO
    use modFACE_header
    use modFACE_allocate
    use modFACE_step
    use modFACE_coupling
    implicit none
    save
    integer,allocatable::ifile_timedata(:)
contains

    subroutine open_timedata_files()
        integer:: ios,k
        character(300)::name
        ! open time file
        allocate(ifile_timedata(nspc))

        do k=1,nspc
            call set_ifile(ifile_timedata(k))
            write (name, '(a,a,a, i2.2, a)') trim(path_folder),trim(casename),'_dat/time_', k, '.dat'
            open (ifile_timedata(k), file=trim(name),status='replace', iostat=ios)
            if (ios.ne.0) then
                write (iout,*) ' *** Cannot open file ', trim(name)
                stop 'Exiting FACE...'
            endif
            write (ifile_timedata(k), '(14a13)')&
                '    time[s]',&
                '    Temp[K]',&
                '      densL',&
                '      densR',&
                '      NsrfL',&
                '      NsrfR',&
                '      GdesL',&
                '      GdesR',&
                '       Qnty',&
                '        Src',&
                '        Rct',&
                '        Rtd',&
                '     inflx1',&
                '       qrad'
        enddo
    end subroutine open_timedata_files

   subroutine close_timedata_files()
        integer:: k
        logical:: test
do k=1,nspc
        inquire(ifile_timedata(k),opened=test)
        write(iout,*) test
        close(ifile_timedata(k))
          enddo
          deallocate(ifile_timedata)

    end subroutine close_timedata_files

    subroutine save_voldata()

        integer j, k,ios
        character*256  name,myfmt1,myfmt2,myfmt3

        !     --- saving snapshot of volume distributions ---

        write(myfmt1,*)    "(a, 1pe13.4e2, a)"
        write(myfmt2,*) "(9a13)"
        write(myfmt3,*) "(i13.4, 8(1pe13.4e2))"

        do k=1,nspc
            write (name, '(a,a,a, i2.2, a, i3.3, a)') trim(path_folder),trim(casename),'dat/vol', k, '_', sfln_voldata, '.dat'
            open (unit=ifile_voldata, file=name, status='replace', iostat=ios)
            if (ios.eq.0) then

                write (ifile_voldata,myfmt1) 'time=', time, ' s'

                write (ifile_voldata,myfmt2) &
                    ' #', &
                    '  x',&
                    ' dens',&
                    '   flx',&
                    '  ero',&
                    '   src',&
                    '  rct',&
                    '   rtd',&
                    '  Cdif'

                do j=0,ngrd
                    write (ifile_voldata, fmt=myfmt3) &
                        j, &
                        x(j), &
                        dens(ndt,j,k), &
                        flx (ndt,j,k), &
                        ero (ndt,j,k), &
                        src (ndt,j,k), &
                        rct (ndt,j,k), &
                        rtd (ndt,j,k), &
                        cdif(ndt,j,k)
                enddo

                close (ifile_voldata)


            else
                write (iout, '(2a)') ' *** error saving ', name

            endif
            sfln_voldata=sfln_voldata+1
        enddo

    end subroutine save_voldata


    subroutine save_heatdata()

        integer j, ios
        character*256  name,myfmt1,myfmt2,myfmt3

        write(myfmt1,*) "(a, 1pe13.4e2, a)"
        write(myfmt2,*) "(6(a))"
        write(myfmt3,*)"(i8.4, 5(1pe13.4e2))"

        write (name, '(a,a,a, i3.3, a)') trim(path_folder),trim(casename),'dat/heat_', sfln_heatdata, '.dat'
        open (ifile_heatdata, file=name, status='replace',iostat=ios)
        if (ios.eq.0) then
            write (ifile_heatdata, myfmt1) 'time=', time, ' s'

            write (ifile_heatdata, myfmt2)'(',&
                '       #',&
                '                  x',&
                '               temp',&
                '               flxt',&
                '                rtt',&
                '               erot)'

            do j=0,ngrd
                write (ifile_heatdata, myfmt3)&
                    j,&
                    x(j),&
                    temp(ndt,j),&
                    flxt(ndt,j),&
                    rtt (ndt,j),&
                    erot(ndt,j)
            enddo

            close (ifile_heatdata)
        else
        write (iout, '(2a)') ' *** error saving ', name
    endif

    sfln_heatdata=sfln_heatdata+1

end subroutine save_heatdata

subroutine save
    real(DP)::tmp
    !     --- saving time file ---
    tmp=2.d0*abs(time-ttm*nint(time/ttm))
    if (tmp .lt. dt_face.or.time.eq.0.0d0) then
        call save_timedata

    endif

    tmp=2.d0*abs(time-tspc*nint(time/tspc))
    if (tmp .lt. dt_face.or.time.eq.0.0d0) then
        call save_voldata
        call save_srfdata
        call save_heatdata
    endif

    !     --- saving restart file ---
    tmp=2.d0*abs(time-tstr*nint(time/tstr))
    if ((tmp .lt. dt_face) .and. (time .ne. 0.d0)) then
    call store_restart(trim(restart_filename))
endif
end subroutine save

subroutine save_timedata
    integer::j,k
    real(DP):: qnty,frmn,rctn
    character*256  myfmt1,myfmt2

    write(myfmt1,*) &
        "('+', ' time=', 1pe13.4e2, ' s; T_l=', 1pe13.4e2,' K; dt=', 1pe13.4e2, ' s; number of iterations ', i3)"
    write(myfmt2,*) "(14(1pe13.4e2))"
    if (verbose_step) write (iout, myfmt1) time, temp(ndt,0), dt_face, cnt

    do k=1,nspc
        qnty=0.d0
        frmn=0.d0
        rctn=0.d0
        do j=0,ngrd-1
            qnty=qnty+0.5d0*(dens(ndt,j,k)+dens(ndt,j+1,k))*dx(j)
            frmn=frmn+0.5d0*(src (ndt,j,k)+src (ndt,j+1,k))*dx(j)
            rctn=rctn+0.5d0*(rct (ndt,j,k)+rct (ndt,j+1,k))*dx(j)
        enddo

        write (ifile_timedata(k),myfmt2)&
            time, &
            temp (ndt,0),&
            dens (ndt,0,   k),&
            dens (ndt,ngrd,k),&
            dsrfl(ndt,     k),&
            dsrfr(ndt,     k),&
            Gdes_l  (ndt,     k),&
            Gdes_r  (ndt,     k),&
            qnty,&
            frmn,&
            rctn,&
            rtd  (ndt,0,k),&
            inflx(1),&
            rad
    enddo

end subroutine save_timedata

subroutine save_srfdata
    character*256  name,myfmt1,myfmt2,myfmt3
    integer k,ios
    !     --- saving snapshot of surface parameters
    write (name, '(a,a,a, i3.3, a)') trim(path_folder),trim(casename),'dat/srf_', sfln_srfdata, '.dat'
    write(myfmt1,*) "(a, 1pe13.4e2, a)"
    write(myfmt2,*) "(12(a13))"
    write(myfmt3,*) "(i13.2, 11(1pe13.4e2))"
    open (ifile_surfdata, file=name, status='replace',iostat=ios)
    if (ios.eq.0) then
        write (ifile_surfdata, myfmt1) 'time=', time, ' s'
        write (ifile_surfdata, myfmt2)&
            'spc#',&
            'NsrfL',&
            'Gabs_l',&
            'Gdes_l',&
            'Gb_l',&
            'Gads_l',&
            'NsrfR',&
            'Gabs_r',&
            'Gdes_r',&
            'Gb_r',&
            'Gads_r',&
            'Jout'

        do k=1,nspc
            write (ifile_surfdata, myfmt3) &
                k, &
                dsrfl(ndt,k),&
                Gabs_l  (ndt,k),&
                Gdes_l  (ndt,k),&
                Gb_l  (ndt,k),&
                Gads_l  (ndt,k),&
                dsrfr(ndt,k),&
                Gabs_r  (ndt,k),&
                Gdes_r  (ndt,k),&
                Gb_r  (ndt,k),&
                Gads_r  (ndt,k),&
                jout (ndt,k)
        enddo

        close (ifile_surfdata)

    else
        write (iout, '(2a)') ' *** error saving ', name
    endif

    sfln_srfdata=sfln_srfdata+1

end subroutine save_srfdata


! store and restore state files
subroutine store_state(filename)
    character(*):: filename
    integer  ifile_store,j, k,ios

    call set_ifile(ifile_store)
    open (unit=ifile_store, file=filename,status='replace', form='formatted', iostat=ios)

    if (ios.eq.0) then
        do k=1,nspc
            do j=0,ngrd
                write (ifile_store) dens(ndt,j,k)
            enddo
        enddo
        do k=1,nspc
            write (ifile_store) dsrfl(ndt,k)
            write (ifile_store) dsrfr(ndt,k)
        enddo

        do j=0,ngrd

            write (ifile_store) temp(ndt,j)

        enddo
        close (ifile_store)

        write(iout,*) 'simulation state stored in '//filename
    else
        write (iout, '(a)') 'ERROR: cannot open state file: ', filename
        stop 'Exiting FACE'
    endif

end subroutine store_state

subroutine restore_state(filename)
    character(*):: filename
    integer ifile_restore,j, k, ios
     call set_ifile(ifile_restore)
    open(ifile_restore, file=trim(filename), iostat=ios,action='read',status='old',form='formatted')

    if ( ios /= 0 ) then
        write(iout,*) 'ERROR: Cannot open history store file :', trim(adjustl(filename))
        stop 'Exiting FACE'
    endif

    write(iout,*) 'Restoring state from file: ', trim(filename)


    do k=1,nspc
        do j=0,ngrd
            read (ifile_restore) dens(ndt,j,k)
        enddo
    enddo
    do k=1,nspc

        read (ifile_restore) dsrfl(ndt,k)
        read (ifile_restore) dsrfr(ndt,k)

    enddo
    do j=0,ngrd
        read (ifile_restore) temp(ndt,j)
    enddo
    close (ifile_restore)

end subroutine restore_state


! ***** store and restore restart files *****

subroutine store_restart(filename)
    character(*):: filename
    integer i, j, k,ios,ifile_restart
    call set_ifile(ifile_restart)
    open (unit=ifile_restart, file=filename,status='replace', form='unformatted', iostat=ios)
    if (ios.eq.0) then
        do k=1,nspc
            do j=0,ngrd
                do i=1,ndt
                    write (ifile_restart) dens(i,j,k), flx (i,j,k), ero (i,j,k),cdif(i,j,k), rct (i,j,k), rtd (i,j,k)
                enddo
            enddo
        enddo
        do j=0,ngrd
            write (ifile_restart) x(j)
        enddo
        do k=1,nspc
            do i=1,ndt
                write (ifile_restart) dsrfl(i,k), rtsl(i,k)
                write (ifile_restart) dsrfr(i,k), rtsr(i,k)
                write (ifile_restart) Gabs_l(i,k), Gdes_l(i,k), Gb_l(i,k), Gads_l(i,k)
                write (ifile_restart) Gabs_r(i,k), Gdes_r(i,k), Gb_r(i,k), Gads_r(i,k)
                write (ifile_restart) jout(i,k)
            enddo
        enddo
        do j=0,ngrd
            do i=1,ndt
                write (ifile_restart) temp(i,j), flxt(i,j), rtt(i,j), erot(i,j)
            enddo
        enddo
        write (ifile_restart) time
        write (ifile_restart) sfln_voldata,sfln_srfdata,sfln_heatdata
        close (ifile_restart)
    else
        write (iout, '(a)') 'ERROR: cannot read restart file ', filename
    endif
end subroutine store_restart

subroutine restore_restart(filename)
    character(*):: filename
    integer i, j, k, ios,ifile_restart
    call set_ifile(ifile_restart)
    open (unit=ifile_restart, file=filename, form='unformatted',iostat=ios,&
        action='read',status='old')
    if ( ios .ne. 0 ) then
        write(iout,*) 'ERROR: cannot open restart file : ',filename
        stop 'Exiting FACE'
    endif
    write(iout,*) 'Restoring state from restart file : ', filename

    do k=1,nspc
        do j=0,ngrd
            do i=1,ndt
                read (ifile_restart) dens(i,j,k), flx (i,j,k), ero(i,j,k),cdif(i,j,k), rct(i,j,k), rtd(i,j,k)
            enddo
        enddo
    enddo

    do j=0,ngrd
        read (ifile_restart) x(j)
    enddo

    do k=1,nspc
        do i=1,ndt
            read (ifile_restart) dsrfl(i,k), rtsl(i,k)
            read (ifile_restart) dsrfr(i,k), rtsr(i,k)
            read (ifile_restart) Gabs_l(i,k), Gdes_l(i,k), Gb_l(i,k), Gads_l(i,k)
            read (ifile_restart) Gabs_r(i,k), Gdes_r(i,k), Gb_r(i,k), Gads_r(i,k)
            read (ifile_restart) jout(i,k)
        enddo
    enddo

    do j=0,ngrd
        do i=1,ndt
            read (ifile_restart) temp(i,j), flxt(i,j), rtt(i,j), erot(i,j)
        enddo
    enddo

    read (ifile_restart) time
    read (ifile_restart) sfln_voldata,sfln_srfdata,sfln_heatdata

    close (ifile_restart)

    call flx_update()

    write (iout,'(a,1pe13.4e2,a)') '  *** simulation restarted with dbl precision file from t=',time, ' s'

end subroutine restore_restart

! ***** restore routine: restore from restart file or state file *****
subroutine restore
    character(string_length) :: restart_file
    character(string_length) :: state_file
    if (verbose_restore) write(iout,*) 'read_restart_file:',read_restart_file
    if (verbose_restore) write(iout,*) 'read_state_file:',read_state_file
    ! cannot restore from restart and state file  at the same time
    if (trim(read_restart_file).ne."no".and.trim(read_state_file).ne."no") then
        write(iout,*) 'ERROR: Cannot restore from restart file and state file simultaneously'
        stop 'Exiting FACE'
    endif

    ! restore from restart file?
    if (trim(read_restart_file).eq."no") then
        write(iout,*) 'no restoration from restart file'
    elseif (read_restart_file.eq."yes") then
        restart_file=trim(path_folder)//'dsave.rst'
        call restore_restart(restart_file)
    else
        restart_file=read_restart_file
        call restore_restart(restart_file)
    endif

    ! restore from state file?
    if (trim(read_state_file).eq."no") then
        write(iout,*) 'no restoration from state file'
    elseif (read_state_file.eq."yes") then
        state_file=trim(path_folder)//'face.state'
        call restore_state(state_file)
    else
        state_file=read_state_file
        call restore_state(state_file)
    endif

end subroutine restore

! ***** *****
subroutine output_final_state
store_state_file=trim(path_folder)//casename//".state"
    call store_state(store_state_file)
    call FACE2fluidcode()

    call print_milestone('dumping file state completed')

end subroutine output_final_state

subroutine close_log()
    if (iout.ne.6) then
        close(iout)
    endif
end subroutine close_log

subroutine print_summary()
write(iout,*) 'Successful execution of FACE !'
write(iout,*) '***************Summary***************'
write(iout,*) '# iteration =', iteration -1
write(iout,'("cpu time of execution= ",f6.3," seconds.")') tcpufinish-tcpustart
write(iout,*) '*************************************'
end subroutine print_summary

 subroutine finalize()
        call close_timedata_files
        call deallocate_variables()
        call close_log
        write(iout,*) 'Quitting FACE after complete normal execution...'
    end subroutine finalize


subroutine write_header_log
character(string_length)::timestamp
call timestring ( timestamp )
write(iout,*) '# created    :', timestamp
write(iout,*) '# casename   :', trim(casename)
write(iout,*) '# pathfolder :', trim(path_folder)
end subroutine write_header_log

subroutine print_milestone(str)
character(*)::str
write(iout,"(a)") " "
write(iout,"('--- ',a,' ---')") str
end subroutine print_milestone

end module modFACE_IO
