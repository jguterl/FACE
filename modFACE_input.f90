module modFACE_input
    use modface_header
    use modFACE_output
    use modFACE_parser
    use modFACE_interface
    use modFACE_allocate
    implicit none
    save

    interface  get_keyword_value
        module procedure get_keyword_value_i
        module procedure get_keyword_value_r
        module procedure get_keyword_value_s
        module procedure get_keyword_value_species_i
        module procedure get_keyword_value_species_r
        module procedure get_keyword_value_species_s
    end interface get_keyword_value

contains




    subroutine read_inputfile(filename)
        implicit none
        character(*),intent(in)::filename
        call print_milestone('Reading input file : ' //trim(filename))
        call open_inputfile(filename)
        ! get the number of line in the file
        call get_nlines()
        ! allocate array of string

        if (.not.allocated(input_lines)) allocate(input_lines(nlines))

        rewind(iuinput)
        ! read the file again and get each lines
        call get_inputlines()

        close(iuinput)
        write(iout,*) 'Reading of input file "', trim(filename) ,'" : DONE '
        ! init input values
        call init_input_single
        call parse_keywords()
        deallocate(input_lines)
         write(iout,*) 'Parsing of input file "', trim(filename) ,'" : DONE '
    end subroutine read_inputfile



    subroutine check_value_input(keyword)
        integer ::k,kk
        character(*) :: keyword


        select case(keyword)

        case('n_species')
            if (nspc<1) then
                write(iout,*) "ERROR: n_species must be at least equal to 1"
                stop
            endif
            call alloc_input_species()

        !check order of the numerical solver
        case('order_solver')
            if (order_solver.ne.1.and.order_solver.ne.2.and.order_solver.ne.5) then
                write(iout,*) "ERROR: order of the solver must be 1 2 or 5. current order : ", order_solver
                stop 'Exiting FACE...'
            else
            ndt=order_solver+1
            write(iout,*) ' order of numerical order: ', order_solver, ' -> ndt =', ndt
            endif


        ! solve heat equation
        case('solve_heat_equation')
            if (solve_heat_eq.ne."no".AND.solve_heat_eq.ne."yes") then
                write(iout,*) "ERROR: solve_heat_equation must be yes or no"
                stop
            endif


        ! check steady-state value
        case('steady_state')
            if (stdst.ne."no".AND.stdst.ne."yes") then
                write(iout,*) "ERROR: steady_state must be yes or no"
                stop
            endif


        ! check that name of species are unique
        case('species_name')
           do k=1,nspc
             do kk=1,nspc
                 if (k.ne.kk.and.namespc(k).eq.namespc(kk)) then
                write(iout,*) "ERROR: two species have the same name: k=",k," and k=",kk
                write(iout,*) "ERROR:  namespc=",namespc(k)," and namespc=",namespc(kk)
                stop 'Exiting FACE'
            endif
       enddo
       enddo
       end select



    end subroutine check_value_input


    subroutine parse_keywords()
        if (verbose_input) write(iout,*) "Parsing keyword"
        call get_keyword_value('order_solver',order_solver)
        call get_keyword_value('read_restart_file',read_restart_file)
        call get_keyword_value('read_state_file',read_state_file)
        call get_keyword_value('wall_thickness',length)
        call get_keyword_value('steady_state', stdst)
        call get_keyword_value('start_time', start_time)
        call get_keyword_value('temp_ramp_start_time', tramp0)
        call get_keyword_value('temp_ramp_stop_time', tramp1)
        call get_keyword_value('end_time', end_time)
        call get_keyword_value('min_dt', dtmin)
        call get_keyword_value('timestep_factor', cdt)
        call get_keyword_value('filter_freq', nucut)
        call get_keyword_value('dump_space_dt', tspc)
        call get_keyword_value('dump_time_dt', ttm)
        call get_keyword_value('dump_restart_dt', tstr)
        call get_keyword_value('temp_ramp_filename', framp)
        call get_keyword_value('solve_heat_equation', solve_heat_eq)
        call get_keyword_value('n_species', nspc)

        call get_keyword_value('n0_max', dens0)
        call get_keyword_value('n0_profile', gprof)
        call get_keyword_value('n0_xmax', gxmax)
        call get_keyword_value('n0_width', gsigm)
        call get_keyword_value('n_max', densm)

        call get_keyword_value('D0', cdif0)
        call get_keyword_value('ED', edif)
        call get_keyword_value('Etr', etr )
        call get_keyword_value('Edt', edtr)

        call get_keyword_value('ns_max', dsrfm)
        call get_keyword_value('ns0_left', dsrfl0)
        call get_keyword_value('ns0_right', dsrfr0)
        call get_keyword_value('Ech_left', echl )
        call get_keyword_value('Ech_right', echr)
        call get_keyword_value('Q_ch_left', qchl)
        call get_keyword_value('Q_ch_right', qchr)
        call get_keyword_value('Eab_left', ebl)
        call get_keyword_value('Eab_right', ebr )
        call get_keyword_value('Es_left', esl)
        call get_keyword_value('Es_right', esr)
        call get_keyword_value('nu0', nu)
        call get_keyword_value('Eimpact_ion', enrg)
        call get_keyword_value('Gamma_min', inflx_min)
        call get_keyword_value('Gamma_max', inflx_max)
        call get_keyword_value('pressure_neutral', prg)
        call get_keyword_value('temp_neutral', tg)
        call get_keyword_value('mass', mass)

        call get_keyword_value('min_ablation_velocity', cero_min)
        call get_keyword_value('max_ablation_velocity', cero_max)
        call get_keyword_value('sputtering_yield', gamero)
        call get_keyword_value('mat_temp_ramp_start', temp0)
        call get_keyword_value('mat_temp_ramp_stop', temp1)
        call get_keyword_value('lattice_constant', lambda)
        call get_keyword_value('cristal_volume_factor', cvlm)
        call get_keyword_value('cristal_surface', csrf)
        call get_keyword_value('lattice_length_factor', clng)
        call get_keyword_value('n_cells', ngrd)
        call get_keyword_value('cell_scaling_factor', alpha)
        call get_keyword_value('thermal_conductivity', thcond)
        call get_keyword_value('heat_capacity', cp)
        call get_keyword_value('density', rho)
        call get_keyword_value('emissivity', emiss)
        call get_keyword_value('heat_formation', qform)
        call get_keyword_value('min_radiation_power', rad_min)
        call get_keyword_value('max_radiation_power', rad_max)
        call get_keyword_value('first_ramp_end_time', t1)
        call get_keyword_value('second_ramp_start_time', t2)
        call get_keyword_value('second_ramp_end_time', t3)
        call get_keyword_value('pulse_period', tp)
        if (verbose_input) write(iout,*) "Parsing keyword: DONE"
    end subroutine parse_keywords

   subroutine write_input_log()
        character(200)::filename
        integer::ios
        call set_ifile(ifile_inputlog)
        write (filename, '(a,a,a, i2.2, a)') trim(path_folder),trim(casename),'_input.face'
        open (ifile_inputlog, file=trim(filename),status='replace', iostat=ios)
        if (ios.ne.0) then
        write (iout,*) ' *** Cannot open file ', trim(filename)
        stop 'Exiting FACE...'
        else
        call write_input_log_keyword('order_solver',order_solver)
        call write_input_log_keyword('read_restart_file',read_restart_file)
        call write_input_log_keyword('read_state_file',read_state_file)
        call write_input_log_keyword('wall_thickness',length)
        call write_input_log_keyword('steady_state', stdst)
        call write_input_log_keyword('start_time', start_time)
        call write_input_log_keyword('end_time', end_time)
        call write_input_log_keyword('temp_ramp_start_time', tramp0)
        call write_input_log_keyword('temp_ramp_stop_time', tramp1)

        call write_input_log_keyword('min_dt', dtmin)
        call write_input_log_keyword('timestep_factor', cdt)
        call write_input_log_keyword('filter_freq', nucut)
        call write_input_log_keyword('dump_space_dt', tspc)
        call write_input_log_keyword('dump_time_dt', ttm)
        call write_input_log_keyword('dump_restart_dt', tstr)
        call write_input_log_keyword('temp_ramp_filename', framp)
        call write_input_log_keyword('solve_heat_equation', solve_heat_eq)
        call write_input_log_keyword('n_species', nspc)

        call write_input_log_keyword('n0_max', dens0)
        call write_input_log_keyword('n0_profile', gprof)
        call write_input_log_keyword('n0_xmax', gxmax)
        call write_input_log_keyword('n0_width', gsigm)
        call write_input_log_keyword('n_max', densm)

        call write_input_log_keyword('D0', cdif0)
        call write_input_log_keyword('ED', edif)
        call write_input_log_keyword('Etr', etr )
        call write_input_log_keyword('Edt', edtr)

        call write_input_log_keyword('ns_max', dsrfm)
        call write_input_log_keyword('ns0_left', dsrfl0)
        call write_input_log_keyword('ns0_right', dsrfr0)
        call write_input_log_keyword('Ech_left', echl )
        call write_input_log_keyword('Ech_right', echr)
        call write_input_log_keyword('Q_ch_left', qchl)
        call write_input_log_keyword('Q_ch_right', qchr)
        call write_input_log_keyword('Eab_left', ebl)
        call write_input_log_keyword('Eab_right', ebr )
        call write_input_log_keyword('Es_left', esl)
        call write_input_log_keyword('Es_right', esr)
        call write_input_log_keyword('nu0', nu)
        call write_input_log_keyword('Eimpact_ion', enrg)
        call write_input_log_keyword('Gamma_min', inflx_min)
        call write_input_log_keyword('Gamma_max', inflx_max)
        call write_input_log_keyword('pressure_neutral', prg)
        call write_input_log_keyword('temp_neutral', tg)
        call write_input_log_keyword('mass', mass)

        call write_input_log_keyword('min_ablation_velocity', cero_min)
        call write_input_log_keyword('max_ablation_velocity', cero_max)
        call write_input_log_keyword('sputtering_yield', gamero)
        call write_input_log_keyword('mat_temp_ramp_start', temp0)
        call write_input_log_keyword('mat_temp_ramp_stop', temp1)
        call write_input_log_keyword('lattice_constant', lambda)
        call write_input_log_keyword('cristal_volume_factor', cvlm)
        call write_input_log_keyword('cristal_surface', csrf)
        call write_input_log_keyword('lattice_length_factor', clng)
        call write_input_log_keyword('n_cells', ngrd)
        call write_input_log_keyword('cell_scaling_factor', alpha)
        call write_input_log_keyword('thermal_conductivity', thcond)
        call write_input_log_keyword('heat_capacity', cp)
        call write_input_log_keyword('density', rho)
        call write_input_log_keyword('emissivity', emiss)
        call write_input_log_keyword('heat_formation', qform)
        call write_input_log_keyword('min_radiation_power', rad_min)
        call write_input_log_keyword('max_radiation_power', rad_max)
        call write_input_log_keyword('first_ramp_end_time', t1)
        call write_input_log_keyword('second_ramp_start_time', t2)
        call write_input_log_keyword('second_ramp_end_time', t3)
        call write_input_log_keyword('pulse_period', tp)
        close(ifile_inputlog)
  endif
    end subroutine write_input_log


    subroutine get_keyword_value_r(keyword,variable)
        character(*)::keyword
        real(DP),intent(inout)::variable
        type(input_val)::inputval
        if (verbose_input) write(iout,*) 'reading real:',keyword
        call assign_keyword_value(keyword,inputval,'r')
        variable=inputval%r
        if (verbose_input) write(iout,*) 'real: ', keyword,'=',variable ,' : ' ,inputval%status
        call check_value_input(keyword)

    end subroutine get_keyword_value_r


    subroutine get_keyword_value_i(keyword,variable)
        character(*)::keyword
        integer,intent(inout)::variable
        type(input_val)::inputval
        if (verbose_input) write(*,*) 'reading int:',keyword
        call assign_keyword_value(keyword,inputval,'i')
        variable=inputval%i
        if (verbose_input) write(*,*) 'int:',keyword,'=',variable ,' : ' ,inputval%status
        call check_value_input(keyword)

    end subroutine get_keyword_value_i

    subroutine get_keyword_value_s(keyword,variable)
        character(*)::keyword
        character(*),intent(inout)::variable
        type(input_val)::inputval
        if (verbose_input) write(*,*) 'reading str:',keyword
        call assign_keyword_value(keyword,inputval,'s')
        variable=adjustl(inputval%s)
        if (verbose_input) write(*,*) 'str:',keyword,'=',variable ,' : ' ,inputval%status
        call check_value_input(keyword)

    end subroutine get_keyword_value_s

    subroutine set_str(str,str0)
    character(*)::str,str0
    integer idx,i
    call StripFrontSpaces(str0)
    idx=index(str0,' ')
    if (idx>1) then
    do i=1,idx
    str(i:i)=str0(i:i)
    enddo
    else
    str=str0
    endif
    end subroutine set_str
    !read keyward values for species
    subroutine get_keyword_value_species_i(keyword,variable)
        character(*)::keyword
        type(input_valsp)::inputval
        integer::k
        integer,intent(inout)::variable(:)
        call check_allocate_species()
        allocate(inputval%s(nspc))
        allocate(inputval%r(nspc))
        allocate(inputval%i(nspc))
        if (verbose_input) write(*,*) 'reading int:',keyword
        call assign_keyword_value_species(keyword,inputval,'i')
        variable=inputval%i
        if (verbose_input) write(iout,*) 'int:',keyword,'=',(variable(k),k=1,nspc) ,' : ' ,inputval%status
        call check_value_input(keyword)

        deallocate(inputval%s)
        deallocate(inputval%r)
        deallocate(inputval%i)
    end subroutine get_keyword_value_species_i



    subroutine get_keyword_value_species_r(keyword,variable)
        character(*)::keyword
        integer::k
        type(input_valsp)::inputval
        real(DP),intent(inout)::variable(:)
        call check_allocate_species()
        allocate(inputval%s(nspc))
        allocate(inputval%r(nspc))
        allocate(inputval%i(nspc))
        if (verbose_input) write(*,*) 'reading real:',keyword
        call assign_keyword_value_species(keyword,inputval,'r')
        variable=inputval%r
        if (verbose_input) write(iout,*) 'real:',keyword,'=',(variable(k),k=1,nspc)  ,' : ' ,inputval%status
        call check_value_input(keyword)
        deallocate(inputval%s)
        deallocate(inputval%r)
        deallocate(inputval%i)
    end subroutine get_keyword_value_species_r

    subroutine get_keyword_value_species_s(keyword,variable)
        character(*)::keyword
        integer ::k
        type(input_valsp)::inputval
        character(lname),intent(inout)::variable(:)

        call check_allocate_species()
        allocate(inputval%s(nspc))
        allocate(inputval%r(nspc))
        allocate(inputval%i(nspc))
        if (verbose_input) write(*,*) 'reading real:',keyword
        call assign_keyword_value_species(keyword,inputval,'s')
        variable=inputval%s
if (verbose_input) write(iout,*) 'str:',keyword,'=',(variable(k),k=1,nspc) ,' : ' ,inputval%status
        deallocate(inputval%s)
        deallocate(inputval%r)
        deallocate(inputval%i)
    end subroutine get_keyword_value_species_s

    subroutine init_input()
        call init_input_single()
    end subroutine

    subroutine init_input_single()

        call init_zero(stdst)
        call init_zero(length)
        call init_zero(start_time)
        call init_zero(tramp0)
        call init_zero(tramp1)
        call init_zero(end_time)
        call init_zero(dtmin)
        call init_zero(cdt)
        call init_zero(nucut)
        call init_zero(tspc)
        call init_zero(ttm)
        call init_zero(tstr)
        call init_zero(framp)
        call init_zero(solve_heat_eq)
        call init_zero(cero_min)
        call init_zero(cero_max)
        call init_zero(gamero)
        call init_zero( temp0)
        call init_zero( temp1)
        call init_zero(lambda)
        call init_zero( cvlm)
        call init_zero(csrf)
        call init_zero(clng)
        call init_zero( ngrd)
        call init_zero(alpha)
        call init_zero( thcond)
        call init_zero( cp)
        call init_zero( rho)
        call init_zero( emiss)
        call init_zero( qform)
        call init_zero(rad_min)
        call init_zero(rad_max)
        call init_zero(t1)
        call init_zero( t2)
        call init_zero( t3)
        call init_zero(tp)
        if (verbose_input) write(iout,*) 'Initialization of single input parameters: OK'
    end subroutine init_input_single



    !      write (6, 1090) cero_min
    !      write (6, 1100) cero_max
    !      write (6, 1110) gamero
    !      write (6, 1120) ngrd
    !      write (6, 1130) alpha
    !      write (6, 1140) thcond
    !      write (6, 1150) cp
    !      write (6, 1160) rho
    !      write (6, 1170) emiss
    !      write (6, 1180) qform
    !      write (6, 1190) rad_min
    !      write (6, 1200) rad_max
    !      write (6, 1210) t1
    !      write (6, 1220) t2
    !      write (6, 1230) t3
    !      write (6, 1240) tp
    !
    !                write (6, 1010) length, start_time, tramp0, tramp1, end_time, dtmin
    !1010  format (' length of the simulation region: ', 1pe12.5, ' m'/
    !     + ' start time: ', 1pe12.5, ' s'/
    !     + ' temperature ramp start: ', 1pe12.5, ' s'/
    !     + ' temperature ramp stop: ', 1pe12.5, ' s'/
    !     + ' simulation time: ', 1pe12.5, ' s'/
    !     + ' minimal time step: ', 1pe12.5, ' s')
    !            write (6, 1030) nspc
    !      write (6, 1040) (i, dens0(i), dsrfl0(i), dsrfr0(i), i=1,nspc)
    !1030  format (' number of species is ', i2)
    !1040  format (' species of sort ', i2,
    !     + '  have initial (max) density ', 1pe12.5, ' m^-3',
    !     + '  initial left  surface density', 1pe12.5, ' m^-2',
    !     + '  initial right surface density', 1pe12.5, ' m^-2')
    !            write (6, 1020) tspc, ttm, tstr
    !1020  format (' spatial parameters saving interval: ', 1pe12.5, ' s'/
    !     + ' temporal parameters saving interval: ', 1pe12.5, ' s'/
    !     + ' restart file saving interval: ', 1pe12.5, ' s')
    !      write (6, 1021) 'Double precision is used'
    !1021  format (a)
    !
    !
    !
    !      write (6, 1050) (i, enrg(i), i=1,nspc)
    !1050  format (' impact energy of species ', i2, ' is ',
    !     + 1pe12.5, ' eV')
    !      write (6, 1060) (i, inflx_min(i), i=1,nspc)
    !1060  format (' min influx of species ', i2, ' is', 1pe12.5,
    !     +        ' m^-2 s^-1')
    !      write (6, 1061) (i, inflx_max(i), i=1,nspc)
    !1061  format (' max influx of species ', i2, ' is', 1pe12.5,
    !     +        ' m^-2 s^-1')
    !      write (6, 1070) (i, prg(i), i=1,nspc)
    !1070  format (' ext. pressure of species ', i2, ' is', 1pe12.5, ' Pa')
    !      write (6, 1080) (i, tg(i), i=1,nspc)
    !1080  format (' ext. temperature of species ', i2, ' is', 1pe12.5,
    !     + ' eV')
    !1090  format (' min ablation speed is ', 1pe12.5, ' m/s')
    !1100  format (' max ablation speed is ', 1pe12.5, ' m/s')
    !1110  format (' sputtering yiled is ', 1pe12.5)
    !1120  format (' number of grid points is ', i6)
    !1130  format (' grid scaling factor is ', 1pe12.5)
    !1140  format (' thermal conductivity is ', 1pe12.5, ' W/(m*K)')
    !1150  format (' heat capacity is ', 1pe12.5, ' J/(kg*K)')
    !1160  format (' density factor is ', 1pe12.5, ' kg/m^3')
    !1170  format (' emissivity is ', 1pe12.5)
    !1180  format (' heat of formation is ', 1pe12.5, ' eV')
    !1190  format (' min radiation power is ', 1pe12.5, ' W/m^2')
    !1200  format (' max radiation power is ', 1pe12.5, ' W/m^2')
    !1210  format (' end time of first ramp is ', 1pe12.5, ' s')
    !1220  format (' start time of second ramp is ', 1pe12.5, ' s')
    !1230  format (' end time of second ramp is ', 1pe12.5, ' s')
    !1240  format (' pulse period is ', 1pe12.5, ' s')
    !c
    !      close (unit=10)
    !      return
    !1999  stop ' *** error occured when reading ''face.ini''!'





    subroutine get_input(face_input)
        type(face_inputs):: face_input
        !if (.not.allocated(input_lines)) then
        !write(iout,*) 'check 0'
        !else
         !write(iout,*) 'check 0: allocated'
        !endif
        call init_input()
        if (face_input%read_input_file) then
            call read_inputfile(face_input%input_filename)
        else
            call set_input_parameters()
        endif
        call write_input_log
    end subroutine get_input

    subroutine set_input_parameters
    write(iout,*) "ERROR: set_input_parameters no implemented yet..."
    stop
    end subroutine
end module modFACE_input