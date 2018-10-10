module modFACE_input
    use modface_header
    use modFACE_output
    use modFACE_parser
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
        if (verbose_input)  write(iout,*) 'Reading of input file "', trim(filename) ,'" : DONE '
        ! init input values
        call init_input_single
        call check_keywords
        call parse_keywords
         call write_input_log
        deallocate(input_lines)
        if (verbose_input)  write(iout,*) 'Parsing of input file "', trim(filename) ,'" : DONE '
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
                call face_error("order of the solver must be 1 2 or 5. current order : ", order_solver)
            else
            ndt=order_solver+1
            if (verbose_init) write(iout,*) ' order of numerical order: ', order_solver, ' -> ndt =', ndt
            endif


        ! solve heat equation
        case('solve_heat_equation')
            if (solve_heat_eq.ne."no".AND.solve_heat_eq.ne."yes") then
                write(iout,*) "ERROR: solve_heat_equation must be yes or no"
                stop
            endif


        ! check steady-state value
        case('steady_state')
            if (steady_state.ne."no".AND.steady_state.ne."yes") then
                write(iout,*) "ERROR: steady_state must be yes or no"
                stop
            endif


        ! check that name of species are unique
        case('species_name')
           do k=1,nspc
             do kk=1,nspc
                 if (k.ne.kk.and.namespc(k).eq.namespc(kk)) then
                call face_error("two species have the same name: k=",k," and k=",kk," namespc=",namespc(k),&
                " and namespc=",namespc(kk))
            endif
       enddo
       enddo

       ! set active_cap (true or false)
        case('active_cap')
         if(active_cap_string.eq."yes") then
         active_cap=.true.
         elseif (active_cap_string.eq."no") then
         active_cap=.false.
         else
         call face_error("unknown option for active_cap (must be yes or no) : ",active_cap_string)
         endif

         case('dump_vol_append')
         if(dump_vol_append_string.eq."yes") then
         dump_vol_append=.true.
         elseif (dump_vol_append_string.eq."no") then
         dump_vol_append=.false.
         else
         call face_error("unknown option for dump_vol_append (must be yes or no) : ",dump_vol_append_string)
         endif

         case('dump_srf_append')
         if(dump_srf_append_string.eq."yes") then
         dump_srf_append=.true.
         elseif (dump_srf_append_string.eq."no") then
         dump_srf_append=.false.
         else
         call face_error("unknown option for dump_srf_append (must be yes or no) : ",dump_srf_append_string)
         endif

! check that the reduction factor is >0 and lt 1
case('reduction_factor_dt')
if (reduction_factor_dt.le.0d0.or.reduction_factor_dt.gt.1d0) then
call face_error("reduction factor dt cannot be <=0 and >1: reduction factor dt=",reduction_factor_dt)
endif
end select
    end subroutine check_value_input

    subroutine check_keywords
    integer :: i
    character(string_length) :: str
    do i=1,nlines
    str=input_lines(i)%keyword
    call StripFrontSpaces(str)
        if (find_keyword_help(trim(str)).eq.-1.and.str(1:1).ne.'#') then
            call face_error('Unknown keyword "', trim(input_lines(i)%keyword) ,'"at line=',i,&
            ' (see list of keywords --list-keywords or -lk)')
        endif
    enddo
    end subroutine check_keywords


    subroutine parse_keywords()
        if (verbose_input) write(iout,*) "Parsing keyword"
        call get_keyword_value('order_solver',order_solver)
        call get_keyword_value('read_restart_file',read_restart_file)
        call get_keyword_value('read_state_file',read_state_file)
        call get_keyword_value('wall_thickness',length)
        call get_keyword_value('steady_state', steady_state)
        call get_keyword_value('start_time', start_time)
        call get_keyword_value('temp_ramp_start_time', tramp0)
        call get_keyword_value('temp_ramp_stop_time', tramp1)
        call get_keyword_value('end_time', end_time)
        call get_keyword_value('dt', dt0_face)
        call get_keyword_value('min_dt', min_dt_face)
        call get_keyword_value('filter_freq', nucut)
        call get_keyword_value('dump_space_dt', dump_space_dt)
        call get_keyword_value('dump_time_dt', dump_time_dt)
         call get_keyword_value('dump_vol_append', dump_vol_append_string)
         call get_keyword_value('dump_srf_append', dump_srf_append_string)
        call get_keyword_value('dump_restart_dt', dump_restart_dt)
        call get_keyword_value('temp_ramp_filename', framp)
        call get_keyword_value('solve_heat_equation', solve_heat_eq)
        call get_keyword_value('n_species', nspc)

        call get_keyword_value('species_name', namespc)
        call get_keyword_value('n0', dens0)
        call get_keyword_value('n0_profile', gprof)
        call get_keyword_value('n0_xmax', gxmax)
        call get_keyword_value('n0_width', gsigm)
        call get_keyword_value('n_max', densm)

        call get_keyword_value('D0', cdif0)
        call get_keyword_value('ED', edif)
        call get_keyword_value('Etr', etr )
        call get_keyword_value('Edt', edtr)
        call get_keyword_value('left_surface_model', left_surface_model)
        call get_keyword_value('right_surface_model', right_surface_model)
        call get_keyword_value('order_desorption_left', order_desorption_left)
        call get_keyword_value('order_desorption_right', order_desorption_right)
        call get_keyword_value('ns0_left', dsrfl0)
        call get_keyword_value('ns0_right', dsrfr0)
        call get_keyword_value('ns_max', dsrfm)
        call get_keyword_value('Eabs_left', Eabs_l )
        call get_keyword_value('Eabs_right', Eabs_r)
        call get_keyword_value('Edes_left', Edes_l)
        call get_keyword_value('Edes_right', Edes_r)
        call get_keyword_value('Eb_left', Eb_l)
        call get_keyword_value('Eb_right', Eb_r )
        call get_keyword_value('Eads_left', Eads_l)
        call get_keyword_value('Eads_right', Eads_r)
        call get_keyword_value('nu0', nu)
        call get_keyword_value('implantation_model', implantation_model)
        call get_keyword_value('implantation_depth', implantation_depth)
        call get_keyword_value('implantation_width', implantation_width)
        call get_keyword_value('Eimpact_ion', enrg)
        call get_keyword_value('Gamma_in', inflx_in)
        call get_keyword_value('Gamma_in_max', inflx_in_max)
        call get_keyword_value('pressure_neutral', gas_pressure)
        call get_keyword_value('temp_neutral', gas_temp)
        call get_keyword_value('mass', mass)
        call get_keyword_value('active_cap', active_cap_string)
        call get_keyword_value('min_ablation_velocity', cero_min)
        call get_keyword_value('max_ablation_velocity', cero_max)
        call get_keyword_value('sputtering_yield', gamero)
        call get_keyword_value('mat_temp_ramp_start', temp0)
        call get_keyword_value('mat_temp_ramp_stop', temp1)
        call get_keyword_value('lattice_constant', lambda)
        call get_keyword_value('volume_factor', cvlm)
        call get_keyword_value('surface_factor', csrf)
        call get_keyword_value('lattice_length_factor', clng)
        call get_keyword_value('n_cells', ngrd)
        call get_keyword_value('cell_scaling_factor', alpha)
        call get_keyword_value('grid_type', grid_type)
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
        call get_keyword_value('pulsed_flux', pulsed_flux)
        call get_keyword_value('pulse_period', tpulse)
        call get_keyword_value('iter_solver_max', iter_solver_max)
        call get_keyword_value('reduction_factor_dt', reduction_factor_dt)
          call get_keyword_value('Nstep_increase_dt',Nstep_increase_dt)

        call get_keyword_value('Nprint_run_info', Nprint_run_info)
        call get_keyword_value('solver_eps',solver_eps)
        call get_keyword_value('solver_udspl',solver_udspl)
        call get_keyword_value('solver_fdspl',solver_fdspl)
        call get_keyword_value('solver_gdspl',solver_gdspl)
        call get_keyword_value('solver_fstp',solver_fstp)

        if (verbose_input) write(iout,*) "Parsing keyword: DONE"
    end subroutine parse_keywords

   subroutine write_input_log()
        character(string_length)::filename
        integer::ios,i
        call set_unit(unit_inputlog)
        filename= trim(path_folder)//trim(casename)//'_input.log'
        open (unit_inputlog, file=trim(filename),status='replace', iostat=ios)
        if (ios.ne.0) then
        call face_error('Cannot open file ', trim(filename))
        else
        do i=1,nlines!
        write(unit_inputlog, '(a50,a3,a120)') adjustl(trim(input_lines(i)%keyword)),' : ',adjustl(trim(input_lines(i)%data))
        enddo
        close(unit_inputlog)
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

    !read keyword values for species
    subroutine get_keyword_value_species_i(keyword,variable)
        character(*)::keyword
        type(input_valsp)::inputval
        integer::k
        integer,intent(inout)::variable(:)
        call check_allocate_species()
        allocate(inputval%s(nspc))
        allocate(inputval%r(nspc))
        allocate(inputval%i(nspc))
        if (verbose_input) write(*,*) 'reading int:',trim(keyword)
        call assign_keyword_value_species(keyword,inputval,'i')
        variable=inputval%i
        if (verbose_input) write(iout,*) 'int:',trim(keyword),'=',(variable(k),k=1,nspc) ,' : ' ,trim(inputval%status)
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
        if (verbose_input) write(*,*) 'reading real:',trim(keyword)
        call assign_keyword_value_species(keyword,inputval,'r')
        variable=inputval%r
        if (verbose_input) write(iout,*) 'real:',trim(keyword),'=',(variable(k),k=1,nspc)  ,' : ' ,trim(inputval%status)
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
        if (verbose_input) write(*,*) 'reading real:',trim(keyword)
        call assign_keyword_value_species(keyword,inputval,'s')
        variable=inputval%s
if (verbose_input) write(iout,*) 'str:',trim(keyword),'=',(variable(k),k=1,nspc) ,' : ' ,trim(inputval%status)
        deallocate(inputval%s)
        deallocate(inputval%r)
        deallocate(inputval%i)
    end subroutine get_keyword_value_species_s

    subroutine init_input()
        call init_input_single()
    end subroutine

    subroutine init_input_single()

        call init_zero(steady_state)
        call init_zero(length)
        call init_zero(start_time)
        call init_zero(tramp0)
        call init_zero(tramp1)
        call init_zero(end_time)
        call init_zero(dt0_face)
!        call init_zero(cdt)
        call init_zero(nucut)
        call init_zero(dump_time_dt)
        call init_zero(dump_space_dt)
        call init_zero(dump_restart_dt)
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
        call init_zero(tpulse)
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
    !      write (6, 1060) (i, inflx_in(i), i=1,nspc)
    !1060  format (' min influx of species ', i2, ' is', 1pe12.5,
    !     +        ' m^-2 s^-1')
    !      write (6, 1061) (i, inflx_in_max(i), i=1,nspc)
    !1061  format (' max influx of species ', i2, ' is', 1pe12.5,
    !     +        ' m^-2 s^-1')
    !      write (6, 1070) (i, gas_pressure(i), i=1,nspc)
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







    subroutine set_input_parameters
    write(iout,*) "ERROR: set_input_parameters no implemented yet..."
    stop
    end subroutine
end module modFACE_input
