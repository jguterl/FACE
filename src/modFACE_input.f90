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
            if (solve_heat_eq_string.ne."no".AND.solve_heat_eq_string.ne."yes") then
                write(iout,*) "ERROR: solve_heat_equation must be yes or no"
                stop
            endif
            if (solve_heat_eq_string.eq."no") then
            solve_heat_eq=.false.
            elseif (solve_heat_eq_string.ne."yes") then
            solve_heat_eq=.true.
            call face_error('heat equation solver must be checked!')
            endif


        ! check steady-state value

case('steady_state')
         if(steady_state_string.eq."yes") then
         steady_state=.true.
         elseif (steady_state_string.eq."no") then
         steady_state=.false.
         else
         call face_error("unknown option for steady_state (must be yes or no) : ",steady_state_string)
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





         case('onthefly_inventory')
         if(print_onthefly_inventory_string.eq."yes") then
         print_onthefly_inventory=.true.
         elseif (print_onthefly_inventory_string.eq."no") then
         print_onthefly_inventory=.false.
         else
         call face_error("unknown option for onthefly_inventory (must be yes or no) : ",print_onthefly_inventory_string)
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
case('variable_timestep')
if(variable_timestep_string.eq."yes") then
         variable_timestep=.true.
         elseif (variable_timestep_string.eq."no") then
         variable_timestep=.false.
         else
         call face_error("unknown option for variable_timestep (must be yes or no) : ",variable_timestep_string)
         endif

         if ((order_solver.gt.1).AND. (variable_timestep)) then
         call face_error("Cannot have variable timestep when order of numerical time scheme >1 (order_solver>1)")
         endif

         case('adjust_reduction_factor')
if(adjust_reduction_factor_string.eq."yes") then
         adjust_reduction_factor=.true.
         elseif (adjust_reduction_factor_string.eq."no") then
         adjust_reduction_factor=.false.
         else
         call face_error("unknown option for adjust_reduction_factor(must be yes or no) : ",adjust_reduction_factor_string)
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
        call get_keyword_value('steady_state', steady_state_string)
        call get_keyword_value('start_time', start_time)
        call get_keyword_value('temp_ramp_start_time', tramp0)
        call get_keyword_value('temp_ramp_stop_time', tramp1)
        call get_keyword_value('end_time', end_time)
        call get_keyword_value('max_iter', max_iter)
        call get_keyword_value('dt', dt0_face)
        call get_keyword_value('min_dt', min_dt_face)
         call get_keyword_value('max_dt', max_dt_face)
        call get_keyword_value('variable_timestep', variable_timestep_string)
        call get_keyword_value('filter_freq', nucut)
        call get_keyword_value('dump_space_dt', dump_space_dt)
        call get_keyword_value('dump_time_dt', dump_time_dt)
         call get_keyword_value('dump_vol_append', dump_vol_append_string)
         call get_keyword_value('dump_srf_append', dump_srf_append_string)
        call get_keyword_value('dump_restart_dt', dump_restart_dt)
        call get_keyword_value('temp_ramp_filename', framp)
        call get_keyword_value('solve_heat_equation', solve_heat_eq_string)
        call get_keyword_value('onthefly_inventory', print_onthefly_inventory_string)
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
        call get_keyword_value('left_surface_model', left_surface_model_string)
        call get_keyword_value('right_surface_model', right_surface_model_string)
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
        call get_keyword_value('diagnostic_depth', diagnostic_depth)
        call get_keyword_value('implantation_width', implantation_width)
        call get_keyword_value('Eimpact_ion', enrg)
        call get_keyword_value('Gamma_in', inflx_in)
        call get_keyword_value('Gamma_in_max', inflx_in_max)
        call get_keyword_value('Gamma_pulse', inflx_in_pulse)
        call get_keyword_value('Gamma_pulse_period', inflx_in_pulse_period)
        call get_keyword_value('Gamma_pulse_starttime', inflx_in_pulse_starttime)
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
        call get_keyword_value('grid_gen_mode', grid_gen_mode)
        call get_keyword_value('grid_dx0', grid_dx0)
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
!        call get_keyword_value('pulsed_flux', pulsed_flux)
!        call get_keyword_value('pulse_period', tpulse)
        call get_keyword_value('iter_solver_max', iter_solver_max)
        call get_keyword_value('reduction_factor_dt', reduction_factor_dt)
        call get_keyword_value('adjust_reduction_factor', adjust_reduction_factor_string)
          call get_keyword_value('Nstep_increase_dt',Nstep_increase_dt)

        call get_keyword_value('Nprint_run_info', Nprint_run_info)
        call get_keyword_value('solver_eps',solver_eps)
        call get_keyword_value('jac_eps',jac_eps)
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
        call init_zero(grid_dx0)
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
!        call init_zero(tpulse)
        if (verbose_input) write(iout,*) 'Initialization of single input parameters: OK'
    end subroutine init_input_single










    subroutine set_input_parameters
    write(iout,*) "ERROR: set_input_parameters no implemented yet..."
    stop
    end subroutine
end module modFACE_input
