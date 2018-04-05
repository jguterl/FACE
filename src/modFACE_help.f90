module modFACE_help
    use modFACE_misc
    use modFACE_IO


    !      use modFACE_parser,only:input_line,nlines
    implicit none
    integer,save:: ihelp,Nhelp
    integer,parameter::nlines_max_help=120

    type helper
        character(string_length)::keyword
        character(string_length)::def
        character(string_length)::units
        character(string_length)::status
        character(string_length)::default
        logical::comment
        logical:: species
    end type helper
    type(helper)::help(nlines_max_help)
    interface set_help
        module procedure set_help_single
        module procedure set_help_species
        module procedure set_help_comment
    end interface
contains

    subroutine print_version
      if (verbose_version) then
        write(iout,*) '###################################'
        call write_short_version
        write(iout,*) '#  Last modified: 23  March 2018  #'
        write(iout,*) '#  R.D. Smirnov & J. Guterl       #'
        write(iout,*) '#  email: rsmirnov@eng.ucsd.edu   #'
        write(iout,*) '#  email: guterlj@fusion.gat.com  #'
        write(iout,*) '###################################'
      endif
    end subroutine print_version

    subroutine write_short_version
        write(iout,*) '#  --->  FACE version 2.1   <----  #'
    end subroutine write_short_version

    subroutine print_help()
        character(30) str1
        character(70) str2
        character(20) str3
        character(15) str4
        character(20) str5
        integer i
        write(iout,*) ' '
        write(iout,*) ' *** List of keywords for FACE input file ***'
        write (str1, '(A30)') 'keyword'
        write (str2, '(A70)') 'def'
        write (str3, '(A20)') 'units'
        write (str4, '(A15)') 'status'
        write (str5, '(A15)') 'default'
        write(iout,'(A5,A30,A70,A20,A15,A20)') '  - ',adjustl(str1),adjustl(str2),adjustl(str3),adjustl(str4),adjustl(str5)
        do i=1,ihelp-1
            write (str1, '(A30)') help(i)%keyword
            write (str2, '(A70)') help(i)%def
            write (str3, '(A20)') help(i)%units
            write (str4, '(A15)') help(i)%status
            write (str5, '(A15)') help(i)%default
            write(iout,'(A5,A30,A70,A20,A15)') '  - ',adjustl(str1),adjustl(str2),adjustl(str3),adjustl(str4),adjustl(str5)
        enddo
    end subroutine print_help

        subroutine print_list_keyword

        integer i
        write(iout,*) ' '
        write(iout,*) ' *** List of keywords for FACE input file ***'
        do i=1,nlines_max_help
             if (help(i)%keyword.ne."none") then
            write (iout, '(a,a)')  ' - ', adjustl(trim(help(i)%keyword))
           endif
           enddo
    end subroutine print_list_keyword

    subroutine init_help()
        integer i
        character(string_length) :: comment_str
         do i=1,nlines_max_help
         help(i)=helper('none','none','none','none','none',.false.,.false.)
        enddo
        ihelp=1
        call set_help('order_solver','Numerical order of solver: 1|2|5','none','non-mandatory',"2")
        call set_help('read_restart_file','read Restart file: yes|no|filename (yes:default "dsave.rst")',&
        'none','non-mandatory',"no")
        call set_help('read_state_file','read state file: yes|no|filename (yes:"face.state").~=read restart file',&
        'none','non-mandatory',"no")
        call set_help('wall_thickness','Wall thickness','[m]','mandatory',"0.01")
        call set_help('start_time','Start time ','[s]','mandatory',"0")
        call set_help('temp_ramp_start_time','Temperature ramp start time','[s]','non-mandatory',"0")
        call set_help('temp_ramp_stop_time ','Temperature ramp stop time','[s]','non-mandatory',"0")
        call set_help('end_time','Simulation time','[s]','mandatory',"1.00e+01")
        call set_help('dt','Nominal time step','[s]','mandatory',"1.00e-03")
        call set_help('min_dt','Minimum time step','[s]','mandatory',"1.00e-06")
        call set_help('iter_solver_max','Max internal iteration for solver','none','non-mandatory',"100")
        call set_help('reduction_factor_dt','Reduction factor for timestep (=1<->no reduction)','none','non-mandatory',"1")
        call set_help('Nstep_increase_dt','increase dt after N quick converged timestp when dt_face<dt','none','non-mandatory',"10")
        call set_help('filter_freq','Low-pass filter cut-off frequency','[s^-1]','non-mandatory',"1e99")
        call set_help('dump_space_dt','Spatial parameters saving time interval','[s]','non-mandatory',"0.00e-00")
        call set_help('dump_time_dt','Temporal parameters saving time interval','[s]','non-mandatory',"0.00e-00")
        call set_help('dump_restart_dt','Restart file saving time interval','[s]','non-mandatory',"1")
        call set_help('steady_state','Steady state yes|no','none','non-mandatory',"no")
        call set_help('temp_ramp_filename','Temperature ramp data file (ramp.dat) used if .ne. 0','none','non-mandatory',"none")
        call set_help('solve_heat_equation','yes|no :solve heat equation','none','non-mandatory',"no")
        comment_str='# ********* Grid parameters ********************************************'
        call set_help(comment_str)
        call set_help('n_cells','Number of cells','none','mandatory',"100")
        call set_help('cell_scaling_factor','Cell width scaling factor','none','mandatory',"1.15305056")
        comment_str='# ********* Parameters for volumetric model and species ***************'
        call set_help(comment_str)
        call set_help('n_species','Number of species','none','mandatory',"3")
        call set_help('species_name','Name of species','none','mandatory',"D H+Tr1 Tr1","species")
        call set_help('n0_profile','Initial density profile is Gaussian (G), step(S), linear(L),peak(P),flat(F)','none'&
            ,'non-mandatory',"G G G","species")
        call set_help('n0','Initial density n0','[m^-3]','non-mandatory',"1.00e+10 1.00e+10 1.00e+10","species")
        call set_help('n0_xmax','Position of maximum of Gaussian profile','[m]','non-mandatory',&
        "0.00e-00 0.00e-00 0.00e-00","species")
        call set_help('n0_width','Standard deviation of Gaussian profile','[m]','non-mandatory',&
        "0.00e-00 0.00e-00 0.00e-00","species")
        call set_help('n_max','Maximum density of species','[m^-3]','non-mandatory',"1.00e+29 1.00e+29 1.00e+29","species")

        call set_help('D0','Diffusion constant of species','[m^2 s^-1]','non-mandatory',"1.00e-07 0.00e-07 0.00e-07","species")
        call set_help('ED','Activation energy of diffusion of species','[eV]','non-mandatory',"0.41 0.00 0.00","species")
        call set_help('Edt', 'Activation energy of detrapping of species', '[eV]' , 'non-mandatory' ,"0.00 0.00 1.40","species")
        call set_help('Etr', 'Activation energy of trapping of species'  , '[eV]' , 'non-mandatory' ,"0.00 1.00 0.00","species")

        comment_str='# ********* Parameters for surface model and species ******************* '
        call set_help(comment_str)
        call set_help('left_surface_model','B: Gamamaout=Kdes*cb^2 S: Gammaout=Kcs^2 N: no flux','[m^-2]',&
        'mandatory',"S N N","species")
        call set_help('right_surface_model','B: Gamamaout=Kdes*cb^2 S: Gammaout=Kcs^2 N: no flux','[m^-2]',&
        'mandatory',"S N N","species")

        call set_help('ns0_left','Initial left surface density of species','[m^-2]','non-mandatory',&
            "1.00e+19 0.00e+19 0.00e+19","species")
        call set_help('ns0_right','Initial right surface density of species','[m^-2]','non-mandatory',&
            "1.00e+19 0.00e+19 0.00e+19","species")
        call set_help('ns_max','Maximum surface density of species','[m^-2]','non-mandatory',"1.00e+19 1.00e+19 1.00e+19","species")
        call set_help('Eabs_left'  ,'Energy of absorption (vacuum->surface)','[eV]','non-mandatory'&
            ,"0.10 0.00 0.00","species")
        call set_help('Eabs_right' ,'Energy of absorption (vacuum->surface)','[eV]','non-mandatory'&
            ,"0.10 0.00 0.00","species")
        call set_help('Edes_left' ,'Energy of desorption (surface->vacuum)','[eV]','mandatory',&
        "1000 0.00 0.00","species")
        call set_help('Edes_right','Energy of desorption (surface->vacuum)','[eV]','mandatory',&
        "1000 0.00 0.00","species")
        call set_help('Eb_left'  ,'Energy of bulk absortion (surface->bulk)','[eV]','mandatory'&
        ,"2.00 0.00 0.00","species")
        call set_help('Eb_right' ,'Energy of bulk absortion (surface->bulk)','[eV]','mandatory'&
            ,"2.00 0.00 0.00","species")
        call set_help('Eads_left','Energy of adsorption (bulk->surface)','[eV]','mandatory',"1.00 0.00 0.00","species")
        call set_help('Eads_right','Energy of adsorption (bulk->surface)','[eV]','mandatory',"1.00 0.00 0.00","species")
        call set_help('nu0','Transition attempt frequency of species','[s^-1]','non-mandatory'&
        ,"1.00e+13 1.00e+13 1.00e+13","species")
        comment_str='# ********* Parameters for implantation model ************************** '
        call set_help(comment_str)
        call set_help('implantation_model','G: gaussian S:Step E: ERFC ','none','non-mandatory',"S S S","species")
        call set_help('implantation_depth','implentation  depth','[m]','non-mandatory',"5.00e-09 5.00e-09 5.00e-09","species")
        call set_help('implantation_width' ,'implentation width' ,'[m]','non-mandatory',"5.00e-09 5.00e-09 5.00e-09","species")

        call set_help('Eimpact_ion','Impact energy of ionized species','[eV]','non-mandatory',&
        "0.00e+00 0.00e+00 0.00e+00","species")
        call set_help('Gamma_in','External flux of ionized species','[m^-2 s^-1]','non-mandatory',&
        "0.00e+00 0.00e+00 0.00e+00","species")
        call set_help('Gamma_in_max','Maximal external flux of ionized species','[m^-2 s^-1]','non-mandatory',&
        "0.00e+00 0.00e+00 0.00e+00","species")
        call set_help('pressure_neutral','External pressure of neutral species','[Pa]','non-mandatory',&
        "0.00e-00 0.00e+00 0.00e+00","species")
        call set_help('temp_neutral','External temperature of neutral species','[eV]','non-mandatory',&
        "0.00e-00 0.00e+00 0.00e+00","species")
        call set_help('mass','Mass of species','[kg]','non-mandatory',"3.343e-27 0.00e+00 0.00e+00","species")
        comment_str='# ********* Parameters for abliation model ****************************** '
        call set_help(comment_str)
        call set_help('min_ablation_velocity','min ablation speed in addition to sputtering','[m s^-1]','non-mandatory'&
            ,"0.00e-00")
        call set_help('max_ablation_velocity','max ablation speed in addition to sputtering','[m s^-1]','non-mandatory'&
            ,"0.00e-00 ")
        call set_help('sputtering_yield','Sputtering yield','none','non-mandatory',"0")

        comment_str='# ********* Parameters for temperature ********************************** '
        call set_help(comment_str)
        call set_help('mat_temp_ramp_start','Material temperature at ramp start (initial at left boundary)','[K]'&
            ,'non-mandatory',"373")
        call set_help('mat_temp_ramp_stop','Material temperature at ramp stop  (initial at right boundary)','[K]',&
           'non-mandatory',"373")
           comment_str='# ********* Material parameters ************************************** '
        call set_help(comment_str)
        call set_help('lattice_constant','Lattice constant of material','[m]','mandatory',"1.00e-10")
        call set_help('cristal_volume_factor','Cristal cell volume factor','none','non-mandatory',"1")
        call set_help('cristal_surface','Surface cell area factor','none','non-mandatory',"1")
        call set_help('lattice_length_factor','Lattice cell length factor','none','non-mandatory',"1")

        call set_help('thermal_conductivity','Thermal conductivity','[W m^-1 K^-1]','non-mandatory',"137.0")
        call set_help('heat_capacity','Heat capacity','[J kg^-1 K^-1]','non-mandatory',"1.400e+02")
        call set_help('density','Density','[kg m^-3]','non-mandatory',"1.930e+04")
        call set_help('emissivity','Emissivity','???','non-mandatory',"0")
        call set_help('heat_formation','Heat of formation per material atom','[eV]','non-mandatory',"0.54")
        call set_help('min_radiation_power','Minimal radiation power','[W m^-2]','non-mandatory',"0.00e-00")
        call set_help('max_radiation_power','Maximal radiation power','[W m^-2]','non-mandatory',"0.00e-00")
          call set_help('#Parameters for material')
        call set_help('first_ramp_end_time','End of first ramp time','[s]','non-mandatory',"1.0e+99")
        call set_help('second_ramp_start_time','Start of second ramp time','[s]','non-mandatory',"1.0e+99")
        call set_help('second_ramp_end_time','End of second ramp time','[s]','non-mandatory',"1.0e+99")
        call set_help('pulsed_flux','pulse incoming plasma flux: yes|no','none','non-mandatory',"no")
        call set_help('pulse_period','Pulse period','[s]','non-mandatory',"1.0e+99")

          comment_str='# ********* Miscelleneaous *********************************************'
        call set_help(comment_str)
        call set_help('active_cap','yes|no','none','non-mandatory',"yes")
        call set_help('verbose','yes|no','none','non-mandatory',"yes")
        call set_help('dump_space_append','yes|no','none','non-mandatory',"no")
        call set_help('dump_time_append','yes|no ','none','non-mandatory',"no")
        call set_help('fluid-interface','yes|no','none','non-mandatory',"no")

        call set_help('Nprint_run_info','print info on current run every Nprint_run_info steps','none','non-mandatory',"100")
        call set_help('solver_eps','solver: precision norm (norm<eps:exit','none','non-mandatory',"3.d-3")
        call set_help('solver_udspl','solver: displacement min vector u','none','non-mandatory',"9.d-1")
        call set_help('solver_fdspl','solver: displacement min function f(u)','none','non-mandatory',"9.d0")
        call set_help('solver_gdspl','solver: displacement min grad u','none','non-mandatory',"1.d-3")
        call set_help('solver_fstp','solver: step reduction for function f(u)','none','non-mandatory',"1.d-1")
   call set_help('couple_factor_dt','dt=dt_fluid*coupling_factor_dt when coupled to fluid code','none','non-mandatory',"1.d-2")
        if (ihelp-1.gt.nlines_max_help) then
        call face_error("Extend size of help array in help module: nlines_max_help=",nlines_max_help)
        endif
        if (verbose_help) call print_milestone('initialization help completed')

    end subroutine init_help

    subroutine  set_help_single(keyword,def,units,status,default)
        character(*):: keyword,def,units,status,default
        help(ihelp)%keyword=trim(keyword)
        if (verbose_help) write(iout,*) "ihelp=",ihelp,"keyword=",keyword
        help(ihelp)%def=trim(def)
        help(ihelp)%units=trim(units)
        help(ihelp)%status=trim(status)
        help(ihelp)%default=trim(default)
        help(ihelp)%species=.false.
        help(ihelp)%comment=.false.
        ihelp=ihelp+1
        Nhelp=ihelp-1
    end subroutine set_help_single

    subroutine  set_help_comment(comment)
        character(*):: comment
        !call StripFrontSpaces(comment)
        help(ihelp)%def=trim(comment)
        help(ihelp)%comment=.true.
        ihelp=ihelp+1
        Nhelp=ihelp-1
    end subroutine set_help_comment

    subroutine  set_help_species(keyword,def,units,status,default,species)
        character(*):: keyword,def,units,status,default,species
        help(ihelp)%keyword=trim(keyword)
        if (verbose_help) write(iout,*) "ihelp=",ihelp,"keyword=",keyword
        help(ihelp)%def=trim(def)
        help(ihelp)%units=trim(units)
        help(ihelp)%status=trim(status)
        help(ihelp)%default=trim(default)
        help(ihelp)%species=.true.
        help(ihelp)%comment=.false.
        ihelp=ihelp+1
        Nhelp=ihelp-1
    end subroutine set_help_species

    subroutine write_default_inputfile(filename,mode)
        character(*)::mode
        integer,parameter::l=15
        character(string_length)::str_keyword
        character(string_length)::str_tmp
        character(8)::str_type
        character(3*(l+1)):: strvaluel
        character(80)::str_def
        character(26)::str_status
        character(100)::str_default
        character(string_length)::str_units
        character(string_length)::fmt,timestamp
        character(*)::filename
        character(l)::str1

        character(3*l)::str_data,str2
        character(l)::strtmp(3)
        integer::i,j,idefault,ios

        call set_unit(idefault)
        open(unit=idefault, file=trim(filename), iostat=ios,action='write')
        if ( ios /= 0 ) then
            call face_error('Cannot write into default input file ', trim(filename))
        endif

        write(iout,*)'Default input file: "', trim(filename) ,'" created'
        call timestring ( timestamp )
        write(idefault,*) '#Default input for FACE. Created: ', timestamp

        do i=1,Nhelp
            if (.not.help(i)%comment) then
                write(fmt,*)'(a',3*l,')'
                ! if keyword is n_species then select if the input file should be with only H (1 spc) or H+Tr (3 species)
                if (help(i)%keyword.eq."n_species") then
                    if (mode.eq."H") then
                        write(str_data,fmt) "1"
                    elseif (mode.eq."H+Tr") then
                        write(str_data,fmt) "3"
                    else
                        call face_error("Unknown mode when writing default input file")
                    endif
               else
                    write(str_data,fmt) adjustl(help(i)%default)
                endif

                ! write keyword
                write(str_keyword,*) trim(adjustl(help(i)%keyword)),' '
                ! write definition
                write(str_def,*) ' ! ',trim(adjustl(help(i)%def)),' '
                ! write units
                write(str_units,*) trim(adjustl(help(i)%units)),' '
                ! write_default
                write(str_default,'(a8,a)') 'deflt: ',trim(adjustl(help(i)%default))
                ! write status
                if (.not.(help(i)%species)) then

                    write(str_type,'(a8)')', single'
                else
                    write(str_type,'(a8)') ',species'

                endif
                write(str_status,'(a1,a13,a8,a1)') '(',adjustr(trim(help(i)%status)) ,adjustl(str_type),') '

                ! write keyword+value
                if (.not.(help(i)%species)) then
                 write(fmt,*)'(a25,a1,a',l+1,')'
                    write(str_tmp,fmt) (adjustl(help(i)%keyword)),' ',adjustl(str_data)
                else
                str2=adjustl(str_data)
                do j=1,3
                        call SplitString(str2,str1,str2,' ')
                        write(fmt,*)'(a14,a1)'
                        write(strtmp(j),fmt) str1,' '

                enddo
                write(fmt,*)'(3a',(l+1),')'
                write(strvaluel,fmt) (strtmp(j),j=1,3)

                    write(fmt,*)'(a25,a1,a)'
                    write(str_tmp,fmt) (adjustl(help(i)%keyword)),' ',adjustl(strvaluel)
                endif


                if (.not.(help(i)%species)) then
                 write(fmt,*)'(a45,a80,a23,a12,a)'
                 write(idefault,fmt) adjustl(str_tmp),str_def,str_status,str_units,trim(str_default)
                else
                write(fmt,*)'(a75,a80,a23,a12,a)'
                write(idefault,fmt) (adjustl(str_tmp)),str_def,str_status,str_units,trim(str_default)
!
!
!                    !BUG INTEL FORT below
!                    !write(fmt,*)'(a30,a1,a',3*(l),'a1,a80,a30,a1,a12,a11,a)'
!                    !write(idefault,fmt) str_keyword,' ',adjustl(trim(strvaluel)), ' ! ',trim(strdef),trim(strsta),' ',adjustl(trim(struni)),&
!                    !' default: ',adjustl(trim(help(i)%default))
!                    write(fmt,*)'(a30,a1,a',3*(l+1),')'!,a80,a30,a1,a12,a11,a)'
!                    write(str_tmp,fmt) trim(adjustl(help(i)%keyword)),' ',adjustl(strvaluel)
!                    write(fmt,*)'(a3,a80,a30,a1,a12,a11,a)'
!                    write(str_tmp2,fmt) " ! ",trim(strdef),trim(str_status)," ",adjustl(trim(struni))," dflt: ",&
!                    adjustl(trim(help(i)%default))
!                    write(idefault,'(a,a)') trim(str_tmp),trim(str_tmp2)
                endif
            else
                write(idefault,*) trim(help(i)%def)
            endif

        enddo
        !
        !       write(*,*,'advance=no') adjustl(str_keyword),' ' strvalue ' ','! ',strtype, help(i)%def, ' (',help(i)%status ,') ',help(i)%units,help(i)%
        !       type(string)
        !       enddo
        close(idefault)

    end subroutine




end module modFACE_help
