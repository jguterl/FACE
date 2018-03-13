module modFACE_help
    use modFACE_misc
    use modFACE_output,only:iout,timestring

    !      use modFACE_parser,only:input_line,nlines
    implicit none
    integer,save:: ihelp=1,Nhelp

    type helper
        character(100),allocatable:: keyword
        character(200),allocatable::def
        character(100),allocatable::units
        character(100),allocatable::status
        character(80),allocatable::default
        logical::comment
        logical:: species
    end type helper
    type(helper)::help(1:100)
    interface set_help
        module procedure set_help_single
        module procedure set_help_species
        module procedure set_help_comment
    end interface
contains

    subroutine display_version()
        write(iout,*) '*********************************'
        write(iout,*) '* FACE version 1.1              *'
        write(iout,*) '* Last modified: 15  Feb 2018   *'
        write(iout,*) '* R.D. Smirnov & J. Guterl      *'
        write(iout,*) '* email: rsmirnov@eng.ucsd.edu  *'
        write(iout,*) '* email: guterlj@fusion.gat.com *'
        write(iout,*) '*********************************'
    end subroutine display_version

    subroutine display_help()
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
    end subroutine display_help

    subroutine init_help()
        call set_help('order_solver','Numerical order of solver: 1|2|5','none','non-mandatory',"2")
        call set_help('read_restart_file','read Restart file: yes|no|filename (yes:default "dsave.rst")',&
        'none','non-mandatory',"no")
        call set_help('read_state_file','read state file: yes|no|filename (yes:"face.state").~=read restart file',&
        'none','non-mandatory',"no")
        call set_help('wall_thickness','Wall thickness','[m]','mandatory',"0.01")
        call set_help('start_time','Start time ','[s]','mandatory',"0")
        call set_help('temp_ramp_start_time','Temperature ramp start time','[s]','non-mandatory',"0")
        call set_help('temp_ramp_stop_time ','Temperature ramp stop time','[s]','non-mandatory',"0")
        call set_help('simulation_end_time','Simulation time','[s]','mandatory',"1.00e+01")
        call set_help('min_dt','Minimal time step','[s]','mandatory',"1.00e-03")
        call set_help('timestep_factor','Time step factor','none','non-mandatory',"1")
        call set_help('filter_freq','Low-pass filter cut-off frequency','[s^-1]','non-mandatory',"1e99")
        call set_help('dump_space_dt','Spatial parameters saving time interval','[s]','non-mandatory',"1.00e-00")
        call set_help('dump_time_dt','Temporal parameters saving time interval','[s]','non-mandatory',"1.00e-00")
        call set_help('dump_restart_dt','Restart file saving time interval','[s]','non-mandatory',"1")
        call set_help('steady_state','Steady state yes|no','none','non-mandatory',"no")
        call set_help('temp_ramp_filename','Temperature ramp data file (ramp.dat) used if .ne. 0','none','non-mandatory',"none")
        call set_help('solve_heat_equation','yes|no :solve heat equation','none','non-mandatory',"no")
        call set_help('# ********* Grid parameters *********')
        call set_help('n_cells','Number of cells','none','mandatory',"100")
        call set_help('cell_scaling_factor','Cell width scaling factor','none','mandatory',"1.15305056")
        call set_help('#')
        call set_help('# ********* Parameters for volumetric model and species ********* ')
        call set_help('#')
        call set_help('n_species','Number of species','none','mandatory',"3")
        call set_help('species_name','Name of species','none','mandatory',"D D D","species")
        call set_help('n0_profile','Initial density profile is Gaussian (G), step(S), linear(L),peak(P),flat(F)','none'&
            ,'non-mandatory',"F F F","species")
        call set_help('n0_max','Initial density n0','[m^-3]','non-mandatory',"1.00e+10 1.00e+10 1.00e+10","species")
        call set_help('n0_xmax','Position of maximum of Gaussian profile','[m]','non-mandatory',&
        "0.00e-00 0.00e-00 0.00e-00","species")
        call set_help('n0_width','Standard deviation of Gaussian profile','[m]','non-mandatory',&
        "0.00e-00 0.00e-00 0.00e-00","species")
        call set_help('n_max','Maximum density of species','[m^-3]','non-mandatory',"1.00e+29 1.00e+29 1.00e+29","species")

        call set_help('D0','Diffusion constant of species','[m^2 s^-1]','non-mandatory',"1.00e-07 0.00e-07 0.00e-07","species")
        call set_help('ED','Activation energy of diffusion of species','[eV]','non-mandatory',"0.41 0.00 0.00","species")
        call set_help('Edt', 'Activation energy of detrapping of species', '[eV]' , 'non-mandatory' ,"0.00 0.00 1.40","species")
        call set_help('Etr', 'Activation energy of trapping of species'  , '[eV]' , 'non-mandatory' ,"0.00 1.00 0.00","species")

        call set_help('# ********* Parameters for surface model and species ********* ')
        call set_help('surface_model','B: Gamamaout=Kdes*cb^2 S: Gammaout=Kcs^2','[m^-2]','mandatory',"S S S","species")
        call set_help('ns0_left','Initial left surface density of species','[m^-2]','non-mandatory',&
            "1.00e+19 0.00e+19 0.00e+19","species")
        call set_help('ns0_right','Initial right surface density of species','[m^-2]','non-mandatory',&
            "1.00e+19 0.00e+19 0.00e+19","species")
        call set_help('ns_max','Maximum surface density of species','[m^-2]','non-mandatory',"1.00e+19 1.00e+19 1.00e+19","species")
        call set_help('Ech_left'  ,'Activation energy of chemisorption of species at left surface','[eV]','non-mandatory'&
            ,"0.10 0.00 0.00","species")
        call set_help('Ech_right' ,'Activation energy of chemisorption of species at right surface','[eV]','non-mandatory'&
            ,"0.10 0.00 0.00","species")
        call set_help('Q_ch_left' ,'Heat of chemisorption of species at left surface','[eV]','mandatory',&
        "1000 0.00 0.00","species")
        call set_help('Q_ch_right','Heat of chemisorption of species ast right surface','[eV]','mandatory',&
        "1000 0.00 0.00","species")
        call set_help('Eab_left'  ,'Activation energy of absorption of species at left surface','[eV]','mandatory'&
        ,"2.00 0.00 0.00","species")
        call set_help('Eab_right' ,'Activation energy of absorption of species at right surface','[eV]','mandatory'&
            ,"2.00 0.00 0.00","species")
        call set_help('Es_left','Energy of solution of species at left surface','[eV]','mandatory',"1.00 0.00 0.00","species")
        call set_help('Es_right','Energy of solution of species at right surface','[eV]','mandatory',"1.00 0.00 0.00","species")
        call set_help('nu0','Transition attempt frequency of species','[s^-1]','non-mandatory'&
        ,"1.00e+13 1.00e+13 1.00e+13","species")
        call set_help('# ********* Parameters for implantation model ********* ')
        call set_help('implantation_model','S: Step E: ERFC T:Trim','none','non-mandatory',"S","species")
        call set_help('Eimpact_ion','Impact energy of ionized species','[eV]','non-mandatory',&
        "0.00e+00 0.00e+00 0.00e+00","species")
        call set_help('Gamma_min','Minimal external flux of ionized species','[m^-2 s^-1]','non-mandatory',&
        "0.00e+00 0.00e+00 0.00e+00","species")
        call set_help('Gamma_max','Maximal external flux of ionized species','[m^-2 s^-1]','non-mandatory',&
        "0.00e+00 0.00e+00 0.00e+00","species")
        call set_help('pressure_neutral','External pressure of neutral species','[Pa]','non-mandatory',&
        "0.00e-00 0.00e+00 0.00e+00","species")
        call set_help('temp_neutral','External temperature of neutral species','[eV]','non-mandatory',&
        "0.00e-00 0.00e+00 0.00e+00","species")
        call set_help('mass','Mass of species','[kg]','non-mandatory',"3.343e-27 0.00e+00 0.00e+00","species")
        call set_help('# ********* Parameters for abliation model ********* ')
        call set_help('min_ablation_velocity','min ablation speed in addition to sputtering','[m s^-1]','non-mandatory'&
            ,"0.00e-00")
        call set_help('max_ablation_velocity','max ablation speed in addition to sputtering','[m s^-1]','non-mandatory'&
            ,"0.00e-00 ")
        call set_help('sputtering_yield','Sputtering yield','none','non-mandatory',"0")
        call set_help('# ********* Parameters for temperature ********* ')
        call set_help('mat_temp_ramp_start','Material temperature at ramp start (initial at left boundary)','[K]'&
            ,'non-mandatory',"373")
        call set_help('mat_temp_ramp_stop','Material temperature at ramp stop  (initial at right boundary)','[K]',&
           'non-mandatory',"373")
           call set_help('# ********* Material parameters ********* ')
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
        call set_help('pulse_period','Pulse period','[s]','non-mandatory',"1.0e+99")
          call set_help('#Miscelleneaous')
        call set_help('verbose','yes|no','none','non-mandatory',"yes")
        call set_help('dump_space_append','yes|no','none','non-mandatory',"no")
        call set_help('dump_time_append','yes|no ','none','non-mandatory',"no")
        call set_help('fluid-interface','yes|no','none','non-mandatory',"no")
        if (verbose_init) call print_milestone('initialization help completed')
    end subroutine init_help

    subroutine  set_help_single(keyword,def,units,status,default)
        character(*):: keyword,def,units,status,default
        help(ihelp)%keyword=trim(keyword)
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
        call StripFrontSpaces(comment)
        help(ihelp)%def=trim(comment)
        help(ihelp)%comment=.true.
        ihelp=ihelp+1
        Nhelp=ihelp-1
    end subroutine set_help_comment

    subroutine  set_help_species(keyword,def,units,status,default,species)
        character(*):: keyword,def,units,status,default,species
        help(ihelp)%keyword=trim(keyword)
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
        character(150)::str_keyword
        character(100)::strdef
        character(14)::strtype
        character(3*l):: strvalue
        character(3*(l+1)):: strvaluel

        character(41):: strsta,struni
        character(60)::fmt,time
        character(*)::filename
        character(l)::str1

        character(3*l)::strdata,str2
        character(l)::strtmp(3)
        integer::i,j,idefault,ios

        idefault=998
        open(unit=idefault, file=trim(filename), iostat=ios,action='write')
        if ( ios /= 0 ) then
            write(iout,*)'ERROR: Cannot write into default input file "', trim(filename)
            stop'Exiting FACE'

        endif

        write(iout,*)'Default input file "', trim(filename) ,'" created'
        call timestring ( time )
        write(idefault,*) '#Default input for FACE. Created: ', time
        do i=1,Nhelp
            if (.not.help(i)%comment) then
                write(fmt,*)'(a',3*l,')'
                ! if keyword is n_species then select if the input file should be with only H (1 spc) or H+Tr (3 species)
                if (help(i)%keyword.eq."n_species") then
                    if (mode.eq."H") then
                        write(strdata,fmt) "1"
                    elseif (mode.eq."H+Tr") then
                        write(strdata,fmt) "3"
                    else
                        write(iout,*) "ERROR: Unknown mode when writing default input file"
                        STOP 'Exiting FACE...'
                    endif
               else
                    write(strdata,fmt) adjustl(help(i)%default)
                endif

                strvalue=''
                write(str_keyword,*) help(i)%keyword
                write(strdef,*) adjustl(help(i)%def)
                write(struni,*) adjustl(help(i)%units)



                if (.not.(help(i)%species)) then
                    write(strtype,'(a)')' , single'
                    write(strsta,'(a1,a13,a8,a1)') '(',adjustr(help(i)%status) ,adjustl(strtype),') '
                         !   write(strsta,*) '(',adjustl(help(i)%status) ,adjustl(strtype),')'

                    write(fmt,*)'(a30,a1,a',l+1,',a1,a80,a23,a1,a12,a11,a)'
                    write(idefault,fmt) str_keyword,' ',adjustl(strdata), ' ! ',strdef,strsta,' ',adjustl(struni),&
                    ' default: ',adjustl(help(i)%default)
                else
                    write(strtype,'(a)') ' , 1:n_species'
                    write(strsta,'(a1,a13,a13,a1)') '(',adjustr(help(i)%status) ,adjustl(strtype),') '
                    str2=adjustl(strdata)

                    do j=1,3
                        call SplitString(str2,str1,str2,' ')

                        write(fmt,*)'(a14,a1)'
                        write(strtmp(j),fmt) (str1),' '
                    enddo
                    write(strvaluel,*) ((strtmp(j)),j=1,3)
                    write(fmt,*)'(a30,a1,a',3*(l),'a1,a80,a30,a1,a12,a11,a)'
                    write(idefault,fmt) str_keyword,' ',adjustl(strvaluel), ' ! ',strdef,trim(strsta),' ',adjustl(struni),&
                    ' default: ',adjustl(help(i)%default)
                endif
            else
                write(idefault,*) help(i)%def
            endif

        enddo
        !
        !       write(*,*,'advance=no') adjustl(str_keyword),' ' strvalue ' ','! ',strtype, help(i)%def, ' (',help(i)%status ,') ',help(i)%units,help(i)%
        !       type(string)
        !       enddo
        close(idefault)

    end subroutine




end module modFACE_help
