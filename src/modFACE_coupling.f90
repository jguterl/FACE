module modFACE_coupling
    use modFACE_header
    use modFACE_error
    implicit none


contains

    subroutine fluidcode2FACE(FACE_input,fluidcode_input)
        type(FACE_inputs), intent(out) :: FACE_input
        type(fluidcode_inputs), intent(in) :: fluidcode_input

        FACE_input%couple_fluidcode=.true.
        FACE_input%casename=fluidcode_input%casename
        FACE_input%logfile=fluidcode_input%log_file
        FACE_input%input_filename=fluidcode_input%input_file
        FACE_input%run_mode='default'
        FACE_input%path=fluidcode_input%path

        face_input%fluidcode_input=fluidcode_input
    end subroutine fluidcode2FACE

    subroutine FACE2fluidcode(FACE_output,fluidcode_output)
        type(FACE_outputs), intent(in) :: FACE_output
        type(fluidcode_outputs), intent(out) :: fluidcode_output
        ! setup output info for fluid code
        fluidcode_output%nspc_face=face_output%nspc
        call alloc_fluidcode_output(fluidcode_output)
        fluidcode_output%final_inventory(1:nspc)=fluidcode_output%final_inventory(1:fluidcode_output%nspc_face)
        fluidcode_output%init_inventory(1:nspc)=fluidcode_output%init_inventory(1:fluidcode_output%nspc_face)
        fluidcode_output%particle_balance=FACE_output%particle_balance
        ! the statement below assumes that only type of species come form the fludi code (species kk=1 <=> species k=1).
        ! But latter we have two populations of H coming from the fluid code
        fluidcode_output%outgassing_flux=FACE_output%outgassing_flux

        fluidcode_output%init_wall_temp=FACE_output%init_wall_temp
        fluidcode_output%final_wall_temp=FACE_output%final_wall_temp
    end subroutine FACE2fluidcode

    subroutine input_fluidcode(fluidcode_input)
        ! overwrite some input parameters from the input file read by FACE with input from the fluid code:
        integer kk,k
        type(fluidcode_inputs), intent(in) :: fluidcode_input
        ! for the moment, only one type of species (free H) is assumed to be implanted in the wall from fluid code.
        ! That can be changed in the future but it requires some udpates of the FACE code.
        if (fluidcode_input%nspc_fluid.gt.1) then
            call face_error(" This current version of FACE does not allow more than one incoming species",&
                "from the fluid code. Update FACE to allow multi-species implantation")
        endif

        if (verbose_couple) write(iout,*) '---- overwriting fluid code input:'

        ! - casename
        casename=fluidcode_input%casename
        if (verbose_couple) write(iout,*) '- casename overwritten :', casename
        ! temperature solver
        solve_heat_eq_string=fluidcode_input%solve_heat_eq_string
        if (verbose_couple) write(iout,*) '- solve_heat_eq overwritten :', solve_heat_eq_string
        ! start time
        start_time=fluidcode_input%time
        if (verbose_couple) write(iout,*) '- start time overwritten : ', start_time
        ! end time
        end_time=start_time+fluidcode_input%dt
        if (verbose_couple) write(iout,*) '- end time overwritten   : ', end_time
        ! end time
        dt0_face=fluidcode_input%dt0_face
        if (verbose_couple) write(iout,*) '- dt_face overwritten   : ', dt_face

        ! ** setting data dumping parameters
        ! space
        if (fluidcode_input%Ndump_space.le.0) then
            dump_space=.false.
        elseif(fluidcode_input%Ndump_space.gt.2) then ! if dumping more than 2 times, then we dump N-2 times since we dump automaticallty at start_time annd end_time
            dump_space_dt=(end_time-start_time)/real(fluidcode_input%Ndump_space-2,DP)
        else ! if dumping 1 or 2 times, then we start at start_time and end_time
            dump_space_dt=1d99
        endif
        ! time
        if (fluidcode_input%Ndump_time.le.0) then
            dump_time=.false.
        elseif(fluidcode_input%Ndump_time.gt.2) then
            dump_time_dt=(end_time-start_time)/real(fluidcode_input%Ndump_time-2,DP)
        else ! if dumping 1 or 2 times, then we start at start_time and end_time
            dump_time_dt=1d99
        endif
        ! restart
        ! no restart file in coupling mode (we assume that we restart from beginning of fluid code timestep which is beginning of a FACE simulation.)
        dump_restart=.false.



        ! state_files
        read_state_file=fluidcode_input%read_state_file
        final_state_file=fluidcode_input%final_state_file
        if (verbose_couple) write(iout,*) '- final_state_file overwritten :', trim(final_state_file)
        if (verbose_couple) write(iout,*) '- read_state_file overwritten :', trim(read_state_file)





        ! check if species name in FACE input file matches impigning species name from fluid code input
        ! but first check if the numbers of species between fluide code and FACe is consistent:

        if (fluidcode_input%nspc_fluid.gt.nspc) then
            call face_error("N species in fluid code input must be <= to N species in FACE")
        endif
        if (verbose_couple) write(iout,*) '- checking if impinging species from fluide codes are the same than', &
            ' the ones set in FACE :'

        do kk=1,fluidcode_input%nspc_fluid
            k=fluidcode_input%indexspc(kk)
            if (fluidcode_input%namespc(kk).ne.namespc(k)) then
                write(iout,*) 'Mismatch in species names between FACE and fluid code'
                write(iout,*) 'species fluide code index=',kk
                write(iout,*) 'species FACE k index=',fluidcode_input%indexspc(kk)
                write(iout,*) 'species fluide code name=',fluidcode_input%namespc(kk)
                write(iout,*) 'species FACE name=',namespc(k)
                call face_error('Issue with coupling')
            endif
        enddo
        if (verbose_couple) write(iout,*) '- checking species name done'

        ! For the moment, we do not allow erosion of the wall and its associated effects of desorption and heat release when running coupled with a fluide code.
        ! This can be changed in the future if needed
        if (cero.ne.0) then
            call face_error('We do not allow FACE to run coupled with fluid code when erosion is on. See coupling module...')
        endif

        ! ** Setup particle flux by overwritting inflx_in
        ! *** first let be sure that pulsed_plasma is not required and thatn ot flux of incomping particles are imposed in the input file
        ! This is not mandatory but we prevent mishandling of input files by user....
        do k=1,nspc
            if (inflx_in(k).ne.0.d0) then
                call face_error(' particles flx in input file must be =0 when running FACE coupled to fluid code: spc k=',&
                    k,'enrg=',inflx_in(k))
            endif
        enddo
        if (pulsed_flux.ne."no") then
            call face_error('Cannot run pulsed fluxwhen running FACE coupled to fluid code')
        endif
         ! ** Overwritting inflx_in with input from fluid code.

        do kk=1,fluidcode_input%nspc_fluid
            k=fluidcode_input%indexspc(kk)
            inflx_in(k)=fluidcode_input%inflx_in(kk)
            ! check if the value of the influx is not negative
            if (fluidcode_input%inflx_in(kk).lt.0) then
                call face_error("negative influx from fluid code: kk=",kk,"flux=",fluidcode_input%inflx_in(kk))
            endif

        enddo

        ! temperature/heat flux input
        ! first check that framp is not called because FACE cannot read temperature from ramp file and get Temp/Q flux input from fluid code at the same time...
        if (framp.ne."none") then
            call face_error("Incompatibility with fluidcode coupling: cannot read external temperature ramp file: framp=",framp)
        endif
        ! also check that no linear temperature ramp are imposed from the input file
        if (tramp1.ne.0d0.or.tramp0.ne.0d0) then
            call face_error("Incompatibility with fluidcode coupling: cannot set temperature  ramp from input file:",&
                "tramp0 and tramp1 must be equal to zero: tramp0=",tramp0," tramp1",tramp1)
        endif

        if (solve_heat_eq_string.eq."yes") then ! if solving heat equation then impose heat flux from fluid code
            ! ** heat flux calculated in FACE as qflx_in=ee*energy*inflx where energy is in [eV].
            ! We thus calculated the energy correspondign to this expression using inflx and qflx_in from the fluid code.
            ! We assume that the entire heat flux is carried by the first species

            ! ** But first, we need to be sure that the user did not impose any particle energy and flux.
            ! This is not mandatory but it may help to avoid some mishandling of the FACE input files with fluid codes
            do k=1,nspc
                if (enrg(k).ne.0) then
                    call face_error('energy of incoming particles must =0 when running FACE coupled to fluid code: spc k='&
                        ,k,'enrg=',enrg(k))
                endif
                if (inflx_in(k).ne.0) then
                    call face_error('flux of incoming particles must be =0 when running FACE coupled to fluid code: spc k='&
                        ,k,'enrg=',inflx_in(k))
                endif
            enddo

            ! ** Calculating enrg=qflx_in/flx_in/ee assuming that all the heat flux is carried by the first species
            do k=1,nspc
                if (k.eq.1) then
                    if (inflx_in(k).eq.0d0) then ! if influx first species =0 then we set it equal to 1 (should not add much particles!-> 1 m^-2 s^-1)
                        inflx_in(k)=1d0
                    endif
                    enrg(k)=fluidcode_input%qflx_in/inflx_in(k)
                else
                    enrg(k)=0
                endif

                enrg(k)=enrg(k)/ee ! adjust units from Watts to eV
            enddo

        else
            ! if not solving the heat eq then we just set the entire wall temperature to the temperature from the fluid code
            temp0=fluidcode_input%tempwall
            temp1=fluidcode_input%tempwall
            restore_state_temp=.false. ! do not restore temperature from state file if temperature imposed by Face code
            ! check that the temperature is positive
            if (fluidcode_input%tempwall.le.0d0) then
                call face_error('temperature from fluid code <=0": T=',fluidcode_input%tempwall)
            endif
        endif

    end subroutine input_fluidcode


    subroutine alloc_fluidcode_input(fluidcode_input)
        type(fluidcode_inputs)::fluidcode_input

        if (.not.allocated(fluidcode_input%namespc))allocate(fluidcode_input%namespc(fluidcode_input%nspc_fluid))
        if (.not.allocated(fluidcode_input%indexspc))allocate(fluidcode_input%indexspc(fluidcode_input%nspc_fluid))
        if (.not.allocated(fluidcode_input%inflx_in))allocate(fluidcode_input%inflx_in(fluidcode_input%nspc_fluid))
        if (.not.allocated(fluidcode_input%Emean))allocate(fluidcode_input%Emean(fluidcode_input%nspc_fluid))
    end subroutine alloc_fluidcode_input

    subroutine dealloc_fluidcode_input(fluidcode_input)
        type(fluidcode_inputs)::fluidcode_input
        if (allocated(fluidcode_input%namespc))deallocate(fluidcode_input%namespc)
        if (allocated(fluidcode_input%indexspc))deallocate(fluidcode_input%indexspc)
        if (allocated(fluidcode_input%inflx_in))deallocate(fluidcode_input%inflx_in)
        if (allocated(fluidcode_input%Emean))deallocate(fluidcode_input%Emean)
    end subroutine dealloc_fluidcode_input

    subroutine alloc_fluidcode_output(fluidcode_output)
        type(fluidcode_outputs)::fluidcode_output
        if (.not.allocated(fluidcode_output%init_inventory))allocate(fluidcode_output%init_inventory(fluidcode_output%nspc_face))
        if (.not.allocated(fluidcode_output%final_inventory))allocate(fluidcode_output%final_inventory(fluidcode_output%nspc_face))
    end subroutine alloc_fluidcode_output

    subroutine alloc_FACE_output(FACE_output)
        type(FACE_outputs)::FACE_output
        if (.not.allocated(FACE_output%init_inventory))allocate(FACE_output%init_inventory(FACE_output%nspc))
        if (.not.allocated(FACE_output%final_inventory))allocate(FACE_output%final_inventory(FACE_output%nspc))
    end subroutine alloc_FACE_output

    subroutine dealloc_fluidcode_output(fluidcode_output)
        type(fluidcode_outputs)::fluidcode_output
        if (allocated(fluidcode_output%init_inventory))  deallocate(fluidcode_output%init_inventory)
        if (allocated(fluidcode_output%final_inventory)) deallocate(fluidcode_output%final_inventory)
    end subroutine dealloc_fluidcode_output



end module modFACE_coupling
