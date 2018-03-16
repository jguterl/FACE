module modFACE_coupling
use modFACE_header
use modFACE_interface
use modFACE_error
contains

!    subroutine write_fluid2wallcode_data(ifile)
!    integer,intent(in)::ifile
!    call header_fluid2wallcode(ifile)
!
!end subroutine write_fluid2wallcode_data

!    subroutine header_fluid2wallcode(ifile)
!    character(16)::str1,str2,str3,str4,str5,str6
!    write(str1,'a16') 'wall index '
!    write(str2,'a16') 'fluid iteration '
!    write(str3,'a16') 'fluid time '
!    write(str4,'a16') 'fuid dt'
!    write(str5,'a16') 'nspecies'
!    write(str6,'a16') 'particle flux'
!    write(str7,'a16') 'ave.particle enrg'
!    write(str8,'a16') 'heat flux'

!    read fluid2wallcode_data()
!    end subroutine
!    write_wallcode2fluid_data()
!    read_FACE2()

!    subroutine FACE2fluidcode()
!
!    end subroutine FACE2fluidcode

 subroutine FACE2fluidcode
 ! send out outgassing flux to fluid code: which one? average over FACE time (tstart to tstop) or final outgassing flux?
 ! should it be steady-state? dt_fluid should be small enough.... to avoid strong variation of jflux out
!if (couple_wallcode) then
!! ?
!endif
!call compute_trace_outgassing(face_output%trace_outgassing)
!call

 end subroutine FACE2fluidcode

subroutine fluidcode2FACE(fluidcode_input)
! overwrite some input parameters from the input file read by FACE with input from the fluid code:
integer kk,k
type(fluidcode_inputs), intent(in) :: fluidcode_input

if (verbose_couple) write(iout,*) '---- overwriting fluid code input:'

! - casename
casename=fluidcode_input%casename
if (verbose_couple) write(iout,*) '- casename overwritten :', casename
! temperature solver
solve_heat_eq=fluidcode_input%solve_heat_eq
if (verbose_couple) write(iout,*) '- solve_heat_eq overwritten :', solve_heat_eq
! start time
start_time=fluidcode_input%time
if (verbose_couple) write(iout,*) '- start time overwritten : ', start_time
! end time
end_time=start_time+fluidcode_input%dt
if (verbose_couple) write(iout,*) '- end time overwritten   : ', end_time
! end time
dt_face=start_time+fluidcode_input%dt_face
if (verbose_couple) write(iout,*) '- dt_face overwritten   : ', dt_face
! state_files
restore_state_file=fluidcode_input%restore_state_file
if (verbose_couple) write(iout,*) '- final_state_file overwritten :', final_state_file
final_state_file=fluidcode_input%final_state_file
if (verbose_couple) write(iout,*) '- restore_state_file overwritten :', restore_state_file





! check if species name in FACE input file matches impigning species name from fluid code input
! but first check if the numbers of species between fluide code and FACe is consistent:

if (fluidcode_input%nspc.gt.nspc) then
call face_error("N species in fluid code input must be <= to N species in FACE")
endif

if (verbose_couple) write(iout,*) '- checking if impinging species from fluide codes are the same than the ones set in FACE :'
do kk=1,fluidcode_input%nspc
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
            call face_error('flux of incoming particles must be =0 when running FACE coupled to fluid code: spc k=',&
            k,'enrg=',inflx_in(k))
        endif
 enddo
 if (pulsed_flux.ne."no") then
 call face_error('Cannot run pulsed fluxwhen running FACE coupled to fluid code')
 endif
  ! ** Overwritting inflx_in with input from fluid code.

 do kk=1,fluidcode_input%nspc
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

if (solve_heat_eq.eq."yes") then ! if solving heat equation then impose heat flux from fluid code
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

elseif  (solve_heat_eq.eq."no") then
! if not solving the heat eq then we just set the entire wall temperature to the temperature from the fluid code
temp0=fluidcode_input%tempwall
temp1=fluidcode_input%tempwall
restore_state_temp=.false. ! do not restore temperature from state file if temperature imposed by Face code
! check that the temperature is positive
if (fluidcode_input%tempwall.le.0d0) then
call face_error('temperature from fluid code <=0": T=',fluidcode_input%tempwall)
endif
else
call face_error('unknown mode for solve_heat_eq fron fluide code')
endif




end subroutine fluidcode2FACE

    end module modFACE_coupling
