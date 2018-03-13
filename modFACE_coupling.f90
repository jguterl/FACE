module modFACE_coupling
use modFACE_header
use modFACE_interface
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

 end subroutine FACE2fluidcode

subroutine fluidcode2FACE_input(fluidcode_input)
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
! state_files
restore_state_file=fluidcode_input%restore_state_file
if (verbose_couple) write(iout,*) '- store_state_file overwritten :', store_state_file
store_state_file=fluidcode_input%store_state_file
if (verbose_couple) write(iout,*) '- restore_state_file overwritten :', restore_state_file

! temperature input mode
if (solve_heat_eq.eq."yes") then
elseif  (solve_heat_eq.eq."no") then
else
write(iout,*) 'ERROR: unknown mode for solve_heat_eq fron fluide code'
stop 'Exiting Face'
endif


! check if species name in FACE input file matches impigning species name from fluid code input
! but first check if the numbers of species between fluide code and FACe is consistent:

if (fluidcode_input%nspc.gt.nspc) then
write(iout,*) "N species in fluid code input must be <= to N species in FACE"
stop 'Exiting FACE'
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
stop 'Exiting FACE'
endif
enddo
if (verbose_couple) write(iout,*) '- checking species name done'
end subroutine

    end module modFACE_coupling
