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
type(fluidcode_inputs), intent(in) :: fluidcode_input
! - casename
casename=fluidcode_input%casename
! temperature solver
solve_heat_eq=fluidcode_input%solve_heat_eq
! start time
start_time=fluidcode_input%time
! end time
end_time=start_time+fluidcode_input%dt
!
! temperature input mode
if (solve_heat_eq.eq."yes") then
elseif  (solve_heat_eq.eq."no") then
else
write(iout,*) 'ERROR: unknown mode for solve_heat_eq fron fluide code'
stop 'Exiting Face'
endif

! history_file
restore_history_file=fluidcode_input%history_file

!fluidcode_input%namespc(1:fluidcode_input%nspc)
! -
if (verbose_couple) then
write(iout,*) '---- fluid code input:'
write(iout,*) '- casename overwritten :', casename
write(iout,*) '- solve_heat_eq overwritten :', solve_heat_eq
write(iout,*) '- history_file overwritten :', restore_history_file
write(iout,*) '- start time overwritten : ', start_time
write(iout,*) '- end time overwritten   : ', end_time
write(iout,*) '- dt_face overwritten    : ', dt_face
write(iout,*) '- checking if impinging species from fluide codes are the same than the ones set in FACE :'

endif
end subroutine

    end module modFACE_coupling
