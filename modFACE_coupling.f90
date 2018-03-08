module modFACE_coupling
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

    end module modFACE_coupling
