module modFACE_error
use modFACE_precision
use modFACE_header
implicit none
interface face_error
module procedure print_error1
module procedure print_error2
module procedure print_error3
module procedure print_error4
module procedure print_error5
module procedure print_error6
end interface face_error
contains
subroutine print_error1(str)
character(*)::str
write(iout,*) "ERROR:", str
write(iout,*) "Abnormal exit of FACE..."
stop
end subroutine print_error1

subroutine print_error2(str,r)
real(DP) :: r
character(*)::str
write(iout,*) "ERROR:", str,r
write(iout,*) "Abnormal exit of FACE..."
stop
end subroutine print_error2

subroutine print_error3(str,i)
integer :: i
character(*)::str
write(iout,*) "ERROR:", str,i
write(iout,*) "Abnormal exit of FACE..."
stop
end subroutine print_error3

subroutine print_error4(str,r,str2,r2)
real(DP) :: r,r2
character(*)::str,str2
write(iout,*) "ERROR:", str,r,str2,r2
write(iout,*) "Abnormal exit of FACE..."
stop
end subroutine print_error4

subroutine print_error5(str,i,str2,i2)
integer :: i,i2
character(*)::str,str2
write(iout,*) "ERROR:", str,i,str2,i2
write(iout,*) "Abnormal exit of FACE..."
stop
end subroutine print_error5

subroutine print_error6(str,str2)
character(*)::str,str2
write(iout,*) "ERROR:", str,str2
write(iout,*) "Abnormal exit of FACE..."
stop
end subroutine print_error6

end module modFACE_error
