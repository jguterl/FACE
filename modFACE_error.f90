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
module procedure print_error7
module procedure print_error8
module procedure print_error9
module procedure print_error10
end interface face_error
interface face_warning
module procedure print_warning1
module procedure print_warning2
module procedure print_warning3
module procedure print_warning4
module procedure print_warning5
module procedure print_warning6
module procedure print_warning7
module procedure print_warning8
module procedure print_warning9
module procedure print_warning10
end interface face_warning
character(string_length):: message_exit_error="Abnormal exit of FACE..."
contains
! errors
subroutine print_error1(str)
character(*)::str
write(iout,*) "ERROR:", str

if (enforce_error) then
write(iout,*) message_exit_error
stop
endif
end subroutine print_error1

subroutine print_error2(str,r)
real(DP) :: r
character(*)::str
write(iout,*) "ERROR:", str,r
if (enforce_error) then
write(iout,*) message_exit_error
stop
endif
end subroutine print_error2

subroutine print_error3(str,i)
integer :: i
character(*)::str
write(iout,*) "ERROR:", str,i
if (enforce_error) then
write(iout,*) message_exit_error
stop
endif
end subroutine print_error3

subroutine print_error4(str,r,str2,r2)
real(DP) :: r,r2
character(*)::str,str2
write(iout,*) "ERROR:", str,r,str2,r2
if (enforce_error) then
write(iout,*) message_exit_error
stop
endif
end subroutine print_error4

subroutine print_error5(str,i,str2,i2)
integer :: i,i2
character(*)::str,str2
write(iout,*) "ERROR:", str,i,str2,i2
if (enforce_error) then
write(iout,*) message_exit_error
stop
endif
end subroutine print_error5

subroutine print_error6(str,str2)
character(*)::str,str2
write(iout,*) "ERROR:", str,str2
if (enforce_error) then
write(iout,*) message_exit_error
stop
endif
end subroutine print_error6

subroutine print_error7(str,i,str2,r,str3,r2)
integer ::i
real(DP) :: r,r2
character(*)::str,str2,str3
write(iout,*) "ERROR:", str,i,str2,r,str3,r2
if (enforce_error) then
write(iout,*) message_exit_error
stop
endif
end subroutine print_error7

subroutine print_error8(str,i,str2,i2,str3,r,str4,r2)
integer ::i,i2
real(DP) :: r,r2
character(*)::str,str2,str3,str4
write(iout,*) "ERROR:", str,i,str2,i2,str3,r,str4,r2
if (enforce_error) then
write(iout,*) message_exit_error
stop
endif
end subroutine print_error8

subroutine print_error9(str,str2,str3)
character(*)::str,str2,str3
write(iout,*) "ERROR:", str,str2,str3
if (enforce_error) then
write(iout,*) message_exit_error
stop
endif
end subroutine print_error9

subroutine print_error10(str,str2,str3,i)
character(*)::str,str2,str3
integer :: i
write(iout,*) "ERROR:", str,str2,str3,i
if (enforce_error) then
write(iout,*) message_exit_error
stop
endif
end subroutine print_error10


! warning
subroutine print_warning1(str)
character(*)::str
write(iout,*) "Warning:", str
stop
end subroutine print_warning1

subroutine print_warning2(str,r)
real(DP) :: r
character(*)::str
write(iout,*) "Warning:", str,r
stop
end subroutine print_warning2

subroutine print_warning3(str,i)
integer :: i
character(*)::str
write(iout,*) "Warning:", str,i
stop
end subroutine print_warning3

subroutine print_warning4(str,r,str2,r2)
real(DP) :: r,r2
character(*)::str,str2
write(iout,*) "Warning:", str,r,str2,r2
stop
end subroutine print_warning4

subroutine print_warning5(str,i,str2,i2)
integer :: i,i2
character(*)::str,str2
write(iout,*) "Warning:", str,i,str2,i2
stop
end subroutine print_warning5

subroutine print_warning6(str,str2)
character(*)::str,str2
write(iout,*) "Warning:", str,str2
stop
end subroutine print_warning6

subroutine print_warning7(str,i,str2,r,str3,r2)
integer ::i
real(DP) :: r,r2
character(*)::str,str2,str3
write(iout,*) "Warning:", str,i,str2,r,str3,r2
stop
end subroutine print_warning7

subroutine print_warning8(str,i,str2,i2,str3,r,str4,r2)
integer ::i,i2
real(DP) :: r,r2
character(*)::str,str2,str3,str4
write(iout,*) "Warning:", str,i,str2,i2,str3,r,str4,r2
stop
end subroutine print_warning8

subroutine print_warning9(str,i,str2,i2,str3,r)
integer ::i,i2
real(DP) :: r
character(*)::str,str2,str3
write(iout,*) "Warning:", str,i,str2,i2,str3,r
stop
end subroutine print_warning9

subroutine print_warning10(str,i,str2,r)
integer ::i
real(DP) :: r
character(*)::str,str2
write(iout,*) "Warning:", str,i,str2,r
stop
end subroutine print_warning10

end module modFACE_error
