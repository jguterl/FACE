module modFACE_misc
use modFACE_output
implicit none

contains
   SUBROUTINE SplitString(instring, string1, string2,delim)
    implicit none
    CHARACTER(*):: instring
    character::delim
    CHARACTER(*),INTENT(OUT):: string1,string2
    INTEGER :: indexx
    indexx = INDEX(instring,delim)
    if (indexx>0) then
    string1 = trim(instring(1:indexx-1))
    string2 = trim(instring(indexx+1:))
    else
    string1=instring
    string2=''
    endif
     END SUBROUTINE SplitString

     subroutine get_single_data(string,delim)

    character::delim
    CHARACTER(*),INTENT(INOUT):: string
    INTEGER :: indexx
    call StripFrontSpaces(string)
    indexx = INDEX(string,delim)
    if (indexx>0) then
    string = trim(string(1:indexx-1))
    endif
     return
     end subroutine

      subroutine get_multiple_data(string,strout,delim)

    character::delim
    CHARACTER(*),INTENT(INOUT):: string,strout
    INTEGER :: indexx

    call StripFrontSpaces(string)

    strout=trim(string)

    indexx = INDEX(trim(string)," ")
    if (indexx>0) then

    strout = trim(string(indexx+1:))
    string = trim(string(1:indexx-1))
    endif
     return
     end subroutine

     subroutine StripFrontSpaces(string)
         character(*) :: string
         integer :: stringLen
         integer :: i,idx=-1
         stringLen= len (string)
         if (stringLen>2) then
             do i=1,stringLen-1
                 if (string(i:i) /= ' ') then
                     idx=i
                     exit
                 endif
             enddo
             if (idx.ge.0) then
             string(1:stringlen-idx)=string(idx:)
             endif
         endif
        ! string=adjustl(string)
     end subroutine StripFrontSpaces

function compare_string(str1,str2) result(match)
  character(*)::str1,str2
      logical::match
            integer::N,lstr1,lstr2,i
      match=.true.
     lstr2=len_trim(str2)
      lstr1=len_trim(str1)
      !write(iout,*) 'L1=',lstr1, 'L2=',lstr2
      N=min(lstr1,lstr2)
      if (N.le.0) then
          match=.false.
      else
          do i=1,N
              if(str1(i:i).ne.str2(i:i)) then
                  match=.false.
              endif

          enddo
      endif

if (lstr1.ne.lstr2) then
match=.false.
endif
! if (len(str1).gt.N) then
!     do i=N+1,len(str1)
!         if (str1(i:i).ne.0) then
!             match=.false.
!         endif
!     enddo
! endif
! !write(iout,*) 'match2=',match
! !
! !
! if (len(str2).gt.N) then
!     do i=N+1,len(str2)
!         if (str2(i:i).eq.) then
!             match=.false.
!         endif
!     enddo
! endif
 !write(iout,*) 'match3=',match
 return

 end function
 end module modFACE_misc
