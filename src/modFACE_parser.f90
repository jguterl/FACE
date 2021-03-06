module modFACE_parser
    use modFACE_output
    use modFACE_header
    use modFACE_help
    use modFACE_misc
    use modFACE_output
    implicit none
    public
    type input_valsp
        real(DP),allocatable:: r(:)
        integer,allocatable :: i(:)
        character(string_length+1),allocatable::s(:)
        character(string_length+1)::status
    end type input_valsp

    type input_val
        real(DP):: r
        integer :: i
        character(string_length+1)::s
        character(string_length+1)::status
    end type input_val

    character, allocatable :: input_line(:)
    integer,parameter::iuinput=785
    integer::nlines,length_line
    character::delimdata=' '
    integer,parameter::maxlength_input=2000
    type string
        character(string_length) :: keyword
        character(string_length) :: data
        character(string_length) :: comment
    end type string
    type(string), allocatable :: input_lines(:)


contains



    function find_keyword(keyword) result(idx)

        character(*)::keyword
        integer:: idx
        integer::i
        logical:: found_str

        found_str=.false.
        idx=-1
        do i=1,nlines
            if (compare_string(keyword,input_lines(i)%keyword)) then
                !    ! check that this keywoard only exists once in the input file
                if (found_str) then
                    write(iout,*) 'ERROR: keyword "', keyword ,'" found twice in the input file,idx=',i
                        call face_error('input_lines(i)%keyword=',input_lines(i)%keyword)
                endif
                idx=i
                found_str=.true.
                if (verbose_parser) write(iout,*) 'find keyword "', keyword ,'" at entry #',idx
            endif
        enddo

    end function




    subroutine check_allocate_species()

        if (nspc<1) then
            write(iout,*) 'ERROR: n_species must be defined in the inputfile before any input for species'
            stop
        endif

    end subroutine check_allocate_species



    subroutine assign_keyword_value(keyword,inputval,typeval)

        character(*)::keyword
        type(input_val) ::inputval
        character::typeval
        integer::idx,idx_help
        character(string_length+1)::str0

        idx=find_keyword(keyword)
        idx_help=find_keyword_help(keyword)
        if (idx_help.eq.-1) then
            call face_error('Unknown keyword "', trim(keyword) ,&
            '" (see module help for list of keywords or run Face with -keywords')
        endif
        if (verbose_parser) write(iout,*) 'keyword:',trim(keyword),' idx=',idx,' idx_help=',idx_help
        write(inputval%status,*) adjustl(help(idx_help)%status)

        if (idx.eq.-1.AND.inputval%status.eq.'mandatory') then
            call face_error('Unknown keyword "', trim(keyword) ,&
            '" (see module help for list of keywords or run Face with -keywords')
        else if (idx.ne.-1) then
            ! seems that an extra blank is added when writing in string, leading to seg fault if size of receiving string is not incremeted by +1
            write(str0,*) input_lines(idx)%data
            if (verbose_parser) write(iout,*) 'getting data from string:',trim(str0)
            call get_single_data(str0,delimdata)

        else
            ! setting default value for keyword
            write(str0,*) help(idx_help)%default
            if (verbose_parser) write(iout,*) 'Setting default value for keyword "' ,trim(keyword),'"'
        endif

        if (verbose_parser) write(iout,*) 'get data for keyword:',trim(keyword)

        if (typeval.eq."r") then
            read(str0,*) inputval%r
        elseif (typeval.eq."i") then
            read(str0,*) inputval%i
        elseif (typeval.eq."s") then
            inputval%s=trim(str0)
        endif
if (verbose_parser) then
if (typeval.eq."r") then
                write(iout,*) "keyword : ",trim(keyword), "; values=",inputval%r
                elseif (typeval.eq."i") then
                write(iout,*) "keyword : ",trim(keyword), "; values=",inputval%i
                elseif (typeval.eq."s") then
                write(iout,*) "keyword : ",trim(keyword), "; values=",trim(inputval%s)
                endif
                 write(iout,*) ' '
endif
        return
    end subroutine assign_keyword_value

    subroutine assign_keyword_value_species(keyword,inputval,typeval)
        character(*)::keyword
        type(input_valsp) ::inputval
        character::typeval
        integer::idx,idx_help,k
        character(string_length+1)::str0,str1

        idx=find_keyword(keyword)
        idx_help=find_keyword_help(keyword)
        if (idx_help.eq.-1) then
            call face_error('Unknown keyword "', keyword ,'" (see module help for list of keywords or run Face with -keywords')
        endif
        write(inputval%status,*) adjustl(help(idx_help)%status)
        if (idx.eq.-1.AND.inputval%status.eq.'mandatory') then
            call face_error('mandatory keyword "', keyword ,'" not found in the inputfile')
        elseif (idx.ne.-1) then
            str0=trim(input_lines(idx)%data)
            !write(iout,*) "str0=",trim(str0)

            do k=1,nspc
                call get_multiple_data(str0,str1,delimdata)
                !write(iout,*) "str0=",trim(str0),"; str1=",str1
                if (typeval.eq."r") then
                    read(str0,*)  inputval%r(k)
                elseif (typeval.eq."i") then
                    read(str0,*) inputval%i(k)
                elseif (typeval.eq."s") then
                inputval%s(k)=trim(str0)
                endif
                str0=trim(str1)
            enddo

if (verbose_parser) then
if (typeval.eq."r") then
                write(iout,*) "keyword : ",trim(keyword), "; values=",(inputval%r(k),k=1,nspc)
                elseif (typeval.eq."i") then
                write(iout,*) "keyword : ",trim(keyword), "; values=",(inputval%i(k),k=1,nspc)
                elseif (typeval.eq."s") then
                write(iout,*) "keyword : ",trim(keyword), "; values=",(trim(inputval%s(k)),k=1,nspc)

                endif
                write(iout,*) ' '
endif
        else
            ! setting default value for keyword
            write(str0,*) help(idx_help)%default
            do k=1,nspc
                call get_multiple_data(str0,str1,delimdata)
                if (typeval.eq."r") then
                    read(str0,*) inputval%r(k)
                elseif (typeval.eq."i") then
                    read(str0,*) inputval%i(k)
                elseif (typeval.eq."s") then
                    inputval%s(k)=str0
                endif
                str0=str1
            enddo
            if (verbose_parser) write(iout,*) 'Setting default value for keyword "' ,keyword,'"'

        endif
        return

    end subroutine assign_keyword_value_species



    subroutine open_inputfile(filename)
        integer ::ios
        character(*),intent(in)::filename

        open(unit=iuinput, file=trim(filename), iostat=ios,action='read')
        if ( ios /= 0 ) then
            write(iout,*) 'Opening of input file "', trim(filename) ,'" : FAIL '
            stop

        endif
        if (verbose_input)  write(iout,*) 'Opening of input file "', trim(filename) ,'" : DONE '
    end subroutine open_inputfile




    subroutine get_nlines()
        integer:: ios
        character(maxlength_input)::line
        nlines = 0
        length_line=0
        do
            read(iuinput, '(A)', iostat=ios) line
            length_line=max(length_line,len_trim(line))
            if (ios /= 0) exit
            nlines = nlines + 1
        enddo


        if (verbose_parser) write(iout,*) "Input file contains " ,nlines," lines"
        if (verbose_parser) write(iout,*) "Max line length " ,length_line," characters"

        if (nlines.lt.1) stop "Error: # lines read in input file <1"

    end subroutine get_nlines

    subroutine get_inputlines()

        integer::i,ios
        character(length_line)::str0,str1,str2,str3,str4

        str0=''
        str1=''
        str2=''
        str3=''
        str4=''
        do i=1,nlines
        if (verbose_parser) write(iout,*) "Reading line: " ,i
            read(iuinput, '(A)', iostat=ios) str0
            if (verbose_parser) write(iout,*) "str0: " ,str0
            call StripFrontSpaces(str0)
            if (verbose_parser) write(iout,*) "str0: " ,str0
            call SplitString(trim(str0), str1, str2," ")
              if (verbose_parser) write(iout,*) "str1: " ,str1
            call SplitString(str2, str3, str4,"!")
            if (verbose_parser) write(iout,*) "str3: " ,str3
            ! BUG REPORTED WITH CLANG compiler on MACOS for this following line
            !input_lines(i)=string(trim(str1),trim(str3),trim(str4))
            input_lines(i)%keyword=trim(str1)
            if (verbose_parser) write(iout,*) "str2: " ,str2
            input_lines(i)%data=trim(str3)

            input_lines(i)%comment=trim(str4)
if (verbose_parser) write(iout,*) "str3: " ,str3
            if (verbose_parser) write(iout,*) "Line # " ,i,' : kw="',trim(input_lines(i)%keyword), &
            '" data="',trim(input_lines(i)%data),'" comment="',trim(input_lines(i)%comment),'"'
        enddo

    end subroutine get_inputlines



    function find_keyword_help(keyword) result(idx_help)
        character(*)::keyword
        integer:: idx_help
        integer::i
        logical:: found_str
        found_str=.false.
        idx_help=-1
        do i=1,Nhelp

            if (compare_string(keyword,help(i)%keyword)) then
                !    ! check that this keywoard only exists once in the input file
                if (found_str) then
                    call face_error('keyword "', keyword ,'" found twice in the input file,idx=',i)
                endif
                idx_help=i
                found_str=.true.
                if (verbose_parser) write(iout,*) 'find keyword in help "', keyword ,'" at entry #',idx_help
            endif
        enddo

    end function find_keyword_help



end module modFACE_parser



