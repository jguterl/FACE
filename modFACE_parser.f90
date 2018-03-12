module modFACE_parser
    use modFACE_output,only:iout,verbose_parser
    use modFACE_header
    use modFACE_help
    use modFACE_misc
    use modFACE_output
    implicit none
    public
    type input_valsp
        real(DP),allocatable:: r(:)
        integer,allocatable :: i(:)
        character(len=lname),allocatable::s(:)
        character(len=20)::status
    end type input_valsp

    type input_val
        real(DP):: r
        integer :: i
        character(150),allocatable::s
        character(len=20)::status
    end type input_val

    character, allocatable :: input_line(:)
    integer,parameter::iuinput=785
    integer::nlines,length_line
    character::delimdata=' '
    integer,parameter::maxlength_input=2000
    type string
        character(100), allocatable :: keyword
        character(100), allocatable :: data
        character(100), allocatable :: comment
    end type string
    type(string), allocatable :: input_lines(:)




    interface write_input_log_keyword
        module procedure write_input_log_keyword_r
        module procedure write_input_log_keyword_s
        module procedure write_input_log_keyword_i
        module procedure write_input_log_keyword_species_r
        module procedure write_input_log_keyword_species_s
        module procedure write_input_log_keyword_species_i
    end interface write_input_log_keyword


contains



    function find_keyword(keyword) result(idx)

        character(*)::keyword
        integer:: idx
        integer::i
        logical:: found_str

        found_str=.false.
        idx=-1
        do i=1,nlines
            if (compare_string(keyword,(input_lines(i)%keyword))) then
                !    ! check that this keywoard only exists once in the input file
                if (found_str) then
                    write(iout,*) 'ERROR: keyword "', keyword ,'" found twice in the input file,idx=',i
                    stop 'Exiting FACE'
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
        character(1000)::str0

        idx=find_keyword(keyword)
        idx_help=find_keyword_help(keyword)
        if (idx_help.eq.-1) then
            write(iout,*) 'ERROR: Unknown keyword "', keyword ,'" (see module help for list of keywords or run Face with -keywords)'
            stop 'Exiting FACE...'
        endif

        write(inputval%status,*) adjustl(help(idx_help)%status)

        if (idx.eq.-1.AND.inputval%status.eq.'mandatory') then
            write(iout,*) 'ERROR: Unknown keyword "', keyword ,'" (see module help for list of keywords or run Face with -keywords)'
            stop 'Exiting FACE...'
        else if (idx.ne.-1) then
            str0=''
            write(str0,*) input_lines(idx)%data
            call get_single_data(str0,delimdata)
        else
            ! setting default value for keyword
            write(str0,*) help(idx_help)%default

            if (verbose_parser) write(iout,*) 'Setting default value for keyword "' ,keyword,'"'
        endif
        if (typeval.eq."r") then
            read(str0,*) inputval%r
        elseif (typeval.eq."i") then
            read(str0,*) inputval%i
        elseif (typeval.eq."s") then
            inputval%s=trim(str0)
        endif

        return
    end subroutine assign_keyword_value

    subroutine assign_keyword_value_species(keyword,inputval,typeval)
        character(*)::keyword
        type(input_valsp) ::inputval
        character::typeval
        integer::idx,idx_help,k
        character(1000)::str0,str1

        idx=find_keyword(keyword)
        idx_help=find_keyword_help(keyword)
        if (idx_help.eq.-1) then
            write(iout,*) 'ERROR: Unknown keyword "', keyword ,'" (see module help for list of keywords or run Face with -keywords)'
            stop 'Exiting FACE...'
        endif
        write(inputval%status,*) adjustl(help(idx_help)%status)
        if (idx.eq.-1.AND.inputval%status.eq.'mandatory') then
            write(iout,*) 'ERROR: mandatory keyword "', keyword ,'" not found in the inputfile'
            stop 'Exiting FACE...'
        else if (idx.ne.-1) then
            str0=''
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

        else
            ! setting default value for keyword
            str0=''
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
        write(iout,*) 'Opening of input file "', trim(filename) ,'" : DONE '
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
            read(iuinput, '(A)', iostat=ios) str0
            call StripFrontSpaces(str0)
            call SplitString(trim(str0), str1, str2," ")

            call SplitString(str2, str3, str4,"!")
            ! BUG REPORTED WITH CLANG compiler on MACOS for this following line
            !input_lines(i)=string(trim(str1),trim(str3),trim(str4))
            input_lines(i)%keyword=trim(str1)
            input_lines(i)%data=trim(str3)
            input_lines(i)%comment=trim(str4)

            if (verbose_parser) write(iout,*) "Line # " ,i,' : kw="',input_lines(i)%keyword, '" data="',input_lines(i)%data,&
                '" comment="',input_lines(i)%comment,'"'
        enddo

    end subroutine get_inputlines



    function find_keyword_help(keyword) result(idx_help)
        character(*)::keyword
        integer:: idx_help
        integer::i
        logical:: found_str
        found_str=.false.
        idx_help=-1
        do i=1,nlines
            if (compare_string(keyword,(help(i)%keyword))) then
                !    ! check that this keywoard only exists once in the input file
                if (found_str) then
                    write(iout,*) 'ERROR: keyword "', keyword ,'" found twice in the input file,idx=',i
                    stop 'Exiting FACE'
                endif
                idx_help=i
                found_str=.true.
                if (verbose_parser) write(iout,*) 'find keyword in help "', keyword ,'" at entry #',idx_help
            endif
        enddo

    end function find_keyword_help


    subroutine write_input_log_keyword_r(keyword,variable)

        character(*)::keyword
        real(DP)::variable

        write(ifile_inputlog,'(a60,f8.3)') adjustl(keyword),variable

    end subroutine write_input_log_keyword_r


    subroutine write_input_log_keyword_i(keyword,variable)

        character(*)::keyword
        integer::variable

        write(ifile_inputlog,'(a60,i6.3)') adjustl(keyword),variable

    end subroutine write_input_log_keyword_i


    subroutine write_input_log_keyword_s(keyword,variable)

        character(*)::keyword
        character(*)::variable
        write(ifile_inputlog,'(a60,a60)') adjustl(keyword),variable

    end subroutine write_input_log_keyword_s


    subroutine write_input_log_keyword_species_r(keyword,variable)

        character(*)::keyword
        character*60::myfmt
        real(DP)::variable(:)
        integer::k

        !write(myfmt,*)'(a60, ,',nspc,'(f6.3))'
        write(ifile_inputlog,*) adjustl(keyword),(variable(k),k=1,nspc)

    end subroutine write_input_log_keyword_species_r


    subroutine write_input_log_keyword_species_i(keyword,variable)

        character(*)::keyword
        integer::variable(:)
        integer::k

        write(ifile_inputlog,*) adjustl(keyword),(variable(k),k=1,nspc)
    end subroutine write_input_log_keyword_species_i


    subroutine write_input_log_keyword_species_s(keyword,variable)

        character(*)::keyword
        character(*)::variable(:)
        integer::k

        write(ifile_inputlog,*)adjustl(keyword),(variable(k),k=1,nspc)

    end subroutine write_input_log_keyword_species_s


end module modFACE_parser



