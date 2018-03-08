!###define ERRORMSG(msg) write(iout,'("Error at ",I4," in file ",/,A,/,"Error message: ",A)') __LINE__,__FILE__,msg
      program FACE
      use modFACE_main
      use modFACE_interface
      implicit none
      type(FACE_inputs)::face_input
      type(fluidcode_inputs)::fluidcode_input
      logical:: couple_test=.false.
!     ******************************************************************
!     * 1-dimensional First wAll simulation CodE "FACE"                *
!     * for modeling of particle transport in the first wall           *
!     * of fusion devices                                              *
!     *                                                                *
!     * Authors: Roman D. Smirnov and Jerome Guterl                                      *
!     * E-mail: rsmirnov@ucsd.edu; rosmirnov@yahoo.com                 *
!     *                                                                *
!     ******************************************************************

      ! selecting input mode:
      ! - couple_test: as if called by fluid code
      ! - otherwise: normal execution with reading of arguments




      if (couple_test) then
      call set_fluidcode_input(fluidcode_input)
      call wrapper_FACE(face_input,fluidcode_input)
      else
      call read_arguments(face_input)
      endif

      ! main call to FACE
      call FACE_main(face_input)

      end program FACE
