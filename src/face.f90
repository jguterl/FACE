#define ERRORMSG(msg) write(iout,'("Error at ",I4," in file ",/,A,/,"Error message: ",A)') __LINE__,__FILE__,msg
      program FACE
      use modFACE_main
      use modFACE_interface
      implicit none

      logical:: couple_fluidcode=.false.

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
      ! - couple_fluid_code: as if called by fluid code
      ! - otherwise: normal execution with reading of arguments




       couple_fluidcode=.true.
      if (couple_fluidcode) then
         ! test of coupling mode'
         call FACE_from_fluidcode

      else
          call FACE_standalone
      endif

      ! main call to FACE



  end program FACE


