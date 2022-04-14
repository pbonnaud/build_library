!   ////////////////////////////////////////////////////////////////////////////////////////////////
!   //                                                                                            //
!   // Copyright (2022) Patrick A. Bonnaud                                                        //
!   //                                                                                            //
!   // This file is part of BUILDLIBRARY (Build Molecular Models for a Library of Molecules).     //
!   //                                                                                            //
!   // BUILDLIBRARY is free software; you can redistribute it and/or modify it under the terms    //
!   // of the GNU General Public License as published by the Free Software Foundation; either     //
!   // version 2 of the License, or (at your option) any later version.                           //
!   //                                                                                            //
!   // BUILDLIBRARY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;  //
!   // without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  //
!   // See the GNU General Public License for more details.                                       //
!   //                                                                                            //
!   // You should have received a copy of the GNU General Public License along with this program. //
!   // If not, see <http://www.gnu.org/licenses/>.                                                //
!   //                                                                                            //
!   ////////////////////////////////////////////////////////////////////////////////////////////////

program BUILDLIBRARY

!   ***********************************************************************************************
!   **              PROGRAM TO BUILD A LIBRARY OF FILES FOR LAMMPS SIMULATIONS                   **
!   ***********************************************************************************************

    use module_physical_constants;

    use module_size_arrays;

    use module_config;

    use module_library;

!   ***********************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4) :: icanal;

    integer (kind=4) :: i, j, k;

    real (kind=8) :: grnd;

    real (kind=8), dimension(1:3) :: CELL_AXIS, CELL_ANGDEG;

    real (kind=8), dimension(1:3,1:3) :: PASSA, PASSB;

    integer (kind=4) :: FLAG_SORT_MOLEC;

    character (len=250) :: CHEXT1, CHEXT2, CHEXT3;

!   ************************************************************************************************

    real (kind=8), dimension(1:6) :: MATA;

!   ### Parametres du substrat #####################################################################
    integer (kind=4) :: ie;

    integer (kind=4) :: ICONF_FINAL;

    character (len=200) :: PATH_DIRECTORY, CH_WORKING_FILE;

!   ### WARNING MESSAGE PARAMETERS #################################################################

    integer (kind=4) :: EOF;

!   ### PARAMETERS THAT WE DON'T CARE ##############################################################

    logical :: PROBE1;

!   ***********************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE, CHPRGM_VERSION, CHDATE;

!   ***********************************************************************************************
                                                                                         !
    icanal = 99;                                                                         !    
                                                                                         !
!   ### Open the output file #######################################################################
                                                                                         !
    open(99,file='OUTPUT');                                                              !
                                                                                         !
!   ### Write the title of the program #############################################################
                                                                                         !
    CHTITLE        = 'Build Files for LAMMPS Library';                                   !
                                                                                         !
    CHPRGM_VERSION = 'v0.07';                                                            !
                                                                                         !
    CHDATE         = '04/05/2021';                                                       !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_PROGRAM_TITLE(99,70,ILENGTH_TITLE,TRIM(CHTITLE),CHPRGM_VERSION,CHDATE);    !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set atom properties as found in the Mendeleiev table #######################################
                                                                                         !
    call SET_ATOM_PROPERTIES(99);                                                        !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read the input file of the current program #################################################
                                                                                         !
    call READ_INPUT(99);                                                                 !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Build the library files ####################################################################
                                                                                         !
    do i = 1, NFILE_LIBRARY;                                                             !
                                                                                         !
        if ( TRIM(CHFILE_FORMAT) == 'lammps' ) then;                                     !
                                                                                         ! 
!           ### Read the input file for lammps #####################################################
                                                                                         !
            call READ_INPUT_LAMMPS(99,                              &                    !
                                   CHNAME_LAMMPS_INPUT_LIBRARY(i));                      !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Read the configuration in the lammps format ########################################
                                                                                         !
            call READ_LAMMPS_CONFIG(99,                      &                           !
                                    CHNAME_FILE_LIBRARY(i),  &                           !
                                    MATA(1:6),               &                           !
                                    FLAG_SORT_MOLEC);                                    !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( TRIM(CHFILE_FORMAT) == 'amber' ) then;                                 !
                                                                                         !
!           ### Read the file containing coordinates of atoms (inpcrd) #############################
                                                                                         !
            call READ_AMBER_INPCRD(99,CHNAME_FILE_LIBRARY(i));                           !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Read the file containing potential parameters (prmtop) #############################
                                                                                         !
            call READ_AMBER_PRMTOP(99,CHNAME_LAMMPS_INPUT_LIBRARY(i));                   !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            call BUILD_PARAMETER_AMBER_TO_LAMMPS(icanal);                                !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else                                                                             !
                                                                                         !
            write(icanal,*);                                                             ! 
                                                                                         !
            write(icanal,'(a23)') 'Stop - Not implemented!';                             !
                                                                                         !
            close(icanal);                                                               !
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
        call BUILD_TEMPLATED_CONFIG(99);                                                 !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Build the template file for lammps #####################################################
                                                                                         !
        CHEXT1 = TRIM(CHNAME_OUTPUT_FILE_LIBRARY(i))//'.template';                       !
                                                                                         !
        call WRITE_LAMMPS_TEMPLATE(99,CHEXT1);                                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Build the parameter file for lammps ####################################################
                                                                                         !
        CHEXT2 = TRIM(CHNAME_OUTPUT_FILE_LIBRARY(i))//'.parameter';                      !
                                                                                         !
        call WRITE_INTERATOMIC_POTENTIALS_TEMPLATE(99,CHEXT2);                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Build the info file for lammps #########################################################
                                                                                         !
        CHEXT3 = TRIM(CHNAME_OUTPUT_FILE_LIBRARY(i))//'.info';                           !
                                                                                         !
        call WRITE_LAMMPS_INFO(99,CHEXT3);                                               !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Move files to the directory where they will be stored in the library ###################
                                                                                         !
        call MOVE_TO_LIBRARY_DIRECTORY(99,CHEXT1,CHEXT2,CHEXT3);                         ! 
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Deallocate arrays ##########################################################################
                                                                                         !
    deallocate(CHNAME_FILE_LIBRARY);                                                     !
                                                                                         !
    deallocate(CHNAME_LAMMPS_INPUT_LIBRARY);                                             !
                                                                                         !
    deallocate(CHNAME_OUTPUT_FILE_LIBRARY);                                              !
                                                                                         !
!   ### Closing the program ########################################################################
                                                                                         !
    write(99,'(a70)') '+'//REPEAT('-',68)//'+';                                          !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(99,'(a70)') '|'//REPEAT(' ',26)//' End of program '//REPEAT(' ',26)//'|';      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(99,'(a70)') '+'//REPEAT('-',68)//'+';                                          !
                                                                                         !
    close(99);                                                                           !
                                                                                         !
end program BUILDLIBRARY

subroutine sgrnd(seed)

      implicit integer(a-z)

!* Period parameters
      parameter(N     =  624)

      dimension mt(0:N-1)
!                     the array for the state vector
      common /block/mti,mt
      save   /block/
!*
!*      setting initial seeds to mt[N] using
!*      the generator Line 25 of Table 1 in
!*      [KNUTH 1981, The Art of Computer Programming
!*         Vol. 2 (2nd Ed.), pp102]
!*
      mt(0)= iand(seed,-1)
      do 1000 mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
 1000 continue

      return

end subroutine sgrnd

!************************************************************************

double precision function grnd()

      implicit integer(a-z)

!* Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
                                 !   constant vector a
      parameter(UMASK = -2147483648)
                                  !  most significant w-r bits
      parameter(LMASK =  2147483647)
                                   ! least significant r bits
! Tempering parameters

      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)

      dimension mt(0:N-1)
!*                     the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/
!*                     mti==N+1 means mt[N] is not initialized
!*
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
!*                        mag01(x) = x * MATA for x=0,1
!*
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
!*
      if(mti.ge.N) then
!*                       generate N words at one time
        if(mti.eq.N+1) then
!*                            if sgrnd() has not been called,
          call sgrnd(4357)
!*                              a default initial seed is used
        endif
!*
        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
!*
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
!*
      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      else
        grnd=dble(y)/(2.0d0**32-1.0d0)
      endif
!*
      return

end function grnd
