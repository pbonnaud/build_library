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

subroutine READ_INPUT(icanal)

!   ************************************************************************************************
!   **                                  READ THE INPUT FILE                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : Canal on which output data are written                                            **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_size_arrays;

    use module_library;

    use module_physical_constants;

    use module_frc_arrays;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

    integer (kind=4) :: EOF;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

!   ************************************************************************************************

    logical :: PROBE1;

    character (len=500) :: CHERROR_MESSAGE;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Read the input file';                                                     !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Check the presence of the input file for the current program ###############################
                                                                                         !
    inquire(FILE='INPUT_BUILDLIBRARY',EXIST=PROBE1);                                     !
                                                                                         !
    if ( PROBE1 .EQV. .FALSE. ) then;                                                    !
                                                                                         !
        CHERROR_MESSAGE = 'No input file - Please check your directory';                 !
                                                                                         !
        call WRITE_ERROR_MESSAGE(icanal,70,CHERROR_MESSAGE);                             !
                                                                                         !
    else                                                                                 !
                                                                                         !
        write(icanal,'(a70)') '| The input file has been found in the directory'// &     !
                              REPEAT(' ',21)//'|';                                       !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read the input file of the current program #################################################
                                                                                         !
    open(1,file='INPUT_BUILDLIBRARY',status='old');                                      !
                                                                                         !
    do i = 1, 6;                                                                         !
                                                                                         !
        read(1,*);                                                                       !
                                                                                         !
    end do                                                                               !
                                                                                         !
    read(1,*) NFILE_LIBRARY, CHFILE_FORMAT;                                              !
                                                                                         !
    write(icanal,'(a42,i8,a20)') '| Number of files to read               : ', &         !
                                 NFILE_LIBRARY,                                &         !
                                 REPEAT(' ',19)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 70 - 42 - 1 - LEN_TRIM(CHFILE_FORMAT);                                      !
                                                                                         !
    write(icanal,'(a70)') '| File format of the atomic coordinates : '// &               !
                          TRIM(CHFILE_FORMAT)//                          &               !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays containing name of input and output files ##################################
                                                                                         !
    allocate(CHNAME_FILE_LIBRARY(1:NFILE_LIBRARY));                                      !
                                                                                         !
    allocate(CHNAME_LAMMPS_INPUT_LIBRARY(1:NFILE_LIBRARY));                              !
                                                                                         !
    allocate(CHNAME_OUTPUT_FILE_LIBRARY(1:NFILE_LIBRARY));                               !
                                                                                         !
!   ### Read names of files to read and the name of output files ###################################
                                                                                         !
    do i = 1, NFILE_LIBRARY;                                                             ! 
                                                                                         !
        read(1,*) CHNAME_FILE_LIBRARY(i),         &                                      !
                  CHNAME_LAMMPS_INPUT_LIBRARY(i), &                                      !
                  CHNAME_OUTPUT_FILE_LIBRARY(i);                                         !
                                                                                         !
    end do                                                                               !
                                                                                         !
    do i = 1, NFILE_LIBRARY;                                                             !
                                                                                         !
        IOSEF1 = 70 - 2 - 1 - LEN_TRIM(CHNAME_FILE_LIBRARY(i));                          !
                                                                                         !
        write(icanal,'(a70)') '| '//TRIM(CHNAME_FILE_LIBRARY(i))// &                     !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        IOSEF1 = 70 - 2 - 1 - LEN_TRIM(CHNAME_LAMMPS_INPUT_LIBRARY(i));                  !
                                                                                         !
        write(icanal,'(a70)') '| '//TRIM(CHNAME_LAMMPS_INPUT_LIBRARY(i))// &             !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        IOSEF1 = 70 - 2 - 1 - LEN_TRIM(CHNAME_OUTPUT_FILE_LIBRARY(i));                   !
                                                                                         !
        write(icanal,'(a70)') '| '//TRIM(CHNAME_OUTPUT_FILE_LIBRARY(i))// &              !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read the flag and the value of the net total charge (adjust) ###############################
                                                                                         !
    read(1,*) iadjust_charges, TOTAL_NET_CHARGE;                                         !
                                                                                         !
    if ( iadjust_charges == 1 ) then;                                                    !
                                                                                         !
        write(icanal,'(a32,f12.5,a26)') '| The (new) total net charge is ', &            !
                                        TOTAL_NET_CHARGE,                   &            !
                                        ' [e]'//REPEAT(' ',21)//'|';                     !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  ! 
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!!!!!!!   ### Name of the library where where the output files will be stored ############################
!   ### Path to the directory (library) where output files will be stored ##########################
                                                                                         !
    read(1,*) CHNAME_LIBRARY_SUBDIRECTORY;                                               !
                                                                                         !
    IOSEF1 = 70 - 2 - 1 - LEN_TRIM(CHNAME_LIBRARY_SUBDIRECTORY);                         !
                                                                                         !
    write(icanal,'(a70)') '| Path to the directory where lammps files will be stored'// &!
                          REPEAT(' ',12)//'|';                                           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHNAME_LIBRARY_SUBDIRECTORY)// &                    !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Close the input file #######################################################################
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   ### Initialization of frc file variables #######################################################
                                                                                         !
!   ibuild_frc = 0;                                                                      !
                                                                                         !
!   FRC_FILE_NAME = 'XXX';                                                               !
                                                                                         !
!   ### Seek the command for building a frc file ###################################################
                                                                                         !
!   open(1,file='INPUT_BUILDLIBRARY',status='old');                                      !
                                                                                         !
!   EOF = 0;                                                                             !
                                                                                         !
!   do                                                                                   !
                                                                                         !
!       read(1,'(a)',iostat=EOF) CHARLINE;                                               !
                                                                                         !
!       if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
!       if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                            !
                                                                                         !
!       if ( INDEX(CHARLINE,'build_frc ') > 0 ) CYCLE;                                   ! If there is still LAMMPS variables in the read line, then 
                                                                                         !
!       read(CHARLINE,*) CHOSEF1, FRC_FILE_NAME;                                         !
                                                                                         !
!   end do                                                                               !
                                                                                         !
!   close(1);                                                                            !
                                                                                         !
!   ### Write an output message if a build frc command was found ###################################
                                                                                         !
!SSSSSSSSSS
!SSSSSSSSSS
!SSSSSSSSSS
!SSSSSSSSSS
!SSSSSSSSSS

!   ### Wrap-up the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| The input file has been read '//REPEAT(' ',38)//'|';        !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine READ_INPUT

