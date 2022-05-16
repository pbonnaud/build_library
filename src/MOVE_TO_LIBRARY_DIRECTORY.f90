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

subroutine MOVE_TO_LIBRARY_DIRECTORY(icanal,CHEXT1,CHEXT2,CHEXT3) 

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal  : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                         **
!   ** CHEXT  : NAME OF THE LAMMPS CONFIGURATION                                                  **
!   **                                                                                            **
!   ************************************************************************************************

    use module_physical_constants;

    use module_size_arrays;

    use module_library;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=250), intent(in) :: CHEXT1, CHEXT2, CHEXT3;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: ilocal_move;

    character (len=250) :: CHARLINE;

    integer (kind=4) :: EOF;

    logical :: PROBE1;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Move files to the library sub-directory';                                 !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the name of the sub-directory ########################################################
                                                                                         !
    IOSEF1 = 70 - 3 - LEN_TRIM(CHNAME_LIBRARY_SUBDIRECTORY);                             !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHNAME_LIBRARY_SUBDIRECTORY)// &                    !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of the flag that tells the program to move generated files ##################
                                                                                         !
    ilocal_move = 0;                                                                     !
                                                                                         !
!   ### Check the presence of the sub-directory ####################################################
                                                                                         !
    if ( TRIM(CHNAME_LIBRARY_SUBDIRECTORY) == '.' ) then;                                !
                                                                                         !
        write(icanal,*) '| Output files will remain in the current directory';           !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    else                                                                                 ! 
                                                                                         !
        ilocal_move = 1;                                                                 !
                                                                                         !
!       ### Write a small script that check the presence of the designated directory ###############
                                                                                         !
        open(22,file='Tempo');                                                           !
                                                                                         !
        write(22,'(a11)') '#!/bin/bash';                                                 !
                                                                                         !
        write(22,*);                                                                     !
                                                                                         !
!       CHOSEF1 = 'PATH_DIR="'//TRIM(CHNAME_LIBRARY_SUBDIRECTORY)//'"';                  !
                                                                                         !
        write(22,*) 'PATH_DIR="'//TRIM(CHNAME_LIBRARY_SUBDIRECTORY)//'"';                !
                                                                                         !
        write(22,*);                                                                     !
                                                                                         !
                                                                                         !
!       write(22,*) 'if [ -d "'//                       &                                !
!                   TRIM(CHNAME_LIBRARY_SUBDIRECTORY)// &                                !
!                   '" ]; then touch EXIST; fi';                                         !

        write(22,*) 'if [ -d $PATH_DIR ]; then touch EXIST; fi';                         !
                                                                                         !
        close(22);                                                                       !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Set executable properties to the script an run it ######################################
                                                                                         !
        call system('chmod +x Tempo');                                                   !
                                                                                         !
        call system('./Tempo');                                                          !
                                                                                         !
!       ### Check for an output generated by the script ############################################
                                                                                         !
        inquire(FILE='EXIST',EXIST=PROBE1);                                              !
                                                                                         !
        write(icanal,'(a11,L1,a58)') '| PROBE1 : ', PROBE1, REPEAT(' ',57)//'|';         !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Remove the script and the generated output (cleaning) ##################################
                                                                                         !
        call system('rm -f EXIST Tempo');                                                !
                                                                                         !
!       ### Create a directory if not found at the given path ###################################### 
                                                                                         !
        if ( PROBE1 .EQV. .FALSE. ) then;                                                !
                                                                                         !
            call system('mkdir ../../'//TRIM(CHNAME_LIBRARY_SUBDIRECTORY));              ! 
                                                                                         !
            IOSEF1 = 70 - 33 - LEN_TRIM(CHNAME_LIBRARY_SUBDIRECTORY);                    !
                                                                                         !
            write(icanal,'(a70)') '| The sub-directory '//                 &             !
                                  TRIM(CHNAME_LIBRARY_SUBDIRECTORY)//      &             !
                                  ' was created'//REPEAT(' ',IOSEF1)//'|';               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else                                                                             !
                                                                                         !
            IOSEF1 = 70 - 51 - LEN_TRIM(CHNAME_LIBRARY_SUBDIRECTORY);                    !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| The sub-directory '//             &                 !
                                  TRIM(CHNAME_LIBRARY_SUBDIRECTORY)//  &                 !
                                  ' already exists in the library'//   &                 !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write the name of the templated files to be moved to the sub-directory #####################
                                                                                         !
    IOSEF1 = 67 - LEN_TRIM(CHEXT1);                                                      !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT1)//REPEAT(' ',IOSEF1)//'|';                   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 67 - LEN_TRIM(CHEXT2);                                                      !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT2)//REPEAT(' ',IOSEF1)//'|';                   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 67 - LEN_TRIM(CHEXT3);                                                      !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT3)//REPEAT(' ',IOSEF1)//'|';                   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Move generated files to the sub-directory ##################################################
                                                                                         !
    if ( ilocal_move == 1 ) then;                                                        !
                                                                                         !
!       call system ('mv '//TRIM(CHEXT1)//' ../../'//         &                              !
!                    TRIM(CHNAME_LIBRARY_SUBDIRECTORY)//'/');                                !
                                                                                         !
!       call system ('mv '//TRIM(CHEXT2)//' ../../'//         &                              !
!                    TRIM(CHNAME_LIBRARY_SUBDIRECTORY)//'/');                                !
                                                                                         !
!       call system ('mv '//TRIM(CHEXT3)//' ../../'//         &                              !
!                    TRIM(CHNAME_LIBRARY_SUBDIRECTORY)//'/');                                !
                                                                                         !
        call system ('mv '//TRIM(CHEXT1)//' '//         &                                !
                     TRIM(CHNAME_LIBRARY_SUBDIRECTORY)//'/');                            !
                                                                                         !
        call system ('mv '//TRIM(CHEXT2)//' '//         &                                !
                     TRIM(CHNAME_LIBRARY_SUBDIRECTORY)//'/');                            !
                                                                                         !
        call system ('mv '//TRIM(CHEXT3)//' '//         &                                !
                     TRIM(CHNAME_LIBRARY_SUBDIRECTORY)//'/');                            !
                                                                                         !
        write(icanal,'(a70)') '| Files were moved to the sub-directory'// &              !
                              REPEAT(' ',30)//'|';                                       !  
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine MOVE_TO_LIBRARY_DIRECTORY
