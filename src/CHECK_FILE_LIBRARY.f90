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

subroutine CHECK_FILE_LIBRARY(icanal,CHFILE_FORMAT,         &
                                     NFILE_LIBRARY,         &
                                     CHNAME_FILE_LIBRARY)

!   ************************************************************************************************
!   **                                  CHECK FILE LIBRARY                                        **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   ** CHFILE_FORMAT         : FORMAT OF FILES THAT ARE INTENDED TO BE READ                       **
!   ** NFILE_LIBRARY         : NUMBER OF FILES COMING FROM THE LIBRARY                            **
!   ** CHNAME_FILE_LIBRARY   : FILE NAME CONTAINING THE MOLECULAR CONFIGURATION OF THE MOLECULE   **


!   ** iprimitive_cell       : THE PRIMITIVE CELL HAS BEEN DEFINED IN THE INPUT FILE              **
!   ** PRIMITIVE_CELL_AXIS   : DIMENSIONS OF THE PRIMITIVE CELL                                   **
!   ** PRIMITIVE_CELL_ANGDEG : ANGLES OF THE PRIMITIVE CELL                                       **
!   ** PRIMITIVE_NATOM       :                                                                    **
!   ** PRIMITIVE_NAT         :                                                                    **
!   ** PRIMITIVE_RI          :                                                                    **
!   ** NBRE_REPEAT_PATTERN   :                                                                    **

!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=20), intent(in) :: CHFILE_FORMAT;

    integer (kind=4), intent(in) :: NFILE_LIBRARY;

    character (len=150), dimension(1:100), intent(inout) :: CHNAME_FILE_LIBRARY;

!   ************************************************************************************************

!   ************************************************************************************************

    integer (kind=4) :: i;

    integer (kind=4), dimension(:), allocatable :: NATOM_FILE_LIBRARY;

    character (len=250) :: CHNAME_LIBRARY_DIRECTORY, CHNAME_TMP1, CHNAME_TMP2;

    logical :: PROBE1, PROBE2, PROBE3, PROBE4;

!   ************************************************************************************************

    integer (kind=4) :: IOSEF1, IOSEF2;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !   
    CHTITLE = 'CHECK FILE LIBRARY';                                                      !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   ### 
                                                                                         ! 
    CHNAME_LIBRARY_DIRECTORY = 'LIB-MOLECULAR-MODELS';                                   !
                                                                                         !
    open(22,file='Tempo');                                                               !
    write(22,'(a11)') '#!/bin/bash';                                                     !
    write(22,*);                                                                         !
    write(22,*) 'if [ -d "'//                    &                                       !
                TRIM(CHNAME_LIBRARY_DIRECTORY)// &                                       !
                '" ]; then touch EXIST; fi';                                             !
                                                                                         !
    close(22);                                                                           !
                                                                                         !
    call system('chmod +x Tempo');                                                       !
    call system('./Tempo');                                                              !
                                                                                         !
    inquire(FILE='EXIST',EXIST=PROBE1);                                                  !         
                                                                                         !
    write(icanal,'(a11,L1,a58)') '| PROBE1 : ', PROBE1, REPEAT(' ',57)//'|';             !
                                                                                         !
    call system('rm -f EXIST Tempo');                                                    !
                                                                                         !
    if ( PROBE1 .EQV. .FALSE. ) then;                                                    ! IF THE LIBRARY DIRECTORY HAS NOT BEEN FOUND IN THE CURRENT DIRECTORY, THEN
                                                                                         !
        do i = 1, 5                                                                      ! LOOP OVER THE PARENT DIRECTORIES
                                                                                         !
            CHNAME_LIBRARY_DIRECTORY = '../'//TRIM(CHNAME_LIBRARY_DIRECTORY);            ! BUILD THE PATH TO THE DIRECTORY
                                                                                         !
            open(22,file='Tempo');                                                       ! CREATE A TEMPORARY SCRIPT TO CHECK FILES
            write(22,*) '#!/bin/bash';                                                   !
            write(22,*);                                                                 !
            write(22,*) 'if [ -d "'//                    &                               !
                        TRIM(CHNAME_LIBRARY_DIRECTORY)// &                               !
                        '" ]; then touch EXIST; fi';                                     !
                                                                                         !
            close(22);                                                                   !
                                                                                         !
            call system('chmod +x Tempo');                                               !
            call system('./Tempo');                                                      !
                                                                                         !
            inquire(FILE='EXIST',EXIST=PROBE2);                                          !
                                                                                         !
            call system('rm -f EXIST Tempo');                                            !
                                                                                         !
            if ( PROBE2 .EQV. .TRUE. ) EXIT;                                             !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    IOSEF1 = LEN_TRIM(CHNAME_LIBRARY_DIRECTORY);                                         !
                                                                                         !
    IOSEF2 = 67 - IOSEF1;                                                                !
                                                                                         !
    if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHNAME_LIBRARY_DIRECTORY)//REPEAT(' ',IOSEF2)//'|'; !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
                                                                                         !
    IOSEF1 = 67 - LEN_TRIM(CHFILE_FORMAT);                                               ! 
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHFILE_FORMAT)//REPEAT(' ',IOSEF1)//'|';            !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    if ( TRIM(CHFILE_FORMAT) == 'xyz' ) then;                                            !
                                                                                         !
        open(22,file='Tempo');                                                           !
        write(22,'(a11)') '#!/bin/bash';                                                 !
        write(22,*);                                                                     !
        write(22,*) 'if [ -d "'//                    &                                   !
                    TRIM(CHNAME_LIBRARY_DIRECTORY)// &                                   !
                    '/XYZ" ]; then touch EXIST; fi';                                     !
        close(22);                                                                       !
                                                                                         !
        call system('chmod +x Tempo');                                                   !
                                                                                         !
        call system('./Tempo');                                                          !
                                                                                         !
        inquire(FILE='EXIST',EXIST=PROBE2);                                              !
                                                                                         !
        call system('rm -f  EXIST Tempo');                                               !
                                                                                         !
        if ( PROBE2 .EQV. .FALSE. ) then;                                                !
            stop;                                                                        !
        else                                                                             !
            CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)//'/XYZ/';          !
        end if                                                                           !
                                                                                         !
    else if ( TRIM(CHFILE_FORMAT) == 'lammps' ) then;                                    !
                                                                                         !
        open(22,file='Tempo');                                                           !
        write(22,'(a11)') '#!/bin/bash';                                                 ! 
        write(22,*);                                                                     !
        write(22,*) 'if [ -d "'//                    &                                   !
                    TRIM(CHNAME_LIBRARY_DIRECTORY)// &                                   !
                    '/LAMMPS" ]; then touch EXIST; fi';                                  !
        close(22);                                                                       ! 
                                                                                         !
        call system('chmod +x Tempo');                                                   ! 
                                                                                         !
        call system('./Tempo');                                                          ! 
                                                                                         !
        inquire(FILE='EXIST',EXIST=PROBE2);                                              !
                                                                                         !
        call system('rm -f  EXIST Tempo');                                               !
                                                                                         !
        if ( PROBE2 .EQV. .FALSE. ) then;                                                !
            stop;                                                                        !
        else                                                                             !
            CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)//'/LAMMPS/';       !
        end if                                                                           !

!       inquire(DIRECTORY=TRIM(CHNAME_LIBRARY_DIRECTORY)//'/LAMMPS',EXIST=PROBE2);

!       write(icanal,*) PROBE2;

!       if ( PROBE2 == .FALSE. ) then
!           stop;
!       else
!           CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)//'/LAMMPS/';
!       end if
    end if

!   write(icanal,*) '| '//TRIM(CHNAME_LIBRARY_DIRECTORY);
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    IOSEF1 = LEN_TRIM(CHNAME_LIBRARY_DIRECTORY);                                         !
                                                                                         ! 
    IOSEF2 = 67 - IOSEF1;                                                                !
                                                                                         !
    if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                        !
                                                                                         ! 
    write(icanal,'(a70)') '| '//TRIM(CHNAME_LIBRARY_DIRECTORY)// &                       !
                              REPEAT(' ',IOSEF2)//'|';                                   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    allocate(NATOM_FILE_LIBRARY(1:NFILE_LIBRARY));                                       !
                                                                                         !
    NATOM_FILE_LIBRARY(1:NFILE_LIBRARY) = 0;                                             !
                                                                                         !
    do i = 1, NFILE_LIBRARY;                                                             !
                                                                                         !
        write(icanal,*) '| '//TRIM(CHNAME_FILE_LIBRARY(i));

        CHNAME_TMP1 = TRIM(CHNAME_LIBRARY_DIRECTORY)//TRIM(CHNAME_FILE_LIBRARY(i));

        if ( TRIM(CHFILE_FORMAT) == 'lammps' ) CHNAME_TMP1 = TRIM(CHNAME_TMP1)//'.template';

        write(icanal,*) TRIM(CHNAME_TMP1);

        inquire(FILE=TRIM(CHNAME_TMP1),EXIST=PROBE3);

!       if ( PROBE3 == .TRUE. ) then;                                                    !
        if ( PROBE3 .EQV. .TRUE. ) then;                                                 !
                                                                                         !
            write(icanal,*) '| THE FILE HAS BEEN FOUND IN THE LIBRARY';

            if ( TRIM(CHFILE_FORMAT) == 'lammps' ) then;                                 !

                CHNAME_TMP2 = TRIM(CHNAME_LIBRARY_DIRECTORY)//TRIM(CHNAME_FILE_LIBRARY(i))//'.parameter';

                inquire(FILE=TRIM(CHNAME_TMP2),EXIST=PROBE4);

                write(icanal,*) TRIM(CHNAME_TMP2);

!               if ( PROBE4 == .FALSE. ) then
                if ( PROBE4 .EQV. .FALSE. ) then
                    write(icanal,*) '| THE FILE HAS NOT BEEN FOUND IN THE LIBRARY';
                    stop;
                else
                    write(icanal,*) '| THE FILE HAS BEEN FOUND IN THE LIBRARY';

                    CHNAME_FILE_LIBRARY(i) = TRIM(CHNAME_LIBRARY_DIRECTORY)//TRIM(CHNAME_FILE_LIBRARY(i));
                end if

            else
                CHNAME_FILE_LIBRARY(i) = TRIM(CHNAME_LIBRARY_DIRECTORY)//TRIM(CHNAME_FILE_LIBRARY(i));

                open(12,file=TRIM(CHNAME_FILE_LIBRARY(i)),status='old');
                read(12,*) NATOM_FILE_LIBRARY(i);
                close(12);

                write(icanal,*) '| NATOM FILE LIBRARY : ', NATOM_FILE_LIBRARY(i);
            end if
        else
            write(icanal,*) '| THE FILE HAS NOT BEEN FOUND IN THE LIBRARY';
            stop;
        end if

    end do
                                                                                         !
    deallocate(NATOM_FILE_LIBRARY);                                                      !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine CHECK_FILE_LIBRARY











