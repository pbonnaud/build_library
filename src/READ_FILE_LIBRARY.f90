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

subroutine READ_XYZ_FILE_LIBRARY(icanal,CHFILE_FORMAT,           &
                                        NFILE_LIBRARY,           &
                                        CHNAME_FILE_LIBRARY,     &
                                        NATOM_FILE_LIBRARY,      &
                                        ATOM_LABEL_LIBRARY,      &
                                        CONFIG_NAT_LIBRARY,      &
                                        CONFIG_RI_LIBRARY,       &
                                        NTYPE_ATOM_LIBRARY,      &
                                        ATOM_MASSES_LIBRARY,     &
                                        CONFIG_ATOMID_LIBRARY,   &
                                        CONFIG_ATOM_TYPE_LIBRARY)
 


!   ************************************************************************************************
!   **                                  CHECK FILE LIBRARY                                        **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                   : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                          **
!   ** CHFILE_FORMAT            : FORMAT OF FILES THAT ARE INTENDED TO BE READ                    **
!   ** NFILE_LIBRARY            : NUMBER OF FILES COMING FROM THE LIBRARY                         **
!   ** CHNAME_FILE_LIBRARY      : FILE NAME CONTAINING THE MOLECULAR CONFIGURATION OF THE         **
!   **                            MOLECULE                                                        **
!   ** NATOM_FILE_LIBRARY       :
!   ** ATOM_LABEL_LIBRARY       :
!   ** CONFIG_NAT_LIBRARY       :
!   ** CONFIG_RI_LIBRARY        :
!   ** NTYPE_ATOM_LIBRARY
!   ** ATOM_MASSES_LIBRARY
!   ** CONFIG_ATOMID_LIBRARY    : ATOM IDs IN THE CONFIGURATION OF THE MOLECULE IN THE LIBRARY    **
!   ** CONFIG_ATOM_TYPE_LIBRARY : 
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=20), intent(in) :: CHFILE_FORMAT;

    integer (kind=4), intent(in) :: NFILE_LIBRARY;

    character (len=150), dimension(1:100), intent(in) :: CHNAME_FILE_LIBRARY;

!   ************************************************************************************************

    integer (kind=4), dimension(1:100), intent(out) :: NATOM_FILE_LIBRARY, NTYPE_ATOM_LIBRARY;

    integer (kind=4), dimension(1:1000,1:100), intent(out) :: CONFIG_ATOMID_LIBRARY, CONFIG_ATOM_TYPE_LIBRARY;

    real (kind=8), dimension(1:3,1:1000,1:100), intent(out) :: CONFIG_RI_LIBRARY;

    real (kind=8), dimension(1:20,1:100), intent(out) :: ATOM_MASSES_LIBRARY; 

    character (len=20), dimension(1:20,1:100) :: ATOM_LABEL_LIBRARY;

    character (len=20), dimension(1:1000,1:100), intent(out) :: CONFIG_NAT_LIBRARY;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: IFOUND, ILENGTH;

    character (len=3) :: CHLENGTH;

!   ************************************************************************************************

    integer (kind=4) :: IOSEF1;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************

    CHTITLE = 'READ XYZ FILE LIBRARY';

    ILENGTH_TITLE = LEN_TRIM(CHTITLE);

    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));

    NTYPE_ATOM_LIBRARY(1:100) = 0;

!   write(icanal,*) 'OK 1';
    NATOM_FILE_LIBRARY(1:100) = 0;

!   write(icanal,*) 'OK 1';
    ATOM_MASSES_LIBRARY(1:20,1:100) = 0.0d0;

!   write(icanal,*) 'OK 1';
    ATOM_LABEL_LIBRARY(1:20,1:100) = 'XXX';    

!   write(icanal,*) 'OK 1';
    CONFIG_RI_LIBRARY(1:3,1:1000,1:100) = 0.0d0;

!   write(icanal,*) 'OK 1';
    CONFIG_NAT_LIBRARY(1:1000,1:100) = 'XXX';

    CONFIG_ATOMID_LIBRARY(1:1000,1:100)    = 0;
    CONFIG_ATOM_TYPE_LIBRARY(1:1000,1:100) = 0;

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    do i = 1, NFILE_LIBRARY

        ILENGTH = LEN_TRIM(CHNAME_FILE_LIBRARY(i));

        IOSEF1 = 70 - 3 - ILENGTH;

        write(icanal,'(a70)') '| '//TRIM(CHNAME_FILE_LIBRARY(i))//REPEAT(' ',IOSEF1)//'|';

        open(12,file=TRIM(CHNAME_FILE_LIBRARY(i)),status='old');
        read(12,*) NATOM_FILE_LIBRARY(i);
        read(12,*);

        do j = 1, NATOM_FILE_LIBRARY(i);
            read(12,*) CONFIG_NAT_LIBRARY(j,i), CONFIG_RI_LIBRARY(1:3,j,i);

            IFOUND = 0;

            if ( NTYPE_ATOM_LIBRARY(i) > 0 ) then
                do k = 1, NTYPE_ATOM_LIBRARY(i);
                    if ( TRIM(CONFIG_NAT_LIBRARY(j,i)) == TRIM(ATOM_LABEL_LIBRARY(k,i)) ) then
                        IFOUND = 1;
                        EXIT;
                    end if
                end do
            end if

            if ( IFOUND == 1 ) CYCLE;

            NTYPE_ATOM_LIBRARY(i) = NTYPE_ATOM_LIBRARY(i) + 1;

            ATOM_LABEL_LIBRARY(NTYPE_ATOM_LIBRARY(i),i) = CONFIG_NAT_LIBRARY(j,i);

        end do

        close(12);

        do j = 1, NATOM_FILE_LIBRARY(i); ! LOOP OVER THE NUMBER OF ATOMS IN THE MOLECULE OF THE LIBRARY
            CONFIG_ATOMID_LIBRARY(j,i) = j;
  
            do k = 1, NTYPE_ATOM_LIBRARY(i);
                if ( TRIM(CONFIG_NAT_LIBRARY(j,i)) == TRIM(ATOM_LABEL_LIBRARY(k,i)) ) then
                    CONFIG_ATOM_TYPE_LIBRARY(j,i) = k;
                    EXIT;
                end if
            end do
        end do

        do j = 1, NTYPE_ATOM_LIBRARY(i);                    !LOOP OVER THE NYBER OF ATOM TYPES IN THE MOLECULE OF THE LIBRARY 
            if ( TRIM(ATOM_LABEL_LIBRARY(j,i)) == 'O' ) ATOM_MASSES_LIBRARY(j,i) = 15.99940d0;
            if ( TRIM(ATOM_LABEL_LIBRARY(j,i)) == 'C' ) ATOM_MASSES_LIBRARY(j,i) = 12.01074d0;
            if ( TRIM(ATOM_LABEL_LIBRARY(j,i)) == 'N' ) ATOM_MASSES_LIBRARY(j,i) = 14.00670d0;
            if ( TRIM(ATOM_LABEL_LIBRARY(j,i)) == 'H' ) ATOM_MASSES_LIBRARY(j,i) =  1.00794d0;
            if ( TRIM(ATOM_LABEL_LIBRARY(j,i)) == 'S' ) ATOM_MASSES_LIBRARY(j,i) = 32.06500d0;
        end do

        do j = 1, NTYPE_ATOM_LIBRARY(i);
            write(icanal,*) '| ', ATOM_LABEL_LIBRARY(j,i), ATOM_MASSES_LIBRARY(j,i);
        end do 

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    end do

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine READ_XYZ_FILE_LIBRARY











