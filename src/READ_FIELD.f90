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

subroutine READ_FIELD(icanal, &
                      PATH_DIRECTORY, &
                      CH_WORKING_FILE, &
                      IFIELD)

!   ************************************************************************************************
!   **                                    READ FIELD FILE                                         **
!   ************************************************************************************************
!   **                                                                                            **
!   ** PATH_DIRECTORY  : PATH TO THE DIRECTORY WHERE THE FILE IS LOCATED                          **
!   ** CH_WORKING_FILE : ROOT NAME OF THE WORKING FILE                                            **
!   **                                                                                            **
!   ** IFIELD          :  PRESENCE OF THE FILE IN THE WORKING DIRECTORY                           **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=200), intent(in) :: PATH_DIRECTORY;

    character (len=5), intent(in) :: CH_WORKING_FILE;

    integer (kind=4), intent(out) :: IFIELD;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: ilabel;

    integer (kind=4) :: SAME_AT, EOF, ICOUNT_FINISH;

    real (kind=8) :: CHARGE;

    character (len=3) :: ATO_CHLABEL;

    character (len=3), dimension(1:10) :: ATOM_LABEL;

    character (len=100) :: CHREAD1;

    integer (kind=4) :: IOSEF1;

    character (len=200) :: CH_DATA_FILE_NAME;

    logical :: PROBE1;

!   ************************************************************************************************

    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    ATOM_LABEL(1:10) = 'XXX';

    IFIELD = 0;

    CH_DATA_FILE_NAME = TRIM(PATH_DIRECTORY)//TRIM(CH_WORKING_FILE)//'_FIELD.dat';

    write(icanal,*) CH_DATA_FILE_NAME;

    inquire(FILE=TRIM(CH_DATA_FILE_NAME),EXIST=PROBE1);

!   if (PROBE1 == .FALSE. ) then
    if (PROBE1 .EQV. .FALSE. ) then;
        write(icanal,*) '| '//TRIM(CH_DATA_FILE_NAME)//' DOESNT EXIST !!!';
!       write(99,*) 'END OF PROGRAM';
!       close(99);
!   end if
    else
        IFIELD = 1;

        ICOUNT_FINISH = 0;

        open(14,file=trim(CH_DATA_FILE_NAME),status='old');

        do        ! READ CHARGES
            read(14,'(a)') CHREAD1;

            if ( index(CHREAD1,'FINISH') > 0 ) then
                ICOUNT_FINISH = ICOUNT_FINISH + 1;
                if (ICOUNT_FINISH == 1 ) EXIT;
            else
                if ( ICOUNT_FINISH == 0 ) then
                    read(CHREAD1,*) ATO_CHLABEL, CHARGE;
                    SAME_AT = 0;
                    do ilabel = 1, NLABEL
                       if ( ATOM_LABEL(ilabel) == ATO_CHLABEL ) then
                           SAME_AT = 1;
                           ATOM_CHAR(ilabel) = CHARGE;
                       end if
                    end do
                    if ( SAME_AT == 0 ) then
                       NLABEL = NLABEL + 1;
                       ATOM_LABEL(NLABEL) = ATO_CHLABEL;
                       ATOM_CHAR(NLABEL) = CHARGE;
                    end if
                end if
            end if
        end do

        do      ! READ CENTERS OF MASS
            EOF = 0;
            read(14,'(a)',iostat=EOF) CHREAD1;

            if ( EOF == 0 ) then
                if ( INDEX(CHREAD1,'FINISH') > 0 ) then
                    EXIT;
                else
                    ICOMA = ICOMA + 1;
                    read(CHREAD1,*) COMA_COORD(1:3,ICOMA);
                end if
            else
                EXIT;
            end if
        end do

        close(14);

!       ### WRITE ATOMIC LABELS ########################################################################
        write(icanal,*);
        write(icanal,*) 'LABELS (FIELD) : ';
        write(icanal,*);

        do i = 1, NLABEL
            write(icanal,*) i, ATOM_LABEL(i);
        end do

    end if

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine READ_FIELD
