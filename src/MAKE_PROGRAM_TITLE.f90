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

subroutine MAKE_PROGRAM_TITLE(icanal,ILENGTH_MAX,ILENGTH_TITLE,CHTITLE,CHPRGM_VERSION,CHDATE)

!   ************************************************************************************************
!   **                               BUILD TITLE FOR THEM PROGRAM                                 **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal         : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                  **
!   ** ILENGTH_MAX    : NUMBER OF COLUMNS ON WHICH THE FILE IS WRITTEN                            ** 
!   ** ILENGTH_TITLE  : LENGTH OF THE STRING CHTITLE                                              **
!   ** CHTITLE        : TITLE OF THE SUBROUTINE                                                   ** 
!   ** CHPRGM_VERSION : VERSION OF THE PROGRAM                                                    **
!   ** CHDATE         : DATE WHEN THE PROGRAM WAS LAST MODIFIED                                   **
!   **                                                                                            **
!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal, ILENGTH_MAX, ILENGTH_TITLE;

    character (len=ILENGTH_TITLE) :: CHTITLE;

    character (len=250), intent(in) :: CHPRGM_VERSION, CHDATE;

!   ************************************************************************************************

    integer (kind=4) :: IREPEAT1, IREPEAT2, IREPEAT3;

    integer (kind=4) :: IOSEF1, IOSEF2;

    character (len=10) :: CHLENGTH_MAX;

!   ************************************************************************************************

    if ( ILENGTH_MAX < 1000 ) write(CHLENGTH_MAX,'(i3)') ILENGTH_MAX;
    if ( ILENGTH_MAX <  100 ) write(CHLENGTH_MAX,'(i2)') ILENGTH_MAX;
    if ( ILENGTH_MAX <   10 ) write(CHLENGTH_MAX,'(i1)') ILENGTH_MAX;

    IREPEAT1 = ILENGTH_MAX - 4 - ILENGTH_TITLE;                         !  NUMBER OF SPACE ON THE LINE WHERE THE TITLE WIL BE WRITTEN
    IREPEAT2 = INT(REAL(IREPEAT1)*0.5d0);

    IOSEF1 = 4 + ILENGTH_TITLE + 2 * IREPEAT2;

    if ( IOSEF1 < ILENGTH_MAX ) then
        IREPEAT3 = IREPEAT2 + ILENGTH_MAX - IOSEF1;
    else
        IREPEAT3 = IREPEAT2;
    end if

    write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') REPEAT('*',ILENGTH_MAX);
    write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') '**'//REPEAT(' ',IREPEAT2)//TRIM(CHTITLE)//REPEAT(' ',IREPEAT3)//'**';
    write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') REPEAT('*',ILENGTH_MAX);
    write(icanal,*);

    IOSEF1 = LEN_TRIM(CHPRGM_VERSION);
    IOSEF2 = LEN_TRIM(CHDATE);

    IREPEAT1 = ILENGTH_MAX - IOSEF1 - IOSEF2;

    write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') TRIM(CHPRGM_VERSION)//REPEAT(' ',IREPEAT1)//TRIM(CHDATE); 
    write(icanal,*);


end subroutine MAKE_PROGRAM_TITLE
