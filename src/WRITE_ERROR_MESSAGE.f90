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

subroutine WRITE_ERROR_MESSAGE(icanal,ILENGTH_MAX,CHERROR_MESSAGE)

!   ************************************************************************************************
!   **                               BUILD TITLES FOR SUBROUTINES                                 **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal        : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                   **
!   ** ILENGTH_MAX   : NUMBER OF COLUMNS ON WHICH THE FILE IS WRITTEN                             ** 
!   ** CHERROR_MESSAGE
!   **                                                                                            **
!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal, ILENGTH_MAX;

    character (len=500) :: CHERROR_MESSAGE;

!   ************************************************************************************************

    integer (kind=4) :: IREPEAT1, IREPEAT2;

    character (len=10) :: CHLENGTH_MAX;

!   ************************************************************************************************

    if ( ILENGTH_MAX < 1000 ) write(CHLENGTH_MAX,'(i3)') ILENGTH_MAX;
    if ( ILENGTH_MAX <  100 ) write(CHLENGTH_MAX,'(i2)') ILENGTH_MAX;
    if ( ILENGTH_MAX <   10 ) write(CHLENGTH_MAX,'(i1)') ILENGTH_MAX;

    IREPEAT1 = ILENGTH_MAX - 8 - LEN_TRIM(CHERROR_MESSAGE);

!   write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') '+ '//TRIM(CHTITLE)//' '//REPEAT('-',IREPEAT1)//'+';
!   write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') '|'//REPEAT(' ',IREPEAT2)//'|';


    write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') REPEAT('!',ILENGTH_MAX);
    write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') '!!! '//TRIM(CHERROR_MESSAGE)//REPEAT(' ',IREPEAT1)//' !!!';
    write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') REPEAT('!',ILENGTH_MAX);
    write(icanal,*);
    write(icanal,'(a22)') '!!! END OF PROGRAM !!!';
    close(icanal);
    stop;

end subroutine WRITE_ERROR_MESSAGE
