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

subroutine MAKE_ROUTINE_TITLE(icanal,ILENGTH_MAX,ILENGTH_TITLE,CHTITLE)

!   ************************************************************************************************
!   **                               BUILD TITLES FOR SUBROUTINES                                 **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal        : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                   **
!   ** ILENGTH_MAX   : NUMBER OF COLUMNS ON WHICH THE FILE IS WRITTEN                             ** 
!   ** ILENGTH_TITLE : LENGTH OF THE STRING CHTITLE                                               **
!   ** CHTITLE       : TITLE OF THE SUBROUTINE                                                    ** 
!   **                                                                                            **
!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal, ILENGTH_MAX, ILENGTH_TITLE;

    character (len=ILENGTH_TITLE) :: CHTITLE;

!   ************************************************************************************************

    integer (kind=4) :: IREPEAT1, IREPEAT2;

    character (len=10) :: CHLENGTH_MAX;

!   ************************************************************************************************

!   write(icanal,*) 'GGGGG';

    if ( ILENGTH_MAX < 1000 ) write(CHLENGTH_MAX,'(i3)') ILENGTH_MAX;
    if ( ILENGTH_MAX <  100 ) write(CHLENGTH_MAX,'(i2)') ILENGTH_MAX;
    if ( ILENGTH_MAX <   10 ) write(CHLENGTH_MAX,'(i1)') ILENGTH_MAX;

    IREPEAT1 = ILENGTH_MAX - 4 - ILENGTH_TITLE;
    IREPEAT2 = ILENGTH_MAX - 2;

    write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') '+ '//TRIM(CHTITLE)//' '//REPEAT('-',IREPEAT1)//'+';
    write(icanal,'(a'//TRIM(CHLENGTH_MAX)//')') '|'//REPEAT(' ',IREPEAT2)//'|';


!   write(icanal,*) 'LLLLL';

end subroutine MAKE_ROUTINE_TITLE
