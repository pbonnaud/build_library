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

subroutine MANAERROR(EOF,NERR,ILINE,JLINE)

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: EOF, NERR, ILINE, JLINE;

!   ************************************************************************************************

    if ( EOF /= 0 ) then
        if ( NERR == 1 ) then
            write(99,*) '!!! ERR 1                                !!!';
            write(99,*) '!!! EOF ', EOF;
            write(99,*) '!!! PROBLEM IN READING THE PARAMETERS    !!!';
            write(99,*) 'PLEASE CHECK YOUR INPUT FILE';
            write(99,*) ' nisp : ', ILINE;
            close(99);
            stop;
        end if
    end if

end subroutine MANAERROR
