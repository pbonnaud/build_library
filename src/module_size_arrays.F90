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

module module_size_arrays

    implicit none;

!   ************************************************************************************************

    integer (kind=4), parameter :: NN1 = 100;     ! MAXIMUM NUMBER OF LIBRARY FILES THAT CAN BE READ
                                                  !
    integer (kind=4), parameter :: NN2 = 2000;    ! MAXIMUM NUMBER OF PARAMETERS READ IN THE LIBRARY FILES 

    integer (kind=4), parameter :: NN3 = 100000; ! MAXIMUM NUMBER IN CONFIGURATION ARRAYS

end module module_size_arrays
