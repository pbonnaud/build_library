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

subroutine SET_TRANSLATION(RIJ_INI,RIJ_FIN,PASSA,PASSB)

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3), intent(in) :: RIJ_INI; 
    real (kind=8), dimension(1:3), intent(out) :: RIJ_FIN; 

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    real (kind=8), dimension(1:3) :: DRIJ, DRIJNO(1:3); 

!   ************************************************************************************************

    RIJ_FIN(1:3) = RIJ_INI(1:3) + DOPLE_TRANS(1:3);
    DRIJNO(1:3)  = MATMUL(PASSB(1:3,1:3),RIJ_FIN(1:3));
    DRIJNO(1:3)  = DRIJNO(1:3) - ANINT( DRIJNO(1:3) );
    RIJ_FIN(1:3) = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));

end subroutine SET_TRANSLATION




