!   //////////////////////////////////////////////////////////////////////////////////////////////////
!   //                                                                                              //
!   //  Copyright (2022) Patrick A. Bonnaud                                                         //
!   //                                                                                              //
!   //  This file is part of BUILDLIBRARY (Build Molecular Models for a Library of Molecules).      //
!   //                                                                                              //
!   //  BUILDLIBRARY is free software; you can redistribute it and/or modify it under the terms of  //
!   //  the GNU General Public License as published by the Free Software Foundation; either         //
!   //  version 2 of the License, or (at your option) any later version.                            //
!   //                                                                                              //
!   //  BUILDLIBRARY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;   //
!   //  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   //
!   //  See the GNU General Public License for more details.                                        //
!   //                                                                                              //
!   //  You should have received a copy of the GNU General Public License along with this program.  //
!   //  If not, see <http://www.gnu.org/licenses/>.                                                 //
!   //                                                                                              //
!   //////////////////////////////////////////////////////////////////////////////////////////////////

subroutine APPLY_PBC(RIJ_IN,RIJ_OUT,PASSA,PASSB)

!   ***********************************************************************************************
!   **                APPLY PERIODIC BOUNDARY CONDITIONS IN THE TRICLINIC BOX                    **
!   ***********************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3), intent(in) :: RIJ_IN;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

    real (kind=8), dimension(1:3), intent(out) :: RIJ_OUT;

!   ************************************************************************************************

    real (kind=8), dimension(1:3) :: RIJNO;

!   ************************************************************************************************

    RIJNO(1:3) = MATMUL(PASSB(1:3,1:3),RIJ_IN(1:3));
    
    RIJNO(1:3) = RIJNO(1:3) - ANINT( RIJNO(1:3) );

    RIJ_OUT(1:3) = MATMUL(PASSA(1:3,1:3),RIJNO(1:3));

end subroutine APPLY_PBC

