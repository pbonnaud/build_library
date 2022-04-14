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

subroutine CHECK_OVERLAPPING_ATOMS(icanal,ISTART,       &
                                          NATOM,        &
                                          CONFIG_RI,    &
                                          PASSA,        &
                                          PASSB,        &
                                          RCUTOFF,      &
                                          IOVERLAPPING)

!   ************************************************************************************************
!   **               THIS ROUTINE IS CHECKING IF THERE IS OVERLAPPING ATOMS                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   ** ISTART : 
!   ** NATOM   : NUMBER OF ATOMS IN THE MOLECULAR CONFIGURATION
!   ** CONFIG_RI : COORDINATES OF ATOMS IN THE MOLECULAR CONFIGURATION 
!   ** PASSA
!   ** PASSB
!   ** RCUTOFF
!   ** IOVERLAPPING : FLAG TO TELL THE PROGRAM IF ATOMS ARE OVERLAPPING 
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: ISTART, NATOM;
 
    real (kind=8), intent(in) :: RCUTOFF;

    real (kind=8), dimension(1:3,1:100000), intent(in) :: CONFIG_RI;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4), intent(out) :: IOVERLAPPING;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    real (kind=8) :: RIJ;

    real (kind=8), dimension(1:3) :: DRIJ;

!   ************************************************************************************************

    IOVERLAPPING = 0;

    do i = ISTART, NATOM
        do j = 1, ISTART-1
            DRIJ(1:3) = CONFIG_RI(1:3,j) - CONFIG_RI(1:3,i);
            
            call APPLY_PBC(DRIJ(1:3),DRIJ(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));

            RIJ = DSQRT( DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)) );

            if ( RIJ <= RCUTOFF ) IOVERLAPPING = IOVERLAPPING + 1;

        end do
    end do

end subroutine CHECK_OVERLAPPING_ATOMS











