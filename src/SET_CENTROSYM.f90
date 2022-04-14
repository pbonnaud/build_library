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

subroutine SET_CENTROSYM(PASSA,PASSB)

!   ************************************************************************************************
!   **                      RENDER THE SIMULATION BOX CENTROSYMMETRIC                             **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

    integer (kind=4) :: i, j;

    integer (kind=4) :: iespece, insp, islab;

    real (kind=8), dimension(1:3) :: RINO;

!   ************************************************************************************************

!   write(99,*) '--> SET CENTROSYM SIM BOX';

!   do j = 1, 2
!       do i = 1, NBINIT(j)
!           RINO(1:3) = MATMUL(PASSB(1:3,1:3),RINIT(i,1:3,j));
!           RINO(1:3) = RINO(1:3) - 0.5d0;
!           RINO(1:3) = RINO(1:3) - ANINT( RINO(1:3) );
!           RINIT(i,1:3,j) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!       end do
!   end do

!   do iespece = 1, espece
!       do insp = 1, NATSP(iespece)
!           RINO(1:3) = MATMUL(PASSB(1:3,1:3),DATA(iespece)%RI(1:3,insp));
!           RINO(1:3) = RINO(1:3) - 0.5d0;
!           RINO(1:3) = RINO(1:3) - ANINT( RINO(1:3) );
!           DATA(iespece)%RI(1:3,insp) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!       end do
!   end do

!   do iespece = 1, espece
!       do insp = 1, NATSP_FINAL(iespece)
!           RINO(1:3) = MATMUL(PASSB(1:3,1:3),DATA_FINAL(iespece)%RI(1:3,insp));
!           RINO(1:3) = RINO(1:3) - 0.5d0;
!           RINO(1:3) = RINO(1:3) - ANINT( RINO(1:3) );
!           DATA_FINAL(iespece)%RI(1:3,insp) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!       end do
!   end do

!   if ( chslab == 1 ) then
!       do islab = 1, NSLAB_FINAL

!           write(*,*) SUB_FINAL(islab)%RI(1:3);
!           RINO(1:3) = MATMUL(PASSB(1:3,1:3),SUB_FINAL(islab)%RI(1:3));
!           RINO(1:3) = RINO(1:3) - 0.5d0;
!           RINO(1:3) = RINO(1:3) - ANINT( RINO(1:3) );
!           SUB_FINAL(islab)%RI(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!       end do
!   end if

end subroutine SET_CENTROSYM
