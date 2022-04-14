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

subroutine APPLY_REPLICA_CREATION(PASSA,PASSB)

!   ************************************************************************************************
!   **                     APPLY REPLICA CREATION OF THE ORIGINAL UNIT CELL                       **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: iespece, insp, islab;

    real (kind=8), dimension(1:3) :: TRANSR, RINO;

!   ************************************************************************************************

!   allocate(DATA_FINAL(1:espece));

!   NATSP_FINAL(1:espece) = NATSP(1:espece) * NMULTI;

!   do i = 1, espece
!       allocate(DATA_FINAL(i)%NAT(1:NATSP_FINAL(i)));
!       allocate(DATA_FINAL(i)%RI(1:3,1:NATSP_FINAL(i)));
!   end do

!   NATSP_TMP(1:10) = 0;

!   TRANSR(1) = - ( REAL(MULTI(1)) - 1.0d0 ) / ( 2.0d0 * REAL(MULTI(1)) );

!   do k = 1, MULTI(1)
!       TRANSR(2) = - ( REAL(MULTI(2)) - 1.0d0 ) / ( 2.0d0 * REAL(MULTI(2)) );
!       do j = 1, MULTI(2)
!          TRANSR(3) = - ( REAL(MULTI(3)) - 1.0d0 ) / ( 2.0d0 * REAL(MULTI(3)) );
!           do i = 1, MULTI(3)
!               do iespece = 1, espece
!                   do insp = 1, NATSP(iespece)
!                       NATSP_TMP(iespece) = NATSP_TMP(iespece) + 1;

!                       DATA_FINAL(iespece)%NAT(NATSP_TMP(iespece))    = DATA(iespece)%NAT(insp); 
!                       DATA_FINAL(iespece)%RI(1:3,NATSP_TMP(iespece)) = DATA(iespece)%RI(1:3,insp);

!                       RINO(1:3) = MATMUL(PASSB(1:3,1:3),DATA_FINAL(iespece)%RI(1:3,NATSP_TMP(iespece)));
!                       RINO(1:3) = RINO(1:3) + TRANSR(1:3);
!                       DATA_FINAL(iespece)%RI(1:3,NATSP_TMP(iespece)) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!                   end do
!               end do

!               if ( MULTI(3) /= 1 ) TRANSR(3) = TRANSR(3) + 1.0d0 / REAL( MULTI(3) );
!           end do
!           TRANSR(2) = TRANSR(2) + 1.0d0 / REAL( MULTI(2) );
!       end do
!       TRANSR(1) = TRANSR(1) + 1.0d0 / REAL( MULTI(1) );
!   end do

!   if ( chslab == 1 )  then
!       write(99,*) '*** MULTI SUB ***********';

!       NSLAB_FINAL = NSLAB * NMULTI;

!       NSLAB_TMP = 0;

!       allocate(SUB_FINAL(1:NSLAB_FINAL));

!       TRANSR(1) = - ( REAL(MULTI(1)) - 1.0d0 ) / ( 2.0d0 * REAL(MULTI(1)) );

!       do k = 1, MULTI(1)
!           TRANSR(2) = - ( REAL(MULTI(2)) - 1.0d0 ) / ( 2.0d0 * REAL(MULTI(2)) );
!           do j = 1, MULTI(2)
!              TRANSR(3) = - ( REAL(MULTI(3)) - 1.0d0 ) / ( 2.0d0 * REAL(MULTI(3)));
!               do i = 1, MULTI(3)

!                   do islab = 1, NSLAB


!                       NSLAB_TMP = NSLAB_TMP + 1;

!                       SUB_FINAL(NSLAB_TMP)%NAT = SUB(islab)%NAT;

!                       RINO(1:3) = MATMUL(PASSB(1:3,1:3),SUB(islab)%RI(1:3));

!                       RINO(1:3) = RINO(1:3) + TRANSR(1:3);

!                       SUB_FINAL(NSLAB_TMP)%RI(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!                   end do

!                   if ( MULTI(3) /= 1 ) TRANSR(3) = TRANSR(3) + 1.0d0 / REAL( MULTI(3) );
!              end do
!              TRANSR(2) = TRANSR(2) + 1.0d0 / REAL( MULTI(2) );
!           end do
!           TRANSR(1) = TRANSR(1) + 1.0d0 / REAL( MULTI(1) );
!        end do

!   end if

!   write(99,*);
!   write(99,*) 'END REPLICA CREATION';
!   write(99,*);

end subroutine APPLY_REPLICA_CREATION
