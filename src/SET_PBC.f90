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

subroutine SET_PBC(PASSA,PASSB)

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

    integer (kind=4) :: h, i;

    real (kind=8), dimension(1:3) :: RINO;

!   ************************************************************************************************

    do h = 1, 2
        do i = 1, NBFINAL(h)
            if ( NATFINAL(i,h) == 'Ow' .OR. NATFINAL(i,h) == 'Hw' ) then
                if ( idople == 0 ) then
                    RINO(1:3)  = MATMUL(PASSB(1:3,1:3),RFINAL(i,1:3,h));
                    RINO(1:3)  = RINO(1:3) - ANINT( RINO(1:3) );
                    RFINAL(i,1:3,h) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
                end if
            end if
        end do
    end do

    do h = 1, 2
        do i = 1, NBFINAL(h)
            if ( NATFINAL(i,h) /= 'Ow' .OR. NATFINAL(i,h) /= 'Hw' ) then
                RINO(1:3) = MATMUL(PASSB(1:3,1:3),RFINAL(i,1:3,h));
                if ( HPORE(1) == 0.0d0 ) RINO(1) = RINO(1) - ANINT( RINO(1) )
                if ( HPORE(2) == 0.0d0 ) RINO(2) = RINO(2) - ANINT( RINO(2) )
                if ( HPORE(3) == 0.0d0 ) RINO(3) = RINO(3) - ANINT( RINO(3) )
                RFINAL(i,1:3,h) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
            end if
        end do
    end do

end subroutine SET_PBC






