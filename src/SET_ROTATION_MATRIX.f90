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

subroutine SET_ROTATION(RIJ_INI,RI_COR,RIJ_FIN,PASSA,PASSB)

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3), intent(in) :: RIJ_INI, RI_COR; 
    real (kind=8), dimension(1:3), intent(out) :: RIJ_FIN; 

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    real (kind=8), dimension(1:3) :: DRIJ, DRIJNO(1:3); 

!   ************************************************************************************************

!   DRIJ(1:3)   = RIJ_INI(1:3) - RI_COR(1:3);
!   DRIJNO(1:3) = MATMUL(PASSB(1:3,1:3),DRIJ(1:3));
!   DRIJNO(1:3) = DRIJNO(1:3) - ANINT( DRIJNO(1:3) );
!   DRIJ(1:3)   = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));

!   DRIJ(1:3) = MATMUL(MATROT(1:3,1:3),DRIJ(1:3));

!   DRIJNO(1:3) = MATMUL(PASSB(1:3,1:3),DRIJ(1:3));
!   DRIJNO(1:3) = DRIJNO(1:3) - ANINT( DRIJNO(1:3) );
!   DRIJ(1:3)   = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));

!   RIJ_FIN(1:3) =  RI_COR(1:3) + DRIJ(1:3);

    RIJ_FIN(1:3) = 0.0d0;

end subroutine SET_ROTATION

subroutine SET_ROTATION_MATRIX(ANGRAD,MATROT)

!   **
!   ** THIS ROUTINE COMPUTE THE ROTATIONAL MATRIX WITH THE QUATERNION TECHNIQUE BASED ON INPUT ANGLES
!   **

!   ** ANGRAD : ANGLES IN [RAD]
!   ** MATROT : ROTATIONAL MATRIX

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3), intent(in) :: ANGRAD;
    
!   ************************************************************************************************

    real (kind=8), dimension(1:3,1:3), intent(out) :: MATROT;

!   ************************************************************************************************

    real (kind=8) :: CSTE1, CSTE2, CSTE3;

    real (kind=8) :: COS_CSTE1, SIN_CSTE1;

    real (kind=8) :: Q0, Q1, Q2, Q3;
    real (kind=8) :: Q0SQ, Q1SQ, Q2SQ, Q3SQ;

!   ************************************************************************************************

    CSTE1 = 0.5d0 * ANGRAD(2);
    CSTE2 = 0.5d0 * ( ANGRAD(3) + ANGRAD(1) );
    CSTE3 = 0.5d0 * ( ANGRAD(3) - ANGRAD(1) );

    COS_CSTE1 = DCOS( CSTE1 );
    SIN_CSTE1 = DSIN( CSTE1 );

    Q0 = COS_CSTE1 * DCOS( CSTE2 );
    Q1 = SIN_CSTE1 * DCOS( CSTE3 );
    Q2 = SIN_CSTE1 * DSIN( CSTE3 );
    Q3 = COS_CSTE1 * DSIN( CSTE2 );

    Q0SQ = Q0 * Q0;
    Q1SQ = Q1 * Q1;
    Q2SQ = Q2 * Q2;
    Q3SQ = Q3 * Q3;

    MATROT(1,1) = Q0SQ + Q1SQ - Q2SQ - Q3SQ;
    MATROT(2,1) = 2.0d0 * ( Q1 * Q2 - Q0 * Q3 );
    MATROT(3,1) = 2.0d0 * ( Q1 * Q3 + Q0 * Q2 );
    MATROT(1,2) = 2.0d0 * ( Q1 * Q2 + Q0 * Q3 );
    MATROT(2,2) = Q0SQ - Q1SQ + Q2SQ - Q3SQ;
    MATROT(3,2) = 2.0d0 * ( Q2 * Q3 - Q0 * Q1 );
    MATROT(1,3) = 2.0d0 * ( Q1 * Q3 - Q0 * Q2 );
    MATROT(2,3) = 2.0d0 * ( Q2 * Q3 + Q0 * Q1 );
    MATROT(3,3) = Q0SQ - Q1SQ - Q2SQ + Q3SQ;

end subroutine SET_ROTATION_MATRIX


