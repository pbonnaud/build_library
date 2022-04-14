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

subroutine APPLY_ION_INSERTION(ADD_NLABEL,  &
                               ADD_LABEL,   &
                               ADD_IONS,    &
                               CELL_AXIS,   &
                               CELL_ANGDEG, &
                               PASSA,       &
                               PASSB)

!
!  SUBROUTINE TO INSERT CALCIUM IONS IN ORDER TO KEEP SYSTEM ELECTRONEUTRALITY
! 

!   ************************************************************************************************
!   **                                                                                            **
!   ** ADD_NLABEL : NUMBER DIFFERENT LABEL TO MODIFY
!   ** ADD_LABEL  : LABEL TO MODIFY
!   ** ADD_NIONS  : NUMBER OF IONS PER LABEL TO INSERT
!   **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************    

    implicit none;

!   ************************************************************************************************    

    real (kind=8), dimension(1:3), intent(in) :: CELL_AXIS, CELL_ANGDEG;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************    

    integer (kind=4), intent(in) :: ADD_NLABEL;

    integer (kind=4), dimension(1:10), intent(in) :: ADD_LABEL, ADD_IONS;

!   ************************************************************************************************    

    integer (kind=4) :: i, j, h;

    integer (kind=4) :: iespece, insp, islab;

    integer (kind=4) :: IMATCH;

    integer (kind=4), dimension(1:10) :: ADD_IESPECE;

    character (len=3), dimension(1:10) :: ATOM_LABEL;

    real (kind=8) :: grnd;

    real (kind=8) :: RIJ;

    real (kind=8), dimension(1:3) :: RINO, RI, DRIJNO, DRIJ; 

!   ************************************************************************************************    

    ADD_IESPECE(1:ADD_NLABEL) = 0;

!   do h = 1, ADD_NLABEL
!       do i = 1, espece
!           write(99,*) i,  nsite(i);
!           do j = 1, nsite(i)
!               write(99,*) j,  mol(i,j)%nat, ATOM_LABEL(ADD_LABEL(h));
!               if (  mol(i,j)%nat == ATOM_LABEL(ADD_LABEL(h)) ) then
!                   ADD_IESPECE(h) = i;
!                   write(99,*) 'OK', ADD_IESPECE(h);
!               end if;
!           end do
!       end do

!   end do

    write(99,*);
    write(99,*) 'INSERTION OF IONS : ';
    write(99,*);

!   do i = 1, ADD_NLABEL
!       write(99,*) i, ADD_LABEL(i), ADD_IONS(i);
!       write(99,*);

!       do j = 1, ADD_IONS(i);

!           IMATCH = 1;

!           do while ( IMATCH == 1 ) 
!               IMATCH = 0;

!               RINO(1) = grnd() - 0.5d0;
!               RINO(2) = grnd() - 0.5d0;
!               RINO(3) = grnd() - 0.5d0;            

!               write(99,*) RINO(1:3);

!               RI(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));

!               write(99,*) RI(1:3);

!               do iespece = 1, espece
!                   do insp = 1, NATSP_TMP(iespece)
!                       if ( DATA_TMP(iespece)%NAT(insp) == 'XXX' ) CYCLE;
!                       DRIJ(1:3) = DATA_TMP(iespece)%RI(1:3,insp) - RI(1:3);
!                       DRIJNO(1:3) = MATMUL(PASSB(1:3,1:3),DRIJ(1:3));
!                       DRIJNO(1:3) = DRIJNO(1:3) - ANINT( DRIJNO(1:3) );
!                       DRIJ(1:3) = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));
!                       RIJ = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));

!                       if ( RIJ < 1.5d0 ) IMATCH = 1;
!                   end do
!               end do

!               write(99,*) 'OK FLUID';

!               do islab = 1, NSLAB_TMP
!                   DRIJ(1:3) = SUB_TMP(islab)%RI(1:3) - RI(1:3);
!                   DRIJNO(1:3) = MATMUL(PASSB(1:3,1:3),DRIJ(1:3));
!                   DRIJNO(1:3) = DRIJNO(1:3) - ANINT( DRIJNO(1:3) );
!                   DRIJ(1:3) = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));
!                   RIJ = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));

!                   if ( RIJ < 1.5d0 ) IMATCH = 1;
!               end do

!               write(99,*) 'OK SOLID';
!           end do

!           NATSP_TMP(ADD_IESPECE(i)) = NATSP_TMP(ADD_IESPECE(i)) + 1;

!           write(99,*) ADD_IESPECE(i), NATSP_TMP(ADD_IESPECE(i));

!           DATA_TMP(ADD_IESPECE(i))%NAT(NATSP_TMP(ADD_IESPECE(i))) = ATOM_LABEL(ADD_LABEL(i));
!           DATA_TMP(ADD_IESPECE(i))%RI(1:3,NATSP_TMP(ADD_IESPECE(i))) = RI(1:3);

!           write(99,'(i4,a3,3f15.6)') j, DATA_TMP(ADD_IESPECE(i))%NAT(NATSP_TMP(ADD_IESPECE(i))), &
!                                         DATA_TMP(ADD_IESPECE(i))%RI(1:3,NATSP_TMP(ADD_IESPECE(i)));
!       end do

!       write(99,*) '********************************************';
!       write(99,*);
!   end do

!   NATOM_TMP = NSLAB_TMP;

!   do i = 1, espece
!       NATOM_TMP = NATOM_TMP + NATSP_TMP(i);
!   end do

!   open(135,file='TMP5.xyz')
!   write(135,*) NATOM_TMP;
!   write(135,'(6f15.6)') CELL_AXIS(1:3), CELL_ANGDEG(1:3);

!   do i = 1, espece
!       do j = 1, NATSP_TMP(i);
!           write(135,'(a3,3f15.6)') DATA_TMP(i)%NAT(j), DATA_TMP(i)%RI(1:3,j);
!       end do
!   end do

!   do islab = 1, NSLAB_TMP
!       write(135,'(a3,3f15.6)') SUB_TMP(islab)%NAT, SUB_TMP(islab)%RI(1:3);
!   end do

!   close(135);

!   allocate(DATA_TMP(1:espece));

!   do i = 1, espece
!       allocate(DATA_TMP(i)%NAT(1:NATSP_PLATE(i)+5000));
!       allocate(DATA_TMP(i)%RI(1:3,1:NATSP_PLATE(i)+5000));

!       DATA_TMP(i)%NAT(1:NATSP_PLATE(i)) = DATA_PLATE(i)%NAT(1:NATSP_PLATE(i));
!       DATA_TMP(i)%RI(1:3,1:NATSP_PLATE(i)) = DATA_PLATE(i)%RI(1:3,1:NATSP_PLATE(i));
!   end do

!   NATSP_TMP(1:espece) = NATSP_PLATE(1:espece);
 
!   INT_ADDED_IONS = INT( ADDED_IONS );

!   write(99,*) 'NUMBER OF IONS TO INSERT : ', INT_ADDED_IONS;

!   do while ( INT_ADDED_IONS /= 0 )
!       TRYRINO(1) = grnd() - 0.5d0;
!       TRYRINO(2) = grnd() - 0.5d0;
!       TRYRINO(3) = grnd() - 0.5d0;

!       write(99,*) TRYRINO(1:3);

!       TRYRI(1:3) = MATMUL(PASSA(1:3,1:3),TRYRINO(1:3));

!       write(99,*) TRYRI(1:3);

!       write(99,*) 'NATSP_PLATE OLD : ', NATSP_PLATE(2);

!       NATSP_PLATE(2) = NATSP_PLATE(2) + 1;

!       write(99,*) 'NATSP_PLATE NEW : ', NATSP_PLATE(2);

!       DATA_TMP(2)%NAT(NATSP_PLATE(2)) = 'Cw';        

!       write(99,*) DATA_TMP(2)%NAT(1), DATA_TMP(2)%NAT(NATSP_PLATE(2));

!       DATA_TMP(2)%RI(1:3,NATSP_PLATE(2)) = TRYRI(1:3);

!       INT_ADDED_IONS = INT_ADDED_IONS - 1;
!   end do

!   deallocate(DATA_TMP);

!   stop;

end subroutine APPLY_ION_INSERTION

