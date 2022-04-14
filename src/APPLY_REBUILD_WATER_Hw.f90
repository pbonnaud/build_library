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

subroutine APPLY_REBUILD_WATER_Hw()

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: f, g, h, i, j, k;

    integer (kind=4) :: iespece, insp, jespece, jnsp, islab;

    real (kind=8) :: grnd;

!   ************************************************************************************************

    integer (kind=4) :: INNCOUNT, IADD, IFOUND1, NXXX, NOw, NHw;

    real (kind=8) :: RIJ;

    real (kind=8), dimension(1:3) :: DRIJ, DRIJNO, RIOW;

    real (kind=8), dimension(1:3,1:4) :: NEIGHBOR_COORD;

!   integer (kind=4) :: NATOM_TMP;

!   integer (kind=4), dimension(1:10) :: NATSP_TMP;

!   type (coord), allocatable, dimension(:) :: DATA_TMP;

!   ************************************************************************************************

!   write(99,*) 'SINGLE Hw ATOMS';
!   write(99,*);

!   write(99,*) 'NSLAB_TMP : ', NSLAB_TMP;

!   do i = 1, espece
!       write(99,*) i, NATSP_TMP(i);
!   end do
!   write(99,*);

!   allocate(DATA_TMP(1:espece));

!   NATSP_TMP(1:espece) = NATSP_PLATE(1:espece) + 2000;

!   do i = 1, espece
!       allocate(DATA_TMP(i)%NAT(1:NATSP_TMP(i)));
!       allocate(DATA_TMP(i)%RI(1:3,1:NATSP_TMP(i)));
!   end do

!   NATSP_TMP(1:espece) = NATSP_PLATE(1:espece);

!   do i = 1, espece
!       DATA_TMP(i)%NAT(1:NATSP_PLATE(i))    = DATA_PLATE(i)%NAT(1:NATSP_PLATE(i));
!       DATA_TMP(i)%RI(1:3,1:NATSP_PLATE(i)) = DATA_PLATE(i)%RI(1:3,1:NATSP_PLATE(i));
!   end do

!   NXXX = 0;

!   do iespece = 1, espece
!       do insp = 1, NATSP_PLATE(iespece); 
!       do insp = 1, NATSP_TMP(iespece); 
!           if ( DATA_PLATE(iespece)%NAT(insp) /= 'Hw' ) CYCLE;    

!           write(99,*) iespece, insp, DATA_TMP(iespece)%NAT(insp);


!           if ( DATA_TMP(iespece)%NAT(insp) /= 'Hw' ) CYCLE;    
!           INNCOUNT = 0;

!           do jespece = 1, espece
!               do jnsp = 1, NATSP_PLATE(jespece)
!               do jnsp = 1, NATSP_TMP(jespece)
!                    if ( DATA_PLATE(jespece)%NAT(jnsp) /= 'Ow' ) CYCLE;
!                    if ( DATA_TMP(jespece)%NAT(jnsp) /= 'Ow' ) CYCLE;
!                    DRIJ(1:3) = DATA_TMP(jespece)%RI(1:3,jnsp) - DATA_TMP(iespece)%RI(1:3,insp);
!                    RIJ = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));

!                    if ( RIJ >= 1.1d0 ) CYCLE;
!                    INNCOUNT = INNCOUNT + 1;
!                    NEIGHBOR_COORD(1:3,INNCOUNT) = DATA_PLATE(jespece)%RI(1:3,jnsp);
!                    NEIGHBOR_COORD(1:3,INNCOUNT) = DATA_TMP(jespece)%RI(1:3,jnsp);
!               end do
!           end do

!           if ( INNCOUNT == 0 ) then
!               write(99,'(2i6,a3,i4)') iespece, insp, DATA_TMP(iespece)%NAT(insp), INNCOUNT;

!               DATA_TMP(iespece)%NAT(insp) = 'XXX';

!               NXXX = NXXX + 1;
!           end if
!       end do
!   end do

!   write(99,*) 'THERE IS ', NXXX, ' Hw TO BE REMOVED';
!   write(99,*);

!   NOw = 0;
!   NHw = 0;

!   do iespece = 1, espece
!       do insp = 1, NATSP_TMP(iespece);
!           if ( DATA_TMP(iespece)%NAT(insp) == 'Hw' ) NHw = NHw + 1;
!           if ( DATA_TMP(iespece)%NAT(insp) == 'Ow' ) NOw = NOw + 1;
!       end do
!   end do

!   write(99,*) 'NOw : ', NOw;
!   write(99,*) 'NHw : ', NHw;
!   write(99,*);


!   NATOM_TMP = NSLAB_PLATE;
!   NATOM_TMP = NSLAB_TMP;

!   do i = 1, espece
!       NATOM_TMP = NATOM_TMP + NATSP_TMP(i);
!   end do

!   open(135,file='TMP4.xyz')
!   write(135,*) NATOM_TMP;
!   write(135,*);

!   do i = 1, espece
!       do j = 1, NATSP_TMP(i);
!           write(135,'(a3,3f15.6)') DATA_TMP(i)%NAT(j), DATA_TMP(i)%RI(1:3,j);
!       end do
!   end do

!   do islab = 1, NSLAB_PLATE
!   do islab = 1, NSLAB_TMP
!       write(135,'(a3,3f15.6)') SUB_PLATE(islab)%NAT, SUB_PLATE(islab)%RI(1:3);
!       write(135,'(a3,3f15.6)') SUB_TMP(islab)%NAT, SUB_TMP(islab)%RI(1:3);
!   end do

!   close(135);

!   deallocate(DATA_PLATE);

!   NATSP_PLATE(1:espece) = NATSP_TMP(1:espece);

!   allocate(DATA_PLATE(1:espece));

!   do i = 1, espece
!       allocate(DATA_PLATE(i)%NAT(1:NATSP_PLATE(i)));
!       allocate(DATA_PLATE(i)%RI(1:3,1:NATSP_PLATE(i)));
!   end do

!   do i = 1, espece
!       DATA_PLATE(i)%NAT(1:NATSP_PLATE(i))    = DATA_TMP(i)%NAT(1:NATSP_PLATE(i));
!       DATA_PLATE(i)%RI(1:3,1:NATSP_PLATE(i)) = DATA_TMP(i)%RI(1:3,1:NATSP_PLATE(i));
!   end do

!   write(99,*) 'OKAY TRANSFER TABLE Hw WATER REBUILD';

!   deallocate(NEIGHBOR_LIST,NEIGHBOR_NBRE,COORD_ADDED,NAT_ADDED);

!   stop;

end subroutine APPLY_REBUILD_WATER_Hw

