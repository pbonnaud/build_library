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

subroutine APPLY_REBUILD_SILICA_CHAINS()

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: f, g, h, i, j, k;

    integer (kind=4) :: islab, jslab;

    real (kind=8) :: grnd;

!   ************************************************************************************************

    integer (kind=4) :: INNCOUNT, IADD, IFOUND1, NSLAB_TMPXXX, NSi, NO, NCa, NOs;

    real (kind=8) :: RIJ, NORM1;

    real (kind=8), dimension(1:3) :: DRIJ, DRIJNO, VECT1;

    real (kind=8), dimension(1:3,1:4) :: NEIGHBOR_COORD, NEIGHBOR_VECT;

!   ************************************************************************************************

!   NSLAB_TMPXXX = NSLAB_TMP;

!   if ( chslab == 1 ) then
!       do islab = 1, NSLAB_PLATE
!       do islab = 1, NSLAB_TMP
!           if ( SUB_PLATE(islab)%NAT /= 'Si' ) CYCLE;    
!           if ( SUB_TMP(islab)%NAT /= 'Si' ) CYCLE;    
!           INNCOUNT = 0;
!           do jslab = 1, NSLAB_PLATE
!           do jslab = 1, NSLAB_TMP
!               if ( SUB_PLATE(jslab)%NAT /= 'O' ) CYCLE;
!               if ( SUB_TMP(jslab)%NAT /= 'O' ) CYCLE;
!               DRIJ(1:3) = SUB_PLATE(jslab)%RI(1:3) - SUB_PLATE(islab)%RI(1:3);
!               DRIJ(1:3) = SUB_TMP(jslab)%RI(1:3) - SUB_TMP(islab)%RI(1:3);
!               RIJ_WOPDB = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));                        
!               DRIJNO(1:3) = MATMUL(PASSB(1:3,1:3),DRIJ(1:3)); 
!               DRIJNO(1:3) = DRIJNO(1:3) - ANINT( DRIJNO(1:3) );
!               DRIJ(1:3) = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));
!               RIJ = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));
!               write(99,*) 'RIJ : ', RIJ;

!               if ( RIJ >= 2.5d0 ) CYCLE;
!               INNCOUNT = INNCOUNT + 1;
!               NEIGHBOR_COORD(1:3,INNCOUNT) = SUB_TMP(jslab)%RI(1:3);
!               write(99,*) 'O', INNCOUNT, RIJ;
!           end do

!           write(99,*);

!           if ( INNCOUNT < 4 ) then
!               write(99,*) islab, SUB_PLATE(islab)%NAT, INNCOUNT;

!               IADD = 4 - INNCOUNT;

!               write(99,*) 'IADD : ', IADD;

!               do jslab = 1, NSLAB_FINAL
!                   if ( SUB_FINAL(jslab)%NAT /= 'O' ) CYCLE;
!                   IFOUND1 = 0;
!                   do i = 1, INNCOUNT
!                       if ( SUB_FINAL(jslab)%RI(1) == NEIGHBOR_COORD(1,i) .AND. &
!                            SUB_FINAL(jslab)%RI(2) == NEIGHBOR_COORD(2,i) .AND. &
!                            SUB_FINAL(jslab)%RI(3) == NEIGHBOR_COORD(3,i) ) then
!                           IFOUND1 = 1;
!                           EXIT;
!                       end if 
!                   end do

!                   if ( IFOUND1 == 1 ) CYCLE;

!                   DRIJ(1:3) = SUB_FINAL(jslab)%RI(1:3) - SUB_TMP(islab)%RI(1:3);
!                   DRIJNO(1:3) = MATMUL(PASSB(1:3,1:3),DRIJ(1:3));
!                   DRIJNO(1:3) = DRIJNO(1:3) - ANINT( DRIJNO(1:3) );
!                   DRIJ(1:3) = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));
!                   RIJ = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));

!                   if ( RIJ >= 2.5d0 ) CYCLE;
!                   write(99,*) IADD, RIJ;
!                   write(99,'(a3,3f15.6)') SUB_FINAL(jslab)%NAT, SUB_FINAL(jslab)%RI(1:3);
!                   write(99,*);

!                   NSLAB_TMPXXX = NSLAB_TMPXXX + 1;

!                   SUB_TMP(NSLAB_TMPXXX)%NAT     = SUB_FINAL(jslab)%NAT;
!                   SUB_TMP(NSLAB_TMP)%RI(1:3) = SUB_FINAL(jslab)%RI(1:3);
!                   SUB_TMP(NSLAB_TMPXXX)%RI(1:3) = SUB_TMP(islab)%RI(1:3) + DRIJ(1:3);

!                   IADD = IADD - 1;    

!                   INNCOUNT = INNCOUNT + 1;
!                   NEIGHBOR_COORD(1:3,INNCOUNT) = SUB_TMP(NSLAB_TMPXXX)%RI(1:3);

!                   if ( IADD == 0 ) EXIT;
!               end do

!               if ( IADD > 0 ) then
!                  write(99,*) 'NEED TO ADD', IADD, ' O';
!                  write(99,*);
!                  do i = 1, INNCOUNT;
!                      NEIGHBOR_VECT(1:3,i) = NEIGHBOR_COORD(1:3,i) - SUB_TMP(islab)%RI(1:3);
!                      VECT1(1:3) = VECT1(1:3) + NEIGHBOR_VECT(1:3,i);
!                  end do

!                  VECT1(1:3) = - VECT1(1:3);
!                  NORM1 = DSQRT(DOT_PRODUCT(VECT1(1:3),VECT1(1:3)));
!                  VECT1(1:3) = 1.65d0 * VECT1(1:3) / NORM1;

!                  NSLAB_TMPXXX = NSLAB_TMPXXX + 1;
!                  SUB_TMP(NSLAB_TMPXXX)%NAT = 'O';
!                  SUB_TMP(NSLAB_TMPXXX)%RI(1:3) = SUB_TMP(islab)%RI(1:3) + VECT1(1:3);

!                  write(99,'(a3,3f15.6)') SUB_TMP(NSLAB_TMPXXX)%NAT, SUB_TMP(NSLAB_TMPXXX)%RI(1:3);
!                  write(99,*);
!               end if

!               write(99,*) '*************************************************';
!           end if
!       end do
!   end if


!   NSLAB_TMP = NSLAB_TMPXXX;

!   NSi = 0;
!   NO  = 0;
!   NCa = 0;
!   NOs = 0;

!   do islab = 1, NSLAB_TMP
!       if ( SUB_TMP(islab)%NAT == 'Si' ) NSi = NSi + 1;
!       if ( SUB_TMP(islab)%NAT == 'O'  ) NO  = NO  + 1;
!       if ( SUB_TMP(islab)%NAT == 'Ca' ) NCa = NCa + 1;
!       if ( SUB_TMP(islab)%NAT == 'Os' ) NOs = NOs + 1;
!   end do

!   write(99,*) 'NSi : ', NSi;
!   write(99,*) 'NO  : ', NO;
!   write(99,*) 'NCa : ', NCa;
!   write(99,*) 'NOs : ', NOs;
!   write(99,*);

!   NATOM_TMP = NSLAB_TMP;

!   do i = 1, espece
!       NATOM_TMP = NATOM_TMP + NATSP_TMP(i);
!   end do

!   open(135,file='TMP2.xyz')
!   write(135,*) NATOM_TMP;
!   write(135,*);

!   do i = 1, espece
!       do j = 1, NATSP_TMP(i);
!           if ( DATA_TMP(i)%NAT(j) == 'Hw' ) CYCLE;
!           if ( DATA_TMP(i)%NAT(j) == 'Ow' ) CYCLE;
!           write(135,'(a3,3f15.6)') DATA_TMP(i)%NAT(j), DATA_TMP(i)%RI(1:3,j);
!       end do
!   end do

!   do islab = 1, NSLAB_TMP
!       write(135,'(a3,3f15.6)') SUB_TMP(islab)%NAT, SUB_TMP(islab)%RI(1:3);
!   end do

!   close(135);

!   deallocate(SUB_PLATE);

!   NSLAB_PLATE = NSLAB_TMP;

!   allocate(SUB_PLATE(1:NSLAB_PLATE));

!   SUB_PLATE(1:NSLAB_PLATE)%NAT = SUB_TMP(1:NSLAB_PLATE)%NAT;

!   do islab = 1, NSLAB_PLATE
!       SUB_PLATE(islab)%RI(1:3) = SUB_TMP(islab)%RI(1:3);
!   end do

!   write(99,*) 'OKAY TRANSFER TABLE SILICA CHAINS REBUILD';

!   deallocate(NEIGHBOR_LIST,NEIGHBOR_NBRE,COORD_ADDED,NAT_ADDED);

!   stop;

end subroutine APPLY_REBUILD_SILICA_CHAINS

