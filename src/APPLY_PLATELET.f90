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

subroutine APPLY_PLATELET()

!   ************************************************************************************************
!   **                                   CREATE A PLATELET PARTICLE                               **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: islab;

    real (kind=8) :: RPLATELET, RXY;

!   ************************************************************************************************

!   RPLATELET = 25.0d0;

!   write(99,*) 'PLATELET RADIUS [A] : ', RPLATELET;
!   write(99,*);
!   write(99,*) 'ALLOCATE DATA_TMP : ';

!   allocate(DATA_TMP(1:espece+1));

!   write(99,*) 'NSLAB_TMP : ',  NSLAB_TMP;


!   do i = 1, espece
!       write(99,*) i, NATSP_FINAL(i);
!       allocate(DATA_TMP(i)%NAT(1:NATSP_FINAL(i)));
!       allocate(DATA_TMP(i)%RI(1:3,1:NATSP_FINAL(i)));
!   end do

!   write(99,*);

!   NATSP_TMP(1:10) = 0;

!   do i = 1, espece
!       do j = 1, NATSP_FINAL(i);
!           RXY = DSQRT( DATA_FINAL(i)%RI(1,j) * DATA_FINAL(i)%RI(1,j) + & 
!                        DATA_FINAL(i)%RI(2,j) * DATA_FINAL(i)%RI(2,j) );
!           if ( RXY <= RPLATELET ) then
!               NATSP_TMP(i) = NATSP_TMP(i) + 1;
!               DATA_TMP(i)%NAT(NATSP_TMP(i)) = DATA_FINAL(i)%NAT(j);
!               DATA_TMP(i)%RI(1:3,NATSP_TMP(i)) = DATA_FINAL(i)%RI(1:3,j);
!           end if 
!       end do
!   end do

!   NATOM_TMP = 0;

!   if ( chslab == 1 ) then
!       allocate(SUB_TMP(1:NSLAB_FINAL));

!       NSLAB_TMP = 0;

!       do i = 1, NSLAB_FINAL
!           write(*,*) SUB_FINAL(i)%RI(1:3)

!           RXY = DSQRT( SUB_FINAL(i)%RI(1) * SUB_FINAL(i)%RI(1) + &
!                        SUB_FINAL(i)%RI(2) * SUB_FINAL(i)%RI(2) );

!           if ( RXY <= RPLATELET ) then
!               NSLAB_TMP = NSLAB_TMP + 1;
!               SUB_TMP(NSLAB_TMP)%NAT     = SUB_FINAL(i)%NAT;
!               SUB_TMP(NSLAB_TMP)%RI(1:3) = SUB_FINAL(i)%RI(1:3);
!           end if

!       end do

!       NATOM_TMP = NSLAB_TMP;

!       write(99,*) 'NSLAB_TMP : ', NSLAB_TMP;

!   end if

!   do i = 1, espece
!       write(99,*) i, NATSP_TMP(i);
!       NATOM_TMP = NATOM_TMP + NATSP_TMP(i);
!   end do

!   open(135,file='PLATELET_BRUT.xyz')
!   write(135,*) NATOM_TMP;
!   write(135,*);

!   do i = 1, espece
!       do j = 1, NATSP_TMP(i);
!           write(135,'(a3,3f15.6)') DATA_TMP(i)%NAT(j), DATA_TMP(i)%RI(1:3,j);
!       end do
!   end do

!   do islab = 1, NSLAB_TMP
!       write(135,'(a3,3f15.6)') SUB_TMP(islab)%NAT, SUB_TMP(islab)%RI(1:3);
!   end do

!   close(135);

!   NATOM_TMP = NSLAB_FINAL;

!   do i = 1, espece
!       NATOM_TMP = NATOM_TMP + NATSP_FINAL(i);
!   end do

!   open(136,file='TMP-FINAL.xyz')
!   write(136,*) NATOM_TMP;
!   write(136,*);

!   do i = 1, espece
!       do j = 1, NATSP_FINAL(i);
!           write(136,'(a3,3f15.6)') DATA_FINAL(i)%NAT(j), DATA_FINAL(i)%RI(1:3,j);
!       end do
!   end do

!   do islab = 1, NSLAB_FINAL
!       write(136,'(a3,3f15.6)') SUB_FINAL(islab)%NAT, SUB_FINAL(islab)%RI(1:3);
!   end do

!   close(136);

!   write(99,*) 'OK 0';

!   allocate(DATA_PLATE(1:espece));

!   write(99,*) 'OK 01';

!   NATSP_PLATE(1:espece) = NATSP_TMP(1:espece);

!   do i = 1, espece
!       write(99,*) 'espece : ', i;
!       allocate(DATA_PLATE(i)%NAT(1:NATSP_TMP(i)));
!       allocate(DATA_PLATE(i)%RI(1:3,1:NATSP_TMP(i)));
!   end do

!   do i = 1, espece
!       DATA_PLATE(i)%NAT(1:NATSP_TMP(i)) = DATA_TMP(i)%NAT(1:NATSP_TMP(i)); 
!       DATA_PLATE(i)%RI(1:3,1:NATSP_TMP(i)) = DATA_TMP(i)%RI(1:3,1:NATSP_TMP(i));
!   end do

!   write(99,*) 'OK 1';

!   if ( chslab == 1 ) then
!       allocate(SUB_PLATE(1:NSLAB_TMP));        

!       NSLAB_PLATE = NSLAB_TMP;

!       SUB_PLATE(1:NSLAB_TMP)%NAT = SUB_TMP(1:NSLAB_TMP)%NAT;

!       do islab = 1, NSLAB_TMP
!           SUB_PLATE(islab)%RI(1:3) = SUB_TMP(islab)%RI(1:3);
!       end do
!   end if

!   write(99,*) 'OK 2';

!   stop;

end subroutine APPLY_PLATELET
