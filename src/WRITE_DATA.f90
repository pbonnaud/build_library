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

subroutine WRITE_DATA()!PATH_DIRECTORY,CH_WORKING_FILE)

!   ************************************************************************************************
!   **                                     READ DATA FILE                                         **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

!   character (len=200), intent(in) :: PATH_DIRECTORY;

!   character (len=5), intent(in) :: CH_WORKING_FILE;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

    integer (kind=4) :: iespece;

    integer (kind=4) :: EOF;

    integer (kind=4) :: SAME_AT;

    character (len=3) :: CHOSEF1;

    character (len=150) :: CHARLINE;

    character (len=200) :: CH_DATA_FILE_NAME;

    logical :: PROBE1;

!   ************************************************************************************************

!   CH_DATA_FILE_NAME = trim(PATH_DIRECTORY)//trim(CH_WORKING_FILE)//'_data.dat';

!   write(99,*) CH_DATA_FILE_NAME;

!   inquire(FILE=trim(CH_DATA_FILE_NAME),EXIST=PROBE1);

!   if (PROBE1 == .FALSE. ) then
!       write(99,*) trim(CH_DATA_FILE_NAME)//' DOESNT EXIST !!!';
!       write(99,*) 'END OF PROGRAM';
!       close(99);
!   end if

!   chslab = 0;

!   open(5,file=trim(CH_DATA_FILE_NAME),status='old');
!   open(150,file='001_data.dat');
!   write(150,*) espece;
!   write(150,*);

!   do i = 1, espece
!       write(150,*) nsite(i), N(i);
!   end do
!   write(150,*);

!   do i = 1, espece
!       do j = 1, nsite(i)
!           write(150,'(a3,1x,a2,i4,7f15.6)') mol(i,j)%nat,       & 
!                                             mol(i,j)%CHDISPREP, &
!                                             mol(i,j)%MULTI,     &
!                                             mol(i,j)%sigma,     &
!                                             mol(i,j)%eps,       &
!                                             mol(i,j)%RI(1:3),   &
!                                             mol(i,j)%q,         &
!                                             mol(i,j)%M;
!       end do
!       write(150,*)
!   end do

!   write(150,*) chslab;
!   write(150,*);
!   write(150,*) NSLAB_FINAL;
!   write(150,*);
!   write(150,*) '           7           4           1';
!   write(150,*);
!   write(150,*);
!   write(150,*);
!   write(150,*);
!   write(150,'(a87)') 'Ow  Si  pn         10.06864       145.36364      4673.83238      1235.58657     2.22355';
!   write(150,'(a87)') 'Ow  O   pn         35.53046       603.07100     20553.64650       618.32431     2.13095';
!   write(150,'(a87)') 'Ow  Os  pn         35.53046       603.07100     20553.64650       618.32431     2.13095';
!   write(150,'(a87)') 'Ow  Ca  pn         17.90063       343.06663     11191.49960       495.36714     2.21464';
!   write(150,*);
!   write(150,'(a87)') 'Hw  Si  pn          2.41564        31.44501       964.78805        90.81095     2.24349';
!   write(150,'(a87)') 'Hw  O   pn          8.42528       133.47966      4167.92019        45.44450     2.14926';
!   write(150,'(a87)') 'Hw  Os  pn          8.42528       133.47966      4167.92019        45.44450     2.14926';
!   write(150,'(a87)') 'Hw  Ca  pn          3.90600        75.51987      2210.72217        36.40762     2.23441';
!   write(150,*);
!   write(150,'(a87)') 'Cw  Si  pn          5.43028        60.54502       985.96590      2471.00924     2.38467';
!   write(150,'(a87)') 'Cw  O   pn         19.36154       260.60148      5650.27847      1236.56659     2.27848';
!   write(150,'(a87)') 'Cw  Os  pn         19.36154       260.60148      5650.27847      1236.56659     2.27848';
!   write(150,'(a87)') 'Cw  Ca  pn         10.57397       157.33140      1636.85374       990.66857     2.37442';
!   write(150,*);
!   write(150,'(a87)') 'Os  Si  pn         10.84292       155.52112      2947.77015      3084.34886     2.28792';
!   write(150,'(a87)') 'Os  O   pn         38.28249       645.71746     14527.40193      1543.50000     2.19000';
!   write(150,'(a87)') 'Os  Os  pn         38.28249       645.71746     14527.40193      1543.50000     2.19000';
!   write(150,'(a87)') 'Os  Ca  pn         19.36154       368.08494      6782.42620      1236.56659     2.27848';
!   write(150,*);
!   write(150,'(a87)') 'Si  Si  pn          3.07514        36.18508       556.68384      6163.40000     2.39500';
!   write(150,'(a87)') 'Si  O   pn         10.84292       155.52112      2947.77015      3084.34886     2.28792';
!   write(150,'(a87)') 'Si  Os  pn         10.84292       155.52112      2947.77015      3084.34886     2.28792';
!   write(150,'(a87)') 'Si  Ca  pn          5.43028        91.13018      1160.39690      2471.00924     2.38467';
!   write(150,*);
!   write(150,'(a87)') 'O   Si  pn         10.84292       155.52112      2947.77015      3084.34886     2.28792';
!   write(150,'(a87)') 'O   O   pn         38.28249       645.71746     14527.40193      1543.50000     2.19000';
!   write(150,'(a87)') 'O   Os  pn         38.28249       645.71746     14527.40193      1543.50000     2.19000';
!   write(150,'(a87)') 'O   Ca  pn         19.36154       368.08494      6782.42620      1236.56659     2.27848';
!   write(150,*);
!   write(150,'(a87)') 'Ca  Si  pn          5.43028        91.13018      1160.39690      2471.00924     2.38467';
!   write(150,'(a87)') 'Ca  O   pn         19.36154       368.08494      6782.42620      1236.56659     2.27848';
!   write(150,'(a87)') 'Ca  Os  pn         19.36154       368.08494      6782.42620      1236.56659     2.27848';
!   write(150,'(a87)') 'Ca  Ca  pn         10.57397       209.44569      2404.48301       990.66857     2.37442';
!   write(150,*);
!   write(150,*) '   53.2337000000000        59.0416140000000        23.6911390000000';
!   close(150);

end subroutine WRITE_DATA

