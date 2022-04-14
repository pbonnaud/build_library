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

subroutine WRITE_FIELD() !PATH_DIRECTORY,CH_WORKING_FILE)

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

!   character (len=200), intent(in) :: PATH_DIRECTORY;

!   character (len=5), intent(in) :: CH_WORKING_FILE;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j, k;

    integer (kind=4) :: imax;

    integer (kind=4) :: IOSEF1;

    integer (kind=4) :: MULTI1, IFLAG1;

    real (kind=8) :: MASSE1, CHARGE1;

    character (len=3) :: LABAT1; 

    character (len=100) :: CHFIELD_NAME;

    character (len=3), dimension(1:10) :: ATOM_LABEL;

    integer (kind=4), dimension(1:10,1:10) :: TESTSPE;

    character (len=200) :: CH_DATA_FILE_NAME;

    logical :: PROBE1;

!   ************************************************************************************************

!   CH_DATA_FILE_NAME = trim(PATH_DIRECTORY)//trim(CH_WORKING_FILE)//'_FIELD.dat';

!   write(99,*) CH_DATA_FILE_NAME;

!   inquire(FILE=trim(CH_DATA_FILE_NAME),EXIST=PROBE1);

!   if (PROBE1 == .TRUE. ) then
!       write(99,*) trim(CH_DATA_FILE_NAME)//' ALREADY EXIST !!!';
!       write(99,*) 'END OF PROGRAM';
!       close(99);
!       stop;
!   end if

!   open(105,file=trim(CH_DATA_FILE_NAME),status='new');
    open(105,file='001_FIELD.dat');

    do j = 1, NLABEL
        write(105,'(a3,1x,f19.13)') ATOM_LABEL(j), ATOM_CHAR(j);
    end do

    write(105,'(a6)') 'FINISH';
    write(105,'(a36)')  '1         1.000000        109.470000';
    write(105,'(a6)') 'FINISH';

!   WRITE POSITION OF THE COM AND EULER ANGLES
    do i = 1, ICOMA
        write(105,'(6f15.6)') COMA_COORD_FINAL(1:3,i), 0.0d0, 0.0d0, 0.0d0;
    end do

    write(105,'(a6)') 'FINISH';

    close(105);

end subroutine WRITE_FIELD
