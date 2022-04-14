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

subroutine WRITE_CONFIG(PATH_DIRECTORY,CH_WORKING_FILE,CELL_AXIS,CELL_ANGDEG);

!   ************************************************************************************************
!   **                           WRITE THE MOLECULAR CONFIGURATION                                **
!   ************************************************************************************************
!   **                                                                                            **
!   ** PATH_DIRECTORY  : PATH TO THE DIRECTORY WHERE THE FILE IS LOCATED                          **
!   ** CH_WORKING_FILE : ROOT NAME OF THE WORKING FILE                                            **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    character (len=200), intent(in) :: PATH_DIRECTORY;

    character (len=5), intent(in) :: CH_WORKING_FILE;

    real (kind=8), dimension(1:3), intent(in) :: CELL_AXIS, CELL_ANGDEG;

!   ************************************************************************************************

    integer (kind=4) :: i;

    integer (kind=4) :: iespece, insp;

    character (len=200) :: CH_DATA_FILE_NAME;

    logical :: PROBE1;

!   ************************************************************************************************

!   CH_DATA_FILE_NAME = trim(PATH_DIRECTORY)//trim(CH_WORKING_FILE)//'.xyz';

!   write(99,*) CH_DATA_FILE_NAME;

!   inquire(FILE=trim(CH_DATA_FILE_NAME),EXIST=PROBE1);

!   if (PROBE1 == .TRUE. ) then
!   if (PROBE1 .EQV. .TRUE. ) then; !
!       write(99,*) trim(CH_DATA_FILE_NAME)//' ALREADY EXIST !!!';
!       write(99,*) 'END OF PROGRAM';
!       close(99);
!       stop;
!   end if

!   open(110,file=trim(CH_DATA_FILE_NAME),status='new');
!   write(110,*) NTOT_ATOMS_FINAL;
!   write(110,'(6f20.13)') CELL_AXIS(1:3), CELL_ANGDEG(1:3);

!   do iespece = 1, espece
!       do insp = 1, NATSP_FINAL(iespece)
!           write(110,'(a3,3f17.8)') DATA_FINAL(iespece)%NAT(insp), DATA_FINAL(iespece)%RI(1:3,insp);
!       end do
!   end do

!   if ( chslab == 1 ) then
!       do i = 1, NSLAB    ! READ SUBSTRATE PARAMETERS
!           write(110,'(a3,3f17.8)') SUB(i)%NAT, SUB(i)%RI(1:3);
!       end do
!   end if

!   close(110);

end subroutine WRITE_CONFIG
