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

subroutine READ_CONFIG(PATH_DIRECTORY,  &
                       CH_WORKING_FILE, &
                       IDATA,           &
                       IFIELD,          &
                       CELL_AXIS,       &
                       CELL_ANGDEG)

!   ************************************************************************************************
!   **                               READ INPUT CONFIGURATION FILE                                **
!   ************************************************************************************************
!   **                                                                                            **
!   ** PATH_DIRECTORY  : PATH TO THE DIRECTORY WHERE THE FILE IS LOCATED                          **
!   ** CH_WORKING_FILE : ROOT NAME OF THE WORKING FILE                                            **
!   ** IDATA           :  PRESENCE OF THE FILE data.dat IN THE WORKING DIRECTORY                  **
!   ** IFIELD          :  PRESENCE OF THE FILE FIELD.dat IN THE WORKING DIRECTORY                 **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: IDATA, IFIELD;

    character (len=200), intent(in) :: PATH_DIRECTORY;

    character (len=5), intent(in) :: CH_WORKING_FILE;

    real (kind=8), dimension(1:3), intent(out) :: CELL_AXIS, CELL_ANGDEG;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

    integer (kind=4) :: istart, iespece, insp, ISLAB;

    integer (kind=4) :: KFOUND;

    real (kind=8) :: SUM_NBATLIST, DELTA_CHARGE;

    real (kind=8) :: NORMRIJ, RIJ;

    real (kind=8), dimension(1:3) :: RI;

    character (len=3) :: CHAIN, CHAIN2;

    real (kind=8) :: ROSEF1;

    character (len=3) :: CHOSEF1, NAT;

    integer (kind=4) :: EOF;

    character (len=200) :: CHREADLINE;

!   ************************************************************************************************

    integer (kind=4) :: INNCOUNT, NTOT_NCT;

    integer (kind=4), allocatable, dimension(:) :: SUBINIT_LABEL;

!   ################################################################################################

    real (kind=8) :: SIGMA1, EPSI1, CHARGE1, MASSE1, MULTI1;
    character (len=3) :: LABAT1, LABAT2;

    character (len=200) :: CH_DATA_FILE_NAME;

    logical :: PROBE1;

!   ************************************************************************************************

!   CH_DATA_FILE_NAME = trim(PATH_DIRECTORY)//trim(CH_WORKING_FILE)//'.xyz';

!   write(99,*) CH_DATA_FILE_NAME;

!   inquire(FILE=trim(CH_DATA_FILE_NAME),EXIST=PROBE1);


!   if (PROBE1 == .FALSE. ) then
!   if (PROBE1 .EQV. .FALSE. ) then;
!       write(99,*) trim(CH_DATA_FILE_NAME)//' DOESNT EXIST !!!';
!       write(99,*) 'END OF PROGRAM';
!       close(99);
!   end if

!   open(11,file=trim(CH_DATA_FILE_NAME),status='old');
!   read(11,*) NTOT_ATOMS;
!   read(11,*) CELL_AXIS(1:3), CELL_ANGDEG(1:3);

!   if ( IDATA == 0 ) then
!       do i = 1, NTOT_ATOMS
!           read(11,*) NAT;
!           if ( espece  == 0 ) then
!               espece = 1;
!               nsite(espece) = 1;
!               NATSP(espece) = 1;
!               mol(espece,nsite(espece))%nat = NAT; 
!           else
!               KFOUND = 0;
!               do h = 1, espece 
!                   do j = 1, nsite(h)
!                       if ( KFOUND == 1 ) EXIT;
!                       if ( mol(h,j)%nat == NAT ) then
!                           NATSP(h) = NATSP(h) + 1;
!                           KFOUND = 1;
!                           EXIT;
!                       end if
!                   end do
!               end do

!               if ( KFOUND == 0 ) then
!                   espece = espece + 1;
!                   nsite(espece) = 1;
!                   NATSP(espece) = 1;
!                   mol(espece,nsite(espece))%nat = NAT;
!               end if
!           end if
!       end do

!       allocate(DATA(1:espece));

!       do i = 1, espece
!           allocate(DATA(i)%NAT(1:NATSP(i)));
!           allocate(DATA(i)%LABEL(1:NATSP(i)));
!           allocate(DATA(i)%INMOL(1:NATSP(i)));
!           allocate(DATA(i)%RI(1:3,1:NATSP(i)));
!           allocate(DATA(i)%RG(1:3,1:NATSP(i)));
!       end do

!       close(11);

!       open(11,file=trim(CH_DATA_FILE_NAME),status='old');
!       read(11,*);
!       read(11,*);
!   end if

!   EOF = 0;
!   NATSP(1:espece) = 0;
!   ISLAB = 0;

!   do while ( EOF == 0 ) 
!       read(11,'(a)',iostat=EOF) CHREADLINE;
!       if ( EOF == 0 ) then
!           read(CHREADLINE,*) NAT, RI(1:3);
!           KFOUND = 0;
!           do i = 1, espece
!               if ( KFOUND == 1 ) EXIT;
!               do j = 1, nsite(i)
!                   if ( mol(i,j)%nat == NAT ) then
!                       NATSP(i) = NATSP(i) + 1;
!                       DATA(i)%NAT(NATSP(i))    = NAT;
!                       DATA(i)%RI(1:3,NATSP(i)) = RI(1:3);
!                       KFOUND = 1;
!                       EXIT;
!                   end if
!               end do
!           end do

!           if ( KFOUND == 0 ) then
!               ISLAB = ISLAB + 1;
!               SUB(ISLAB)%NAT = NAT;
!               SUB(ISLAB)%RI(1:3) = RI(1:3);
!           end if
!       end if
!   end do

!   write(99,*) 'ISLAB : ', ISLAB;

!   NSLAB = ISLAB;

!   write(99,*);
!   write(99,*) 'READ CONFIGURATION : OK';
!   write(99,*);

end subroutine READ_CONFIG
