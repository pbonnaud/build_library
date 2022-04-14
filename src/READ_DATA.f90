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

subroutine READ_DATA(icanal,          &
                     PATH_DIRECTORY,  &
                     CH_WORKING_FILE, &
                     IDATA)

!   ************************************************************************************************
!   **                                     READ DATA FILE                                         **
!   ************************************************************************************************
!   **                                                                                            **
!   ** PATH_DIRECTORY  : PATH TO THE DIRECTORY WHERE THE FILE IS LOCATED                          **
!   ** CH_WORKING_FILE : ROOT NAME OF THE WORKING FILE                                            **
!   **                                                                                            **
!   ** IDATA           : PRESENCE OF THE FILE IN THE WORKING DIRECTORY                            **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=200), intent(in) :: PATH_DIRECTORY;

    character (len=5), intent(in) :: CH_WORKING_FILE;

    integer (kind=4), intent(out) :: IDATA;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

    integer (kind=4) :: iespece;

    integer (kind=4) :: EOF;

    integer (kind=4) :: SAME_AT;

    character (len=3) :: CHOSEF1;

    character (len=3), dimension(1:10) :: ATOM_LABEL;

    character (len=150) :: CHARLINE;

    character (len=200) :: CH_DATA_FILE_NAME;

    logical :: PROBE1;

!   ************************************************************************************************

!   write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!   IDATA = 0;

!   CH_DATA_FILE_NAME = TRIM(PATH_DIRECTORY)//TRIM(CH_WORKING_FILE)//'_data.dat';

!   write(icanal,*) CH_DATA_FILE_NAME;

!   inquire(FILE=TRIM(CH_DATA_FILE_NAME),EXIST=PROBE1);

!   if (PROBE1 == .FALSE. ) then
!   if (PROBE1 .EQV. .FALSE. ) then;
!       write(icanal,*) '| '//TRIM(CH_DATA_FILE_NAME)//' DOESNT EXIST !!!';
!       write(icanal,*);
!       espece = 0;
!       write(icanal,*) '| SET espece : ', espece;
!       nsite(1) = 0;
!       write(icanal,*) '| SET INITIAL nsite(1) : ', nsite(1);

!       chslab = 0;
!       write(99,*) 'END OF PROGRAM';
!       close(99);
!   end if
!       write(icanal,*);
!   else
!       IDATA = 1;

!       chslab = 0;

!       open(5,file=trim(CH_DATA_FILE_NAME),status='old');
!       read(5,*) espece;
!       read(5,*);

!       write(icanal,*) '| NUMBER OF SPECIES IN THE SYSTEM : ', espece;

!       if ( espece > 0 ) then
!           do i = 1, espece
!               read(5,*) nsite(i), N(i);
!               write(icanal,*) nsite(i), N(i);
!           end do
!           read(5,*);

!           do i = 1, espece
!               do j = 1, nsite(i)

!                   read(5,'(a)') CHARLINE;

!                   EOF = 0;

!                   read(CHARLINE,*,iostat=EOF) mol(i,j)%nat, mol(i,j)%sigma, mol(i,j)%eps, mol(i,j)%RI(1:3), &
!                                               mol(i,j)%q, mol(i,j)%M;

!                   if ( EOF /= 0 ) then
!                       EOF = 0;
!                       read(CHARLINE,*,iostat=EOF) mol(i,j)%nat, mol(i,j)%CHDISPREP, mol(i,j)%MULTI, mol(i,j)%sigma, &
!                                                   mol(i,j)%eps, mol(i,j)%RI(1:3), mol(i,j)%q, mol(i,j)%M;

!                       if ( EOF /= 0 ) then
!                           write(icanal,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!';
!                           write(icanal,*) '!!! PROBLEM IN READING THE SUBSTRATE PARAMETER    !!!';
!                           write(icanal,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!';
!                           close(icanal);
!                           stop;
!                       end if
!                   else
!                       mol(i,j)%CHDISPREP = 'lj';
!                       mol(i,j)%MULTI = 1;
!                   end if
!               end do
!               read(5,*)
!           end do
!       end if

!       read(5,*) chslab;

!       if ( chslab == 1 )  then
!           read(5,*)
!           read(5,*) Nslab;
!           read(5,*)
!       end if

!       close(5);

!       ### COMPUTE MAXIMUM NUMBER OF SITES PER MOLECULES AND PER SPECIES ##############################
!       NATSP(1:10) = 0;
!       NSITE_TOT(1:10) = 0;

!       do i = 1, espece
!           do j = 1, nsite(i)
!               NSITE_TOT(i) = NSITE_TOT(i) + mol(i,j)%MULTI;
!           end do
!           NATSP(i) = NSITE_TOT(i) * N(i);
!       end do

!       ### FIND ATOMIC LABEL IN DATA FILE #############################################################
!       NLABEL = 0;

!       ATOM_LABEL(1:10) = 'XXX';
 
!       do iespece = 1, espece
!           do j = 1, nsite(iespece)

!               if ( NLABEL == 0 ) then
!                   NLABEL = NLABEL + 1;
!                   ATOM_LABEL(NLABEL) = mol(iespece,j)%nat;
!               else
!                   SAME_AT = 0;
!                   do i = 1, NLABEL
!                       if ( ATOM_LABEL(i) == mol(iespece,j)%nat ) SAME_AT = 1;
!                   end do

!                   if ( SAME_AT == 0 ) then
!                       NLABEL = NLABEL + 1;
!                       ATOM_LABEL(NLABEL) = mol(iespece,j)%nat;
!                   end if
!               end if
!           end do
!       end do

!       ### WRITE ATOMIC LABELS ########################################################################
!      write(icanal,*);
!      write(icanal,*) '| LABELS : ';
!      write(icanal,*);
!   
!       do i = 1, NLABEL
!           write(icanal,*) i, ATOM_LABEL(i);
!       end do

!   end if

!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
!   write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine READ_DATA

