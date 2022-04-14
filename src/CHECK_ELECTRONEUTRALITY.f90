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

subroutine CHECK_ELECTRONEUTRALITY(ADD_NLABEL,ADD_LABEL,ADD_IONS)

!
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

    integer (kind=4), intent(out) :: ADD_NLABEL;

    integer (kind=4), dimension(1:10), intent(out) :: ADD_LABEL, ADD_IONS;

!    real (kind=8), intent(out) :: ADDED_IONS;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j; 

    integer (kind=4) :: JLABEL, NNEWSPEC;

    character (len=3), dimension(1:10) :: ATOM_LABEL;

    real (kind=8) :: SUM_CONFIG, SUM_CONFIG_NEW, SUM_ADDED_IONS, ADDED_IONS;

    real (kind=8) :: QNEWSPEC, QMAX, DIFF;

    real (kind=8), dimension(1:NLABEL) :: SUM_QLABEL;

!   ************************************************************************************************

!   write(99,*) '>>> CHECK ELECTRONEUTRALITY';
!   write(99,*);

!   do i = 1, NLABEL
!       write(99,'(a3,1x,f15.6)') ATOM_LABEL(i), ATOM_CHAR(i);
!   end do

!   write(99,*);

!   SUM_CONFIG = 0.0d0;

!   SUM_QLABEL(1:NLABEL) = 0.0d0;

!   do i = 1, espece
!       do j = 1, NATSP_TMP(i)
!           if ( DATA_TMP(i)%NAT(j) == 'XXX' ) CYCLE;
!           do h = 1, NLABEL
!               if ( DATA_TMP(i)%NAT(j) == ATOM_LABEL(h) ) then
!                   SUM_CONFIG = SUM_CONFIG + ATOM_CHAR(h);
!                   SUM_QLABEL(h) = SUM_QLABEL(h) + ATOM_CHAR(h);
!                   EXIT;
!               end if
!           end do
!       end do
!   end do

!   do i = 1, NSLAB_TMP
!       do h = 1, NLABEL
!           if ( SUB_TMP(i)%NAT == ATOM_LABEL(h) ) then
!               SUM_CONFIG = SUM_CONFIG + ATOM_CHAR(h);
!               SUM_QLABEL(h) = SUM_QLABEL(h) + ATOM_CHAR(h);
!               EXIT;
!           end if
!       end do
!   end do

!   do i = 1, NLABEL
!        write(99,'(a3,1x,f15.6)') ATOM_LABEL(i), SUM_QLABEL(i);
!   end do

!   write(99,*);
!   write(99,*) 'SUM_CONFIG [e] : ', SUM_CONFIG;

!   if ( ABS(SUM_CONFIG) > 1E-6 ) then 
!       write(99,*) 'COMPUTE THE NUMBER OF CALCIUM IONS TO BE ADDED : ';

!       ADDED_IONS = 0.0d0;
!       do j = 1, NLABEL
!           if ( ATOM_LABEL(j) == 'Cw' ) then
!               ADDED_IONS = ABS(SUM_CONFIG) / ABS(ATOM_CHAR(j));
!               JLABEL = j;
!               EXIT;
!           end if
!       end do
!       ADDED_IONS = ANINT( ADDED_IONS );

!       SUM_ADDED_IONS = ADDED_IONS * ATOM_CHAR(JLABEL)
 
!       SUM_CONFIG_NEW = SUM_ADDED_IONS + SUM_CONFIG;

!       write(99,*) ADDED_IONS, INT(ADDED_IONS) * ATOM_CHAR(JLABEL), SUM_CONFIG + INT(ADDED_IONS) * ATOM_CHAR(JLABEL);
!       write(99,*) ADDED_IONS, SUM_ADDED_IONS, SUM_CONFIG_NEW;

!       if ( ABS(SUM_CONFIG_NEW) > 1E-6 ) then
    
!           QMAX = 2.0d0;
!           QNEWSPEC = 5.0d0;
!           NNEWSPEC = 0;    

!           do while ( SUM_CONFIG_NEW > 1E-6 ) 
!           do while ( QNEWSPEC > QMAX ) 
!               ADDED_IONS = ADDED_IONS - 1;
!               NNEWSPEC = NNEWSPEC + 1;
!               SUM_ADDED_IONS = ADDED_IONS * ATOM_CHAR(JLABEL);
!               DIFF = -(SUM_CONFIG + SUM_ADDED_IONS);            
!               write(99,*) 'DIFF : ', DIFF;

!               QNEWSPEC = DIFF / REAL(NNEWSPEC);
!               write(99,*) 'QNEWSPEC : ', QNEWSPEC;
!           end do

!           SUM_ADDED_IONS = ADDED_IONS * ATOM_CHAR(JLABEL);

!           write(99,*) ATOM_LABEL(JLABEL), ATOM_CHAR(JLABEL), ADDED_IONS, SUM_ADDED_IONS;
!           write(99,*) 'Cw2', QNEWSPEC, NNEWSPEC, QNEWSPEC * REAL(NNEWSPEC);
            
!           SUM_CONFIG_NEW = SUM_ADDED_IONS + SUM_CONFIG + QNEWSPEC * REAL(NNEWSPEC);

!           write(99,*) SUM_CONFIG_NEW;

!           NLABEL = NLABEL + 1;
!           ATOM_LABEL(NLABEL) = 'Cw2';
!           ATOM_CHAR(NLABEL) = QNEWSPEC;

!           do i = 1, NLABEL
!               write(99,'(a3,1x,f15.6)') ATOM_LABEL(i), ATOM_CHAR(i);
!           end do

!           write(99,*);

!           ADD_NLABEL = 2;

!           ADD_LABEL(1) = JLABEL; 
!           ADD_IONS(1) = ADDED_IONS;

!           N(2) = N(2) + ADD_IONS(1);

!           ADD_LABEL(2) = NLABEL; 
!           ADD_IONS(2) = NNEWSPEC;

!           espece = espece + 1;

!           allocate(DATA_TMP(espece)%NAT(1:NNEWSPEC));
!           allocate(DATA_TMP(espece)%RI(1:3,1:NNEWSPEC));

!           DATA_TMP(espece)%NAT(1:NNEWSPEC) = 'XXX';            
!           DATA_TMP(espece)%RI(1:3,1:NNEWSPEC) = 0.0d0;

!           nsite(espece) = 1;
!           N(espece) = NNEWSPEC;

!           mol(espece,nsite(espece))%nat       = 'Cw2';
!           mol(espece,nsite(espece))%q         = QNEWSPEC;
!           mol(espece,nsite(espece))%M         = 40.078d0;
!           mol(espece,nsite(espece))%CHDISPREP = 'lj';
!           mol(espece,nsite(espece))%sigma     = 3.019d0;
!           mol(espece,nsite(espece))%eps       = 62.74d0;
!           mol(espece,nsite(espece))%RI(1:3)   = 0.0d0;
!           mol(espece,nsite(espece))%MULTI     = 1;
!       end if

!   end if 

!   stop;

end subroutine CHECK_ELECTRONEUTRALITY

