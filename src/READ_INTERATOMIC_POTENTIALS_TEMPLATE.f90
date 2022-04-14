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

subroutine READ_INTERATOMIC_POTENTIALS_TEMPLATE(icanal,CHEXT,                    &
                                                       NTYPE_TEMPLATE,           &
                                                       NTYPE_BOND_TPLTE,         &
                                                       NTYPE_ANGLE_TPLTE,        &
                                                       NTYPE_IMPROPER_TPLTE,     &
                                                       TEMPLATE_MASSES,          &
                                                       ATOM_LABEL_TPLTE,         &
                                                       POTENTIAL_CLASS2_TPLTE,   &
                                                       BOND_COEFFS_TPLTE,        &
                                                       ANGLE_COEFFS_TPLTE,       &
                                                       BONDBOND_COEFFS_TPLTE,    &
                                                       BONDANGLE_COEFFS_TPLTE,   &
                                                       IMPROPER_COEFFS_TPLTE,    &
                                                       ANGLEANGLE_COEFFS_TPLTE,  &
                                                       IPOTENTIAL_CLASS2_TPLTE)

!NTYPE_DIHEDRAL_TPLTE,)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                  : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                         **
!   ** CHEXT                   : NAME OF THE LAMMPS CONFIGURATION                                 **
!   ** NTYPE_TEMPLATE          : NUMBER OF ATOM TYPE IN THE TEMPLATED MOLECULE                    **
!   ** NTYPE_BOND_TPLTE        : NUMBER OF BOND TYPES IN THE TEMPLATED MOLECULE                   **
!   ** NTYPE_ANGLE_TPLTE       : NUMBER OF ANGLE TYPES IN THE TEMPLATED MOLECULE                  ** 
!   ** NTYPE_IMPROPER_TPLTE    : NUMBER OF IMPROPER TYPES IN THE TEMPLATED MOLECULE                  **
!   **                                                                                            **
!   ** TEMPLATE_MASSES         : MASS OF ATOMS IN THE TEMPLATED MOLECULE                                  **
!   ** ATOM_LABEL_TPLTE        : LABEL OF ATOMS IN THE TEMPLATED MOLECULE                                 **
!   ** POTENTIAL_CLASS2_TPLTE  : PAIR POTENTIAL PARAMETERS FOR INTERMOLECULAR INTERACTIONS OF   **  
!   **                           THE ATOMS BELONGING TO THE TEMPLATED MOLECULE                  **
!   ** BOND_COEFFS_TPLTE       : BOND COEFFICIENTS OF THE TEMPLATED MOLECULE                         **    
!   ** IPOTENTIAL_CLASS2_TPLTE : FLAG TO TELL THE PROGRAM THAT THERE IS CLASS2 POTENTIAL 
!   **                           PARAMETERS IN THE TEMPLATED MOLECULE
!   ** ANGLE_COEFFS_TPLTE      : ANGLE COEFFICIENTS OF THE TEMPLATED MOLECULE                        **
!   ** BONDBOND_COEFFS_TPLTE   :                                                                    **
!   ** BONDANGLE_COEFFS_TPLTE  :                                                                  **
!   ** IMPROPER_COEFFS_TPLTE   : IMPROPER COEFFICIENTS OF THE TEMPLATED MOLECULE                    **
!   ** ANGLEANGLE_COEFFS_TPLTE :


!   ** NTYPE_DIHEDRAL_TPLTE : NUMBER OF DIHEDRAL TYPES IN THE TEMPLATED MOLECULE                  **



!   **                                                                                            **
!   ************************************************************************************************

!   use data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: NTYPE_TEMPLATE, NTYPE_BOND_TPLTE, NTYPE_ANGLE_TPLTE;

    integer (kind=4), intent(in) :: NTYPE_IMPROPER_TPLTE;

    character (len=250), intent(in) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4), intent(out) :: IPOTENTIAL_CLASS2_TPLTE;

!   integer (kind=4), intent(out) :: NTYPE_BOND_TPLTE, NTYPE_ANGLE_TPLTE;

!   integer (kind=4), intent(out) :: NTYPE_DIHEDRAL_TPLTE, NTYPE_IMPROPER_TPLTE;

    real (kind=8), dimension(1:20), intent(out) :: TEMPLATE_MASSES;

    character (len=20), dimension(1:20), intent(out) :: ATOM_LABEL_TPLTE;

    real (kind=8), dimension(1:2,1:100), intent(out) :: POTENTIAL_CLASS2_TPLTE;

    real (kind=8), dimension(1:4,1:10), intent(out) :: BOND_COEFFS_TPLTE, ANGLE_COEFFS_TPLTE;

    real (kind=8), dimension(1:4,1:10), intent(out) :: IMPROPER_COEFFS_TPLTE, BONDANGLE_COEFFS_TPLTE;
!   real (kind=8), dimension(1:4,1:10), intent(out) :: BONDANGLE_COEFFS_TPLTE;

    real (kind=8), dimension(1:3,1:10), intent(out) :: BONDBOND_COEFFS_TPLTE;

    real (kind=8), dimension(1:6,1:10), intent(out) :: ANGLEANGLE_COEFFS_TPLTE;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: EOF, EOF2;

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4;

    real (kind=8) :: ROSEF1, ROSEF2;

    character (len=15) :: CHOSEF1, CHOSEF2, CHOSEF3;

    character (len=250) :: CHARLINE, CHARLINE2;

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************

    CHTITLE = 'READ INTERATOMIC POTENTIALS TEMPLATE';

    ILENGTH_TITLE = LEN_TRIM(CHTITLE);

    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));

!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    IOSEF1 = 67 - LEN_TRIM(CHEXT);

    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF1)//'|';
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!   NTYPE_BOND_TPLTE     = 0;
!   NTYPE_ANGLE_TPLTE    = 0; 
!   NTYPE_DIHEDRAL_TPLTE = 0; 
!   NTYPE_IMPROPER_TPLTE = 0; 

    IMPROPER_COEFFS_TPLTE(1:4,1:10) = 0.0d0;
    BOND_COEFFS_TPLTE(1:4,1:10)     = 0.0d0;
    ANGLE_COEFFS_TPLTE(1:4,1:10)    = 0.0d0;

    IPOTENTIAL_CLASS2_TPLTE = 0;

    ATOM_LABEL_TPLTE(1:20) = 'XXX';

!   open(5,file='000_'//TRIM(CHEXT));
    open(5,file=TRIM(CHEXT));

    EOF = 0;

    do
        read(5,'(a)',iostat=EOF) CHARLINE;

        if ( EOF /= 0 ) EXIT;

        if ( INDEX(CHARLINE,'Masses') > 0 ) then
            read(5,*);
            do i = 1, NTYPE_TEMPLATE;
                read(5,*) IOSEF1, TEMPLATE_MASSES(IOSEF1);
            end do

            do i = 1, NTYPE_TEMPLATE
                if ( ABS( TEMPLATE_MASSES(i) -  1.00794d0 ) < 0.1d0 ) ATOM_LABEL_TPLTE(i) = 'H';
                if ( ABS( TEMPLATE_MASSES(i) - 12.01070d0 ) < 0.1d0 ) ATOM_LABEL_TPLTE(i) = 'C';
                if ( ABS( TEMPLATE_MASSES(i) - 14.00670d0 ) < 0.1d0 ) ATOM_LABEL_TPLTE(i) = 'N';
                if ( ABS( TEMPLATE_MASSES(i) - 15.99940d0 ) < 0.1d0 ) ATOM_LABEL_TPLTE(i) = 'O';
                if ( ABS( TEMPLATE_MASSES(i) - 28.08550d0 ) < 0.1d0 ) ATOM_LABEL_TPLTE(i) = 'Si';
                if ( ABS( TEMPLATE_MASSES(i) - 32.06400d0 ) < 0.1d0 ) ATOM_LABEL_TPLTE(i) = 'S';
                if ( ABS( TEMPLATE_MASSES(i) - 40.07800d0 ) < 0.1d0 ) ATOM_LABEL_TPLTE(i) = 'Ca';
                if ( ABS( TEMPLATE_MASSES(i) - 55.84500d0 ) < 0.1d0 ) ATOM_LABEL_TPLTE(i) = 'Fe';
            end do

            write(icanal,'(a70)') '| TEMPLATE MASSES WERE READ'//REPEAT(' ',42)//'|';
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHARLINE,'Pair Coeffs # lj/class2/coul/long') > 0 ) then
            IPOTENTIAL_CLASS2_TPLTE = 1;

            read(5,*);
            do i = 1, NTYPE_TEMPLATE;
                read(5,*) IOSEF1, POTENTIAL_CLASS2_TPLTE(1:2,IOSEF1);
            end do

            write(icanal,'(a70)') '| PAIR COEFFICIENTS WERE READ'//REPEAT(' ',40) //'|';
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHARLINE,' Bond Coeffs # class2') > 0 ) then
            EOF2 = 0;

            read(5,*);
            do i = 1, NTYPE_BOND_TPLTE;                                 ! LOOP OVER THE NUMBER OF BOND TYPES IN THE MOLECULAR CONFIGURATION
                read(5,*,iostat=EOF2) IOSEF1, BOND_COEFFS_TPLTE(1:4,IOSEF1);
                if ( EOF2 /= 0 ) EXIT;
            end do

            write(icanal,'(a70)') '| BOND COEFFICIENTS WERE READ';
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHARLINE,' Angle Coeffs # class2') > 0 ) then
            EOF2 = 0;

            read(5,*);
            do i = 1, NTYPE_ANGLE_TPLTE;                          ! LOOP OVER THE NUMBER OF ANGLE TYPES IN THE MOLECULAR CONFIGURATION
                read(5,*,iostat=EOF2) IOSEF1, ANGLE_COEFFS_TPLTE(1:4,IOSEF1);
                if ( EOF2 /= 0 ) EXIT;
            end do

            write(icanal,'(a70)') '| ANGLE COEFFICIENTS WERE READ';
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHARLINE,' BondBond Coeffs') > 0 ) then
            read(5,*);
            do i = 1, NTYPE_ANGLE_TPLTE
                read(5,*) IOSEF1, BONDBOND_COEFFS_TPLTE(1:3,IOSEF1);
            end do

            write(icanal,'(a70)') '| BONDBOND COEFFICIENTS WERE READ'//REPEAT(' ',36)//'|';
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHARLINE,' BondAngle Coeffs') > 0 ) then
            read(5,*);
            do i = 1, NTYPE_ANGLE_TPLTE
                read(5,*) IOSEF1, BONDANGLE_COEFFS_TPLTE(1:4,IOSEF1);
            end do

            write(icanal,'(a70)') '| BONDANGLE COEFFICIENTS WERE READ'//REPEAT(' ',35)//'|';
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
 
        else if ( INDEX(CHARLINE,' Improper Coeffs') > 0 ) then
            EOF2 = 0;

            read(5,*);
            do i = 1, NTYPE_IMPROPER_TPLTE;
                read(5,*,iostat=EOF2) IOSEF1, IMPROPER_COEFFS_TPLTE(1:2,IOSEF1);
                if ( EOF2 /= 0 ) EXIT;
            end do
!           do
!               read(5,'(a)',iostat=EOF2) CHARLINE2;

!               if ( EOF2 /= 0 ) then
!                   EXIT;
!               else
!                   read(CHARLINE2,*,iostat=EOF2) IOSEF1, IMPROPER_COEFFS_TPLTE(1:2,IOSEF1);
!                   if ( EOF2 /= 0 ) EXIT;

!                   NTYPE_IMPROPER_TPLTE = NTYPE_IMPROPER_TPLTE + 1;
!               end if
!           end do

!           write(icanal,'(a25,i4,a41)') '| NTYPE_IMPROPER_TPLTE : ', NTYPE_IMPROPER_TPLTE, REPEAT(' ',40)//'|';
            write(icanal,'(a70)') '| IMPROPER COEFFICIENTS WERE READ'; 
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHARLINE,' AngleAngle Coeffs') > 0 ) then
            read(5,*);
            do i = 1, NTYPE_IMPROPER_TPLTE
                read(5,*) IOSEF1, ANGLEANGLE_COEFFS_TPLTE(1:6,IOSEF1);
            end do

            write(icanal,'(a70)') '| ANGLEANGLE COEFFICIENTS WERE READ'//REPEAT(' ',34)//'|';
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
        end if
    end do

    close(5);

    do i = 1, NTYPE_TEMPLATE
        write(icanal,'(a2,i8,f12.6,a3)') '| ', i, TEMPLATE_MASSES(i), TRIM(ATOM_LABEL_TPLTE(i));
    end do

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine READ_INTERATOMIC_POTENTIALS_TEMPLATE
