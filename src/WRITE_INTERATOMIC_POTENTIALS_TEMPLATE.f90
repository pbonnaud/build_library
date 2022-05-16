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

subroutine WRITE_INTERATOMIC_POTENTIALS_TEMPLATE(icanal,CHEXT)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                  : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                         **
!   **                                                                                            **
!   ** CHEXT                   : NAME OF THE LAMMPS CONFIGURATION                                 **
!   **                                                                                            **
!   ************************************************************************************************

    use module_physical_constants;

    use module_size_arrays;

    use module_config;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=250), intent(in) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: EOF, EOF2;

    character (len=250) :: CHARLINE, CHARLINE2;

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Write interatomic potentials template';                                   !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   ! 
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !  
!   ### Write the name of the interatomic potentials templated file ################################
                                                                                         !
    IOSEF1 = 67 - LEN_TRIM(CHEXT);                                                       !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Write the interatomic potentials templated file ############################################
                                                                                         !
    open(5,file=TRIM(CHEXT));                                                            !
                                                                                         !
    write(5,'(a7)') ' Masses';                                                           !
                                                                                         !
    write(5,*);                                                                          !
                                                                                         !
    do i = 1, NTYPE_ATOM;                                                                !
                                                                                         !
        write(5,'(i8,f12.6)') i, ATOM_MASSE(i);                                          !
                                                                                         !
    end do                                                                               !
                                                                                         !
    if ( IPOTENTIAL_CLASS2 == 1 ) then;                                                  !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
        CHOSEF1 = 'Pair Coeffs # '//TRIM(POTENTIAL_CLASS2_CHTYPE);                       !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1) + 1;                                                  !
                                                                                         !
        if ( IOSEF1 < 100 ) write(CHOSEF2,'(i2)') IOSEF1;                                !
                                                                                         !
        write(5,'(a'//TRIM(CHOSEF2)//')') ' '//TRIM(CHOSEF1);                            !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
        do i = 1, NTYPE_ATOM;                                                            !
                                                                                         !
            write(5,'(i8,2f15.6)') i, POTENTIAL_CLASS2(1:2,i);                           !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(5,*);                                                                          !
                                                                                         !
!   ### Write bond coefficients ####################################################################
                                                                                         !
    if ( NTYPE_BOND > 0 ) then;                                                          !
                                                                                         !
        CHOSEF1 = 'Bond Coeffs # '//TRIM(CH_BOND_STYLE);                                 !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1) + 1;                                                  !
                                                                                         !
        if ( IOSEF1 < 100 ) write(CHOSEF2,'(i2)') IOSEF1;                                !
                                                                                         !
        write(5,'(a'//TRIM(CHOSEF2)//')') ' '//TRIM(CHOSEF1);                            !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
        write(CHOSEF3,'(i1)') NPARAM_BONDS;                                              !
                                                                                         !
        do i = 1, NTYPE_BOND;                                                            ! LOOP OVER THE NUMBER OF BOND TYPES IN THE MOLECULAR CONFIGURATION
                                                                                         !
            write(5,'(i8,'//TRIM(CHOSEF3)//'f15.6)') i, BOND_COEFFS(1:NPARAM_BONDS,i);   !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write angle coefficients ###################################################################
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        CHOSEF1 = 'Angle Coeffs # '//TRIM(CH_ANGLE_STYLE);                               !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1) + 1;                                                  !
                                                                                         !
        if ( IOSEF1 < 100 ) write(CHOSEF2,'(i2)') IOSEF1;                                !
                                                                                         !
        write(5,'(a'//TRIM(CHOSEF2)//')') ' '//TRIM(CHOSEF1);                            !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
        write(CHOSEF3,'(i1)') NPARAM_ANGLES;                                             !
                                                                                         !
        do i = 1, NTYPE_ANGLE;                                                           ! LOOP OVER THE NUMBER OF ANGLE TYPES IN THE MOLECULAR CONFIGURATION
                                                                                         !
            write(5,'(i8,'//TRIM(CHOSEF3)//'f12.6)') i,                              &   !
                                                     ANGLE_COEFFS(1:NPARAM_ANGLES,i);    !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
        if ( NPARAM_ANGLES == 4 ) then;                                                  !
                                                                                         !
            write(5,'(a16)') ' BondBond Coeffs';                                         !
                                                                                         !
            write(5,*);                                                                  !
                                                                                         !
            do i = 1, NTYPE_ANGLE;                                                       !
                                                                                         !
                write(5,'(i8,3f12.6)') i, BONDBOND_COEFFS(1:3,i);                        !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(5,*);                                                                  !
                                                                                         !
            write(5,'(a17)') ' BondAngle Coeffs';                                        !
                                                                                         !
            write(5,*);                                                                  !
                                                                                         !
            do i = 1, NTYPE_ANGLE;                                                       !
                                                                                         !
                write(5,'(i8,4f12.6)') i, BONDANGLE_COEFFS(1:4,i);                       !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(5,*);                                                                  !  
                                                                                         !
        end if                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write dihedral coefficients ################################################################
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        CHOSEF1 = 'Dihedral Coeffs # '//TRIM(CH_DIHEDRAL_STYLE);                         !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1) + 1;                                                  !
                                                                                         !
        if ( IOSEF1 < 100 ) write(CHOSEF2,'(i2)') IOSEF1;                                !
                                                                                         !
        write(5,'(a'//TRIM(CHOSEF2)//')') ' '//TRIM(CHOSEF1);                            !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
        IOSEF2 = 0;                                                                      !
                                                                                         !
        if ( TRIM(CH_DIHEDRAL_STYLE) == 'charmmfsw' ) then;                              !
                                                                                         !
            CHOSEF3 = 'i8,f12.6,2i8,f12.6';                                              !
                                                                                         !
            IOSEF2 = 1;                                                                  ! 
                                                                                         !
        else                                                                             !
                                                                                         !
            write(CHOSEF3,'(i1)') NPARAM_DIHEDRALS;                                      !
                                                                                         !
            CHOSEF3 = 'i8,'//TRIM(CHOSEF3)//'f12.6';                                     !
                                                                                         !
        end if                                                                           !
                                                                                         !
        do i = 1, NTYPE_DIHEDRAL;                                                        !
                                                                                         !
            if ( IOSEF2 == 1 ) then;                                                     !
                                                                                         !
                write(5,'('//TRIM(CHOSEF3)//')') i,                         &            !
                                                 DIHEDRAL_COEFFS(1,i),      &            !
                                                 INT(DIHEDRAL_COEFFS(2,i)), &            !
                                                 INT(DIHEDRAL_COEFFS(3,i)), &            !
                                                 DIHEDRAL_COEFFS(4,i);                   !
                                                                                         !
            else                                                                         !
                                                                                         !
                write(5,'('//TRIM(CHOSEF3)//')') i,  &                                   !
                                                 DIHEDRAL_COEFFS(1:NPARAM_DIHEDRALS,i);  !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(5,*);                                                                      !

        if ( NPARAM_DIHEDRALS == 6 ) then;

            write(5,'(a25)') ' MiddleBondTorsion Coeffs';

            write(5,*);

            do i = 1, NTYPE_DIHEDRAL;

                write(5,'(i8,4f12.6)') i, MIDDLEBONDTORSION_COEFFS(1:4,i);

            end do

            write(5,*);

            write(5,*) ' EndBondTorsion Coeffs';

            write(5,*);

            do i = 1, NTYPE_DIHEDRAL;

                write(5,'(i8,8f12.6)') i, ENDBONDTORSION_COEFFS(1:8,i);

            end do

            write(5,*);

            write(5,'(a29)') ' AngleTorsion Coeffs # class2';

            write(5,*);

            do i = 1, NTYPE_DIHEDRAL;

                write(5,'(i8,8f12.6)') i, ANGLETORSION_COEFFS(1:8,i);

            end do

            write(5,*);

            write(5,'(a25)') ' AngleAngleTorsion Coeffs';

            write(5,*);

            do i = 1, NTYPE_DIHEDRAL;

                write(5,'(i8,3f12.6)') i, ANGLEANGLETORSION_COEFFS(1:3,i);

            end do

            write(5,*);

            write(5,'(a18)') ' BondBond13 Coeffs';

            write(5,*);

            do i = 1, NTYPE_DIHEDRAL;

                write(5,'(i8,3f12.6)') i, BONDBOND13_COEFFS(1:3,i);

            end do

            write(5,*);

        end if

    end if

!   ### Write improper coefficients ################################################################
                                                                                         !
    if ( NTYPE_IMPROPER > 0 ) then;                                                      !
                                                                                         !
        write(5,'(a25)') ' Improper Coeffs # class2';                                    !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
        do i = 1, NTYPE_IMPROPER;                                                        !
                                                                                         !
            write(5,'(i8,2f12.6)') i, IMPROPER_COEFFS(1:2,i);                            !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
        write(5,'(a18)') ' AngleAngle Coeffs';                                           !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
        do i = 1, NTYPE_IMPROPER;                                                        !
                                                                                         !
            write(5,'(i8,6f12.6)') i, ANGLEANGLE_COEFFS(1:6,i);                          !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    close(5);                                                                            !
                                                                                         ! 
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine WRITE_INTERATOMIC_POTENTIALS_TEMPLATE
