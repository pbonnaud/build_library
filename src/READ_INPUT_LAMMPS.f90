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

subroutine READ_INPUT_LAMMPS(icanal,           &
                             CHEXT)

!   ************************************************************************************************
!   **                                 READ LAMMPS FILES                                          **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal        : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                   **
!   ** CHEXT         : NAME OF THE LAMMPS CONFIGURATION                                           **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_physical_constants;

    use module_library;

    use module_config;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=150), intent(in) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: ixy, ixz, iyz;

    integer (kind=4) :: ibonds, ibond_coeffs;

    integer (kind=4) :: iangles, iangle_coeffs;

    integer (kind=4) :: iimpropers, iimproper_coeffs, iangletorsion_coeffs, idihedrals;

    integer (kind=4) :: imcon, ICONNECT, IATOM_TYPE, ATOMS_FLAG;

    integer (kind=4) :: IFOUND;

    real (kind=8) :: XLO, XHI, YLO, YHI, ZLO, ZHI;

    real (kind=8), allocatable, dimension(:) :: BOND_PCFF_VERSION;

    character (len=4), dimension(1:6,1:134) :: PCFF_EQUIVALENCE;

    character (len=150) :: CHARLINE;

    integer (kind=4) :: EOF, EOF2;

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5;

    real (kind=8) :: ROSEF1, ROSEF2, ROSEF3;

    real (kind=8), dimension(1:4) :: TAB_ROSEF; 

    character (len=250) :: CHOSEF1, CHOSEF2, CHOSEF3, CHOSEF4, CHOSEF5, CHOSEF6;

    character (len=250) :: CHAIN_LENGTH;

    logical :: PROBE1;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE, ILENGTH_CHEXT;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Read input LAMMPS';                                                       !
                                                                                         ! 
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the name of the configuration file to be considered ##################################
                                                                                         !
    ILENGTH_CHEXT = LEN_TRIM(CHEXT);                                                     !
                                                                                         !
    IOSEF1 = 67 - ILENGTH_CHEXT;                                                         !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of variables and arrays for readin the configuration ########################
                                                                                         !
    NINPUT_PAIR_POTENTIAL = 0;                                                           !
                                                                                         !
    PAIR_POTENTIAL_SPECIES(1:2,1:1000) = 0;                                              !
                                                                                         !
    INPUT_PAIR_POTENTIAL(1:2,1:1000) = 0.0d0;                                            !
                                                                                         !
    PAIR_POTENTIAL_CHTYPE(1:1000) = 'XXX';                                               !
                                                                                         !
!   ##### Start reading the configuration file #####################################################
                                                                                         !
    open(2,file=TRIM(CHEXT));                                                            !
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
                                                                                         !
        read(2,'(a)',iostat=EOF) CHAIN_LENGTH;                                           !
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHAIN_LENGTH,'atom_style ') > 0 ) then;                               !
                                                                                         !
            read(CHAIN_LENGTH,*) CHOSEF1, CH_ATOM_STYLE;                                 !
                                                                                         !
            IOSEF1 = 70 - 23 - LEN_TRIM(CH_ATOM_STYLE);                                  !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| The atom style is : '//TRIM(CH_ATOM_STYLE)// &      !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'bond_style ') > 0 ) then;                          !
                                                                                         !
            read(CHAIN_LENGTH,*) CHOSEF1, CH_BOND_STYLE;                                 !
                                                                                         !
            IOSEF1 = 70 - 23 - LEN_TRIM(CH_BOND_STYLE);                                  !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| The bond style is : '//TRIM(CH_BOND_STYLE)// &      !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'angle_style ') > 0 ) then;                         !
                                                                                         !
            read(CHAIN_LENGTH,*) CHOSEF1, CH_ANGLE_STYLE;                                !
                                                                                         !
            IOSEF1 = 70 - 24 - LEN_TRIM(CH_ANGLE_STYLE);                                 !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| The angle style is : '//TRIM(CH_ANGLE_STYLE)// &    !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'dihedral_style ') > 0 ) then;                      !
                                                                                         !
            read(CHAIN_LENGTH,*) CHOSEF1, CH_DIHEDRAL_STYLE;                             !
                                                                                         !
            IOSEF1 = 70 - 27 - LEN_TRIM(CH_DIHEDRAL_STYLE);                              !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| The dihedral style is : '// &                       !
                                  TRIM(CH_DIHEDRAL_STYLE)//      &                       !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !

         else if ( INDEX(CHAIN_LENGTH,'improper_style ') > 0 ) then;

             write(icanal,*) 'Improper style - not implemented - stop!';
   
             stop;
                                               
                                                                                         !
         else if ( INDEX(CHAIN_LENGTH,'pair_coeff ') > 0 ) then;                         !
                                                                                         ! 
            NINPUT_PAIR_POTENTIAL = NINPUT_PAIR_POTENTIAL + 1;                           !
                                                                                         !
            read(CHAIN_LENGTH,*) CHOSEF1,                                           &    !
                                 PAIR_POTENTIAL_SPECIES(1:2,NINPUT_PAIR_POTENTIAL);      !
                                                                                         !
            IOSEF1 = 1;                                                                  !
                                                                                         !
            if ( INDEX(CHAIN_LENGTH,'lj/cut/coul/long') > 0 ) then;                      !
                                                                                         !
                IOSEF1 = INDEX(CHAIN_LENGTH,'lj/cut/coul/long');                         !
                                                                                         !
                IOSEF2 = IOSEF1 + 16;                                                    !
                                                                                         !
                PAIR_POTENTIAL_CHTYPE(NINPUT_PAIR_POTENTIAL) = &                         !
                CHAIN_LENGTH(IOSEF1:IOSEF2);                                             !
                                                                                         !
            end if                                                                       !
                                                                                         !
            IOSEF2 = IOSEF2 + 1;                                                         !
                                                                                         !
            IOSEF3 = LEN_TRIM(CHAIN_LENGTH);                                             !

            CHOSEF2 = CHAIN_LENGTH(IOSEF2:IOSEF3); 

            read(CHOSEF2,*) INPUT_PAIR_POTENTIAL(1:2,NINPUT_PAIR_POTENTIAL);             !
                                                                                         !
            write(icanal,'(a2,2i4,a60)') '| ',                                              &
                                         PAIR_POTENTIAL_SPECIES(1:2,NINPUT_PAIR_POTENTIAL), &
                                         ' '//REPEAT('-',58)//'|';                       !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 70 - 3 - LEN_TRIM(PAIR_POTENTIAL_CHTYPE(NINPUT_PAIR_POTENTIAL));    !
                                                                                         !
            write(icanal,'(a70)') '| '//                                               & !
                                  TRIM(PAIR_POTENTIAL_CHTYPE(NINPUT_PAIR_POTENTIAL))// & !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              ! 
                                                                                         !
            write(icanal,'(a2,2f15.6,a38)') '| ',                                            &
                                            INPUT_PAIR_POTENTIAL(1:2,NINPUT_PAIR_POTENTIAL), &
                                            REPEAT(' ',37)//'|';                         !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(2);                                                                            !
                                                                                         !
    write(icanal,'(a2,i4,a64)') '| ',                          &                         !
                                NINPUT_PAIR_POTENTIAL,         &                         !
                                ' pair potentials were read'// &                         !
                                REPEAT(' ',37)//'|';                                     !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine READ_INPUT_LAMMPS
