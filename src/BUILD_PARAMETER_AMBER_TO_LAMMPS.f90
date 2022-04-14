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

subroutine BUILD_PARAMETER_AMBER_TO_LAMMPS(icanal)

!   ************************************************************************************************
!   **                                 READ LAMMPS FILES                                          **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                          **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_physical_constants;

    use module_size_arrays;

    use module_config;

    use module_library;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=150) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: iindex, AMBER_ICO;

    integer (kind=4) :: IBOND, IANGLE, IDIHEDRAL;

    integer (kind=4) :: ibonds, ibond_coeffs;

    integer (kind=4) :: iangles, iangle_coeffs;

    integer (kind=4) :: iimpropers, iimproper_coeffs, iangletorsion_coeffs, idihedrals;

    integer (kind=4) :: IFOUND;

    integer (kind=4) :: IMAX_ATOM_TYPE;

    integer (kind=4), allocatable, dimension(:) :: TMP_ICO;

    real (kind=8) :: SUM_CHARGES, DELTA_CHARGE, CHARGE_CORRECTION;

    real (kind=8) :: LAMMPS_EPSILONIJ, LAMMPS_SIGMAIJ, AMBER_AIJ, AMBER_BIJ;

    integer (kind=4) :: EOF, EOF2;

!   ************************************************************************************************

!   integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5;

!   real (kind=8) :: ROSEF1, ROSEF2, ROSEF3;

    real (kind=8), dimension(1:4) :: TAB_ROSEF; 

!   character (len=10) :: CHOSEF1, CHOSEF2, CHOSEF3, CHOSEF4, CHOSEF5, CHOSEF6;

    character (len=250) :: CHAIN_LENGTH;

    logical :: PROBE1;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE, ILENGTH_CHEXT;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Build parameters (from AMBER to LAMMPS)';                                 !
                                                                                         ! 
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the total number of atoms in the configuration  ######################################
                                                                                         !
    write(icanal,'(a10,i8,a52)') '| NATOM : ', NATOM, REPEAT(' ',51)//'|';               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Build and write the total number of bonds in the configuration #############################
                                                                                         ! 
    NBOND = AMBER_NBONH + AMBER_NBONA;                                                   !
                                                                                         !
    write(icanal,'(a10,i8,a52)') '| NBOND : ', NBOND, REPEAT(' ',51)//'|';               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Build and write the total number of angles in the configuration ############################
                                                                                         !
    NANGLE = AMBER_NTHETH + AMBER_NTHETA;                                                !
                                                                                         !
    write(icanal,'(a11,i8,a51)') '| NANGLE : ',       &                                  !
                                 NANGLE,              &                                  !
                                 REPEAT(' ',50)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Build and write the total number of dihedrals in the configuration #########################
                                                                                         !
    NDIHEDRAL = AMBER_NPHIH + AMBER_NPHIA;                                               !
                                                                                         !
    write(icanal,'(a14,i8,a48)') '| NDIHEDRAL : ',    &                                  !
                                 NDIHEDRAL,           &                                  !
                                 REPEAT(' ',47)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the total number of impropers in the configuration to zero #############################
                                                                                         !
    NIMPROPER = 0;                                                                       !
                                                                                         !
!   ### Write the total number of atom types in the molecular configuration ########################
                                                                                         !
    write(icanal,'(a15,i8,a47)') '| NTYPE_ATOM : ',   &                                  !
                                 NTYPE_ATOM,          &                                  !
                                 REPEAT(' ',46)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the total number of bond types in the molecular configuration ########################
                                                                                         !
    write(icanal,'(a15,i8,a47)') '| NTYPE_BOND : ',   &                                  !
                                 NTYPE_BOND,          &                                  !
                                 REPEAT(' ',46)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the total number of angle types in the molecular configuration #######################
                                                                                         !
    write(icanal,'(a16,i8,a46)') '| NTYPE_ANGLE : ',    &                                !
                                  NTYPE_ANGLE,          &                                !
                                  REPEAT(' ',45)//'|';                                   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the total number of dihedral types in the molecular configuration ####################
                                                                                         !
    write(icanal,'(a19,i8,a43)') '| NTYPE_DIHEDRAL : ', &                                !
                                 NTYPE_DIHEDRAL,        &                                !
                                 REPEAT(' ',42)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the total number of improper types to zero #############################################
                                                                                         !
    NTYPE_IMPROPER = 0;                                                                  !
                                                                                         !
!   ### Conversion of partial charges (AMBER units to electron charge units) #######################
                                                                                         !
    CONFIG_QI(1:NATOM) = CONFIG_QI(1:NATOM) / 18.2223d0;                                 ! IN [e]
                                                                                         !
    SUM_CHARGES = SUM(CONFIG_QI(1:NATOM));                                               ! IN [e]
                                                                                         !
    write(icanal,'(a20,f12.4,a38)') '| SUM_CHARGES [e] : ', &                            !
                                    SUM_CHARGES,            &                            !
                                    REPEAT(' ',37)//'|';                                 !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Readjust charges ###########################################################################
                                                                                         !
    if ( iadjust_charges == 1 ) then;                                                    !
                                                                                         !
        DELTA_CHARGE = TOTAL_NET_CHARGE - SUM_CHARGES;                                   !
                                                                                         !
        write(icanal,'(a21,f20.6,a29)') '| DELTA_CHARGE [e] : ', &                       !
                                        DELTA_CHARGE,            &                       !
                                        REPEAT(' ',28)//'|';                             !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        CHARGE_CORRECTION = DELTA_CHARGE / REAL( NATOM );                                !
                                                                                         ! 
        write(icanal,'(a26,f20.6,a24)') '| CHARGE_CORRECTION [e] : ', &                  !
                                        CHARGE_CORRECTION,            &                  !
                                        REPEAT(' ',23)//'|';                             !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        CONFIG_QI(1:NATOM) = CONFIG_QI(1:NATOM) + CHARGE_CORRECTION;                     ! IN [e]
                                                                                         !
        SUM_CHARGES = SUM(CONFIG_QI(1:NATOM));                                           ! IN [e]
                                                                                         !
        write(icanal,'(a25,f12.4,a33)') '| (new) SUM_CHARGES [e] : ', &                  !
                                    SUM_CHARGES,                      &                  !
                                    REPEAT(' ',32)//'|';                                 !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set readable atom names ####################################################################
                                                                                         !
    do i = 1, NATOM;                                                                     !
                                                                                         !
        do j = 1, 200;                                                                   !
                                                                                         !
            ROSEF1 = ABS( CONFIG_MASS(i) - MOLAR_MASS_ELEMENT(j) );                      !
                                                                                         !
            if ( ROSEF1 < 0.1d0 ) CONFIG_NAT(i) = ELEMENT_SYMBOL_NAME(j);                !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               ! 
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set atom IDs ###############################################################################
                                                                                         !
    do i = 1, NATOM;                                                                     !
                                                                                         !
        CONFIG_ATOMID(i) = i;                                                            !
                                                                                         !
    end do                                                                               ! 
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set molecule ID ############################################################################
                                                                                         !
    CONFIG_MOLECULEID(1:NATOM) = 1;                                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Check atom types ###########################################################################
                                                                                         !
    IMAX_ATOM_TYPE = MAXVAL(CONFIG_ATOM_TYPE(1:NATOM));                                  !
                                                                                         !
    write(icanal,'(a32,2i8,a22)') '| IMAX_ATOM_TYPE | NTYPE_ATOM : ', &                  !
                                  IMAX_ATOM_TYPE,                     &                  !
                                  NTYPE_ATOM,                         &                  !
                                  REPEAT(' ',21)//'|';                                   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate array containing the list of atom labels ##########################################
                                                                                         !
    allocate(ATOM_LABEL(1:NTYPE_ATOM));                                                  !
                                                                                         !
!   ### Initialization of the array containing the list of atom labels #############################
                                                                                         !
    ATOM_LABEL(1:NTYPE_ATOM) = 'XXX';                                                    !
                                                                                         !
!   ### Set the list of atom labels ################################################################
                                                                                         !
    do i = 1, NTYPE_ATOM;                                                                !
                                                                                         ! 
        do j = 1, NATOM;                                                                 !
                                                                                         !
            if ( CONFIG_ATOM_TYPE(j) == i ) then;                                        !
                                                                                         !
                ATOM_LABEL(i) = CONFIG_NAT(j);                                           !
                                                                                         !
                EXIT;                                                                    !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
        IOSEF1 = 70 - 7 - 4 - 4 - 1 - LEN_TRIM(ATOM_LABEL(i));                           !
                                                                                         !
        write(icanal,'(a7,i4,a59)') '| Type ', i, ' is '//TRIM(ATOM_LABEL(i))// &        !
                                    REPEAT(' ',IOSEF1)//'|';                             !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate array containing the list of atom masses ##########################################
                                                                                         !
    allocate(ATOM_MASSE(1:NTYPE_ATOM));                                                  !
                                                                                         !
!   ### Initialization of the array containing the list of atom masses #############################
                                                                                         !
    ATOM_MASSE(1:NTYPE_ATOM) = 0.0d0;                                                    !
                                                                                         !
!   ### Set the list of atom masses ################################################################
                                                                                         !
    do i = 1, NTYPE_ATOM;                                                                !
                                                                                         ! 
        do j = 1, NATOM;                                                                 !
                                                                                         !
            if ( CONFIG_ATOM_TYPE(j) == i ) then;                                        !
                                                                                         !
                ATOM_MASSE(i) = CONFIG_MASS(j);                                          !
                                                                                         !
                EXIT;                                                                    !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(icanal,'(a21,i4,a4,f12.4,a29)') '| Molar mass of type ', i, &              !
                                              ' is ', ATOM_MASSE(i),      &              !
                                              REPEAT(' ',28)//'|';                       !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Deallocate array containing masses for amber ###############################################
                                                                                         !
    deallocate(CONFIG_MASS);                                                             ! 
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays containing nonbonded parameters for lammps #################################
                                                                                         !
    allocate(POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM));                                        !
                                                                                         !
    allocate(POTENTIAL_CLASS2_CROSS_VALUES(1:2,1:AMBER_NLJ_INTER));                      !
                                                                                         !
    allocate(POTENTIAL_CLASS2_CROSS_ATOMID(1:2,1:AMBER_NLJ_INTER));                      !
                                                                                         !
    allocate(TMP_ICO(1:AMBER_NLJ_INTER));                                                !
                                                                                         !
!   ### Initialization of arrays containing nonbonded parameters for lammps ########################
                                                                                         !
    POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM) = 0.0d0;                                          !
                                                                                         !
    POTENTIAL_CLASS2_CROSS_VALUES(1:2,1:AMBER_NLJ_INTER) = 0.0d0;                        !
                                                                                         !
    POTENTIAL_CLASS2_CROSS_ATOMID(1:2,1:AMBER_NLJ_INTER) = 0;                            !
                                                                                         !
    TMP_ICO(1:AMBER_NLJ_INTER) = 0;                                                      !
                                                                                         !
    POTENTIAL_CLASS2_CHTYPE = 'XXXX';                                                    !
                                                                                         !
    NPOTENTIAL_CLASS2_CROSS = 0;                                                         !
                                                                                         !
!   ### Set nonbonded parameters (Lennard-Jones 12-6) ##############################################
                                                                                         !
    IPOTENTIAL_CLASS2 = 1;                                                               ! 
                                                                                         !
    POTENTIAL_CLASS2_CHTYPE = 'lj/charmmfsw/coul/long';                                  !
                                                                                         !
    do i = 1, NTYPE_ATOM;                                                                !
                                                                                         !
        do j = 1, NTYPE_ATOM;                                                            !
                                                                                         !
            iindex = NTYPE_ATOM * ( i - 1 ) + j;                                         !
                                                                                         !
            AMBER_ICO = NONBONDED_PARM_INDEX(iindex);                                    !
                                                                                         !
            if ( NONBONDED_PARM_INDEX(iindex) < 0 ) then;                                !
                                                                                         !
                write(icanal,'(a70)') '| /!\ Warning /!\'//REPEAT(' ',52)//'|';          !  
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a70)') '| A Lennard-Jones 12-10 '// &                     !
                                      'interaction was found'//    &                     !
                                      REPEAT(' ',24)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a11,4i6,a35)') '| iindex : ', iindex, i, j, &             !
                                              AMBER_ICO, REPEAT(' ',34)//'|';            !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            else                                                                         !
                                                                                         !
                AMBER_AIJ = LENNARD_JONES_ACOEF(AMBER_ICO);                              ! IN [kcal/mol.A^12]
                                                                                         !
                AMBER_BIJ = LENNARD_JONES_BCOEF(AMBER_ICO);                              ! IN [kcal/mol.A^6]
                                                                                         !
                LAMMPS_EPSILONIJ = 0.0d0;                                                !
                                                                                         !
                LAMMPS_SIGMAIJ = 0.0d0;                                                  !
                                                                                         !
                if ( ABS(AMBER_AIJ) > 0 ) then;                                          !
                                                                                         ! 
                    LAMMPS_EPSILONIJ = 0.25d0 * AMBER_BIJ * AMBER_BIJ / AMBER_AIJ;       ! IN [kcal/mol]
                                                                                         !
                end if                                                                   !
                                                                                         !
                if ( ABS(AMBER_BIJ) > 0 ) then;                                          !
                                                                                         !
                    ROSEF1 = 1.0d0 / 6.0d0;                                              ! 
                                                                                         !
                    LAMMPS_SIGMAIJ = ( AMBER_AIJ / AMBER_BIJ )**ROSEF1;                  ! IN [A]
                                                                                         !
                end if                                                                   !
                                                                                         !
                if ( i == j ) then;                                                      !
                                                                                         !
                    POTENTIAL_CLASS2(1,i) = LAMMPS_EPSILONIJ;                            ! IN [kcal/mol]
                                                                                         !
                    POTENTIAL_CLASS2(2,i) = LAMMPS_SIGMAIJ;                              ! IN [A]
                                                                                         !
                else                                                                     !
                                                                                         !
                    IFOUND = 0;                                                          !
                                                                                         !
                    if ( NPOTENTIAL_CLASS2_CROSS > 0 ) then;                             !
                                                                                         !
!                       IFOUND = FINDLOC(TMP_ICO(1:NPOTENTIAL_CLASS2_CROSS),AMBER_ICO);  !
                                                                                         !
                        do k = 1, NPOTENTIAL_CLASS2_CROSS;                               !
                                                                                         !
                            if ( TMP_ICO(k) == AMBER_ICO ) then;                         !
                                                                                         !
                                IFOUND = 1;                                              !
                                                                                         !
                                EXIT;                                                    !
                                                                                         !
                            end if                                                       !
                                                                                         !
                        end do                                                           !
                                                                                         !
                    end if                                                               !
                                                                                         !
                    if ( IFOUND > 0 ) CYCLE;                                             !
                                                                                         !
                    NPOTENTIAL_CLASS2_CROSS = NPOTENTIAL_CLASS2_CROSS + 1;               !
                                                                                         !
                    IOSEF1 = NPOTENTIAL_CLASS2_CROSS;                                    !
                                                                                         !
                    TMP_ICO(IOSEF1) = AMBER_ICO;                                         !
                                                                                         !
                    POTENTIAL_CLASS2_CROSS_VALUES(1,IOSEF1) = LAMMPS_EPSILONIJ;          ! IN [kcal/mol]
                                                                                         !
                    POTENTIAL_CLASS2_CROSS_VALUES(2,IOSEF1) = LAMMPS_SIGMAIJ;            ! IN [A]
                                                                                         !
                    POTENTIAL_CLASS2_CROSS_ATOMID(1,IOSEF1) = i;                         !
                                                                                         !
                    POTENTIAL_CLASS2_CROSS_ATOMID(2,IOSEF1) = j;                         !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Nonbonded interactions were converted '// &                 !
                          'in the lammps format'//REPEAT(' ',9)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Deallocate arrays containing nonbonded parameters for amber ################################
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        deallocate(LENNARD_JONES_ACOEF);                                                 !
                                                                                         !
        deallocate(LENNARD_JONES_BCOEF);                                                 !
                                                                                         !
        deallocate(TMP_ICO);                                                             !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays containing bond coefficients for lammps ####################################
                                                                                         !
    allocate(BOND_COEFFS(1:4,1:NTYPE_BOND));                                             !
                                                                                         !
!   ### Initialization of arrays containing bond coefficients for lammps ###########################
                                                                                         !
    BOND_COEFFS(1:4,1:NTYPE_BOND) = 0.0d0;                                               !
                                                                                         !
!   ### Set bond parameters ########################################################################
                                                                                         !
    ibond_coeffs = 1;                                                                    !
                                                                                         !
    CH_BOND_STYLE = 'harmonic';                                                          !
                                                                                         !
    NPARAM_BONDS = 2;                                                                    !
                                                                                         !
    do i = 1, NTYPE_BOND;                                                                !
                                                                                         !
        BOND_COEFFS(1,i) = BOND_FORCE_CONSTANT(i);                                       ! IN [kcal/mol.A^-2] 
                                                                                         !
        BOND_COEFFS(2,i) = BOND_EQUIL_VALUE(i);                                          !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Deallocate arrays containing bond coefficients for amber ###################################
                                                                                         !
    deallocate(BOND_FORCE_CONSTANT);                                                     !
                                                                                         !
    deallocate(BOND_EQUIL_VALUE);                                                        !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays containing angle coefficients for lammps ###################################
                                                                                         !
    allocate(ANGLE_COEFFS(1:4,1:NTYPE_ANGLE));                                           !
                                                                                         !
!   ### Initialization of arrays containing angle coefficients for lammps ##########################
                                                                                         !
    ANGLE_COEFFS(1:4,1:NTYPE_ANGLE) = 0.0d0;                                             !
                                                                                         !
!   ### Set angle parameters for lammps ############################################################
                                                                                         !
    iangle_coeffs = 1;                                                                   !
                                                                                         !
    CH_ANGLE_STYLE = 'harmonic';                                                         !
                                                                                         !
    NPARAM_ANGLES = 2;                                                                   !
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         ! 
                                                                                         ! 
        do i = 1, NTYPE_ANGLE;                                                           !
                                                                                         !
            ANGLE_COEFFS(1,i) = ANGLE_FORCE_CONSTANT(i);                                 ! IN [kcal/mol.rad^-2]  
                                                                                         !
            ANGLE_COEFFS(2,i) = ANGLE_EQUIL_VALUE(i) * 360.d0 / TWOPI;                   ! IN [rad.*deg./rad.] = [deg.]
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Deallocate arrays containing angle coefficients for amber ##################################
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        deallocate(ANGLE_FORCE_CONSTANT);                                                !
                                                                                         !
        deallocate(ANGLE_EQUIL_VALUE);                                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays containing dihedral coefficients for lammps ################################
                                                                                         !
    allocate(DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL));                                     !
                                                                                         !
!   ### Initialization of arrays containing dihedral coefficients for lammps #######################
                                                                                         !
    DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL) = 0.0d0;                                       !
                                                                                         !
!   ### Set dihedral parameters for lammps #########################################################
                                                                                         !
    CH_DIHEDRAL_STYLE = 'charmmfsw';                                                     ! 
                                                                                         !
    NPARAM_DIHEDRALS = 4;                                                                !
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        do i = 1, NTYPE_DIHEDRAL;                                                        !
                                                                                         !
            DIHEDRAL_COEFFS(1,i) = DIHEDRAL_FORCE_CONSTANT(i) * 0.5d0;                   ! IN [kcal/mol]
                                                                                         !
            DIHEDRAL_COEFFS(2,i) = DIHEDRAL_PERIODICITY(i);                              ! IN [1]
                                                                                         !
            DIHEDRAL_COEFFS(3,i) = DIHEDRAL_PHASE(i) * 360.d0 / TWOPI;                   ! IN [rad.*deg./rad.] = [deg.]
                                                                                         !
            DIHEDRAL_COEFFS(4,i) = 0.0d0;                                                ! Scaling factor set to 0 for amber ff (see dihedral_style charmm manual page of lammps)
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Deallocate arrays containing dihedral coefficients for amber ###############################
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0) then;                                                       !
                                                                                         !
        deallocate(DIHEDRAL_FORCE_CONSTANT);                                             !
                                                                                         !
        deallocate(DIHEDRAL_PERIODICITY);                                                !
                                                                                         !
        deallocate(DIHEDRAL_PHASE);                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set improper parameters for lammps #########################################################
                                                                                         !
    iimproper_coeffs = 0;                                                                !
                                                                                         !
    NTYPE_IMPROPER = 0;                                                                  !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays containing bond properties for lammps ######################################
                                                                                         !
    allocate(BOND_ATOMID(1:2,1:NBOND));                                                  !
                                                                                         !
    allocate(BOND_TYPE(1:NBOND));                                                        !
                                                                                         !
!   ### Initialization of arrays containing bond properties for lammps #############################
                                                                                         !
    BOND_ATOMID(1:2,1:NBOND) = 0;                                                        !
                                                                                         !
    BOND_TYPE(1:NBOND) = 0;                                                              !
                                                                                         !
!   ### Set bond properties for lammps #############################################################
                                                                                         !
    ibonds = 1;                                                                          !
                                                                                         !
    IBOND = 0;                                                                           !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, AMBER_NBONH;                                                               !
                                                                                         !
        IBOND = IBOND + 1;                                                               !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF2 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF3 = IOSEF2 + 1;                                                             !
                                                                                         !
        BOND_ATOMID(1,IBOND) = BONDS_INC_HYDROGEN(IOSEF1) / 3 + 1;                       !
                                                                                         !
        BOND_ATOMID(2,IBOND) = BONDS_INC_HYDROGEN(IOSEF2) / 3 + 1;                       !
                                                                                         !
        BOND_TYPE(IBOND) = BONDS_INC_HYDROGEN(IOSEF3);                                   !
                                                                                         !
        IOSEF1 = IOSEF3;                                                                 !
                                                                                         !
    end do                                                                               !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         ! 
    do i = 1, AMBER_NBONA;                                                               !
                                                                                         !
        IBOND = IBOND + 1;                                                               !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF2 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF3 = IOSEF2 + 1;                                                             !
                                                                                         !
        BOND_ATOMID(1,IBOND) = BONDS_WITHOUT_HYDROGEN(IOSEF1) / 3 + 1;                   !
                                                                                         !
        BOND_ATOMID(2,IBOND) = BONDS_WITHOUT_HYDROGEN(IOSEF2) / 3 + 1;                   !
                                                                                         !
        BOND_TYPE(IBOND) = BONDS_WITHOUT_HYDROGEN(IOSEF3);                               !
                                                                                         !
        IOSEF1 = IOSEF3;                                                                 !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Bonds and bond types were converted '// &                   !
                          'in the lammps format'//REPEAT(' ',11)//'|';                   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Deallocate arrays containing bond properties for amber #####################################
                                                                                         !
    deallocate(BONDS_INC_HYDROGEN);                                                      !
                                                                                         !
    if ( AMBER_NBONA > 0 ) then;                                                         !
                                                                                         !
        deallocate(BONDS_WITHOUT_HYDROGEN);                                              !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays containing angle properties for lammps #####################################
                                                                                         !
    allocate(ANGLE_ATOMID(1:3,1:NANGLE));                                                !
                                                                                         !
    allocate(ANGLE_TYPE(1:NANGLE));                                                      !
                                                                                         !
!   ### Initialization of arrays containing angle properties for lammps ############################
                                                                                         !
    ANGLE_ATOMID(1:3,1:NANGLE) = 0;                                                      !
                                                                                         !
    ANGLE_TYPE(1:NANGLE) = 0;                                                            !
                                                                                         !
!   ### Set angle properties for lammps ############################################################
                                                                                         !
    iangles = 1;                                                                         !
                                                                                         !
    IANGLE = 0;                                                                          !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, AMBER_NTHETH;                                                              !
                                                                                         !
        IANGLE = IANGLE + 1;                                                             !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF2 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF3 = IOSEF2 + 1;                                                             !
                                                                                         !
        IOSEF4 = IOSEF3 + 1;                                                             !
                                                                                         !
        ANGLE_ATOMID(1,IANGLE) = ANGLES_INC_HYDROGEN(IOSEF1) / 3 + 1;                    !
                                                                                         !
        ANGLE_ATOMID(2,IANGLE) = ANGLES_INC_HYDROGEN(IOSEF2) / 3 + 1;                    !
                                                                                         !
        ANGLE_ATOMID(3,IANGLE) = ANGLES_INC_HYDROGEN(IOSEF3) / 3 + 1;                    !
                                                                                         !
        ANGLE_TYPE(IANGLE) = ANGLES_INC_HYDROGEN(IOSEF4);                                !
                                                                                         !
        IOSEF1 = IOSEF4;                                                                 !
                                                                                         !
    end do                                                                               !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, AMBER_NTHETA;                                                              !
                                                                                         !
        IANGLE = IANGLE + 1;                                                             !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF2 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF3 = IOSEF2 + 1;                                                             !
                                                                                         !
        IOSEF4 = IOSEF3 + 1;                                                             !
                                                                                         !
        ANGLE_ATOMID(1,IANGLE) = ANGLES_WITHOUT_HYDROGEN(IOSEF1) / 3 + 1;                !
                                                                                         !
        ANGLE_ATOMID(2,IANGLE) = ANGLES_WITHOUT_HYDROGEN(IOSEF2) / 3 + 1;                !
                                                                                         !
        ANGLE_ATOMID(3,IANGLE) = ANGLES_WITHOUT_HYDROGEN(IOSEF3) / 3 + 1;                !
                                                                                         !
        ANGLE_TYPE(IANGLE) = ANGLES_WITHOUT_HYDROGEN(IOSEF4);                            !
                                                                                         !
        IOSEF1 = IOSEF4;                                                                 !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Angles and angle types were converted '// &                 !
                          'in the lammps format'//REPEAT(' ',9)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Deallocate arrays containing angle properties for amber ####################################
                                                                                         !
    if ( AMBER_NTHETH > 0 ) then;                                                        !
                                                                                         !
        deallocate(ANGLES_INC_HYDROGEN);                                                 !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( AMBER_NTHETA > 0 ) then;                                                        !
                                                                                         !
        deallocate(ANGLES_WITHOUT_HYDROGEN);                                             !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays containing dihedral properties for lammps ##################################
                                                                                         !
    allocate(DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL));                                          !
                                                                                         !
    allocate(DIHEDRAL_TYPE(1:NDIHEDRAL));                                                !
                                                                                         !
!   ### Initialization of arrays containing dihedral properties for lammps #########################
                                                                                         !
    DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL) = 0;                                                !
                                                                                         !
    DIHEDRAL_TYPE(1:NDIHEDRAL) = 0;                                                      !
                                                                                         !
!   ### Set dihedral properties for lammps #########################################################
                                                                                         !
    idihedrals = 1;                                                                      !
                                                                                         !
    IDIHEDRAL = 0;                                                                       !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, AMBER_NPHIH;                                                               !
                                                                                         !
        IDIHEDRAL = IDIHEDRAL + 1;                                                       !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF2 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF3 = IOSEF2 + 1;                                                             !
                                                                                         !
        IOSEF4 = IOSEF3 + 1;                                                             !
                                                                                         !
        IOSEF5 = IOSEF4 + 1;                                                             !
                                                                                         !
        DIHEDRAL_ATOMID(1,IDIHEDRAL) = ABS( DIHEDRALS_INC_HYDROGEN(IOSEF1) ) / 3 + 1;    !
                                                                                         !
        DIHEDRAL_ATOMID(2,IDIHEDRAL) = ABS( DIHEDRALS_INC_HYDROGEN(IOSEF2) ) / 3 + 1;    !
                                                                                         !
        DIHEDRAL_ATOMID(3,IDIHEDRAL) = ABS( DIHEDRALS_INC_HYDROGEN(IOSEF3) ) / 3 + 1;    !
                                                                                         !
        DIHEDRAL_ATOMID(4,IDIHEDRAL) = ABS( DIHEDRALS_INC_HYDROGEN(IOSEF4) ) / 3 + 1;    !
                                                                                         !
        DIHEDRAL_TYPE(IDIHEDRAL) = DIHEDRALS_INC_HYDROGEN(IOSEF5);                       !
                                                                                         !
        IOSEF1 = IOSEF5;                                                                 !
                                                                                         !
    end do                                                                               !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, AMBER_NPHIA;                                                               !
                                                                                         !
        IDIHEDRAL = IDIHEDRAL + 1;                                                       !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF2 = IOSEF1 + 1;                                                             !
                                                                                         !
        IOSEF3 = IOSEF2 + 1;                                                             !
                                                                                         !
        IOSEF4 = IOSEF3 + 1;                                                             !
                                                                                         !
        IOSEF5 = IOSEF4 + 1;                                                             !
                                                                                         !
        DIHEDRAL_ATOMID(1,IDIHEDRAL) = ABS( DIHEDRALS_WITHOUT_HYDROGEN(IOSEF1) ) / 3 + 1;!
                                                                                         !
        DIHEDRAL_ATOMID(2,IDIHEDRAL) = ABS( DIHEDRALS_WITHOUT_HYDROGEN(IOSEF2) ) / 3 + 1;!
                                                                                         !
        DIHEDRAL_ATOMID(3,IDIHEDRAL) = ABS( DIHEDRALS_WITHOUT_HYDROGEN(IOSEF3) ) / 3 + 1;!
                                                                                         !
        DIHEDRAL_ATOMID(4,IDIHEDRAL) = ABS( DIHEDRALS_WITHOUT_HYDROGEN(IOSEF4) ) / 3 + 1;!
                                                                                         !
        DIHEDRAL_TYPE(IDIHEDRAL) = DIHEDRALS_WITHOUT_HYDROGEN(IOSEF5);                   !
                                                                                         !
        IOSEF1 = IOSEF5;                                                                 !
                                                                                         !
    end do                                                                               !
                                                                                         !
                                                                                         !
    write(icanal,'(a70)') '| Dihedrals and dihedral types were converted '// &           !
                          'in the lammps format'//REPEAT(' ',3)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Deallocate arrays containing dihedral properties for amber #################################
                                                                                         !
    if ( AMBER_NPHIH > 0 ) then;                                                         !
                                                                                         !
        deallocate(DIHEDRALS_INC_HYDROGEN);                                              !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( AMBER_NPHIA > 0 ) then;                                                         !
                                                                                         !
        deallocate(DIHEDRALS_WITHOUT_HYDROGEN);                                          !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set improper properties for lammps #########################################################
                                                                                         !
    iimpropers = 0;                                                                      !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine BUILD_PARAMETER_AMBER_TO_LAMMPS
