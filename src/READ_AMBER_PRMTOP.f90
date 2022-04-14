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

subroutine READ_AMBER_PRMTOP(icanal,CHEXT) 

!   ************************************************************************************************
!   **                                 READ LAMMPS FILES                                          **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal        : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                   **
!   ** CHEXT         : NAME OF THE LAMMPS CONFIGURATION                                           **
!   **                                                                                            **
!   ** FLAG_SORT_MOLEC   : FLAG TO CHECK IF MOLECULES AND ATOMS NEED TO BE SORTED                 ** 
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

    character (len=150), intent(in) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: NREAD_LINE;

    integer (kind=4) :: ixy, ixz, iyz;

    integer (kind=4) :: ibonds, ibond_coeffs;

    integer (kind=4) :: iangles, iangle_coeffs;

    integer (kind=4) :: iimpropers, iimproper_coeffs, iangletorsion_coeffs, idihedrals;

    integer (kind=4) :: imcon, ICONNECT, IATOM_TYPE, ATOMS_FLAG;

    integer (kind=4) :: iread_charge;

    integer (kind=4) :: IFOUND;

    character (len=150) :: CHARLINE;

    integer (kind=4) :: EOF, EOF2;

!   integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5;

    integer (kind=4), dimension(1:31) :: TAB_IOSEF;

!   real (kind=8) :: ROSEF1, ROSEF2, ROSEF3;

    real (kind=8), dimension(1:4) :: TAB_ROSEF; 

!   character (len=150) :: CHOSEF1, CHOSEF2, CHOSEF3, CHOSEF4, CHOSEF5, CHOSEF6;

    character (len=250) :: CHAIN_LENGTH;

    logical :: PROBE1;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE, ILENGTH_CHEXT;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Read AMBER prmtop file';                                                  !
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
!   ### Check the presence of the file in the current directory #################################### 
                                                                                         !
    inquire(file=TRIM(CHEXT),exist=PROBE1);                                              !
                                                                                         !
    if ( PROBE1 .EQV. .FALSE. ) then;                                                    !
                                                                                         !
        write(icanal,*) 'Error: the file was not found in the current directory - stop'; !
                                                                                         !
        stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of reading flags ############################################################
                                                                                         !
    iread_charge = 0;                                                                    !
                                                                                         !
!   ### Start reading the file #####################################################################
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
        if ( INDEX(CHAIN_LENGTH,'%FLAG') == 0 ) CYCLE;                                   !
                                                                                         !
        if ( INDEX(CHAIN_LENGTH,'POINTERS') > 0 ) then;                                  !
                                                                                         !
            read(2,*);                                                                   !
                                                                                         !
!           ### Read size of arrays ################################################################
                                                                                         !
            read(2,*,iostat=EOF) NATOM, NTYPE_ATOM, AMBER_NBONH,                       & !
                                 TAB_IOSEF(1), AMBER_NTHETH, TAB_IOSEF(2),             & !
                                 AMBER_NPHIH,  TAB_IOSEF(3:6),                         & !
                                 NRESIDUES, AMBER_NBONA, AMBER_NTHETA, AMBER_NPHIA,    & !
                                 NTYPE_BOND, NTYPE_ANGLE, NTYPE_DIHEDRAL, AMBER_NATYP;   !
                                                                                         !
!           ### Reading error management ###########################################################
                                                                                         !
            if ( EOF /= 0 ) then;                                                        !
                                                                                         !
                write(icanal,*) 'Error: there was a problem in reading POINTERS - stop'; !
                                                                                         !
                stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### 

            NTYPE_ATOM_SQ = NTYPE_ATOM * NTYPE_ATOM;                                     !
                                                                                         !
!           ### Determine the number of L-J interactions to read in the amber file #################
                                                                                         !
            AMBER_NLJ_INTER = ( NTYPE_ATOM_SQ + NTYPE_ATOM ) / 2;                        !
                                                                                         !
!           ### Write size of arrays ###############################################################
                                                                                         !
            write(icanal,'(a10,i8,a52)') '| NATOM : ', NATOM, REPEAT(' ',51)//'|';       !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            if ( NRESIDUES > 0 ) then;                                                   !
                                                                                         !
                write(icanal,'(a14,i8,a48)') '| NRESIDUES : ',    &                      !
                                             NRESIDUES,           &                      !
                                             REPEAT(' ',47)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            write(icanal,'(a19,i8,a43)') '| NTYPE_ATOM     : ', &                        !
                                         NTYPE_ATOM,            &                        !
                                         REPEAT(' ',42)//'|';                            !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            if ( NTYPE_BOND > 0 ) then;                                                  !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NTYPE_BOND     : ', &                    !
                                             NTYPE_BOND,            &                    !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( NTYPE_ANGLE > 0 ) then;                                                 !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NTYPE_ANGLE    : ', &                    !
                                             NTYPE_ANGLE,           &                    !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( NTYPE_DIHEDRAL > 0 ) then;                                              !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NTYPE_DIHEDRAL : ', &                    !
                                             NTYPE_DIHEDRAL,        &                    !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NATYP > 0 ) then;                                                 !
                                                                                         !
                write(icanal,'(a16,i8,a46)') '| AMBER_NATYP : ', &                       !
                                             AMBER_NATYP,        &                       !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NBONH > 0 ) then;                                                 !
                                                                                         !
                write(icanal,'(a16,i8,a46)') '| AMBER_NBONH : ', &                       !
                                             AMBER_NBONH,        &                       !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NBONA > 0 ) then;                                                 !
                                                                                         !
                write(icanal,'(a16,i8,a46)') '| AMBER_NBONA : ', &                       !
                                             AMBER_NBONA,        &                       !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NTHETH > 0 ) then;                                                !
                                                                                         !
                write(icanal,'(a17,i8,a45)') '| AMBER_NTHETH : ', &                      !
                                             AMBER_NTHETH,        &                      !
                                             REPEAT(' ',44)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NTHETA > 0 ) then;                                                !
                                                                                         !
                write(icanal,'(a17,i8,a45)') '| AMBER_NTHETA : ', &                      !
                                             AMBER_NTHETA,        &                      !
                                             REPEAT(' ',44)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NPHIH > 0 ) then;                                                 !
                                                                                         !
                write(icanal,'(a15,i8,a47)') '| AMBER_NPHIH : ',  &                      !
                                             AMBER_NPHIH,         &                      !
                                             REPEAT(' ',46)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NPHIA > 0 ) then;                                                 !
                                                                                         !
                write(icanal,'(a15,i8,a47)') '| AMBER_NPHIA : ',  &                      !
                                             AMBER_NPHIA,         &                      !
                                             REPEAT(' ',46)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Allocate arrays with the size of the the number of atoms ###########################
                                                                                         !
            allocate(CONFIG_QI(1:NATOM));                                                !
                                                                                         !
            allocate(CONFIG_MASS(1:NATOM));                                              !
                                                                                         !
            allocate(CONFIG_VI(1:3,1:NATOM));                                            !
                                                                                         !
            allocate(CONFIG_NAT(1:NATOM));                                               !
                                                                                         !
            allocate(CONFIG_ATOMID(1:NATOM));                                            !
                                                                                         !
            allocate(CONFIG_MOLECULEID(1:NATOM));                                        !
                                                                                         !
            allocate(CONFIG_ATOMIC_NUMBER(1:NATOM));                                     !
                                                                                         !
            allocate(CONFIG_NUMBER_EXCLUDED_ATOMS(1:NATOM));                             !
                                                                                         !
            allocate(CONFIG_ATOM_TYPE(1:NATOM));                                         !
                                                                                         !
            allocate(CONFIG_PCFF_TYPE(1:NATOM));                                         !
                                                                                         !
!           ### Initialization of arrays with the size of the number of atoms ######################
                                                                                         !
            CONFIG_QI(1:NATOM)         = 0.0d0;                                          !
                                                                                         !
            CONFIG_MASS(1:NATOM)       = 0.0d0;                                          !
                                                                                         !
            CONFIG_VI(1:3,1:NATOM)     = 0.0d0;                                          !
                                                                                         !
            CONFIG_NAT(1:NATOM)        = 'XXX';                                          !
                                                                                         !
            CONFIG_ATOMID(1:NATOM)     = 0;                                              !
                                                                                         !
            CONFIG_MOLECULEID(1:NATOM) = 0;                                              !
                                                                                         !
            CONFIG_ATOMIC_NUMBER(1:NATOM) = 0;                                           !
                                                                                         !
            CONFIG_NUMBER_EXCLUDED_ATOMS(1:NATOM) = 0;                                   !
                                                                                         !
            CONFIG_ATOM_TYPE(1:NATOM)  = 0;                                              !
                                                                                         !
            CONFIG_PCFF_TYPE(1:NATOM) = 'XXX';                                           !
                                                                                         !
!           ### Allocate arrays with the size of the number of atom types ##########################
                                                                                         !
            allocate(NONBONDED_PARM_INDEX(1:NTYPE_ATOM_SQ));                             !
                                                                                         !
!           ### Initialization of arrays with the size of the number of atom types #################
                                                                                         !
            NONBONDED_PARM_INDEX(1:NTYPE_ATOM_SQ) = 0;                                   !
                                                                                         !
            if ( NRESIDUES > 0 ) then;                                                   !
                                                                                         !
!               ### Allocate arrays with the size of the number of residues ########################
                                                                                         ! 
                allocate(RESIDUE_LABEL(1:NRESIDUES));                                    !
                                                                                         !
                allocate(RESIDUE_POINTER(1:NRESIDUES));                                  !
                                                                                         !  
!               ### Initialization of arrays with the size of the number of residues ###############
                                                                                         !
                RESIDUE_LABEL(1:NRESIDUES) = 'XXX';                                      !
                                                                                         !
                RESIDUE_POINTER(1:NRESIDUES) = 0;                                        !
                                                                                         !
             end if                                                                      !
                                                                                         !
             if ( NTYPE_BOND > 0 ) then;                                                 !
                                                                                         !
!                ### Allocate arrays with the size of the number of bond types #####################
                                                                                         !
                 allocate(BOND_FORCE_CONSTANT(1:NTYPE_BOND));                            !
                                                                                         !
                 allocate(BOND_EQUIL_VALUE(1:NTYPE_BOND));                               !
                                                                                         !
!                ### Initialization of arrays with the size of the number of bond types ############
                                                                                         !
                 BOND_FORCE_CONSTANT(1:NTYPE_BOND) = 0.0d0;                              !
                                                                                         !
                 BOND_EQUIL_VALUE(1:NTYPE_BOND) = 0.0d0;                                 !
                                                                                         !
             end if                                                                      !
                                                                                         !
             if ( NTYPE_ANGLE > 0 ) then;                                                ! 
                                                                                         !
!               ### Allocate arrays with the size of the number of angle types #####################
                                                                                         !
                allocate(ANGLE_FORCE_CONSTANT(1:NTYPE_ANGLE));                           !
                                                                                         !
                allocate(ANGLE_EQUIL_VALUE(1:NTYPE_ANGLE));                              !
                                                                                         !
!               ### Initialization of arrays with the size of the number of angle types ###########
                                                                                         !
                ANGLE_FORCE_CONSTANT(1:NTYPE_ANGLE) = 0.0d0;                             !
                                                                                         !
                ANGLE_EQUIL_VALUE(1:NTYPE_ANGLE) = 0.0d0;                                !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( NTYPE_DIHEDRAL > 0 ) then;                                              !
                                                                                         !
!               ### Allocate arrays with the size of the number of dihedral types ##################
                                                                                         !
                allocate(DIHEDRAL_FORCE_CONSTANT(1:NTYPE_DIHEDRAL));                     !
                                                                                         !
                allocate(DIHEDRAL_PERIODICITY(1:NTYPE_DIHEDRAL));                        !
                                                                                         !
                allocate(DIHEDRAL_PHASE(1:NTYPE_DIHEDRAL));                              !   
                                                                                         !
                allocate(SCEE_SCALE_FACTOR(1:NTYPE_DIHEDRAL));                           !
                                                                                         !
                allocate(SCNB_SCALE_FACTOR(1:NTYPE_DIHEDRAL));                           !
                                                                                         !
!               ### Initialization of arrays with the size of the number of dihedral types #########
                                                                                         !
                DIHEDRAL_FORCE_CONSTANT(1:NTYPE_DIHEDRAL) = 0.0d0;                       !
                                                                                         !
                DIHEDRAL_PERIODICITY(1:NTYPE_DIHEDRAL) = 0.0d0;                          !
                                                                                         !
                DIHEDRAL_PHASE(1:NTYPE_DIHEDRAL) = 0.0d0;                                !
                                                                                         !
                SCEE_SCALE_FACTOR(1:NTYPE_DIHEDRAL) = 0.0d0;                             !
                                                                                         !
                SCNB_SCALE_FACTOR(1:NTYPE_DIHEDRAL) = 0.0d0;                             ! 
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NATYP > 0 ) then;                                                 !
                                                                                         !
!               ### Allocate arrays with the size of the number of atom types in parameter file ####
                                                                                         !
                allocate(SOLTY(1:AMBER_NATYP));                                          !
                                                                                         !
!               ### Initialization of arrays with the size of the nbre of atom types ###############
                                                                                         !
                SOLTY(1:AMBER_NATYP) = 0.0d0;                                            !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NLJ_INTER > 0 ) then;                                             !
                                                                                         !
                write(icanal,'(a20,i8,a42)') '| AMBER_NLJ_INTER : ', &                   !
                                             AMBER_NLJ_INTER,        &                   !
                                             REPEAT(' ',41)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               ### Allocate arrays for reading Lennard-Jones interaction parameters ###############
                                                                                         !
                allocate(LENNARD_JONES_ACOEF(1:AMBER_NLJ_INTER));                        !
                                                                                         !
                allocate(LENNARD_JONES_BCOEF(1:AMBER_NLJ_INTER));                        !
                                                                                         !
!               ### Initialization of arrays for reading Lennard-Jones interaction parameters ######
                                                                                         !  
                LENNARD_JONES_ACOEF(1:AMBER_NLJ_INTER) = 0.0d0;                          !
                                                                                         !
                LENNARD_JONES_BCOEF(1:AMBER_NLJ_INTER) = 0.0d0;                          !
                                                                                         ! 
            end if                                                                       !
                                                                                         !
            if ( AMBER_NBONH > 0 ) then;                                                 !
                                                                                         !
                IOSEF1 = 3 * AMBER_NBONH;                                                !
                                                                                         !
!               ### Allocate arrays for bonds containing hydrogen atoms ############################
                                                                                         !
                allocate(BONDS_INC_HYDROGEN(1:IOSEF1));                                  !
                                                                                         !
!               ### Initialization of arrays for bonds containing hydrogen atoms ###################
                                                                                         !
                BONDS_INC_HYDROGEN(1:IOSEF1) = 0;                                        !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NBONA > 0 ) then;                                                 !
                                                                                         !
                IOSEF1 = 3 * AMBER_NBONA;                                                !
                                                                                         !
!               ### Allocate arrays for bonds that does not contain hydrogen atoms #################
                                                                                         !
                allocate(BONDS_WITHOUT_HYDROGEN(1:IOSEF1));                              !
                                                                                         !
!               ### Initialization of arrays for bonds that does not contain hydrogen atoms ########
                                                                                         !
                BONDS_WITHOUT_HYDROGEN(1:IOSEF1) = 0;                                    !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NTHETH > 0 ) then;                                                !
                                                                                         !
                IOSEF1 = 4 * AMBER_NTHETH;                                               !
                                                                                         !
!               ### Allocate arrays for angles containing hydrogen atoms ###########################
                                                                                         !
                allocate(ANGLES_INC_HYDROGEN(1:IOSEF1));                                 !
                                                                                         !
!               ### Initialization of arrays for angles containing hydrogen atoms ##################
                                                                                         !
                ANGLES_INC_HYDROGEN(1:IOSEF1) = 0;                                       !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NTHETA > 0 ) then;                                                !
                                                                                         !
                IOSEF1 = 4 * AMBER_NTHETA;                                               !
                                                                                         !
!               ### Allocate arrays for angles that does not contain hydrogen atoms ################
                                                                                         !
                allocate(ANGLES_WITHOUT_HYDROGEN(1:IOSEF1));                             !
                                                                                         !
!               ### Initialization of arrays for angles that does not contain hydrogen atoms #######
                                                                                         !
                ANGLES_WITHOUT_HYDROGEN(1:IOSEF1) = 0;                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NPHIH > 0 ) then;                                                 !
                                                                                         !
                IOSEF1 = 5 * AMBER_NPHIH;                                                !
                                                                                         !
!               ### Allocate arrays for dihedrals containing hydrogen atoms ########################
                                                                                         !
                allocate(DIHEDRALS_INC_HYDROGEN(1:IOSEF1));                              !
                                                                                         !
!               ### Initialization of arrays for dihedrals containing hydrogen atoms ###############
                                                                                         !
                DIHEDRALS_INC_HYDROGEN(1:IOSEF1) = 0;                                    !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( AMBER_NPHIA > 0 ) then;                                                 ! 
                                                                                         !
                IOSEF1 = 5 * AMBER_NPHIA;                                                !
                                                                                         !
!               ### Allocate arrays for dihedrals that does not contain hydrogen atoms #############
                                                                                         !
                allocate(DIHEDRALS_WITHOUT_HYDROGEN(1:IOSEF1));                          !
                                                                                         !
!               ### Initialization of arrays for dihedrals that does not contain hydrogen atoms ####
                                                                                         !
                DIHEDRALS_WITHOUT_HYDROGEN(1:IOSEF1) = 0;                                !
                                                                                         !
            end if                                                                       !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'ATOM_NAME') > 0 ) then;                            !
                                                                                         !
            write(icanal,'(a70)') '| Entering ATOM_NAME section '//REPEAT('>',40)//'|';  !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            NREAD_LINE = CEILING( REAL( NATOM ) / 20.0d0 );                              !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           ### Read the format of ATOM_NAME #######################################################
                                                                                         !
            read(2,'(a)',iostat=EOF) CHAIN_LENGTH;                                       !
                                                                                         !
            if ( EOF /= 0 ) then;                                                        !
                                                                                         !
                write(icanal,*) 'There was a problem in reading the format of ATOM_NAME - stop';
                                                                                         !
                stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
            IOSEF1 = INDEX(CHAIN_LENGTH,'%FORMAT');                                      !                
                                                                                         !
            if ( INDEX(CHAIN_LENGTH,'%FORMAT') > 0 ) then;                               !
                                                                                         !
                IOSEF1 = IOSEF1 + 8;                                                     !
                                                                                         !
                IOSEF2 = LEN_TRIM(CHAIN_LENGTH) - 1;                                     !
                                                                                         ! 
                CHOSEF1 = CHAIN_LENGTH(IOSEF1:IOSEF2);                                   !
                                                                                         !
                IOSEF3 = 70 - 16 - 1 - LEN_TRIM(CHOSEF1);                                !
                                                                                         !
                write(icanal,'(a70)') '| The format is '//     &                         !
                                      TRIM(CHOSEF1)//          &                         !
                                      REPEAT(' ',IOSEF3)//'|';                           !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            else                                                                         !
                                                                                         !
                write(icanal,*) 'Error in reading the format for ATOM_NAME - stop';      !
                                                                                         !
                stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NATOM ) then;                                             !
                                                                                         !
                    read(2,'('//TRIM(CHOSEF1)//')',iostat=EOF) CONFIG_NAT(IOSEF1);       !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 19;                                                !
                                                                                         !
                    if ( IOSEF2 > NATOM ) IOSEF2 = NATOM;                                !
                                                                                         !
                    read(2,'('//TRIM(CHOSEF1)//')',iostat=EOF) CONFIG_NAT(IOSEF1:IOSEF2);!
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
                if ( EOF /= 0 ) then; 

                    write(icanal,*) 'Error: there was a problem in reading ATOM_NAME - stop';

                    stop; !//////////////////////////////////////////////////////////////!

                end if

            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Atom names were read'//REPEAT(' ',47)//'|';         !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              ! 
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'CHARGE') > 0 ) then;                               !
                                                                                         !
            write(icanal,'(a70)') '| Entering CHARGE section '//REPEAT('>',43)//'|';     !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NATOM ) / 5.0d0 );                               !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NATOM ) then;                                             !
                                                                                         !
                    read(2,*) CONFIG_QI(IOSEF1);                                         !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NATOM ) IOSEF2 = NATOM;                                !
                                                                                         !
                    read(2,*) CONFIG_QI(IOSEF1:IOSEF2);                                  !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Charges were read'//REPEAT(' ',50)//'|';            !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            iread_charge = 1;                                                            !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'ATOMIC_NUMBER') > 0 ) then;                        !
                                                                                         !
            write(icanal,'(a70)') '| Entering ATOMIC_NUMBER section '// &                !
                                  REPEAT('>',36)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NATOM ) / 10.0d0 );                              !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NATOM ) then;                                             !
                                                                                         !
                    read(2,*) CONFIG_ATOMIC_NUMBER(IOSEF1);                              !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > NATOM ) IOSEF2 = NATOM;                                !
                                                                                         !
                    read(2,*) CONFIG_ATOMIC_NUMBER(IOSEF1:IOSEF2);                       !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Atomic numbers were read'//REPEAT(' ',43)//'|';     !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'MASS') > 0 ) then;                                 !
                                                                                         !
            write(icanal,'(a70)') '| Entering MASS section '// &                         !
                                  REPEAT('>',45)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NATOM ) / 5.0d0 );                               !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NATOM ) then;                                             !
                                                                                         !
                    read(2,*) CONFIG_MASS(IOSEF1);                                       !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NATOM ) IOSEF2 = NATOM;                                !
                                                                                         !
                    read(2,*) CONFIG_MASS(IOSEF1:IOSEF2);                                !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Masses were read'//REPEAT(' ',51)//'|';             !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'ATOM_TYPE_INDEX') > 0 ) then;                      !
                                                                                         !
            write(icanal,'(a70)') '| Entering ATOM_TYPE_INDEX section '// &              !
                                  REPEAT('>',34)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NATOM ) / 10.0d0 );                              !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NATOM ) then;                                             !
                                                                                         !
                    read(2,*) CONFIG_ATOM_TYPE(IOSEF1);                                  !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > NATOM ) IOSEF2 = NATOM;                                !
                                                                                         !
                    read(2,*) CONFIG_ATOM_TYPE(IOSEF1:IOSEF2);                           !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Atom types were read'//REPEAT(' ',47)//'|';         !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'NUMBER_EXCLUDED_ATOMS') > 0 ) then;                !
                                                                                         !
            write(icanal,'(a70)') '| Entering NUMBER_EXCLUDED_ATOMS section '// &        !
                                  REPEAT('>',28)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NATOM ) / 10.0d0 );                              !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NATOM ) then;                                             !
                                                                                         !
                    read(2,*) CONFIG_NUMBER_EXCLUDED_ATOMS(IOSEF1);                      !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > NATOM ) IOSEF2 = NATOM;                                !
                                                                                         !
                    read(2,*) CONFIG_NUMBER_EXCLUDED_ATOMS(IOSEF1:IOSEF2);               !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Number of excluded atoms were read'// &             !
                                  REPEAT(' ',33)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'NBONDED_PARM_INDEX') > 0 ) then;                   !
                                                                                         !
            write(icanal,'(a70)') '| Entering NBONDED_PARM_INDEX section '// &           !
                                  REPEAT('>',31)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NTYPE_ATOM_SQ ) / 10.0d0 );                      !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NTYPE_ATOM_SQ ) then;                                     !
                                                                                         !
                    read(2,*) NONBONDED_PARM_INDEX(IOSEF1);                              !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > NTYPE_ATOM_SQ ) IOSEF2 = NTYPE_ATOM_SQ;                !
                                                                                         !
                    read(2,*) NONBONDED_PARM_INDEX(IOSEF1:IOSEF2);                       !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Nonbonded parameter index was read'// &             !
                                  REPEAT(' ',33)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'RESIDUE_LABEL') > 0 ) then;                        !
                                                                                         !
            write(icanal,'(a70)') '| Entering RESIDUE_LABEL section '// &                !
                                  REPEAT('>',36)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NRESIDUES ) / 20.0d0 );                          !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NRESIDUES ) then;                                         !
                                                                                         !
                    read(2,*) RESIDUE_LABEL(IOSEF1);                                     !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 19;                                                !
                                                                                         !
                    if ( IOSEF2 > NRESIDUES ) IOSEF2 = NRESIDUES;                        !
                                                                                         !
                    read(2,*) RESIDUE_LABEL(IOSEF1:IOSEF2);                              !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Residue labels were read'// &                       !
                                  REPEAT(' ',43)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'RESIDUE_POINTER') > 0 ) then;                      !
                                                                                         !
            write(icanal,'(a70)') '| Entering RESIDUE_POINTER section '// &              !
                                  REPEAT('>',34)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NRESIDUES ) / 10.0d0 );                          !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NRESIDUES ) then;                                         !
                                                                                         !
                    read(2,*) RESIDUE_POINTER(IOSEF1);                                   !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > NRESIDUES ) IOSEF2 = NRESIDUES;                        !
                                                                                         !
                    read(2,*) RESIDUE_POINTER(IOSEF1:IOSEF2);                            !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Residue pointers were read'// &                     !
                                  REPEAT(' ',41)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'BOND_FORCE_CONSTANT') > 0 ) then;                  !
                                                                                         !
            write(icanal,'(a70)') '| Entering BOND_FORCE_CONSTANT section '// &          !
                                  REPEAT('>',30)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NTYPE_BOND ) / 5.0d0 );                          !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NTYPE_BOND ) then;                                        !
                                                                                         !
                    read(2,*) BOND_FORCE_CONSTANT(IOSEF1);                               !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NTYPE_BOND ) IOSEF2 = NTYPE_BOND;                      !
                                                                                         !
                    read(2,*) BOND_FORCE_CONSTANT(IOSEF1:IOSEF2);                        !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Bond force constants were read'// &                 !
                                  REPEAT(' ',37)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'BOND_EQUIL_VALUE') > 0 ) then;                     !
                                                                                         !
            write(icanal,'(a70)') '| Entering BOND_EQUIL_VALUE section '// &             !
                                  REPEAT('>',33)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NTYPE_BOND ) / 5.0d0 );                          !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NTYPE_BOND ) then;                                        !
                                                                                         !
                    read(2,*) BOND_EQUIL_VALUE(IOSEF1);                                  !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NTYPE_BOND ) IOSEF2 = NTYPE_BOND;                      !
                                                                                         !
                    read(2,*) BOND_EQUIL_VALUE(IOSEF1:IOSEF2);                           !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Bond equilibrium values were read'// &              !
                                  REPEAT(' ',34)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'ANGLE_FORCE_CONSTANT') > 0 ) then;                 !
                                                                                         !
            write(icanal,'(a70)') '| Entering ANGLE_FORCE_CONSTANT section '// &         !
                                  REPEAT('>',29)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NTYPE_ANGLE ) / 5.0d0 );                         !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NTYPE_ANGLE ) then;                                       !
                                                                                         !
                    read(2,*) ANGLE_FORCE_CONSTANT(IOSEF1);                              !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NTYPE_ANGLE ) IOSEF2 = NTYPE_ANGLE;                    !
                                                                                         !
                    read(2,*) ANGLE_FORCE_CONSTANT(IOSEF1:IOSEF2);                       !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Angle force constants were read'// &                !
                                  REPEAT(' ',36)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'ANGLE_EQUIL_VALUE') > 0 ) then;                    !
                                                                                         !
            write(icanal,'(a70)') '| Entering ANGLE_EQUIL_VALUE section '// &            !
                                  REPEAT('>',32)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NTYPE_ANGLE ) / 5.0d0 );                         !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NTYPE_ANGLE ) then;                                       !
                                                                                         !
                    read(2,*) ANGLE_EQUIL_VALUE(IOSEF1);                                 !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NTYPE_ANGLE ) IOSEF2 = NTYPE_ANGLE;                    !
                                                                                         !
                    read(2,*) ANGLE_EQUIL_VALUE(IOSEF1:IOSEF2);                          !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Angle equilibrium values were read'// &             !
                                  REPEAT(' ',33)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'DIHEDRAL_FORCE_CONSTANT') > 0 ) then;              !
                                                                                         !
            write(icanal,'(a70)') '| Entering DIHEDRAL_FORCE_CONSTANT section '// &      !
                                  REPEAT('>',26)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NTYPE_DIHEDRAL ) / 5.0d0 );                      !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NTYPE_DIHEDRAL ) then;                                    !
                                                                                         !
                    read(2,*) DIHEDRAL_FORCE_CONSTANT(IOSEF1);                           !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NTYPE_DIHEDRAL ) IOSEF2 = NTYPE_DIHEDRAL;              !
                                                                                         !
                    read(2,*) DIHEDRAL_FORCE_CONSTANT(IOSEF1:IOSEF2);                    !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Dihedral force constants were read'// &             !
                                  REPEAT(' ',33)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'DIHEDRAL_PERIODICITY') > 0 ) then;                 !
                                                                                         !
            write(icanal,'(a70)') '| Entering DIHEDRAL_PERIODICITY section '// &         !
                                  REPEAT('>',29)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NTYPE_DIHEDRAL ) / 5.0d0 );                      !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NTYPE_DIHEDRAL ) then;                                    !
                                                                                         !
                    read(2,*) DIHEDRAL_PERIODICITY(IOSEF1);                              !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NTYPE_DIHEDRAL ) IOSEF2 = NTYPE_DIHEDRAL;              !
                                                                                         !
                    read(2,*) DIHEDRAL_PERIODICITY(IOSEF1:IOSEF2);                       !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Dihedral periodicities were read'// &               !
                                  REPEAT(' ',35)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              ! 
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'DIHEDRAL_PHASE') > 0 ) then;                       !
                                                                                         !
            write(icanal,'(a70)') '| Entering DIHEDRAL_PHASE section '// &               !
                                  REPEAT('>',35)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NTYPE_DIHEDRAL ) / 5.0d0 );                      !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NTYPE_DIHEDRAL ) then;                                    !
                                                                                         !
                    read(2,*) DIHEDRAL_PHASE(IOSEF1);                                    !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NTYPE_DIHEDRAL ) IOSEF2 = NTYPE_DIHEDRAL;              !
                                                                                         !
                    read(2,*) DIHEDRAL_PHASE(IOSEF1:IOSEF2);                             !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Dihedral phases were read'// &                      !
                                  REPEAT(' ',42)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              ! 
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'SCEE_SCALE_FACTOR') > 0 ) then;                    !
                                                                                         !
            write(icanal,'(a70)') '| Entering SCEE_SCALE_FACTOR section '// &            !
                                  REPEAT('>',32)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NTYPE_DIHEDRAL ) / 5.0d0 );                      !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NTYPE_DIHEDRAL ) then;                                    !
                                                                                         !
                    read(2,*) SCEE_SCALE_FACTOR(IOSEF1);                                 !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NTYPE_DIHEDRAL ) IOSEF2 = NTYPE_DIHEDRAL;              !
                                                                                         !
                    read(2,*) SCEE_SCALE_FACTOR(IOSEF1:IOSEF2);                          !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| SCEE scale factors were read'// &                   !
                                  REPEAT(' ',39)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              ! 
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'SCNB_SCALE_FACTOR') > 0 ) then;                    !
                                                                                         !
            write(icanal,'(a70)') '| Entering SCNB_SCALE_FACTOR section '// &            !
                                  REPEAT('>',32)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( NTYPE_DIHEDRAL ) / 5.0d0 );                      !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == NTYPE_DIHEDRAL ) then;                                    !
                                                                                         !
                    read(2,*) SCNB_SCALE_FACTOR(IOSEF1);                                 !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > NTYPE_DIHEDRAL ) IOSEF2 = NTYPE_DIHEDRAL;              !
                                                                                         !
                    read(2,*) SCNB_SCALE_FACTOR(IOSEF1:IOSEF2);                          !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| SCNB scale factors were read'// &                   !
                                  REPEAT(' ',39)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              ! 
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'SOLTY') > 0 ) then;                                !
                                                                                         !
            write(icanal,'(a70)') '| Entering SOLTY section '// &                        !
                                  REPEAT('>',44)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( AMBER_NATYP ) / 5.0d0 );                         !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == AMBER_NATYP ) then;                                       !
                                                                                         !
                    read(2,*) SOLTY(IOSEF1);                                             !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > AMBER_NATYP ) IOSEF2 = AMBER_NATYP;                    !
                                                                                         !
                    read(2,*) SOLTY(IOSEF1:IOSEF2);                                      !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| SOLTY values were read'// &                         !
                                  REPEAT(' ',45)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !  
                                                                                         !                              
        else if ( INDEX(CHAIN_LENGTH,'LENNARD_JONES_ACOEF') > 0 ) then;                  !
                                                                                         !
            write(icanal,'(a70)') '| Entering LENNARD_JONES_ACOEF section '// &          !
                                  REPEAT('>',30)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( AMBER_NLJ_INTER ) / 5.0d0 );                     !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == AMBER_NLJ_INTER ) then;                                   !
                                                                                         !
                    read(2,*) LENNARD_JONES_ACOEF(IOSEF1);                               !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > AMBER_NLJ_INTER ) IOSEF2 = AMBER_NLJ_INTER;            !
                                                                                         !
                    read(2,*) LENNARD_JONES_ACOEF(IOSEF1:IOSEF2);                        !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Lennard-Jones A coefficients were read'// &         !
                                  REPEAT(' ',29)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !  
                                                                                         ! 
        else if ( INDEX(CHAIN_LENGTH,'LENNARD_JONES_BCOEF') > 0 ) then;                  !
                                                                                         !
            write(icanal,'(a70)') '| Entering LENNARD_JONES_BCOEF section '// &          !
                                  REPEAT('>',30)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( AMBER_NLJ_INTER ) / 5.0d0 );                     !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == AMBER_NLJ_INTER ) then;                                   !
                                                                                         !
                    read(2,*) LENNARD_JONES_BCOEF(IOSEF1);                               !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 4;                                                 !
                                                                                         !
                    if ( IOSEF2 > AMBER_NLJ_INTER ) IOSEF2 = AMBER_NLJ_INTER;            !
                                                                                         !
                    read(2,*) LENNARD_JONES_BCOEF(IOSEF1:IOSEF2);                        !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Lennard-Jones B coefficients were read'// &         !
                                  REPEAT(' ',29)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !  
                                                                                         ! 
        else if ( INDEX(CHAIN_LENGTH,'BONDS_INC_HYDROGEN') > 0 ) then;                   !
                                                                                         !
            IOSEF5 = 3 * AMBER_NBONH;                                                    !
                                                                                         !
            write(icanal,'(a70)') '| Entering BONDS_INC_HYDROGEN section '// &           !
                                  REPEAT('>',31)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( IOSEF5 ) / 10.0d0 );                             !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == IOSEF5 ) then;                                            !
                                                                                         !
                    read(2,*) BONDS_INC_HYDROGEN(IOSEF1);                                !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > IOSEF5 ) IOSEF2 = IOSEF5;                              !
                                                                                         !
                    read(2,*) BONDS_INC_HYDROGEN(IOSEF1:IOSEF2);                         !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Bonds containing hydrogen atoms were read'// &      !
                                  REPEAT(' ',26)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'BONDS_WITHOUT_HYDROGEN') > 0 ) then;               !
                                                                                         !
            IOSEF5 = 3 * AMBER_NBONA;                                                    !
                                                                                         !
            write(icanal,'(a70)') '| Entering BONDS_WITHOUT_HYDROGEN section '// &       !
                                  REPEAT('>',27)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( IOSEF5 ) / 10.0d0 );                             !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == IOSEF5 ) then;                                            !
                                                                                         !
                    read(2,*) BONDS_WITHOUT_HYDROGEN(IOSEF1);                            !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > IOSEF5 ) IOSEF2 = IOSEF5;                              !
                                                                                         !
                    read(2,*) BONDS_WITHOUT_HYDROGEN(IOSEF1:IOSEF2);                     !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Bonds that does not contain '// &                   !
                                  'hydrogen atoms were read'//       &                   !
                                  REPEAT(' ',15)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'ANGLES_INC_HYDROGEN') > 0 ) then;                  !
                                                                                         !
            IOSEF5 = 4 * AMBER_NTHETH;                                                   !
                                                                                         !
            write(icanal,'(a70)') '| Entering ANGLES_INC_HYDROGEN section '// &          !
                                  REPEAT('>',30)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( IOSEF5 ) / 10.0d0 );                             !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == IOSEF5 ) then;                                            !
                                                                                         !
                    read(2,*) ANGLES_INC_HYDROGEN(IOSEF1);                               !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > IOSEF5 ) IOSEF2 = IOSEF5;                              !
                                                                                         !
                    read(2,*) ANGLES_INC_HYDROGEN(IOSEF1:IOSEF2);                        !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Angles containing hydrogen atoms were read'// &     !
                                  REPEAT(' ',25)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'ANGLES_WITHOUT_HYDROGEN') > 0 ) then;              !
                                                                                         !
            IOSEF5 = 4 * AMBER_NTHETA;                                                   !
                                                                                         !
            write(icanal,'(a70)') '| Entering ANGLES_WITHOUT_HYDROGEN section '// &      !
                                  REPEAT('>',26)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( IOSEF5 ) / 10.0d0 );                             !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == IOSEF5 ) then;                                            !
                                                                                         !
                    read(2,*) ANGLES_WITHOUT_HYDROGEN(IOSEF1);                           !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > IOSEF5 ) IOSEF2 = IOSEF5;                              !
                                                                                         !
                    read(2,*) ANGLES_WITHOUT_HYDROGEN(IOSEF1:IOSEF2);                    !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Angles that does not contain '// &                  !
                                  'hydrogen atoms were read'//        &                  !  
                                  REPEAT(' ',14)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'DIHEDRALS_INC_HYDROGEN') > 0 ) then;               !
                                                                                         !
            IOSEF5 = 5 * AMBER_NPHIH;                                                    !
                                                                                         !
            write(icanal,'(a70)') '| Entering DIHEDRALS_INC_HYDROGEN section '// &       !
                                  REPEAT('>',27)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( IOSEF5 ) / 10.0d0 );                             !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == IOSEF5 ) then;                                            !
                                                                                         !
                    read(2,*) DIHEDRALS_INC_HYDROGEN(IOSEF1);                            !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > IOSEF5 ) IOSEF2 = IOSEF5;                              !
                                                                                         !
                    read(2,*) DIHEDRALS_INC_HYDROGEN(IOSEF1:IOSEF2);                     !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Dihedrals containing hydrogen atoms were read'// &  !  
                                  REPEAT(' ',22)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'DIHEDRALS_WITHOUT_HYDROGEN') > 0 ) then;           !
                                                                                         !
            IOSEF5 = 5 * AMBER_NPHIA;                                                    !
                                                                                         !
            write(icanal,'(a70)') '| Entering DIHEDRALS_WITHOUT_HYDROGEN section '// &   !
                                  REPEAT('>',23)//'|';                                   !
                                                                                         !            
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            read(2,*);                                                                   ! 
                                                                                         !
            NREAD_LINE = CEILING( REAL( IOSEF5 ) / 10.0d0 );                             !
                                                                                         !
            write(icanal,'(a11,i8,a51)') '| There is ',     &                            !
                                         NREAD_LINE,        &                            !
                                         ' lines to read'// &                            !
                                         REPEAT(' ',36)//'|';                            ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do i = 1, NREAD_LINE;                                                        !
                                                                                         ! 
                IOSEF1 = IOSEF1 + 1;                                                     !
                                                                                         !
                if ( IOSEF1 == IOSEF5 ) then;                                            !
                                                                                         !
                    read(2,*) DIHEDRALS_WITHOUT_HYDROGEN(IOSEF1);                        !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 9;                                                 !
                                                                                         !
                    if ( IOSEF2 > IOSEF5 ) IOSEF2 = IOSEF5;                              !
                                                                                         !
                    read(2,*) DIHEDRALS_WITHOUT_HYDROGEN(IOSEF1:IOSEF2);                 !
                                                                                         !
                    IOSEF1 = IOSEF2;                                                     ! 
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Dihedrals that does not contain '// &               !
                                  'hydrogen atoms were read'//           &               !  
                                  REPEAT(' ',11)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Closing the input data file ################################################################
                                                                                         !
    close(2);                                                                            !
                                                                                         !
!   ### Check that necessary data were read in the prmtop file #####################################
                                                                                         !
    if ( iread_charge == 0 ) then;                                                       !
                                                                                         !
        write(icanal,*) 'Error: charges were not read in the prmtop file - stop';        !
                                                                                         !
        stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine READ_AMBER_PRMTOP
