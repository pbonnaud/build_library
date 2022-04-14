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

subroutine READ_LAMMPS_CONFIG(icanal,CHEXT,                      &
                                     MATA,                       &
                                     FLAG_SORT_MOLEC)

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

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=150), intent(in) :: CHEXT;

    real (kind=8), dimension(1:6), intent(out) :: MATA;

    integer (kind=4), intent(out) :: FLAG_SORT_MOLEC;

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

    character (len=10) :: CHOSEF1, CHOSEF2, CHOSEF3, CHOSEF4, CHOSEF5, CHOSEF6;

    character (len=250) :: CHAIN_LENGTH;

    logical :: PROBE1;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE, ILENGTH_CHEXT;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Read LAMMPS configuration';                                               !
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
    NATOM     = 0;

    NBOND     = 0;

    NANGLE    = 0;

    NDIHEDRAL = 0;

    NIMPROPER = 0;

    NTYPE_ATOM     = 0;

    NTYPE_BOND     = 0;

    NTYPE_ANGLE    = 0; 

    NTYPE_DIHEDRAL = 0;

    NTYPE_IMPROPER = 0; 

    XLO = 0.0d0;

    XHI = 0.0d0;

    YLO = 0.0d0;

    YHI = 0.0d0;

    ZLO = 0.0d0;

    ZHI = 0.0d0;

    MATA(1:6) = 0.0d0;

    FLAG_SORT_MOLEC = 0;

    IPOTENTIAL_CLASS2 = 0;

    ibonds     = 0;

    iangles    = 0;

    iimpropers = 0;

    idihedrals = 0;

    ibond_coeffs         = 0;

    iangle_coeffs        = 0;

    iimproper_coeffs     = 0;

    iangletorsion_coeffs = 0;

    NPARAM_BONDS = 0;

    NPARAM_ANGLES = 0;

    NPARAM_DIHEDRALS = 0;

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
        if ( INDEX(CHAIN_LENGTH,'atoms') > 0 ) then;                                     !
                                                                                         ! 
!           ### Read the number of atoms in the current molecular configuration ####################
                                                                                         !
            read(CHAIN_LENGTH,*) NATOM;                                                  !
                                                                                         !
            write(icanal,'(a19,i8,a43)') '| NATOM          : ', &                        !
                                         NATOM,                 &                        !
                                         REPEAT(' ',42)//'|';                            !
                                                                                         ! 
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Allocate arrays with the size of the the number of atoms ###########################
                                                                                         !
            allocate(CONFIG_QI(1:NATOM));

            allocate(CONFIG_RI(1:3,1:NATOM));

            allocate(CONFIG_VI(1:3,1:NATOM));

            allocate(CONFIG_NAT(1:NATOM));

            allocate(CONFIG_ATOMID(1:NATOM));

            allocate(CONFIG_MOLECULEID(1:NATOM));

            allocate(CONFIG_ATOM_TYPE(1:NATOM));

            allocate(CONFIG_PCFF_TYPE(1:NATOM));

            CONFIG_QI(1:NATOM)         = 0.0d0;

            CONFIG_RI(1:3,1:NATOM)     = 0.0d0;

            CONFIG_VI(1:3,1:NATOM)     = 0.0d0;

            CONFIG_NAT(1:NATOM)        = 'XXX';

            CONFIG_ATOMID(1:NATOM)     = 0;

            CONFIG_MOLECULEID(1:NATOM) = 0;

            CONFIG_ATOM_TYPE(1:NATOM)  = 0;

            CONFIG_PCFF_TYPE(1:NATOM) = 'XXX';

!           stop; !//////////////////////////////////////////////////////////////////////! 
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'atom types') > 0 ) then;                           !
                                                                                         !
            read(CHAIN_LENGTH,*) NTYPE_ATOM;                                             !
                                                                                         !
            write(icanal,'(a19,i8,a43)') '| NTYPE_ATOM     : ', &
                                         NTYPE_ATOM,            &
                                         REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            allocate(ATOM_LABEL(1:NTYPE_ATOM));

            allocate(ATOM_MASSE(1:NTYPE_ATOM));

            allocate(POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM));

            ATOM_LABEL(1:NTYPE_ATOM)           = 'XXX';

            ATOM_MASSE(1:NTYPE_ATOM)           = 0.0d0;

            POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM) = 0.0d0;

            POTENTIAL_CLASS2_CHTYPE = 'XXXX';                                            !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!

        else if ( INDEX(CHAIN_LENGTH,'bonds') > 0 ) then;

            read(CHAIN_LENGTH,*) NBOND;

            write(icanal,'(a19,i8,a43)') '| NBOND          : ', &
                                         NBOND,                 &
                                         REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            allocate(BOND_ATOMID(1:2,1:NBOND));

            allocate(BOND_TYPE(1:NBOND));

            allocate(BOND_PCFF_VERSION(1:NBOND));

            allocate(BOND_PCFF_TYPE(1:2,1:NBOND));

            allocate(BOND_PCFF_EQUIVALENCE(1:2,1:NBOND));

            BOND_ATOMID(1:2,1:NBOND) = 0;

            BOND_TYPE(1:NBOND)       = 0;

            BOND_PCFF_VERSION(1:NBOND) = 0.0d0;

            BOND_PCFF_TYPE(1:2,1:NBOND) = 'XXX';

            BOND_PCFF_EQUIVALENCE(1:2,1:NBOND) = 'XXX';

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'bond types') > 0 ) then;

            read(CHAIN_LENGTH,*) NTYPE_BOND;

            write(icanal,'(a19,i8,a43)') '| NTYPE_BOND     : ', &
                                         NTYPE_BOND,            &
                                         REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            allocate(BOND_COEFFS(1:4,1:NTYPE_BOND));

            BOND_COEFFS(1:4,1:NTYPE_BOND) = 0.0d0;

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'angles') > 0 ) then;

            read(CHAIN_LENGTH,*) NANGLE;

            write(icanal,'(a19,i8,a43)') '| NANGLE         : ', &
                                          NANGLE,               &
                                          REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            allocate(ANGLE_ATOMID(1:3,1:NANGLE));

            allocate(ANGLE_TYPE(1:NANGLE));

            ANGLE_ATOMID(1:3,1:NANGLE) = 0;

            ANGLE_TYPE(1:NANGLE)       = 0;

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'angle types') > 0 ) then;

            read(CHAIN_LENGTH,*) NTYPE_ANGLE;

            write(icanal,'(a19,i8,a43)') '| NTYPE_ANGLE    : ', &
                                         NTYPE_ANGLE,           &
                                         REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            allocate(ANGLE_COEFFS(1:4,1:NTYPE_ANGLE));

            allocate(BONDBOND_COEFFS(1:3,1:NTYPE_ANGLE));

            allocate(BONDANGLE_COEFFS(1:4,1:NTYPE_ANGLE));

            ANGLE_COEFFS(1:4,1:NTYPE_ANGLE)     = 0.0d0;

            BONDBOND_COEFFS(1:3,1:NTYPE_ANGLE)  = 0.0d0;

            BONDANGLE_COEFFS(1:4,1:NTYPE_ANGLE) = 0.0d0;

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'dihedrals') > 0 ) then;

            read(CHAIN_LENGTH,*) NDIHEDRAL;

            write(icanal,'(a19,i8,a43)') '| NDIHEDRAL      : ', &
                                         NDIHEDRAL,             &
                                         REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            allocate(DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL));

            allocate(DIHEDRAL_TYPE(1:NDIHEDRAL));

            DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL) = 0;

            DIHEDRAL_TYPE(1:NDIHEDRAL) = 0;

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'dihedral types') > 0 ) then;

            read(CHAIN_LENGTH,*) NTYPE_DIHEDRAL;

            write(icanal,'(a19,i8,a43)') '| NTYPE_DIHEDRAL : ', &
                                         NTYPE_DIHEDRAL,        &
                                         REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            allocate(DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL));

            allocate(MIDDLEBONDTORSION_COEFFS(1:4,1:NTYPE_DIHEDRAL));

            allocate(ENDBONDTORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL));

            allocate(ANGLETORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL));

            allocate(ANGLEANGLETORSION_COEFFS(1:3,1:NTYPE_DIHEDRAL));

            allocate(BONDBOND13_COEFFS(1:3,1:NTYPE_DIHEDRAL));

            DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL)          = 0.0d0;

            MIDDLEBONDTORSION_COEFFS(1:4,1:NTYPE_DIHEDRAL) = 0.0d0;

            ENDBONDTORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL)    = 0.0d0;

            ANGLETORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL)      = 0.0d0;

            ANGLEANGLETORSION_COEFFS(1:3,1:NTYPE_DIHEDRAL) = 0.0d0;

            BONDBOND13_COEFFS(1:3,1:NTYPE_DIHEDRAL)        = 0.0d0;

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'impropers') > 0 ) then;

            read(CHAIN_LENGTH,*) NIMPROPER;

            write(icanal,'(a19,i8,a43)') '| NIMPROPER      : ', &
                                         NIMPROPER,             &
                                         REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            allocate(IMPROPER_ATOMID(1:4,1:NIMPROPER));

            allocate(IMPROPER_TYPE(1:NIMPROPER));

            IMPROPER_ATOMID(1:4,1:NIMPROPER) = 0;

            IMPROPER_TYPE(1:NIMPROPER)       = 0;

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'improper types') > 0 ) then;

            read(CHAIN_LENGTH,*) NTYPE_IMPROPER;

            write(icanal,'(a19,i8,a43)') '| NTYPE_IMPROPER : ', &
                                         NTYPE_IMPROPER,        &
                                         REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            allocate(IMPROPER_COEFFS(1:2,1:NTYPE_IMPROPER));

            allocate(ANGLEANGLE_COEFFS(1:6,1:NTYPE_IMPROPER));

            IMPROPER_COEFFS(1:2,1:NTYPE_IMPROPER)   = 0.0d0;

            ANGLEANGLE_COEFFS(1:6,1:NTYPE_IMPROPER) = 0.0d0;

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'xlo xhi') > 0 ) then;

            read(CHAIN_LENGTH,*) XLO, XHI;

            MATA(1) = XHI - XLO;

            write(icanal,'(a19,f15.6,a36)') '| MATA(1)    [A] : ', &
                                            MATA(1),               &
                                            REPEAT(' ',35)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'ylo yhi') > 0 ) then;

            read(CHAIN_LENGTH,*) YLO, YHI;

            MATA(2) = YHI - YLO;

            write(icanal,'(a19,f15.6,a36)') '| MATA(2)    [A] : ', &
                                            MATA(2),               &
                                            REPEAT(' ',35)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'zlo zhi') > 0 ) then;

            read(CHAIN_LENGTH,*) ZLO, ZHI;

            MATA(3) = ZHI - ZLO;

            write(icanal,'(a19,f15.6,a36)') '| MATA(3)    [A] : ', &
                                            MATA(3),               &
                                            REPEAT(' ',35)//'|';
 
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'xy xz yz') > 0 ) then; 

            read(CHAIN_LENGTH,*) MATA(6), MATA(5), MATA(4);

            write(icanal,'(a19,f15.6,a36)') '| XY TILT [A]    : ', &
                                            MATA(6),               &
                                            REPEAT(' ',35)//'|';

            write(icanal,'(a19,f15.6,a36)') '| XZ TILT [A]    : ', &
                                            MATA(5),               &
                                            REPEAT(' ',35)//'|';

            write(icanal,'(a19,f15.6,a36)') '| YZ TILT [A]    : ', &
                                            MATA(4),               &
                                            REPEAT(' ',35)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'Masses') > 0 ) then;

            read(2,*);

            do i = 1, NTYPE_ATOM;

               read(2,*) IOSEF1, ATOM_MASSE(IOSEF1);

            end do

            do i = 1, NTYPE_ATOM

!               if ( ABS( ATOM_MASSE(i) -  1.00794d0 ) < 0.1d0 ) ATOM_LABEL(i) = 'H';
!               if ( ABS( ATOM_MASSE(i) - 12.01070d0 ) < 0.1d0 ) ATOM_LABEL(i) = 'C';
!               if ( ABS( ATOM_MASSE(i) - 14.0067d0  ) < 0.1d0 ) ATOM_LABEL(i) = 'N';
!               if ( ABS( ATOM_MASSE(i) - 15.9994d0  ) < 0.1d0 ) ATOM_LABEL(i) = 'O';
!               if ( ABS( ATOM_MASSE(i) - 28.0855d0  ) < 0.1d0 ) ATOM_LABEL(i) = 'Si';
!               if ( ABS( ATOM_MASSE(i) - 32.064d0   ) < 0.1d0 ) ATOM_LABEL(i) = 'S';
!               if ( ABS( ATOM_MASSE(i) - 40.078d0   ) < 0.1d0 ) ATOM_LABEL(i) = 'Ca';
!               if ( ABS( ATOM_MASSE(i) - 55.8450d0  ) < 0.1d0 ) ATOM_LABEL(i) = 'Fe';

                do j = 1, 200;                                                           !
                                                                                         !
                    ROSEF1 = ABS( ATOM_MASSE(i) - MOLAR_MASS_ELEMENT(j) );               !
                                                                                         !
                    if ( ROSEF1 < 0.1d0 ) ATOM_LABEL(i) = ELEMENT_SYMBOL_NAME(j);        !
                                                                                         !
                end do                                                                   !
                                                                                         !
            end do                                                                       ! 
                                                                                         !
            write(icanal,'(a70)') '| Masses were read'//REPEAT(' ',51)//'|';             !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!       else if ( ( INDEX(CHAIN_LENGTH,'Pair Coeffs') > 0 ) .OR. &                       !
!               ( INDEX(CHAIN_LENGTH,'Bond Coeffs') > 0 ) ) then;                        !
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'Pair Coeffs') > 0 ) then;                          !
                                                                                         !
            IPOTENTIAL_CLASS2 = 1;                                                       !
                                                                                         !
            if ( INDEX(CHAIN_LENGTH,'# lj/class2/coul/long') > 0 ) then;                 !
                                                                                         !
                POTENTIAL_CLASS2_CHTYPE = 'lj/class2/coul/long';                         !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           if ( INDEX(CHAIN_LENGTH,'Bond Coeffs') > 0 ) IPOTENTIAL_CLASS2 = 2;          !
                                                                                         !
            read(2,*);                                                                   !
                                                                                         !
            do i = 1, NTYPE_ATOM;                                                        !
                                                                                         ! 
                read(2,*) IOSEF1, POTENTIAL_CLASS2(1:2,IOSEF1);                          !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Pair Coeffs were read'//REPEAT(' ',46)//'|';        !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!       else if ( INDEX(CHAIN_LENGTH,'Bond Coeffs # class2') > 0 ) then;
        else if ( INDEX(CHAIN_LENGTH,'Bond Coeffs') > 0 ) then;                          !
                                                                                         !
            NPARAM_BONDS = 1;                                                            !
                                                                                         !
            if ( TRIM(CH_BOND_STYLE) == 'class2' ) NPARAM_BONDS = 4;                     !
                                                                                         !
            if ( TRIM(CH_BOND_STYLE) == 'harmonic' ) NPARAM_BONDS = 2;                   !

!           if ( INDEX(CHAIN_LENGTH,'# class2') > 0 ) then;                              !
                                                                                         !
!               NPARAM_BONDS = 4                                                         !
                                                                                         !
!           else                                                                         !
                                                                                         !
!               NPARAM_BONDS = 2                                                         !
                                                                                         !
!           end if

            read(2,*);

            do i = 1, NTYPE_BOND;

                read(2,*) IOSEF1, BOND_COEFFS(1:NPARAM_BONDS,IOSEF1);

            end do

            ibond_coeffs = 1;

            write(icanal,'(a70)') '| Bond Coeffs were read'//REPEAT(' ',46)//'|';        !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !  
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!       else if ( INDEX(CHAIN_LENGTH,'Angle Coeffs # class2') > 0 ) then;
        else if ( INDEX(CHAIN_LENGTH,'Angle Coeffs') > 0 ) then;

            NPARAM_ANGLES = 2;

            if ( INDEX(CHAIN_LENGTH,'# class2') > 0 ) NPARAM_ANGLES = 4; 

            read(2,*);

            do i = 1, NTYPE_ANGLE;

                read(2,*) IOSEF1, ANGLE_COEFFS(1:NPARAM_ANGLES,IOSEF1);

            end do

            iangle_coeffs = 1;

            write(icanal,'(a70)') '| Angle Coeffs were read'//REPEAT(' ',45)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !   
        else if ( INDEX(CHAIN_LENGTH,'BondBond Coeffs') > 0 ) then;

            read(2,*);

            do i = 1, NTYPE_ANGLE;

                read(2,*) IOSEF1, BONDBOND_COEFFS(1:3,IOSEF1);

            end do

            write(icanal,'(a70)') '| BondBond Coeffs were read'//REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'BondAngle Coeffs') > 0 ) then;

            read(2,*);

            do i = 1, NTYPE_ANGLE;

                read(2,*) IOSEF1, BONDANGLE_COEFFS(1:4,IOSEF1);

            end do

            write(icanal,'(a70)') '| BondAngle Coeffs were read'//REPEAT(' ',41)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!           stop;

        else if ( INDEX(CHAIN_LENGTH,'Dihedral Coeffs') > 0 ) then;

            NPARAM_DIHEDRALS = 6;

!           if ( INDEX(CHAIN_LENGTH,'# opls') > 0 ) NPARAM_DIHEDRALS = 4;

            if ( TRIM(CH_DIHEDRAL_STYLE) == 'opls' ) NPARAM_DIHEDRALS = 4;

            read(2,*);

            do i = 1, NTYPE_DIHEDRAL;

                read(2,*) IOSEF1, DIHEDRAL_COEFFS(1:NPARAM_DIHEDRALS,IOSEF1);

            end do

            write(icanal,'(a70)') '| Dihedral Coeffs were read'//REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!           stop; !//////////////////////////////////////////////////////////////////////!

        else if ( INDEX(CHAIN_LENGTH,'MiddleBondTorsion Coeffs') > 0 ) then;

            read(2,*);

            do i = 1, NTYPE_DIHEDRAL;

                read(2,*) IOSEF1, MIDDLEBONDTORSION_COEFFS(1:4,IOSEF1);

            end do

            write(icanal,'(a70)') '| MiddleBondTorsion Coeffs were read'//REPEAT(' ',33)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHAIN_LENGTH,'EndBondTorsion Coeffs') > 0 ) then;

            read(2,*);

            do i = 1, NTYPE_DIHEDRAL;

                read(2,*) IOSEF1, ENDBONDTORSION_COEFFS(1:8,IOSEF1);

            end do

            write(icanal,'(a70)') '| EndBondTorsion Coeffs were read'//REPEAT(' ',36)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHAIN_LENGTH,'AngleTorsion Coeffs # class2') > 0 ) then;

            read(2,*);

            do i = 1, NTYPE_DIHEDRAL;

                read(2,*) IOSEF1, ANGLETORSION_COEFFS(1:8,IOSEF1);

            end do

            iangletorsion_coeffs = 1;

            write(icanal,'(a70)') '| AngleTorsion Coeffs were read'//REPEAT(' ',38)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHAIN_LENGTH,'AngleAngleTorsion Coeffs') > 0 ) then;

            read(2,*);

            do i = 1, NTYPE_DIHEDRAL;

                read(2,*) IOSEF1, ANGLEANGLETORSION_COEFFS(1:3,IOSEF1);

            end do

            write(icanal,'(a70)') '| AngleAngleTorsion Coeffs were read'//REPEAT(' ',33)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHAIN_LENGTH,'BondBond13 Coeffs') > 0 ) then;

            read(2,*);

            do i = 1, NTYPE_DIHEDRAL;

                read(2,*) IOSEF1, BONDBOND13_COEFFS(1:3,IOSEF1);

            end do

            write(icanal,'(a70)') '| BondBond13 Coeffs were read'//REPEAT(' ',40)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHAIN_LENGTH,'Improper Coeffs # class2') > 0 ) then;

            read(2,*);

            do i = 1, NTYPE_IMPROPER;

                read(2,*) IOSEF1, IMPROPER_COEFFS(1:2,IOSEF1);

            end do

            iimproper_coeffs = 1;

            write(icanal,'(a70)') '| Improper Coeffs were read'//REPEAT(' ',42)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHAIN_LENGTH,'AngleAngle Coeffs') > 0 ) then;

            read(2,*);

            do i = 1, NTYPE_IMPROPER;

                read(2,*) IOSEF1, ANGLEANGLE_COEFFS(1:6,IOSEF1);

            end do

            write(icanal,'(a70)') '| AngleAngle Coeffs were read'//REPEAT(' ',40)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHAIN_LENGTH,'Atoms') > 0 ) then;

            ATOMS_FLAG = 0;

            if ( TRIM(CH_ATOM_STYLE) == 'charge' ) ATOMS_FLAG = 5;

            if ( TRIM(CH_ATOM_STYLE) == 'full' ) ATOMS_FLAG = 13;

!           if ( INDEX(CHAIN_LENGTH,'# charge') > 0 ) ATOMS_FLAG = 5;

!           if ( INDEX(CHAIN_LENGTH,'# full')   > 0 ) ATOMS_FLAG = 13;

            write(icanal,'(a15,i8,a47)') '| ATOMS_FLAG : ', ATOMS_FLAG, REPEAT(' ',46)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            read(2,*);

            do i = 1, NATOM;

                read(2,'(a)') CHAIN_LENGTH;

                EOF2 = 0;

                if ( ATOMS_FLAG == 5 ) then;                                             !
                                                                                         !
                    read(CHAIN_LENGTH,*,iostat=EOF2) CONFIG_ATOMID(i),    &              !
                                                     CONFIG_ATOM_TYPE(i), &              !
                                                     CONFIG_QI(i),        &              !
                                                     CONFIG_RI(1:3,i);                   !
                                                                                         !
                else if ( ATOMS_FLAG == 13 ) then;                                       !
                                                                                         ! 
                    read(CHAIN_LENGTH,*,iostat=EOF2) CONFIG_ATOMID(i),     &
                                                     CONFIG_MOLECULEID(i), &
                                                     CONFIG_ATOM_TYPE(i),  &
                                                     CONFIG_QI(i),         &
                                                     CONFIG_RI(1:3,i);

                    if ( CONFIG_ATOMID(i) > FLAG_SORT_MOLEC ) FLAG_SORT_MOLEC = CONFIG_ATOMID(i);

                else

                    read(CHAIN_LENGTH,*,iostat=EOF2) CONFIG_ATOMID(i), &
                                                     CONFIG_MOLECULEID(i), &
                                                     CONFIG_ATOM_TYPE(i), &
                                                     CONFIG_QI(i), &
                                                     CONFIG_RI(1:3,i);
 
                    if ( EOF2 /= 0 ) then;

                         read(CHAIN_LENGTH,*) CONFIG_ATOMID(i),    &
                                              CONFIG_ATOM_TYPE(i), &
                                              CONFIG_QI(i),        &
                                              CONFIG_RI(1:3,i);

                    end if

                end if

                do j = 1, NTYPE_ATOM
                    if ( IATOM_TYPE == j ) then
                        CONFIG_NAT(i) = ATOM_LABEL(j);
                        EXIT;
                    end if
                end do
            end do

            write(icanal,'(a70)') '| Atoms were read'//REPEAT(' ',52)//'|';              ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'Velocities') > 0 ) then;                           !

           read(2,*);

           do i = 1, NATOM;

               read(2,*) IOSEF1, CONFIG_VI(1:3,i);

           end do

           write(icanal,'(a70)') '| Velocities were read'//REPEAT(' ',47)//'|';

           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        else if ( INDEX(CHAIN_LENGTH,'Bonds') > 0 ) then;

            read(2,*);

            do i = 1, NBOND

                read(2,'(a)') CHAIN_LENGTH;

                read(CHAIN_LENGTH,*) IOSEF1, BOND_TYPE(IOSEF1), BOND_ATOMID(1:2,IOSEF1);

!               if ( INDEX(CHAIN_LENGTH,'#') > 0 ) then;

!                   read(CHAIN_LENGTH,*) IOSEF2, IOSEF3, IOSEF4, IOSEF5, CHOSEF1, BOND_PCFF_TYPE(1:2,IOSEF1) !CHOSEF2, CHOSEF3;

!                   do j = 1, NATOM;

!                       if ( CONFIG_ATOMID(j) == BOND_ATOMID(1,IOSEF1) ) then;

!                           CONFIG_PCFF_TYPE(j) = TRIM(BOND_PCFF_TYPE(1,IOSEF1));

!                       end if

!                       if ( CONFIG_ATOMID(j) == BOND_ATOMID(2,IOSEF1) ) then;

!                           CONFIG_PCFF_TYPE(j) = TRIM(BOND_PCFF_TYPE(2,IOSEF1));

!                       end if

!                   end do

!               end if

            end do

            ibonds = 1;

            write(icanal,'(a70)') '| Bonds were read'//REPEAT(' ',52)//'|';              !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'Angles') > 0 ) then;                               !
                                                                                         ! 
            do i = 1, NANGLE;                                                            !
                                                                                         !
                read(2,*) IOSEF1, ANGLE_TYPE(IOSEF1), ANGLE_ATOMID(1:3,IOSEF1);          !
                                                                                         !
            end do                                                                       !
                                                                                         !
            iangles = 1;

            write(icanal,'(a70)') '| Angles were read'//REPEAT(' ',51)//'|';             !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'Dihedrals') > 0 ) then;                            !
                                                                                         !
            do i = 1, NDIHEDRAL;

                read(2,*) IOSEF1, DIHEDRAL_TYPE(IOSEF1), DIHEDRAL_ATOMID(1:4,IOSEF1);

            end do

            idihedrals = 1;

            write(icanal,'(a70)') '| Dihedrals were read'//REPEAT(' ',48)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
         else if ( INDEX(CHAIN_LENGTH,'Impropers') > 0 ) then;                           !
                                                                                         ! 
            IOSEF2 = 0;

            do i = 1, NIMPROPER;

                read(2,*) IOSEF1, IMPROPER_TYPE(IOSEF1), IMPROPER_ATOMID(1:4,IOSEF1);

                if ( IMPROPER_TYPE(IOSEF1) > IOSEF2 ) IOSEF2 = IMPROPER_TYPE(IOSEF1);

            end do

            iimpropers = 1;

            write(icanal,'(a11,i8,a51)') '| IOSEF2 : ', IOSEF2, REPEAT(' ',50)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            write(icanal,'(a70)') '| Impropers were read'//REPEAT(' ',48)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!           stop; !//////////////////////////////////////////////////////////////////////!

        end if

    end do

    close(2);

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| Unread parameters'//REPEAT(' ',50)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Pair potentials ############################################################################
                                                                                         !
    if ( ( IPOTENTIAL_CLASS2 == 0 ) .AND. ( NINPUT_PAIR_POTENTIAL > 0 ) ) then;          !
                                                                                         !
        write(icanal,'(a70)') '| Pair potential parameters were read in the '// &        !
                              'lammps input file'//                             &        !
                              REPEAT(' ',7)//'|';                                        !
                                                                                         !
        IPOTENTIAL_CLASS2 = 1;                                                           !
                                                                                         !
        do i = 1, NTYPE_ATOM;                                                            !
                                                                                         !
            do j = 1, NINPUT_PAIR_POTENTIAL;                                             !
                                                                                         !
                if ( PAIR_POTENTIAL_SPECIES(1,j) /= i ) CYCLE;                           !
                                                                                         !
                if ( PAIR_POTENTIAL_SPECIES(2,j) /= i ) CYCLE;                           !
                                                                                         !
                POTENTIAL_CLASS2(1:2,i) = INPUT_PAIR_POTENTIAL(1:2,j);                   !
                                                                                         !
                POTENTIAL_CLASS2_CHTYPE = TRIM(PAIR_POTENTIAL_CHTYPE(j));                !
                                                                                         !
            end do                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ##### Bond coefficients ########################################################################
                                                                                         !
    if ( ibonds == 1 .AND. ibond_coeffs == 0 ) then;

        write(icanal,'(a70)') '| Bond Coefficients were not found '//REPEAT(' ',34)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        write(icanal,'(a70)') '| CHECK #1 : READ WITHOUT # class2 '//REPEAT(' ',34)//'|';

        ibond_coeffs = 0;

        open(2,file=TRIM(CHEXT));

        EOF = 0;

        do
            read(2,'(a)',iostat=EOF) CHAIN_LENGTH;

            if ( EOF /= 0 ) EXIT;

            if ( INDEX(CHAIN_LENGTH,'Bond Coeffs') > 0 ) then;

                if ( INDEX(CHAIN_LENGTH,'BondBond Coeffs') == 0 ) then;

                    read(2,*);

                    do i = 1, NTYPE_BOND;

                        read(2,*) IOSEF1, BOND_COEFFS(1:4,IOSEF1);

                    end do

                    ibond_coeffs = 1;

                end if

            end if

        end do

        close(2);

        if ( ibond_coeffs == 0 ) then;
 
            write(icanal,'(a70)') '| CHECK #1 HAS FAILED '//REPEAT(' ',47)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            write(icanal,'(a70)') '| CHECK #2 : READ FROM THE pcff.frc FILE '//REPEAT(' ',38)//'|';

            inquire(FILE='../pcff.frc',EXIST=PROBE1);

            if ( PROBE1 .EQV. .TRUE. ) then;

                PCFF_EQUIVALENCE(1:6,1:134) = 'XXX';

                open(2,file='../pcff.frc');

                EOF = 0;

                do 
                    read(2,'(a)',iostat=EOF) CHAIN_LENGTH;

                    if ( EOF /= 0 ) EXIT;

                    if ( INDEX(CHAIN_LENGTH,'#equivalence          cff91') > 0 ) then;
                        read(2,*);
                        read(2,*);
                        read(2,*);
                        read(2,*);
                        read(2,*);

                        do k = 1, 134
                            read(2,*) ROSEF1, IOSEF1, PCFF_EQUIVALENCE(1:6,k);
                        end do

                    else if ( INDEX(CHAIN_LENGTH,'#quartic_bond         cff91') > 0 ) then;
                        read(2,*);
                        read(2,*);
                        read(2,*);
                        read(2,*);
                        read(2,*);

                        do k = 1, 127
                            read(2,'(a)',iostat=EOF) CHAIN_LENGTH;

                            if ( EOF /= 0 ) EXIT;

                            if ( INDEX(CHAIN_LENGTH,'#') > 0 ) then;
                                EXIT;
                            else
                                read(CHAIN_LENGTH,*) ROSEF1, IOSEF1, &
                                                     CHOSEF1, CHOSEF2, TAB_ROSEF(1:4); !  ROSEF1, ROSEF2, ROSEF3, ROSEF4; 

                                do i = 1, NTYPE_BOND;


                                    do j = 1, 134
                                        if ( TRIM(BOND_PCFF_TYPE(1,i)) == TRIM(PCFF_EQUIVALENCE(1,j)) ) then
                                            BOND_PCFF_EQUIVALENCE(1,i) = PCFF_EQUIVALENCE(3,j);
                                        end if
                                    end do

                                    do j = 1, 134
                                        if ( TRIM(BOND_PCFF_TYPE(2,i)) == TRIM(PCFF_EQUIVALENCE(1,j)) ) then
                                            BOND_PCFF_EQUIVALENCE(2,i) = PCFF_EQUIVALENCE(3,j);
                                        end if
                                    end do

                                    if ( TRIM(CHOSEF1) == TRIM(BOND_PCFF_EQUIVALENCE(1,i)) .AND. &
                                         TRIM(CHOSEF2) == TRIM(BOND_PCFF_EQUIVALENCE(2,i)) ) then
                                      
                                        if ( BOND_PCFF_VERSION(i) < ROSEF1 ) then
                                            write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';
                                            BOND_PCFF_VERSION(i) = ROSEF1;

                                            BOND_COEFFS(1:4,i) = TAB_ROSEF(1:4);
                                           
                                            write(icanal,'(a70)') '| 1 : '//REPEAT(' ',63)//'|';
                                            write(icanal,'(4a20)') BOND_PCFF_TYPE(1:2,i), CHOSEF1, CHOSEF2;
                                            write(icanal,*) BOND_PCFF_EQUIVALENCE(1:2,i);
                                            write(icanal,'(4f15.6)') BOND_COEFFS(1:4,i); 
                                        end if
                                    end if

                                    if ( TRIM(CHOSEF2) == TRIM(BOND_PCFF_EQUIVALENCE(1,i)) .AND. &
                                         TRIM(CHOSEF1) == TRIM(BOND_PCFF_EQUIVALENCE(2,i)) ) then

!                                       write(icanal,*) 'YYYYYYYYY';

                                        if ( BOND_PCFF_VERSION(i) < ROSEF1 ) then
                                            write(icanal,*) '----------------------------------Y';
                                            BOND_PCFF_VERSION(i) = ROSEF1;

                                            BOND_COEFFS(1:4,i) = TAB_ROSEF(1:4);

                                            write(icanal,*) '2';
                                            write(icanal,'(4a20)') BOND_PCFF_TYPE(1:2,i), CHOSEF1, CHOSEF2;
                                            write(icanal,*) BOND_PCFF_EQUIVALENCE(1:2,i);
                                            write(icanal,'(4f15.6)') BOND_COEFFS(1:4,i);
                                        end if
                                    end if

                                end do
                            end if
                        end do
                    end if
                end do

                close(2);
            end if

        end if

    end if

!   ##### Angle coefficients ########################################################################
                                                                                         !
    if ( iangles == 1 .AND. iangle_coeffs == 0 ) then;                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        write(icanal,'(a70)') '| Angle Coefficients were not found '//REPEAT(' ',33)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        write(icanal,'(a70)') '| CHECK #1 : READ WITHOUT # class2 '//REPEAT(' ',34)//'|';

        iangle_coeffs = 0;

        open(2,file=TRIM(CHEXT));

        EOF = 0;

        do
            read(2,'(a)',iostat=EOF) CHAIN_LENGTH;

            if ( EOF /= 0 ) EXIT;

            if ( INDEX(CHAIN_LENGTH,'Angle Coeffs') > 0 ) then;
                if ( INDEX(CHAIN_LENGTH,'BondAngle Coeffs') == 0 .AND. &
                     INDEX(CHAIN_LENGTH,'AngleAngle Coeffs') == 0 ) then;
                    read(2,*);
                    do i = 1, NTYPE_ANGLE;
                        read(2,*) IOSEF1, ANGLE_COEFFS(1:4,IOSEF1);
                    end do

                    iangle_coeffs = 1;
                end if
            end if
        end do

        close(2);

        if ( iangle_coeffs == 0 ) then;

            write(icanal,'(a70)') '| CHECK #1 HAS FAILED '//REPEAT(' ',47)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            write(icanal,'(a70)') '| CHECK #2 : READ FROM THE pcff.frc FILE '//REPEAT(' ',38)//'|';

            write(icanal,*) 'NOT IMPLEMENTED: STOP!';

!           stop; !//////////////////////////////////////////////////////////////////////!

        end if

    end if

!   ##### Improper coefficients ####################################################################

    if ( iimpropers == 1 .AND. iimproper_coeffs == 0 ) then;

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        write(icanal,'(a70)') '| Improper Coefficients were not found '//REPEAT(' ',30)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        write(icanal,'(a70)') '| CHECK #1 : READ WITHOUT # class2 '//REPEAT(' ',34)//'|';

        iimproper_coeffs = 0;

        open(2,file=TRIM(CHEXT));

        EOF = 0;

        do
            read(2,'(a)',iostat=EOF) CHAIN_LENGTH;

            if ( EOF /= 0 ) EXIT;

            if ( INDEX(CHAIN_LENGTH,'Improper Coeffs') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_IMPROPER;

                    read(2,*) IOSEF1, IMPROPER_COEFFS(1:2,IOSEF1);

                end do

                iimproper_coeffs = 1;

            end if

        end do

        close(2);

        if ( iimproper_coeffs == 0 ) then;

            write(icanal,'(a70)') '| CHECK #1 HAS FAILED '//REPEAT(' ',47)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            write(icanal,'(a70)') '| CHECK #2 : READ FROM THE pcff.frc FILE '//REPEAT(' ',38)//'|';

            write(icanal,*) 'NOT IMPLEMENTED: STOP!';

            stop;

        end if

    end if

!   ##### AngleTorsion coefficients ################################################################
                                                                                         !
    if ( ( idihedrals           == 1 ) .AND.     &                                       ! 
         ( iangletorsion_coeffs == 0 ) .AND.     &                                       !
         ( NPARAM_ANGLES        == 4 )   ) then;                                         !

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        write(icanal,'(a70)') '| AngleTorsion Coefficients were not found '// &
                              REPEAT(' ',26)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        write(icanal,'(a70)') '| CHECK #1 : READ WITHOUT # class2 '// &
                              REPEAT(' ',34)//'|';

        iangletorsion_coeffs = 0;

        open(2,file=TRIM(CHEXT));

        EOF = 0;

        do
            read(2,'(a)',iostat=EOF) CHAIN_LENGTH;

            if ( EOF /= 0 ) EXIT;

            if ( INDEX(CHAIN_LENGTH,'AngleTorsion Coeffs') > 0 ) then;

                if ( INDEX(CHAIN_LENGTH,'AngleAngleTorsion Coeffs') == 0 ) then;

                    read(2,*);

                    do i = 1, NTYPE_DIHEDRAL;

                        read(2,*) IOSEF1, ANGLETORSION_COEFFS(1:8,IOSEF1);

                    end do

                    iangletorsion_coeffs = 1;

                end if

            end if

        end do

        close(2);

        if ( iangletorsion_coeffs == 0 ) then;

            write(icanal,'(a70)') '| CHECK #1 HAS FAILED '//REPEAT(' ',47)//'|';

            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            write(icanal,'(a70)') '| CHECK #2 : READ FROM THE pcff.frc FILE '//REPEAT(' ',38)//'|';

            write(icanal,*) 'NOT IMPLEMENTED: STOP!';

            stop;

        end if

    end if

!   ### Closing the routine ########################################################################

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

!   stop; !//////////////////////////////////////////////////////////////////////////////!

end subroutine READ_LAMMPS_CONFIG
