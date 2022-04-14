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

subroutine RANDOM_INSERTION_MOLECULE(icanal,NFILE_LIBRARY,             &
                                            NMOLECULE_INSERTION,       &
                                            CHNAME_FILE_LIBRARY,       &
                                            NATOM_FILE_LIBRARY,        &
                                            NBOND_LIBRARY,             &
                                            NANGLE_LIBRARY,            &
                                            NDIHEDRAL_LIBRARY,         &
                                            NIMPROPER_LIBRARY,         &
                                            ATOM_LABEL_LIBRARY,        &
                                            CONFIG_NAT_LIBRARY,        &
                                            CONFIG_QI_LIBRARY,         &
                                            CONFIG_RI_LIBRARY,         &
                                            POTENTIAL_CLASS2_LIBRARY,  &
                                            BOND_COEFFS_LIBRARY,       &
                                            ANGLE_COEFFS_LIBRARY,      &
                                            BONDBOND_COEFFS_LIBRARY,   &
                                            BONDANGLE_COEFFS_LIBRARY,  &
                                            IMPROPER_COEFFS_LIBRARY,   &
                                            ANGLEANGLE_COEFFS_LIBRARY, &
                                            PASSA,                     &
                                            PASSB,                     &
                                            NTYPE_ATOM_LIBRARY,        &
                                            NTYPE_BOND_LIBRARY,        &
                                            NTYPE_ANGLE_LIBRARY,       &
                                            NTYPE_IMPROPER_LIBRARY,    &
                                            ATOM_MASSES_LIBRARY,       &
                                            CONFIG_ATOMID_LIBRARY,     &
                                            CONFIG_ATOM_TYPE_LIBRARY,  &
                                            BOND_TYPE_LIBRARY,         &
                                            ANGLE_TYPE_LIBRARY,        &
                                            IMPROPER_TYPE_LIBRARY,     &
                                            BOND_ATOMID_LIBRARY,       &
                                            ANGLE_ATOMID_LIBRARY,      &
                                            IMPROPER_ATOMID_LIBRARY,   &
                                            IPOTENTIAL_CLASS2_LIBRARY, &
                                            NATOM,                     &
                                            NBOND,                     &
                                            NANGLE,                    &
                                            NDIHEDRAL,                 &                        
                                            NIMPROPER,                 &
                                            NTYPE_ATOM,                &
                                            NTYPE_BOND,                &
                                            NTYPE_ANGLE,               &
                                            NTYPE_DIHEDRAL,            &
                                            NTYPE_IMPROPER,            &
                                            ATOM_MASSE,                &
                                            ATOM_LABEL,                &
                                            CONFIG_NAT,                &
                                            CONFIG_QI,                 &
                                            CONFIG_RI,                 &
                                            POTENTIAL_CLASS2,          &
                                            BOND_COEFFS,               &
                                            ANGLE_COEFFS,              &
                                            BONDBOND_COEFFS,           &
                                            BONDANGLE_COEFFS,          &
                                            IMPROPER_COEFFS,           &
                                            ANGLEANGLE_COEFFS,         &
                                            CONFIG_ATOMID,             &
                                            CONFIG_MOLECULEID,         &
                                            CONFIG_ATOM_TYPE,          &
                                            BOND_TYPE,                 &
                                            ANGLE_TYPE,                &
                                            IMPROPER_TYPE,             &
                                            BOND_ATOMID,               &
                                            ANGLE_ATOMID,              &
                                            IMPROPER_ATOMID,           &
                                            IPOTENTIAL_CLASS2)

!   ************************************************************************************************
!   **                                  READ THE INPUT FILE                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                    : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   ** NFILE_LIBRARY             : NUMBER OF FILES COMING FROM THE LIBRARY                            **
!   ** NMOLECULE_INSERTION       : NUMBER OF INSERTED MOLECULES FOR A GIVE FILE                       **
!   ** CHNAME_FILE_LIBRARY       : FILE NAME CONTAINING THE MOLECULAR CONFIGURATION OF THE MOLECULE   **
!   ** NATOM_FILE_LIBRARY        : NUMBER OF ATOMS IN THE FILES FROM THE LIBRARY              **
!   ** NBOND_LIBRARY             : NUMBER OF BONDS IN THE FILES FROM THE LIBRARY              **
!   ** NANGLE_LIBRARY            : NUMBER OF ANGLES IN THE FILES FROM THE LIBRARY             **
!   ** NDIHEDRAL_LIBRARY         : NUMBER OF DIHEDRAL ANGLES IN THE FILES FROM THE LIBRARY    **
!   ** NIMPROPER_LIBRARY         : NUMBER OF IMPROPERS IN THE FILES FROM THE LIBRARY          **
!   ** ATOM_LABEL_LIBRARY        :  **
!   ** CONFIG_NAT_LIBRARY        :  **
!   ** CONFIG_QI_LIBRARY         : CHARGES OF ATOMS IN THE FILES FROM THE LIBRARY  **
!   ** CONFIG_RI_LIBRARY         : COORDINATES OF ATOMS IN THE FILES FROM THE LIBRARY **
!   ** POTENTIAL_CLASS2_LIBRARY  : LENNARD-JONES POTENTIAL PARAMETERS IN THE FILES FROM THE LIBRARY **
!   ** BOND_COEFFS_LIBRARY       :
!   ** ANGLE_COEFFS_LIBRARY      :
!   ** BONDBOND_COEFFS_LIBRARY   :
!   ** BONDANGLE_COEFFS_LIBRARY  : 
!   ** IMPROPER_COEFFS_LIBRARY   :
!   ** ANGLEANGLE_COEFFS_LIBRARY :
!   ** PASSA                     :  **
!   ** PASSB                     :  **
!   ** NTYPE_ATOM_LIBRARY        : NUMBER OF ATOM TYPES IN THE LIBRARY FILE **
!   ** NTYPE_BOND_LIBRARY        : NUMBER OF BOND TYPES IN THE LIBRARY FILE **
!   ** NTYPE_ANGLE_LIBRARY       : NUMBER OF ANGLE TYPES IN THE LIBRARY FILE **
!   ** NTYPE_IMPROPER_LIBRARY    : NUMBER OF IMPROPER TYPES IN THE LIBRARY FILE **
!   ** ATOM_MASSES_LIBRARY       :
!   ** CONFIG_ATOMID_LIBRARY     :
!   ** CONFIG_ATOM_TYPE_LIBRARY  : 
!   ** BOND_TYPE_LIBRARY         :
!   ** ANGLE_TYPE_LIBRARY        :
!   ** IMPROPER_TYPE_LIBRARY     :
!   ** BOND_ATOMID_LIBRARY       :
!   ** ANGLE_ATOMID_LIBRARY      :
!   ** IMPROPER_ATOMID_LIBRARY   :
!   ** IPOTENTIAL_CLASS2_LIBRARY : FLAG
!   **
!   ** NATOM                    : NUMBER OF ATOMS IN THE FINAL CONFIGURATION  ** 
!   ** NBOND                    : NUMBER OF BONDS IN THE FINAL CONFIGURATION   **
!   ** NANGLE                   : NUMBER OF ANGLES IN THE FINAL CONFIGURATION   **
!   ** NDIHEDRAL                : NUMBER OF DIHEDRALS IN THE FINAL CONFIGURATION **                       
!   ** NIMPROPER                : NUMBER OF IMPROPERS IN THE FINAL CONFIGURATION **
!   ** NTYPE_ATOM               : NUMBER OF ATOM TYPES IN THE FINAL MOLECULAR CONFIGURATION          **
!   ** NTYPE_BOND               : NUMBER OF BOND TYPES IN THE FINAL MOLECULAR CONFIGURATION          **
!   ** NTYPE_ANGLE              : NUMBER OF ANGLE TYPES IN THE FINAL MOLECULAR CONFIGURATION **
!   ** NTYPE_DIHEDRAL           : NUMBER OF DIHEDRAL TYPES IN THE FINAL MOLECULAR CONFIGURATION **
!   ** NTYPE_IMPROPER           : NUMBER OF IMPROPER TYPES IN THE FINAL MOLECULAR CONFIGURATION **
!   ** ATOM_MASSE               : MASS OF ATOMS IN THE MOLECULAR CONFIGURATION                       **
!   ** ATOM_LABEL               : LIST OF ATOM LABEL IN THE FINAL MOLECULAR CONFIGURATION            **
!   ** CONFIG_NAT               :
!   ** CONFIG_QI                : CHARGES OFR ATOMS IN THE FINAL MOLECULAR CONFIGURATION
!   ** CONFIG_RI                : COORDINATES OF ATOMS IN THE FINAL MOLECULAR CONFIGURATION
!   ** POTENTIAL_CLASS2         : LENNARD-JONES POTENTIAL PARAMETERS IN THE FINAL MOLECULAR CONFIGURATION **
!   ** BOND_COEFFS              :
!   ** ANGLE_COEFFS             : 
!   ** BONDBOND_COEFFS          :
!   ** BONDANGLE_COEFFS         :
!   ** IMPROPER_COEFFS          :
!   ** ANGLEANGLE_COEFFS        :
!   ** CONFIG_ATOMID            :
!   ** CONFIG_MOLECULEID        :
!   ** CONFIG_ATOM_TYPE         : 
!   ** BOND_TYPE                :
!   ** ANGLE_TYPE               :
!   ** IMPROPER_TYPE            :
!   ** BOND_ATOMID              :
!   ** ANGLE_ATOMID             :
!   ** IMPROPER_ATOMID          :
!   ** IPOTENTIAL_CLASS2        : FLAG
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_physical_constants;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: NFILE_LIBRARY;

    integer (kind=4), dimension(1:100), intent(in) :: NMOLECULE_INSERTION;

    character (len=150), dimension(1:100), intent(in) :: CHNAME_FILE_LIBRARY;

    integer (kind=4), dimension(1:100), intent(in) :: NATOM_FILE_LIBRARY, NTYPE_ATOM_LIBRARY;

    integer (kind=4), dimension(1:100), intent(in) :: NBOND_LIBRARY, NTYPE_BOND_LIBRARY;

    integer (kind=4), dimension(1:100), intent(in) :: NANGLE_LIBRARY, NTYPE_ANGLE_LIBRARY;

    integer (kind=4), dimension(1:100), intent(in) :: NDIHEDRAL_LIBRARY;

    integer (kind=4), dimension(1:100), intent(in) :: NIMPROPER_LIBRARY, NTYPE_IMPROPER_LIBRARY;

    integer (kind=4), dimension(1:100), intent(in) :: IPOTENTIAL_CLASS2_LIBRARY;

    integer (kind=4), dimension(1:1000,1:100), intent(in) :: CONFIG_ATOMID_LIBRARY, CONFIG_ATOM_TYPE_LIBRARY;

    integer (kind=4), dimension(1:1000,1:100), intent(in) :: BOND_TYPE_LIBRARY, ANGLE_TYPE_LIBRARY, IMPROPER_TYPE_LIBRARY;

    integer (kind=4), dimension(1:2,1:1000,1:100), intent(in) :: BOND_ATOMID_LIBRARY;

    integer (kind=4), dimension(1:3,1:1000,1:100), intent(in) :: ANGLE_ATOMID_LIBRARY;

    integer (kind=4), dimension(1:4,1:1000,1:100), intent(in) :: IMPROPER_ATOMID_LIBRARY;

    real (kind=8), dimension(1:20,1:100), intent(in) :: ATOM_MASSES_LIBRARY;

    real (kind=8), dimension(1:1000,1:100) ::  CONFIG_QI_LIBRARY;

    real (kind=8), dimension(1:3,1:1000,1:100), intent(in) :: CONFIG_RI_LIBRARY;

    real (kind=8), dimension(1:2,1:100,1:100), intent(in) :: POTENTIAL_CLASS2_LIBRARY;

    real (kind=8), dimension(1:2,1:10,1:100), intent(in) :: IMPROPER_COEFFS_LIBRARY;

    real (kind=8), dimension(1:3,1:10,1:100), intent(in) :: BONDBOND_COEFFS_LIBRARY;

    real (kind=8), dimension(1:4,1:10,1:100), intent(in) :: BOND_COEFFS_LIBRARY, ANGLE_COEFFS_LIBRARY, BONDANGLE_COEFFS_LIBRARY;

    real (kind=8), dimension(1:6,1:10,1:100), intent(in) :: ANGLEANGLE_COEFFS_LIBRARY;

    character (len=20), dimension(1:20,1:100), intent(in) :: ATOM_LABEL_LIBRARY;

    character (len=20), dimension(1:1000,1:100), intent(in) :: CONFIG_NAT_LIBRARY;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4), intent(out) :: IPOTENTIAL_CLASS2;

    integer (kind=4), intent(out) :: NATOM, NBOND, NANGLE, NDIHEDRAL, NIMPROPER;

    integer (kind=4), intent(out) :: NTYPE_ATOM, NTYPE_BOND, NTYPE_ANGLE, NTYPE_DIHEDRAL, NTYPE_IMPROPER;

    integer (kind=4), dimension(1:100000), intent(out) :: CONFIG_ATOMID, CONFIG_ATOM_TYPE;

    integer (kind=4), dimension(1:100000), intent(out) :: CONFIG_MOLECULEID;

    integer (kind=4), dimension(1:100000), intent(out) :: BOND_TYPE, ANGLE_TYPE, IMPROPER_TYPE;

    integer (kind=4), dimension(1:2,1:100000), intent(out) :: BOND_ATOMID;

    integer (kind=4), dimension(1:3,1:100000), intent(out) :: ANGLE_ATOMID;

    integer (kind=4), dimension(1:4,1:100000), intent(out) :: IMPROPER_ATOMID;

    real (kind=8), dimension(1:20), intent(out) :: ATOM_MASSE;

    real (kind=8), dimension(1:100000), intent(out) :: CONFIG_QI;

    real (kind=8), dimension(1:3,1:100000), intent(out) :: CONFIG_RI;

    character (len=20), dimension(1:20), intent(out) :: ATOM_LABEL;

    character (len=20), dimension(1:100000), intent(out) :: CONFIG_NAT;

    real (kind=8), dimension(1:2,1:100), intent(out) :: POTENTIAL_CLASS2;

    real (kind=8), dimension(1:2,1:10), intent(out) :: IMPROPER_COEFFS;

    real (kind=8), dimension(1:3,1:10), intent(out) :: BONDBOND_COEFFS;

    real (kind=8), dimension(1:4,1:10), intent(out) :: BOND_COEFFS, ANGLE_COEFFS, BONDANGLE_COEFFS;

    real (kind=8), dimension(1:6,1:10), intent(out) :: ANGLEANGLE_COEFFS;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m;

    integer (kind=4) :: IFOUND;

    integer (kind=4) :: IINSERTED, IATOM, IOVERLAPPING;

    integer (kind=4) :: IATOM_TYPE_TRANSLATE, IMOLECULEID;

    integer (kind=4) :: IBOND_TYPE_TRANSLATE, IANGLE_TYPE_TRANSLATE, IIMPROPER_TYPE_TRANSLATE;

    real (kind=8) :: grnd;

    real (kind=8) :: SUM_MASSES, INV_SUM_MASSES;

    real (kind=8), dimension(1:3) :: RINO, RI, RAND_ANGLES, RG, DRIJ;

    real (kind=8), dimension(1:3,1:3) :: MATROT;

!   ************************************************************************************************

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************

    CHTITLE = 'RANDOM INSERTION MOLECULE';

    ILENGTH_TITLE = LEN_TRIM(CHTITLE);

    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));

    CONFIG_QI(1:100000)     = 0.0d0;
    CONFIG_RI(1:3,1:100000) = 0.0d0;

    CONFIG_NAT(1:100000) = 'XXXX';

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

    ATOM_MASSE(1:20) = 0.0d0;

    ATOM_LABEL(1:20) = 'XXX';

    CONFIG_ATOMID(1:100000)     = 0;
    CONFIG_MOLECULEID(1:100000) = 0; 

    BOND_ATOMID(1:2,1:100000)     = 0;
    ANGLE_ATOMID(1:3,1:100000)    = 0;
    IMPROPER_ATOMID(1:4,1:100000) = 0;

    CONFIG_ATOM_TYPE(1:100000) = 0;
    BOND_TYPE(1:100000)        = 0;
    ANGLE_TYPE(1:100000)       = 0;
    IMPROPER_TYPE(1:100000)    = 0;

    POTENTIAL_CLASS2(1:2,1:100) = 0.0d0;

    BOND_COEFFS(1:4,1:10)       = 0.0d0;
    ANGLE_COEFFS(1:4,1:10)      = 0.0d0;
    BONDBOND_COEFFS(1:3,1:10)   = 0.0d0;
    BONDANGLE_COEFFS(1:4,1:10)  = 0.0d0;
    IMPROPER_COEFFS(1:2,1:10)   = 0.0d0;
    ANGLEANGLE_COEFFS(1:6,1:10) = 0.0d0;

    IPOTENTIAL_CLASS2 = 0;

    IATOM_TYPE_TRANSLATE     = 0;
    IBOND_TYPE_TRANSLATE     = 0;
    IANGLE_TYPE_TRANSLATE    = 0;
    IIMPROPER_TYPE_TRANSLATE = 0;

    IMOLECULEID = 0;

    do i = 1, NFILE_LIBRARY;                  ! LOOP OVER THE MOLECULE TYPE READ IN THE LIBRARY 

        SUM_MASSES = 0.0d0;                   ! INITIALIZATION OF THE SUM OF MASSES

        RG(1:3) = 0.0d0;                      ! INITIALIZATION OF THE CENTER OF MASS

        write(icanal,*) '| NATOM_FILE_LIBRARY : ', NATOM_FILE_LIBRARY(i);
        write(icanal,*) '| NTYPE_ATOM_LIBRARY : ', NTYPE_ATOM_LIBRARY(i);

        do j = 1, NATOM_FILE_LIBRARY(i);      ! LOOP OVER THE NUMBER OF ATOMS IN THE MOLECULE TYPE I

            do k = 1, NTYPE_ATOM_LIBRARY(i);  ! LOOP OVER THE NUMBER OF ATOM TYPES IN THE MOLECULE I

                if ( TRIM(CONFIG_NAT_LIBRARY(j,i)) == TRIM(ATOM_LABEL_LIBRARY(k,i)) ) then
                    RG(1:3) = RG(1:3) + ATOM_MASSES_LIBRARY(k,i) * CONFIG_RI_LIBRARY(1:3,j,i);
                    SUM_MASSES = SUM_MASSES + ATOM_MASSES_LIBRARY(k,i);
                    EXIT;
                end if
            end do

        end do

        INV_SUM_MASSES = 1.0d0 / SUM_MASSES;

        write(icanal,*) '| SUM MASSES : ', SUM_MASSES, INV_SUM_MASSES;

        RG(1:3) = RG(1:3) * INV_SUM_MASSES;

        write(icanal,*) '| RG : ', RG(1:3);
        write(icanal,*);

        IINSERTED = 0;

        do while ( IINSERTED < NMOLECULE_INSERTION(i) )
            RINO(1) = grnd() - 0.5d0;
            RINO(2) = grnd() - 0.5d0;
            RINO(3) = grnd() - 0.5d0;

            RI(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));

            do k = 1, 1000
                RAND_ANGLES(1) = TWOPI * grnd();
                RAND_ANGLES(2) = TWOPI * grnd();
                RAND_ANGLES(3) = TWOPI * grnd();

                write(icanal,*) RAND_ANGLES(1:3);

                call SET_ROTATION_MATRIX(RAND_ANGLES(1:3),MATROT(1:3,1:3));

                IATOM = NATOM;

                do j = 1, NATOM_FILE_LIBRARY(i);      ! LOOP OVER THE NUMBER OF ATOMS IN THE MOLECULE TYPE I
                    IATOM = IATOM + 1;
 
                    DRIJ(1:3) = CONFIG_RI_LIBRARY(1:3,j,i) - RG(1:3);

                    call APPLY_PBC(DRIJ(1:3),DRIJ(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));

                    DRIJ(1:3) = MATMUL(MATROT(1:3,1:3),DRIJ(1:3));

                    CONFIG_RI(1:3,IATOM) = RI(1:3) + DRIJ(1:3);    
                    CONFIG_NAT(IATOM)    = CONFIG_NAT_LIBRARY(j,i);
                    CONFIG_QI(IATOM)     = CONFIG_QI_LIBRARY(j,i);

                    CONFIG_ATOMID(IATOM)    = NATOM + CONFIG_ATOMID_LIBRARY(j,i);
                    CONFIG_ATOM_TYPE(IATOM) = IATOM_TYPE_TRANSLATE + CONFIG_ATOM_TYPE_LIBRARY(j,i);
                end do

                call CHECK_OVERLAPPING_ATOMS(icanal,NATOM+1,                  &
                                                    IATOM,                    &
                                                    CONFIG_RI(1:3,1:100000),  &
                                                    PASSA(1:3,1:3),           &
                                                    PASSB(1:3,1:3),           &
                                                    1.5d0,                    &
                                                    IOVERLAPPING);

                if ( IOVERLAPPING == 0 ) EXIT; 
            end do

            if ( IOVERLAPPING == 0 ) then

                IMOLECULEID = IMOLECULEID + 1;                   ! UPDATE THE MOLECULE ID

                CONFIG_MOLECULEID(NATOM+1:IATOM) = IMOLECULEID; ! APPLY THE NEW MOLECULE ID 

                IOSEF1 = NATOM;
                IOSEF2 = NBOND     + 1;                     !
                IOSEF3 = NANGLE    + 1;

                NATOM     = IATOM;                            ! UPDATE THE NUMBER OF ATOMS IN THE FINAL CONFIGURATION
                NBOND     = NBOND     + NBOND_LIBRARY(i);     ! UPDATE THE NUMBER OF BONDS IN THE FINAL CONFIGURATION 
                NANGLE    = NANGLE    + NANGLE_LIBRARY(i);    ! UPDATE THE NUMBER OF ANGLES IN THE FINAL CONFIGURATION
                NDIHEDRAL = NDIHEDRAL + NDIHEDRAL_LIBRARY(i); ! UPDATE THE NUMBER OF DIHEDRALS IN THE FINAL CONFIGURATION

                IINSERTED = IINSERTED + 1;

                BOND_ATOMID(1:2,IOSEF2:NBOND)         = IOSEF1 + BOND_ATOMID_LIBRARY(1:2,1:NBOND_LIBRARY(i),i);
                ANGLE_ATOMID(1:3,IOSEF3:NANGLE)       = IOSEF1 + ANGLE_ATOMID_LIBRARY(1:3,1:NANGLE_LIBRARY(i),i);

                BOND_TYPE(IOSEF2:NANGLE)        = IBOND_TYPE_TRANSLATE     + BOND_TYPE_LIBRARY(1:NBOND_LIBRARY(i),i);
                ANGLE_TYPE(IOSEF3:NANGLE)       = IANGLE_TYPE_TRANSLATE    + ANGLE_TYPE_LIBRARY(1:NANGLE_LIBRARY(i),i);

!               do j = 1, NANGLE_LIBRARY(i)
!                   write(icanal,'(4i8)') ANGLE_TYPE_LIBRARY(j,i), ANGLE_ATOMID_LIBRARY(1:3,j,i);
!               end do

                if ( NIMPROPER_LIBRARY(i) > 0 ) then
                    IOSEF4 = NIMPROPER + 1;
                    NIMPROPER = NIMPROPER + NIMPROPER_LIBRARY(i); ! UPDATE THE NUMBER OF IMPROPERS IN THE FINAL CONFIGURATION
                    IMPROPER_ATOMID(1:4,IOSEF4:NIMPROPER) = IOSEF1 + IMPROPER_ATOMID_LIBRARY(1:4,1:NIMPROPER_LIBRARY(i),i);
                    IMPROPER_TYPE(IOSEF4:NIMPROPER) = IIMPROPER_TYPE_TRANSLATE + IMPROPER_TYPE_LIBRARY(1:NIMPROPER_LIBRARY(i),i);
                end if 

                write(icanal,*) '| IINSERTED : ', IINSERTED, 100.0d0 * REAL(IINSERTED) / REAL(NMOLECULE_INSERTION(i));
                write(icanal,*);
            end if
        end do

        IOSEF1 = NTYPE_ATOM     + 1;                                                     ! UPDATE THE NUMBER OF ATOMS IN THE FINAL MOLECULAR CONFIGURATION
        IOSEF2 = NTYPE_BOND     + 1;                     !
        IOSEF3 = NTYPE_ANGLE    + 1;
        IOSEF4 = NTYPE_IMPROPER + 1;

        NTYPE_ATOM     = NTYPE_ATOM     + NTYPE_ATOM_LIBRARY(i);                         ! UPDATE THE NUMBER OF ATOM TYPES     
        NTYPE_BOND     = NTYPE_BOND     + NTYPE_BOND_LIBRARY(i);                         ! UPDATE THE NUMBER OF BOND TYPES 
        NTYPE_ANGLE    = NTYPE_ANGLE    + NTYPE_ANGLE_LIBRARY(i);                        ! UPDATE THE NUMBER OF ANGLE TYPES  
        NTYPE_IMPROPER = NTYPE_IMPROPER + NTYPE_IMPROPER_LIBRARY(i);                     ! UPDATE THE NUMBER OF IMPROPER TYPES

        ATOM_LABEL(IOSEF1:NTYPE_ATOM) = ATOM_LABEL_LIBRARY(1:NTYPE_ATOM_LIBRARY(i),i);   ! UPDATE THE LIST OF LABEL OF ATOMS 

        ATOM_MASSE(IOSEF1:NTYPE_ATOM) = ATOM_MASSES_LIBRARY(1:NTYPE_ATOM_LIBRARY(i),i);  ! UPDATE THE LIST OF MASSES 

        POTENTIAL_CLASS2(1:2,IOSEF1:NTYPE_ATOM) = POTENTIAL_CLASS2_LIBRARY(1:2,1:NTYPE_ATOM_LIBRARY(i),i);           ! UPDATE LENNARD-JONES-POTENTIAL PARAMETERS

        BOND_COEFFS(1:4,IOSEF2:NTYPE_BOND)           = BOND_COEFFS_LIBRARY(1:4,1:NTYPE_BOND_LIBRARY(i),i);           ! 
        ANGLE_COEFFS(1:4,IOSEF3:NTYPE_ANGLE)         = ANGLE_COEFFS_LIBRARY(1:4,1:NTYPE_ANGLE_LIBRARY(i),i);         !
        BONDBOND_COEFFS(1:3,IOSEF3:NTYPE_ANGLE)      = BONDBOND_COEFFS_LIBRARY(1:3,1:NTYPE_ANGLE_LIBRARY(i),i);      !
        BONDANGLE_COEFFS(1:4,IOSEF3:NTYPE_ANGLE)     = BONDANGLE_COEFFS_LIBRARY(1:4,1:NTYPE_ANGLE_LIBRARY(i),i);     !
        IMPROPER_COEFFS(1:2,IOSEF4:NTYPE_IMPROPER)   = IMPROPER_COEFFS_LIBRARY(1:2,1:NTYPE_IMPROPER_LIBRARY(i),i);   !
        ANGLEANGLE_COEFFS(1:6,IOSEF4:NTYPE_IMPROPER) = ANGLEANGLE_COEFFS_LIBRARY(1:6,1:NTYPE_IMPROPER_LIBRARY(i),i); !

        if ( IPOTENTIAL_CLASS2_LIBRARY(i) == 1 ) IPOTENTIAL_CLASS2 = 1; 

        IATOM_TYPE_TRANSLATE     = IATOM_TYPE_TRANSLATE     + NTYPE_ATOM_LIBRARY(i);
        IBOND_TYPE_TRANSLATE     = IBOND_TYPE_TRANSLATE     + NTYPE_BOND_LIBRARY(i);
        IANGLE_TYPE_TRANSLATE    = IANGLE_TYPE_TRANSLATE    + NTYPE_ANGLE_LIBRARY(i);
        IIMPROPER_TYPE_TRANSLATE = IIMPROPER_TYPE_TRANSLATE + NTYPE_IMPROPER_LIBRARY(i);
    end do

    write(icanal,*) '| NATOM     : ', NATOM;
    write(icanal,*) '| NBOND     : ', NBOND;
    write(icanal,*) '| NANGLE    : ', NANGLE;
    write(icanal,*) '| NDIHEDRAL : ', NDIHEDRAL;
    write(icanal,*) '| NIMPROPER : ', NIMPROPER;
    write(icanal,*);
    write(icanal,*) '| NTYPE_ATOM     : ', NTYPE_ATOM;
    write(icanal,*) '| NTYPE_BOND     : ', NTYPE_BOND;
    write(icanal,*) '| NTYPE_ANGLE    : ', NTYPE_ANGLE;
    write(icanal,*) '| NTYPE_DIHEDRAL : ', NTYPE_DIHEDRAL;
    write(icanal,*) '| NTYPE_IMPROPER : ', NTYPE_IMPROPER;

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine RANDOM_INSERTION_MOLECULE











