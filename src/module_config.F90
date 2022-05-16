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

module module_config

    use module_size_arrays;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: IPOTENTIAL_CLASS2;

    integer (kind=4) :: NATOM,                   &
                        NBOND,                   &
                        NANGLE,                  &
                        NDIHEDRAL,               &
                        NIMPROPER,               &
                        NRESIDUES,               &
                        NTYPE_ATOM,              &
                        NTYPE_ATOM_SQ,           &
                        NTYPE_BOND,              &
                        NTYPE_ANGLE,             &
                        NTYPE_DIHEDRAL,          &
                        NTYPE_IMPROPER,          &
                        NPARAM_BONDS,            &
                        NPARAM_ANGLES,           &
                        NPARAM_DIHEDRALS,        &
                        AMBER_NATYP,             &
                        AMBER_NLJ_INTER,         &
                        AMBER_NBONH,             &
                        AMBER_NBONA,             &
                        AMBER_NTHETH,            &
                        AMBER_NTHETA,            &
                        AMBER_NPHIH,             &
                        AMBER_NPHIA,             &
                        NPOTENTIAL_CLASS2_CROSS;

    integer (kind=4), allocatable, dimension(:) :: CONFIG_ATOM_TYPE,             &
                                                   BOND_TYPE,                    &
                                                   ANGLE_TYPE,                   &
                                                   IMPROPER_TYPE,                &
                                                   DIHEDRAL_TYPE,                &
                                                   CONFIG_ATOMID,                &
                                                   CONFIG_MOLECULEID,            &
                                                   CONFIG_ATOMIC_NUMBER,         &
                                                   CONFIG_NUMBER_EXCLUDED_ATOMS, &
                                                   NONBONDED_PARM_INDEX,         &
                                                   RESIDUE_POINTER,              &
                                                   BONDS_INC_HYDROGEN,           &
                                                   BONDS_WITHOUT_HYDROGEN,       &
                                                   ANGLES_INC_HYDROGEN,          &
                                                   ANGLES_WITHOUT_HYDROGEN,      &
                                                   DIHEDRALS_INC_HYDROGEN,       &
                                                   DIHEDRALS_WITHOUT_HYDROGEN;

    integer (kind=4), allocatable, dimension(:,:) :: BOND_ATOMID,                   &
                                                     POTENTIAL_CLASS2_CROSS_ATOMID;

    integer (kind=4), allocatable, dimension(:,:) :: ANGLE_ATOMID;

    integer (kind=4), allocatable, dimension(:,:) :: IMPROPER_ATOMID, &
                                                     DIHEDRAL_ATOMID;

    character (len=4), allocatable, dimension(:) :: CONFIG_PCFF_TYPE;

    character (len=4), allocatable, dimension(:,:) :: BOND_PCFF_TYPE, BOND_PCFF_EQUIVALENCE;

!   ************************************************************************************************

    real (kind=8), allocatable, dimension(:) :: ATOM_MASSE;

    real (kind=8), allocatable, dimension(:) :: CONFIG_QI,               &
                                                CONFIG_MASS,             &
                                                BOND_FORCE_CONSTANT,     &
                                                BOND_EQUIL_VALUE,        &
                                                ANGLE_FORCE_CONSTANT,    &
                                                ANGLE_EQUIL_VALUE,       &
                                                DIHEDRAL_FORCE_CONSTANT, &
                                                DIHEDRAL_PERIODICITY,    &
                                                DIHEDRAL_PHASE,          &
                                                SCEE_SCALE_FACTOR,       &
                                                SCNB_SCALE_FACTOR,       &
                                                SOLTY,                   &
                                                LENNARD_JONES_ACOEF,     &
                                                LENNARD_JONES_BCOEF;

    real (kind=8), allocatable, dimension(:,:) :: POTENTIAL_CLASS2, &
                                                  POTENTIAL_CLASS2_CROSS_VALUES;

    real (kind=8), allocatable, dimension(:,:) :: CONFIG_RI,  &
                                                  CONFIG_VI;

    real (kind=8), allocatable, dimension(:,:) :: IMPROPER_COEFFS;

    real (kind=8), allocatable, dimension(:,:) :: BONDBOND_COEFFS,           &
                                                  ANGLEANGLETORSION_COEFFS,  &
                                                  BONDBOND13_COEFFS;

    real (kind=8), allocatable, dimension(:,:) :: BOND_COEFFS,               &
                                                  ANGLE_COEFFS,              &
                                                  BONDANGLE_COEFFS,          &
                                                  MIDDLEBONDTORSION_COEFFS;

    real (kind=8), allocatable, dimension(:,:) :: ANGLEANGLE_COEFFS,         &
                                                  DIHEDRAL_COEFFS;

    real (kind=8), allocatable, dimension(:,:) :: ENDBONDTORSION_COEFFS,     &
                                                  ANGLETORSION_COEFFS;

!   ************************************************************************************************

    character (len=250) :: CH_ATOM_STYLE,           &
                           CH_BOND_STYLE,           &
                           CH_ANGLE_STYLE,          &
                           CH_DIHEDRAL_STYLE,       &
                           CH_IMPROPER_STYLE,       &
                           POTENTIAL_CLASS2_CHTYPE;

    character (len=20), allocatable, dimension(:) :: ATOM_LABEL;

    character (len=20), allocatable, dimension(:) :: CONFIG_NAT,    &
                                                     RESIDUE_LABEL;

end module module_config
