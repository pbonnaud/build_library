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

module module_library

    use module_size_arrays;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: NFILE_LIBRARY;

    integer (kind=4) :: IMAX_NATOM, IMAX_NBOND, IMAX_NANGLE, IMAX_NDIHEDRAL, IMAX_NIMPROPER;

    integer (kind=4) :: IMAX_NTYPE_NATOM,    &
                        IMAX_NTYPE_BOND,     &
                        IMAX_NTYPE_ANGLE,    &
                        IMAX_NTYPE_DIHEDRAL, &
                        IMAX_NTYPE_IMPROPER;

    integer (kind=4) :: NINPUT_PAIR_POTENTIAL;

    integer (kind=4) :: NMOLECULE_INSERTION,        &
                                                   IPOTENTIAL_CLASS2_LIBRARY,  &
                                                   NATOM_LIBRARY,         &
                                                   NBOND_LIBRARY,              &
                                                   NANGLE_LIBRARY,             &
                                                   NDIHEDRAL_LIBRARY,          &
                                                   NIMPROPER_LIBRARY,          &
                                                   NTYPE_ATOM_LIBRARY,         &
                                                   NTYPE_BOND_LIBRARY,         &
                                                   NTYPE_ANGLE_LIBRARY,        &
                                                   NTYPE_DIHEDRAL_LIBRARY,     &
                                                   NTYPE_IMPROPER_LIBRARY;

    integer (kind=4), dimension(1:2,1:1000) :: PAIR_POTENTIAL_SPECIES; 

    integer (kind=4), allocatable, dimension(:) :: CONFIG_ATOM_TYPE_LIBRARY,    &
                                                     BOND_TYPE_LIBRARY,           & 
                                                     ANGLE_TYPE_LIBRARY,          &
                                                     DIHEDRAL_TYPE_LIBRARY,       &
                                                     IMPROPER_TYPE_LIBRARY;

    integer (kind=4), allocatable, dimension(:) :: CONFIG_ATOMID_LIBRARY;

    integer (kind=4), allocatable, dimension(:,:) :: BOND_ATOMID_LIBRARY,     &
                                                       ANGLE_ATOMID_LIBRARY,    &
                                                       DIHEDRAL_ATOMID_LIBRARY, &
                                                       IMPROPER_ATOMID_LIBRARY;

!   ************************************************************************************************

    real (kind=8), dimension(1:2,1:1000) :: INPUT_PAIR_POTENTIAL;

    real (kind=8), allocatable, dimension(:) :: ATOM_MASSES_LIBRARY;

    real (kind=8), allocatable, dimension(:,:) :: POTENTIAL_CLASS2_LIBRARY;

    real (kind=8), allocatable, dimension(:) :: CONFIG_QI_LIBRARY;

    real (kind=8), allocatable, dimension(:,:) :: IMPROPER_COEFFS_LIBRARY;

    real (kind=8), allocatable, dimension(:,:) :: BONDBOND_COEFFS_LIBRARY,           &
                                                  ANGLEANGLETORSION_COEFFS_LIBRARY,  &
                                                  BONDBOND13_COEFFS_LIBRARY;

    real (kind=8), allocatable, dimension(:,:) :: CONFIG_RI_LIBRARY;

    real (kind=8), allocatable, dimension(:,:) :: BOND_COEFFS_LIBRARY,               &
                                                  ANGLE_COEFFS_LIBRARY,              &
                                                  BONDANGLE_COEFFS_LIBRARY,          &
                                                  MIDDLEBONDTORSION_COEFFS_LIBRARY;

    real (kind=8), allocatable, dimension(:,:) :: ANGLEANGLE_COEFFS_LIBRARY,         &
                                                  DIHEDRAL_COEFFS_LIBRARY;

    real (kind=8), allocatable, dimension(:,:) :: ENDBONDTORSION_COEFFS_LIBRARY,     &
                                                  ANGLETORSION_COEFFS_LIBRARY;

!   ************************************************************************************************

    character (len=20) :: CHFILE_FORMAT;

    character (len=250) :: CHNAME_LIBRARY_SUBDIRECTORY;

    character (len=250), allocatable, dimension(:) :: CHNAME_FILE_LIBRARY,         &
                                                      CHNAME_LAMMPS_INPUT_LIBRARY, &
                                                      CHNAME_OUTPUT_FILE_LIBRARY;

    character (len=250), dimension(1:1000) :: PAIR_POTENTIAL_CHTYPE;

    character (len=20), allocatable, dimension(:) :: ATOM_LABEL_LIBRARY;

    character (len=20), allocatable, dimension(:) :: CONFIG_NAT_LIBRARY;

end module module_library
