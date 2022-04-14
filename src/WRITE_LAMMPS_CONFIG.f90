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

subroutine WRITE_LAMMPS_CONFIG(icanal,CHEXT,NATOM,              &
                                            NBOND,              &
                                            NANGLE,             &
                                            NDIHEDRAL,          &
                                            NIMPROPER,          &
                                            NTYPE_ATOM,         &
                                            NTYPE_BOND,         &
                                            NTYPE_ANGLE,        &
                                            NTYPE_DIHEDRAL,     &
                                            NTYPE_IMPROPER,     &
                                            ATOM_MASSES,        &
                                            ATOM_LABEL,         & 
                                            CELL_AXIS,          &
                                            MATA,               &
                                            CONFIG_NAT,         &
                                            CONFIG_QIJ,         &
                                            CONFIG_RIJ,         &
                                            CONFIG_VIJ,         &
                                            POTENTIAL_CLASS2,   &
                                            BOND_COEFFS,        &
                                            ANGLE_COEFFS,       &
                                            BONDBOND_COEFFS,    &
                                            BONDANGLE_COEFFS,   &
                                            IMPROPER_COEFFS,    &
                                            ANGLEANGLE_COEFFS,  &
                                            CONFIG_ATOMID,      &
                                            CONFIG_MOLECULEID,  &
                                            CONFIG_ATOM_TYPE,   &
                                            BOND_TYPE,          &
                                            ANGLE_TYPE,         &
                                            IMPROPER_TYPE,      &
                                            BOND_ATOMID,        &
                                            ANGLE_ATOMID,       &
                                            IMPROPER_ATOMID,    &
                                            IPOTENTIAL_CLASS2)

!   ************************************************************************************************
!   **                                 READ LAMMPS FILES                                          **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                          **
!   ** CHEXT  : NAME OF THE LAMMPS CONFIGURATION                                                  **
!   **                                                                                            **
!   ** NATOM             : NUMBER OF ATOMS                                                        **
!   ** NBOND             : NUMBER OF BONDS IN THE MOLECULAR CONFIGURATION                         **
!   ** NANGLE            : NUMBER OF ANGLES IN THE MOLECULAR CONFIGURATION                        **
!   ** NDIHEDRAL         : NUMBER OF DIHEDRAL ANGLES IN THE MOLECULAR CONFIGURATION               **
!   ** NIMPROPER         : NUMBER OF IMPROPER ANGLES IN THE MOLECULAR CONFIGURATION               **
!   ** NTYPE_ATOM        : NUMBER OF ATOM TYPE IN THE SIMULATION BOX                              **
!   ** NTYPE_BOND        : NUMBER OF BOND TYPE IN THE MOLECULAR CONFIGURATION                     **
!   ** NTYPE_ANGLE       : NUMBER OF ANGLE TYPE IN THE MOLECULAR CONFIGURATION                    **
!   ** NTYPE_DIHEDRAL    : NUMBER OF DIHEDRAL ANGLE TYPE IN THE MOLECULAR CONFIGURATION           **
!   ** NTYPE_IMPROPER    : NUMBER OF IMPROPER ANGLE TYPE IN THE MOLECULAR CONFIGURATION           **
!   ** ATOM_MASSES       : MASSES OF ATOMS IN THE SIMULATION BOX                                  **
!   ** ATOM_LABEL        : LABELS OF ATOMS IN THE SIMULATION BOX                                  **
!   ** CELL_AXIS         : CELL DIMENSIONS                                                        **
!   ** MATA              : PASSAGE MATRIX                                                         **
!   ** CONFIG_NAT        : NATURE OF ATOMS IN THE SIMULATION BOX                                  **
!   ** CONFIG_QIJ        : CHARGE OF ATOMS IN THE SIMULATION BOX                                  **
!   ** CONFIG_RIJ        : POSITION OF ATOMS IN THE SIMULATION BOX                                **
!   ** CONFIG_VIJ        : VWLOCITIES OF ATOMS IN THE SIMULATION BOX                              **
!   ** POTENTIAL_CLASS2  : PAIR POTENTIAL COEFFICIENT FOR THE SIMULATION                          **
!   ** BOND_COEFFS       : PARAMETERS **
!   ** ANGLE_COEFFS      : PARAMETERS **
!   ** BONDBOND_COEFFS   : PARAMETERS **
!   ** BONDANGLE_COEFFS  : PARAMETERS **
!   ** IMPROPER_COEFFS   : PARAMETERS **
!   ** ANGLEANGLE_COEFFS : PARAMETERS **
!   ** CONFIG_ATOMID     : ATOM ID OF THE CONFIGURATION THAT HAS BEEN READ                        **  
!   ** CONFIG_MOLECULEID : MOLECULE ID OF THE CONFIGURATION THAT IS READ                          **
!   ** BOND_TYPE         : TYPE OF BONDS IN THE MOLECULAR CONFIGURATION                           **
!   ** ANGLE_TYPE        : TYPE OF ANGLES IN THE MOLECULAR CONFIGURATION                          **
!   ** IMPROPER_TYPE     : TYPE OF IMPROPER IN THE MOLECULAR CONFIGURATION                        **
!   ** BOND_ATOMID       : ATOM IDs INVOLVED IN BONDS                                             **
!   ** ANGLE_ATOMID      : ATOM IDs INVOLVED IN ANGLES                                            **
!   ** IMPROPER_ATOMID   : ATOM IDs INVOLVED IN IMPROPERS                                         **
!   ** IPOTENTIAL_CLASS2 :
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=20), intent(in) :: CHEXT;

    integer (kind=4), intent(in) :: NATOM, NBOND, NTYPE_ATOM, NTYPE_BOND, IPOTENTIAL_CLASS2;

    integer (kind=4), intent(in) :: NANGLE, NDIHEDRAL, NIMPROPER;

    integer (kind=4), intent(in) :: NTYPE_ANGLE, NTYPE_DIHEDRAL, NTYPE_IMPROPER;

    real (kind=8), dimension(1:20), intent(in) :: ATOM_MASSES;

    real (kind=8), dimension(1:3), intent(in) :: CELL_AXIS;

    real (kind=8), dimension(1:6), intent(in) :: MATA;

    character (len=20), dimension(1:20), intent(in) :: ATOM_LABEL;

    character (len=3), dimension(1:100000), intent(in) :: CONFIG_NAT;

    real (kind=8), dimension(1:100000), intent(in) :: CONFIG_QIJ;

    real (kind=8), dimension(1:3,1:100000), intent(in) :: CONFIG_RIJ, CONFIG_VIJ;

    integer (kind=4), dimension(1:100000), intent(in) :: CONFIG_ATOMID, CONFIG_MOLECULEID, CONFIG_ATOM_TYPE;

    real (kind=8), dimension(1:2,1:100), intent(in) :: POTENTIAL_CLASS2;

    real (kind=8), dimension(1:2,1:10), intent(in) :: IMPROPER_COEFFS;

    real (kind=8), dimension(1:4,1:10), intent(in) :: BOND_COEFFS, ANGLE_COEFFS;

    real (kind=8), dimension(1:4,1:10), intent(in) :: BONDANGLE_COEFFS;

    real (kind=8), dimension(1:6,1:10), intent(in) :: ANGLEANGLE_COEFFS;

    real (kind=8), dimension(1:3,1:10), intent(in) :: BONDBOND_COEFFS;

    integer (kind=4), dimension(1:100000), intent(in) :: BOND_TYPE, ANGLE_TYPE, IMPROPER_TYPE;

    integer (kind=4), dimension(1:2,1:100000), intent(in) :: BOND_ATOMID;

    integer (kind=4), dimension(1:3,1:100000), intent(in) :: ANGLE_ATOMID;

    integer (kind=4), dimension(1:4,1:100000), intent(in) :: IMPROPER_ATOMID;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: IATOM_TYPE;

!   ************************************************************************************************

    open(103,file=TRIM(CHEXT));
    write(103,'(a52)') '# System description ###############################';
    write(103,'(a6,3f15.6)') '#     ', CELL_AXIS(1:3);
    write(103,'(i8,1x,a5)') NATOM,     'atoms';
    if ( NBOND     > 0 ) write(103,'(i8,1x,a5)') NBOND,     'bonds';
    if ( NANGLE    > 0 ) write(103,'(i8,1x,a6)') NANGLE,    'angles';
    if ( NDIHEDRAL > 0 ) write(103,'(i8,1x,a9)') NDIHEDRAL, 'dihedrals';
    if ( NIMPROPER > 0 ) write(103,'(i8,1x,a9)') NIMPROPER, 'impropers';
    write(103,*);
    write(103,'(i8,1x,a10)') NTYPE_ATOM,     'atom types';
    if ( NBOND     > 0 ) write(103,'(i8,1x,a10)') NTYPE_BOND,     'bond types';
    if ( NANGLE    > 0 ) write(103,'(i8,1x,a11)') NTYPE_ANGLE,    'angle types';
    if ( NDIHEDRAL > 0 ) write(103,'(i8,1x,a14)') NTYPE_DIHEDRAL, 'dihedral types';
    if ( NIMPROPER > 0 ) write(103,'(i8,1x,a14)') NTYPE_IMPROPER, 'improper types';
    write(103,*);
    write(103,'(2f15.6,a11)') -0.5d0 * MATA(1), 0.5d0 * MATA(1), '    xlo xhi';
    write(103,'(2f15.6,a11)') -0.5d0 * MATA(2), 0.5d0 * MATA(2), '    ylo yhi';
    write(103,'(2f15.6,a11)') -0.5d0 * MATA(3), 0.5d0 * MATA(3), '    zlo zhi';

    if ( ( ABS(MATA(6)) > 1.0E-4 ) .OR. &
         ( ABS(MATA(5)) > 1.0E-4 ) .OR. &
         ( ABS(MATA(4)) > 1.0E-4 ) ) write(103,'(3f15.6,a9)') MATA(6), &
                                                              MATA(5), &
                                                              MATA(4), &
                                                              ' xy xz yz';
    write(103,*);
    write(103,'(a6)') 'Masses';
    write(103,*);
    do i = 1, NTYPE_ATOM
        write(103,'(i4,f12.6)') i, ATOM_MASSES(i);
    end do

    if ( IPOTENTIAL_CLASS2 == 1 ) then
        write(103,*);
        write(103,'(a33)') 'Pair Coeffs # lj/class2/coul/long';
        write(103,*);

        do i = 1, NTYPE_ATOM;
            write(103,'(i4,2f12.4)') i, POTENTIAL_CLASS2(1:2,i);
        end do
    end if

    if ( NTYPE_BOND /= 0 ) then
        write(103,*);
        write(103,'(a20)') 'Bond Coeffs # class2';
        write(103,*);
        do i = 1, NTYPE_BOND
            write(103,'(i4,4f12.4)') i, BOND_COEFFS(1:4,i);
        end do
    end if

    if ( NTYPE_ANGLE /= 0 ) then
        write(103,*);
        write(103,'(a21)') 'Angle Coeffs # class2';
        write(103,*);
        do i = 1, NTYPE_ANGLE
            write(103,'(i4,4f12.4)') i, ANGLE_COEFFS(1:4,i);
        end do

        write(103,*);
        write(103,'(a15)') 'BondBond Coeffs';
        write(103,*);

        do i = 1, NTYPE_ANGLE
            write(103,'(i4,3f12.4)') i, BONDBOND_COEFFS(1:3,i);
        end do

        write(103,*);
        write(103,'(a16)') 'BondAngle Coeffs';
        write(103,*);

        do i = 1, NTYPE_ANGLE
            write(103,'(i4,4f12.4)') i, BONDANGLE_COEFFS(1:4,i);
        end do
    end if

    if ( NTYPE_IMPROPER /= 0 ) then
        write(103,*);
        write(103,'(a24)') 'Improper Coeffs # class2';
        write(103,*);

        do i = 1, NTYPE_IMPROPER
            write(103,'(i4,2f12.4)') i, IMPROPER_COEFFS(1:2,i);
        end do

        write(103,*);
        write(103,'(a17)') 'AngleAngle Coeffs';
        write(103,*);

        do i = 1, NTYPE_IMPROPER
            write(103,'(i4,6f12.4)') i, ANGLEANGLE_COEFFS(1:6,i);
        end do
    end if

    write(103,*);
    write(103,'(a5)') 'Atoms';
    write(103,*);

    do i = 1, NATOM
        if ( CONFIG_MOLECULEID(i) == 0 ) then
            write(103,'(2i8,4f21.15)') CONFIG_ATOMID(i),    &
                                       CONFIG_ATOM_TYPE(i), &
                                       CONFIG_QIJ(i),       &
                                       CONFIG_RIJ(1:3,i);
        else
            write(103,'(3i8,4f21.15)') CONFIG_ATOMID(i),     &
                                       CONFIG_MOLECULEID(i), &
                                       CONFIG_ATOM_TYPE(i),  &
                                       CONFIG_QIJ(i),        &
                                       CONFIG_RIJ(1:3,i);
        end if
    end do

    write(103,*);
    write(103,'(a10)') 'Velocities';
    write(103,*);

    do i = 1, NATOM
        write(103,'(i8,3f21.15)') CONFIG_ATOMID(i), CONFIG_VIJ(1:3,i);
    end do

    if ( NBOND > 0 ) then
        write(103,*);
        write(103,'(a5)') 'Bonds';
        write(103,*);

        do i = 1, NBOND
            write(103,'(4i8)') i, BOND_TYPE(i), BOND_ATOMID(1:2,i);
        end do
    end if

    if ( NANGLE > 0 ) then
        write(103,*);
        write(103,'(a6)') 'Angles';
        write(103,*);

        do i = 1, NANGLE
            write(103,'(5i8)') i, ANGLE_TYPE(i), ANGLE_ATOMID(1:3,i);
        end do
    end if

    if ( NIMPROPER > 0 ) then
        write(103,*);
        write(103,'(a9)') 'Impropers';
        write(103,*);

        do i = 1, NIMPROPER
            write(103,'(6i8)') i, IMPROPER_TYPE(i), IMPROPER_ATOMID(1:4,i);
        end do
    end if

    close(103);

end subroutine WRITE_LAMMPS_CONFIG

