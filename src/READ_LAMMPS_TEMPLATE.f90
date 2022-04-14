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

subroutine READ_LAMMPS_TEMPLATE(icanal,CHEXT,                    &
                                       NATOM_TPLTE,              &
                                       NBOND_TPLTE,              &
                                       NANGLE_TPLTE,             &
                                       NDIHEDRAL_TPLTE,          &
                                       NIMPROPER_TPLTE,          &
                                       TEMPLATE_TYPE,            &
                                       NTYPE_TEMPLATE,           &
                                       NTYPE_BOND_TPLTE,         &
                                       NTYPE_ANGLE_TPLTE,        &
                                       NTYPE_IMPROPER_TPLTE,     &
                                       CONFIG_QI_TPLTE,          &
                                       CONFIG_RI_TPLTE,          &
                                       CONFIG_ATOMID_TPLTE,      &
                                       BOND_TYPE_TPLTE,          &
                                       ANGLE_TYPE_TPLTE,         &
                                       IMPROPER_TYPE_TPLTE,      &
                                       BOND_ATOMID_TPLTE,        &
                                       ANGLE_ATOMID_TPLTE,       &
                                       IMPROPER_ATOMID_TPLTE)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal  : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                         **
!   ** ILENGTH : LENGTH OF THE CHARACTER CHAIN CHEXT                                              **
!   ** CHEXT  : NAME OF THE LAMMPS CONFIGURATION                                                  **
!   **                                                                                            **
!   ** NATOM_TPLTE          : NUMBER OF ATOMS IN THE TEMPLATE                                        **
!   ** NBOND_TPLTE          : NUMBER OF COVALENT BONDS IN THE MOLECULE DESRCIBED IN THE TEMPLATE     **
!   ** NANGLE_TPLTE         : NUMBER OF ANGLES IN THE TEMPLATED MOLECULE  **
!   ** NDIHEDRAL_TPLTE      : NUMBER OF DIHEDRAL ANGLES IN THE TEMPLATED MOLECULE **
!   ** NIMPROPER_TPLTE      : NUMBER OF IMPROPER ANGLES IN THE TEMPLATED MOLECULE **

!   ** TEMPLATE_TYPE        : TYPE OF ATOMS IN THE TEMPLATED MOLECULE                                **
!   ** NTYPE_TEMPLATE       : NUMBER OF ATOM TYPES IN THE TEMPLATED MOLECULE                         **
!   ** NTYPE_BOND_TPLTE     : NUMBER OF BOND TYPES IN THE TEMPLATED MOLECULE                         **
!   ** NTYPE_ANGLE_TPLTE    : NUMBER OF ANGLE TYPES IN THE TEMPLATED MOLECULE                      **
!   ** NTYPE_IMPROPER_TPLTE : 


!   ** CONFIG_QI_TPLTE     : PARTIAL CHARGES OF ATOMS IN THE TEMPLATED MOLECULE                   **

!   ** CONFIG_RI_TPLTE     : COORDINATES OF ATOMS IN THE TEMPLATED MOLECULE                       **
!   ** CONFIG_ATOMID_TPLTE : ATOM ID IN THE TEMPLATED MOLECULE                                    **
!   ** BOND_TYPE_TPLTE     : TYPE OF BONDS IN THE MOLECULAR CONFIGURATION                         **
!   ** ANGLE_TYPE_TPLTE    : TYPE OF ANGLES IN THE MOLECULAR CONFIGURATION                        **
!   ** IMPROPER_TYPE_TPLTE : TYPE OF IMPROPER ANGLES IN THE MOLECULAR CONFIGURATION               **


!   ** BOND_ATOMID_TPLTE     : ATOM ID RELATED TO THE BONDS IN THE MOLECULAR CONFIGURATION            **
!   ** ANGLE_ATOMID_TPLTE    : ATOM ID RELATED TO THE ANGLES IN THE MOLECULAR CONFIGURATION
!   ** IMPROPER_ATOMID_TPLTE : ATOM ID RELATED TO THE IMPROPER ANGLES IN THE MOLECULAR CONFIGURATION

!   ** NANGLE_TPLTE      : NUMBER OF BENDING ANGLES IN THE MOLECULE DESCRIBED IN THE TEMPLATE     **
!   ** NDIHEDRAL_TPLTE   : NUMBER OF DIHEDRAL ANGLES IN THE MOLECULES DESCRIBED IN THE TEMPLATE   **
!   ** NIMPROPER_TPLTE   : NUMBER OF IMPROPER ANGLES IN THE MOLECULES DESCRIBED IN THE TEMPLATE   **



!   ** TEMPLATE_BOND     : DEFINE BONDS AMONG ATOMS OF THE TEMPLATED MOLECULE                     **
!   ** TEMPLATE_ANGLE    : DEFINE ANGLES AMONG ATOMS OF THE TEMPLATED MOLECULES                   **
!   ** TEMPLATE_IMPROPER : DEFINE IMPROPER ANGLES AMONG ATOMS OF THE TEMPLATED MOLECULE           **

!   **                                                                                            **
!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=250), intent(in) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4), intent(out) :: NTYPE_TEMPLATE;

    integer (kind=4), intent(out) :: NATOM_TPLTE, NBOND_TPLTE, NANGLE_TPLTE, NDIHEDRAL_TPLTE, NIMPROPER_TPLTE;

    integer (kind=4), intent(out) :: NTYPE_BOND_TPLTE, NTYPE_ANGLE_TPLTE, NTYPE_IMPROPER_TPLTE;

    integer (kind=4), dimension(1:1000), intent(out) :: TEMPLATE_TYPE, BOND_TYPE_TPLTE;

    integer (kind=4), dimension(1:1000), intent(out) :: ANGLE_TYPE_TPLTE, IMPROPER_TYPE_TPLTE; 

    integer (kind=4), dimension(1:1000), intent(out) :: CONFIG_ATOMID_TPLTE;

    integer (kind=4), dimension (1:2,1:1000), intent(out) :: BOND_ATOMID_TPLTE;

    integer (kind=4), dimension (1:3,1:1000), intent(out) :: ANGLE_ATOMID_TPLTE;

    integer (kind=4), dimension (1:4,1:1000), intent(out) :: IMPROPER_ATOMID_TPLTE;

    real (kind=8), dimension (1:1000), intent(out) :: CONFIG_QI_TPLTE;

    real (kind=8), dimension (1:3,1:1000), intent(out) :: CONFIG_RI_TPLTE;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    character (len=250) :: CHARLINE;

    integer (kind=4) :: EOF;

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4;

    character (len=10) :: CHOSEF1, CHOSEF2, CHOSEF3;

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************

    CHTITLE = 'READ LAMMPS TEMPLATE';

    ILENGTH_TITLE = LEN_TRIM(CHTITLE);

    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));

    write(icanal,*) TRIM(CHEXT);

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    NATOM_TPLTE     = 0;
    NBOND_TPLTE     = 0;
    NANGLE_TPLTE    = 0;
    NDIHEDRAL_TPLTE = 0;
    NIMPROPER_TPLTE = 0;

    NTYPE_TEMPLATE       = 0;
    NTYPE_BOND_TPLTE     = 0;
    NTYPE_ANGLE_TPLTE    = 0;
    NTYPE_IMPROPER_TPLTE = 0;

    TEMPLATE_TYPE(1:1000)       = 0;
    BOND_TYPE_TPLTE(1:1000)     = 0;
    ANGLE_TYPE_TPLTE(1:1000)    = 0;
    IMPROPER_TYPE_TPLTE(1:1000) = 0;

    CONFIG_ATOMID_TPLTE(1:1000) = 0;

    BOND_ATOMID_TPLTE(1:2,1:1000)     = 0;
    ANGLE_ATOMID_TPLTE(1:3,1:1000)    = 0;
    IMPROPER_ATOMID_TPLTE(1:4,1:1000) = 0;

    CONFIG_QI_TPLTE(1:1000) = 0.0d0;

    CONFIG_RI_TPLTE(1:3,1:1000) = 0.0d0;

    open(2,file=TRIM(CHEXT));

    EOF = 0;

    do
        read(2,'(a)',iostat=EOF) CHARLINE;

        if ( EOF /= 0 ) EXIT;

        if ( INDEX(CHARLINE,'atoms') > 0 ) then
            read(CHARLINE,*) NATOM_TPLTE;

        else if ( INDEX(CHARLINE,'bonds') > 0 ) then
            read(CHARLINE,*) NBOND_TPLTE;

        else if ( INDEX(CHARLINE,'angles') > 0 ) then
            read(CHARLINE,*) NANGLE_TPLTE;

        else if ( INDEX(CHARLINE,'dihedrals') > 0 ) then
            read(CHARLINE,*) NDIHEDRAL_TPLTE;

        else if ( INDEX(CHARLINE,'impropers') > 0 ) then
            read(CHARLINE,*) NIMPROPER_TPLTE;

        else if ( INDEX(CHARLINE,'Charges') > 0 ) then
            read(2,*);
            do i = 1, NATOM_TPLTE;
                read(2,*) IOSEF1, CONFIG_QI_TPLTE(IOSEF1);  !TEMPLATE_QIJ(IOSEF1);
            end do

        else if ( INDEX(CHARLINE,'Coords') > 0 ) then
            read(2,*);
            do i = 1, NATOM_TPLTE;
                read(2,*) CONFIG_ATOMID_TPLTE(i), CONFIG_RI_TPLTE(1:3,CONFIG_ATOMID_TPLTE(i));
            end do

        else if ( INDEX(CHARLINE,'Types') > 0 ) then
            read(2,*);
            do i = 1, NATOM_TPLTE;
                read(2,*) IOSEF1, TEMPLATE_TYPE(IOSEF1);
                 if ( TEMPLATE_TYPE(IOSEF1) > NTYPE_TEMPLATE ) NTYPE_TEMPLATE = TEMPLATE_TYPE(IOSEF1);
            end do

        else if ( INDEX(CHARLINE,'Bonds') > 0 ) then
            read(2,*);
            do i = 1, NBOND_TPLTE;
                read(2,*) IOSEF1, BOND_TYPE_TPLTE(IOSEF1), BOND_ATOMID_TPLTE(1:2,IOSEF1);      !TEMPLATE_BOND(1:3,IOSEF1);
            end do

            if ( BOND_TYPE_TPLTE(IOSEF1) > NTYPE_BOND_TPLTE ) NTYPE_BOND_TPLTE = BOND_TYPE_TPLTE(IOSEF1);

        else if ( INDEX(CHARLINE,'Angles') > 0 ) then
            read(2,*);
            do i = 1, NANGLE_TPLTE;
                read(2,*) IOSEF1, ANGLE_TYPE_TPLTE(IOSEF1), ANGLE_ATOMID_TPLTE(1:3,IOSEF1);
            end do

            if ( ANGLE_TYPE_TPLTE(IOSEF1) > NTYPE_ANGLE_TPLTE ) NTYPE_ANGLE_TPLTE = ANGLE_TYPE_TPLTE(IOSEF1);

!           do i = 1, NANGLE_TPLTE;
!               write(icanal,'(5i8)') i, ANGLE_TYPE_TPLTE(i), ANGLE_ATOMID_TPLTE(1:3,i);
!           end do

         else if ( INDEX(CHARLINE,'Impropers') > 0 ) then
            read(2,*);
            do i = 1, NIMPROPER_TPLTE;
                read(2,*) IOSEF1, IMPROPER_TYPE_TPLTE(IOSEF1), IMPROPER_ATOMID_TPLTE(1:4,IOSEF1);
            end do

            if ( IMPROPER_TYPE_TPLTE(IOSEF1) > NTYPE_IMPROPER_TPLTE ) NTYPE_IMPROPER_TPLTE = IMPROPER_TYPE_TPLTE(IOSEF1);
        end if
    end do

    close(2);

    write(icanal,'(a16,i8,a46)') '| NATOM_TPLTE : ', NATOM_TPLTE, REPEAT(' ',45)//'|';
    write(icanal,*) '| NBOND_TPLTE          : ', NBOND_TPLTE;
    write(icanal,*) '| NANGLE_TPLTE         : ', NANGLE_TPLTE;
    write(icanal,*) '| NDIHEDRAL_TPLTE      : ', NDIHEDRAL_TPLTE;
    write(icanal,*) '| NIMPROPER_TPLTE      : ', NIMPROPER_TPLTE;
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,*) '| NTYPE_TEMPLATE       : ', NTYPE_TEMPLATE;
    write(icanal,*) '| NTYPE_BOND_TPLTE     : ', NTYPE_BOND_TPLTE; 
    write(icanal,*) '| NTYPE_ANGLE_TPLTE    : ', NTYPE_ANGLE_TPLTE;
    write(icanal,*) '| NTYPE_IMPROPER_TPLTE : ', NTYPE_IMPROPER_TPLTE;
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine READ_LAMMPS_TEMPLATE
