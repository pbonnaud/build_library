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

subroutine WRITE_LAMMPS_TEMPLATE(icanal,CHEXT) 

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal  : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                         **
!   **                                                                                            **
!   ** CHEXT  : NAME OF THE LAMMPS CONFIGURATION                                                  **
!   **                                                                                            **
!   ************************************************************************************************

    use module_physical_constants;

    use module_size_arrays;

    use module_library;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=250), intent(in) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    character (len=250) :: CHARLINE;

    integer (kind=4) :: EOF;

!   integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4;

!   character (len=10) :: CHOSEF1, CHOSEF2, CHOSEF3;

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Write LAMMPS template';                                                   !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   ### Write the name of the templated file #######################################################
                                                                                         !
    IOSEF1 = 67 - LEN_TRIM(CHEXT);                                                       !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Write the templated file ###################################################################
                                                                                         !
    open(2,file=TRIM(CHEXT));                                                            !
                                                                                         !
    write(2,'(a14)') '# Molecule of ';                                                   !
                                                                                         !
    write(2,*);                                                                          !
                                                                                         !
    if ( NATOM_LIBRARY > 0 ) write(2,'(i8,a6)') NATOM_LIBRARY, ' atoms';                 !
                                                                                         !
    if ( NBOND_LIBRARY > 0 ) write(2,'(i8,a6)') NBOND_LIBRARY, ' bonds';                 !
                                                                                         !
    if ( NANGLE_LIBRARY > 0 ) write(2,'(i8,a7)') NANGLE_LIBRARY, ' angles';              !
                                                                                         !
    if ( NDIHEDRAL_LIBRARY > 0 ) write(2,'(i8,a10)') NDIHEDRAL_LIBRARY, ' dihedrals';    !
                                                                                         !
    if ( NIMPROPER_LIBRARY > 0 ) write(2,'(i8,a10)') NIMPROPER_LIBRARY, ' impropers';    !
                                                                                         !
    write(2,*);                                                                          !
                                                                                         !
    if ( NATOM_LIBRARY > 0 ) then; 

        write(2,'(a7)') 'Charges';

        write(2,*);

        do i = 1, NATOM_LIBRARY;

            write(2,'(i8,f15.6)') CONFIG_ATOMID_LIBRARY(i), &
                                  CONFIG_QI_LIBRARY(CONFIG_ATOMID_LIBRARY(i));
        end do

        write(2,*);

        write(2,'(a6)') 'Coords';

        write(2,*);

        do i = 1, NATOM_LIBRARY;

            write(2,'(i8,3f15.6)') CONFIG_ATOMID_LIBRARY(i),  &
                                   CONFIG_RI_LIBRARY(1:3,CONFIG_ATOMID_LIBRARY(i));

        end do

        write(2,*);

        write(2,'(a5)') 'Types';

        write(2,*);

        do i = 1, NATOM_LIBRARY;

            write(2,'(2i8)') CONFIG_ATOMID_LIBRARY(i), &
                             CONFIG_ATOM_TYPE_LIBRARY(CONFIG_ATOMID_LIBRARY(i));

        end do

        write(2,*);

    end if

    if ( NBOND_LIBRARY > 0 ) then;

        write(2,'(a5)') 'Bonds';

        write(2,*);

        do i = 1, NBOND_LIBRARY;

            write(2,'(4i8)') i, BOND_TYPE_LIBRARY(i), BOND_ATOMID_LIBRARY(1:2,i);

        end do

        write(2,*);

    end if

    if ( NANGLE_LIBRARY > 0 ) then;

        write(2,'(a6)') 'Angles';

        write(2,*);

        do i = 1, NANGLE_LIBRARY;

            write(2,'(5i8)') i, ANGLE_TYPE_LIBRARY(i), ANGLE_ATOMID_LIBRARY(1:3,i);

        end do

        write(2,*);

    end if

    if ( NDIHEDRAL_LIBRARY > 0 ) then;

        write(2,'(a9)') 'Dihedrals';

        write(2,*);

        do i = 1, NDIHEDRAL_LIBRARY;

            write(2,'(6i8)') i, DIHEDRAL_TYPE_LIBRARY(i), DIHEDRAL_ATOMID_LIBRARY(1:4,i);

        end do

        write(2,*);

    end if 

!   ### Write atom IDs involved in impropers #######################################################
                                                                                         !
    if ( NIMPROPER_LIBRARY > 0 ) then;                                                   !
                                                                                         !
        write(2,'(a9)') 'Impropers';                                                     !
                                                                                         !
        write(2,*);                                                                      !
                                                                                         !
        do i = 1, NIMPROPER_LIBRARY;                                                     !
                                                                                         !
            write(2,'(6i8)') i,                              &                           !
                             IMPROPER_TYPE_LIBRARY(i),       &                           !
                             IMPROPER_ATOMID_LIBRARY(1:4,i);                             !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    close(2);                                                                            !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         ! 
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine WRITE_LAMMPS_TEMPLATE
