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

subroutine SET_FORCE_FIELD_PROPERTIES(icanal)

!   ************************************************************************************************
!   **                           Set force field properties for lammps                            **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : Canal on which the output file is written                                         **
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

!   ************************************************************************************************

    integer (kind=4) :: IOSEF1; 

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE, ILENGTH_CHEXT;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Set force field properties for lammps';                                   !
                                                                                         ! 
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set and write the atom style ###############################################################
                                                                                         !
    CH_ATOM_STYLE = 'full';                                                              !
                                                                                         !
    IOSEF1 = 70 - 23 - LEN_TRIM(CH_ATOM_STYLE);                                          !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| The atom style is : '//TRIM(CH_ATOM_STYLE)// &              !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Set and write the style of pair interactions ###############################################
                                                                                         !
    POTENTIAL_CLASS2_CHTYPE = 'none';                                                    !
                                                                                         !
    if ( TRIM(CHFILE_FORMAT) == 'lammps-ligpargen' ) then;                               !
                                                                                         !
        IPOTENTIAL_CLASS2 = 1;                                                           !
                                                                                         !
        POTENTIAL_CLASS2_CHTYPE = 'lj/cut/coul/long';                                    !
                                                                                         !
        IOSEF1 = 70 - 22 - 1 - LEN_TRIM(POTENTIAL_CLASS2_CHTYPE);                        !
                                                                                         !
        write(icanal,'(a70)') '| The pair style is : '//      &                          !
                              TRIM(POTENTIAL_CLASS2_CHTYPE)// &                          !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !  
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set and write the style of bond interactions ###############################################
                                                                                         !
    if ( TRIM(CHFILE_FORMAT) == 'lammps-ligpargen' ) then;                               !
                                                                                         !
        CH_BOND_STYLE = 'harmonic';                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
    IOSEF1 = 70 - 23 - LEN_TRIM(CH_BOND_STYLE);                                          !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| The bond style is : '//TRIM(CH_BOND_STYLE)// &              !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set and write the style of angle interactions ############################################## 
                                                                                         !
    if ( TRIM(CHFILE_FORMAT) == 'lammps-ligpargen' ) then;                               !
                                                                                         !
        CH_ANGLE_STYLE = 'harmonic';                                                     !
                                                                                         !
    end if                                                                               !
                                                                                         !
    IOSEF1 = 70 - 24 - LEN_TRIM(CH_ANGLE_STYLE);                                         !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| The angle style is : '//TRIM(CH_ANGLE_STYLE)// &            !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set and write the style of dihedral interactions ###########################################
                                                                                         !
    if ( TRIM(CHFILE_FORMAT) == 'lammps-ligpargen' ) then;                               !
                                                                                         !
        CH_DIHEDRAL_STYLE = 'opls';                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
    IOSEF1 = 70 - 27 - LEN_TRIM(CH_DIHEDRAL_STYLE);                                      !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| The dihedral style is : '// &                               !
                          TRIM(CH_DIHEDRAL_STYLE)//      &                               !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set and write the style of improper interactions ###########################################
                                                                                         !
   if ( TRIM(CHFILE_FORMAT) == 'lammps-ligpargen' ) then;                                !
                                                                                         !
        CH_IMPROPER_STYLE = 'none';                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
    IOSEF1 = 70 - 27 - LEN_TRIM(CH_IMPROPER_STYLE);                                      !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| The improper style is : '// &                               !
                          TRIM(CH_IMPROPER_STYLE)//      &                               !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine SET_FORCE_FIELD_PROPERTIES
