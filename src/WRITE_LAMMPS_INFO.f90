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

subroutine WRITE_LAMMPS_INFO(icanal,CHEXT)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                  : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                         **
!   **                                                                                            **
!   ** CHEXT                   : NAME OF THE LAMMPS CONFIGURATION                                 **
!   **                                                                                            **
!   ************************************************************************************************

    use module_physical_constants;

    use module_size_arrays;

    use module_config;

    use module_library;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=250), intent(in) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: EOF, EOF2;

!   integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4;

!   real (kind=8) :: ROSEF1, ROSEF2;

!   character (len=250) :: CHOSEF1, CHOSEF2, CHOSEF3;

    character (len=250) :: CHARLINE, CHARLINE2;

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Write additional information in a *.info file';                           !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   ! 
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !  
!   ### Write the name of the interatomic potentials templated file ################################
                                                                                         !
    IOSEF1 = 67 - LEN_TRIM(CHEXT);                                                       !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the interatomic potentials templated file ############################################
                                                                                         !
    open(5,file=TRIM(CHEXT));                                                            !
                                                                                         !
    write(5,'(a24)') 'units               real';                                         !
                                                                                         !
    write(5,*);                                                                          !
                                                                                         !
    write(5,'(a24)') 'atom_style          full';                                         !
                                                                                         !
    write(5,*);                                                                          !
                                                                                         !
!   ### Write the style for pair interactions ######################################################
                                                                                         !
    if ( IPOTENTIAL_CLASS2 == 1 ) then;                                                  !
                                                                                         !
        IOSEF1 = 20 + LEN_TRIM(POTENTIAL_CLASS2_CHTYPE);                                 !
                                                                                         !
        if ( IOSEF1 < 100 ) write(CHOSEF1,'(i2)') IOSEF1;                                !
                                                                                         !
        if ( IOSEF1 < 10 ) write(CHOSEF1,'(i1)') IOSEF1;                                 !
                                                                                         !
        CHOSEF2 = 'a'//TRIM(CHOSEF1)//',3f8.1';                                          !
                                                                                         !
!        write(icanal,*) TRIM(CHOSEF2); 

        write(5,'('//TRIM(CHOSEF2)//')') 'pair_style          '//        &               !
                                         TRIM(POTENTIAL_CLASS2_CHTYPE),  &               !
                                         11.0,  12.0, 12.0;                              !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write the style for bond interactions ######################################################
                                                                                         !
    if ( NTYPE_BOND > 0 ) then;                                                          !
                                                                                         !
        IOSEF1 = 20 + LEN_TRIM(CH_BOND_STYLE);                                           !
                                                                                         !
        if ( IOSEF1 < 100 ) write(CHOSEF1,'(i2)') IOSEF1;                                !
                                                                                         !
        if ( IOSEF1 < 10 ) write(CHOSEF1,'(i1)') IOSEF1;                                 !
                                                                                         !
        CHOSEF2 = 'a'//TRIM(CHOSEF1);                                                    !
                                                                                         !
        write(5,'('//TRIM(CHOSEF2)//')') 'bond_style          '//TRIM(CH_BOND_STYLE);    !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write the style for angle interactions #####################################################
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        IOSEF1 = 20 + LEN_TRIM(CH_ANGLE_STYLE);                                          !
                                                                                         !
        if ( IOSEF1 < 100 ) write(CHOSEF1,'(i2)') IOSEF1;                                !
                                                                                         !
        if ( IOSEF1 < 10 ) write(CHOSEF1,'(i1)') IOSEF1;                                 !
                                                                                         !
        CHOSEF2 = 'a'//TRIM(CHOSEF1);                                                    !
                                                                                         !
        write(5,'('//TRIM(CHOSEF2)//')') 'angle_style         '//TRIM(CH_ANGLE_STYLE);   !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write the style for dihedral interactions ##################################################
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        IOSEF1 = 20 + LEN_TRIM(CH_DIHEDRAL_STYLE);                                       !
                                                                                         !
        if ( IOSEF1 < 100 ) write(CHOSEF1,'(i2)') IOSEF1;                                !
                                                                                         !
        if ( IOSEF1 < 10 ) write(CHOSEF1,'(i1)') IOSEF1;                                 !
                                                                                         !
        CHOSEF2 = 'a'//TRIM(CHOSEF1);                                                    !
                                                                                         !
        write(5,'('//TRIM(CHOSEF2)//')') 'dihedral_style      '// &                      !
                                         TRIM(CH_DIHEDRAL_STYLE);                        !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Additional information about bonds #########################################################
                                                                                         !
    if ( TRIM(CHFILE_FORMAT) == 'amber' ) then;                                          !
                                                                                         !
        write(5,'(a25)') 'special_bonds       amber';                                    !
                                                                                         !
        write(5,*);                                                                      !

!special_bonds       lj/coul  0.0  0.0  0.5
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Set mixing rules for cross interactions ####################################################
                                                                                         !
    if ( TRIM(CHFILE_FORMAT) == 'amber' ) then;                                          !
                                                                                         !
        write(5,'(a35)') 'pair_modify         mix  arithmetic';                          !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         !
!pair_modify         mix  geometric  tail  yes

    end if                                                                               !
                                                                                         !
!   ### Set lammps command for coulomb interactions ################################################
                                                                                         !
    write(5,'(a32)') 'kspace_style        pppm  1.0e-5';                                 !
                                                                                         !
    write(5,*);                                                                          !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Add cross parameters for nonbonded interations #############################################
                                                                                         !
    if ( NPOTENTIAL_CLASS2_CROSS > 0 ) then;                                             !
                                                                                         !
        do i = 1, NPOTENTIAL_CLASS2_CROSS;                                               !
                                                                                         !
            write(5,'(a15,2i6,2f12.6)') 'pair_coeff     ',                  &            !
                                        POTENTIAL_CLASS2_CROSS_ATOMID(1,i), &            !
                                        POTENTIAL_CLASS2_CROSS_ATOMID(2,i), &            !
                                        POTENTIAL_CLASS2_CROSS_VALUES(1,i), &            !
                                        POTENTIAL_CLASS2_CROSS_VALUES(2,i);              !
                                                                                         !   
        end do                                                                           !
                                                                                         !
        write(5,*);                                                                      !
                                                                                         ! 
    end if                                                                               !
                                                                                         !
    close(5);                                                                            !
                                                                                         ! 
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine WRITE_LAMMPS_INFO
