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

subroutine READ_AMBER_INPCRD(icanal,CHEXT)

!   ************************************************************************************************
!   **                                 READ LAMMPS FILES                                          **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal        : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                   **
!   **                                                                                            **
!   ** CHEXT         : NAME OF THE LAMMPS CONFIGURATION                                           **
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

    integer (kind=4) :: NREAD_LINE, IATOM, JATOM;

!   ************************************************************************************************

    logical :: PROBE1;
 
!   ************************************************************************************************

!   integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE, ILENGTH_CHEXT;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Read AMBER inpcrd file';                                                  !
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
!   ##### Start reading the head of the configuration file #########################################
                                                                                         !
    open(2,file=TRIM(CHEXT));                                                            !
                                                                                         !
    read(2,*);                                                                           !    
                                                                                         !
    read(2,*) NATOM;                                                                     !
                                                                                         !
    write(icanal,'(a11,i8,a51)') '| There is ',                  &                       !
                                 NATOM,                          &                       !
                                 ' atoms in the configuration'// &                       !
                                 REPEAT(' ',23)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate configuration arrays ##############################################################
                                                                                         !
    allocate(CONFIG_RI(1:3,1:NATOM));                                                    !
                                                                                         !
!   ### Initialization of configuration arrays #####################################################
                                                                                         !
    CONFIG_RI(1:3,1:NATOM) = 0.0d0;                                                      !
                                                                                         !
!   ### Read coordinates ###########################################################################
                                                                                         !
    NREAD_LINE = CEILING( REAL( NATOM ) / 2.0d0 );                                       !
                                                                                         !
    write(icanal,'(a11,i8,a51)') '| There is ',     &                                    !
                                 NREAD_LINE,        &                                    !
                                 ' lines to read'// &                                    !
                                 REPEAT(' ',36)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IATOM = 0;                                                                           !
                                                                                         !
    do i = 1, NREAD_LINE;                                                                !
                                                                                         !
        IATOM = IATOM + 1;                                                               !
                                                                                         !
        if ( IATOM == NATOM ) then;                                                      !
                                                                                         !
            read(2,*) CONFIG_RI(1:3,IATOM);                                              !
                                                                                         !
        else                                                                             !
                                                                                         ! 
            JATOM = IATOM + 1;                                                           !
                                                                                         !
            read(2,*) CONFIG_RI(1:3,IATOM), CONFIG_RI(1:3,JATOM);                        !
                                                                                         !
            IATOM = JATOM;                                                               !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
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
end subroutine READ_AMBER_INPCRD
