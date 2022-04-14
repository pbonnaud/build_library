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

subroutine SET_ATOM_PROPERTIES(icanal)

!   ************************************************************************************************
!   **                                    READ XYZ FILES                                          **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                          **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

!   real (kind=8) :: ROSEF1;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Set atom properties';                                                     !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   ### Initialization of arrays ###################################################################
                                                                                         !
    MOLAR_MASS_ELEMENT(1:200) = 0.0d0;                                                   !
                                                                                         !
    ELEMENT_RADIUS(1:200) = 0.0d0;                                                       !
                                                                                         !
    ELEMENT_FULL_NAME(1:200) = 'xxx';                                                    !
                                                                                         !
    ELEMENT_SYMBOL_NAME(1:200) = 'XXX';                                                  !
                                                                                         !
!   ### Set the molar mass of elements #############################################################
                                                                                         !
    MOLAR_MASS_ELEMENT(1)  = 1.00794d0;                                                  !
    MOLAR_MASS_ELEMENT(2)  = 4.002602d0;                                                 !
    MOLAR_MASS_ELEMENT(3)  = 6.941d0;                                                    !
    MOLAR_MASS_ELEMENT(4)  = 9.0121831d0;                                                !
    MOLAR_MASS_ELEMENT(5)  = 10.811d0;                                                   !
    MOLAR_MASS_ELEMENT(6)  = 12.0107d0;
    MOLAR_MASS_ELEMENT(7)  = 14.0067d0;
    MOLAR_MASS_ELEMENT(8)  = 15.9994d0;
    MOLAR_MASS_ELEMENT(9)  = 18.998403163d0;
    MOLAR_MASS_ELEMENT(11) = 22.98976928d0;
    MOLAR_MASS_ELEMENT(14) = 28.0855d0;
    MOLAR_MASS_ELEMENT(15) = 30.973761998d0;
    MOLAR_MASS_ELEMENT(16) = 32.064d0;
    MOLAR_MASS_ELEMENT(17) = 35.453d0;
    MOLAR_MASS_ELEMENT(20) = 40.078d0;
    MOLAR_MASS_ELEMENT(26) = 55.8450d0;
    MOLAR_MASS_ELEMENT(35) = 79.904d0;

!   ### Set the full name of elements ##############################################################
                                                                                         !
    ELEMENT_FULL_NAME(1)  = 'hydrogen';                                                  !
    ELEMENT_FULL_NAME(2)  = 'helium';                                                    !
    ELEMENT_FULL_NAME(3)  = 'lithium';                                                   !
    ELEMENT_FULL_NAME(4)  = 'beryllium';                                                 !
    ELEMENT_FULL_NAME(5)  = 'bore';                                                      !
    ELEMENT_FULL_NAME(6)  = 'carbon';                                                    !
    ELEMENT_FULL_NAME(7)  = 'nitrogen';
    ELEMENT_FULL_NAME(8)  = 'oxygen';
    ELEMENT_FULL_NAME(9)  = 'fluor';
    ELEMENT_FULL_NAME(11) = 'sodium';
    ELEMENT_FULL_NAME(14) = 'silicon';
    ELEMENT_FULL_NAME(15) = 'phosphore';
    ELEMENT_FULL_NAME(16) = 'sulfur';
    ELEMENT_FULL_NAME(17) = 'chlore';
    ELEMENT_FULL_NAME(20) = 'calcium';
    ELEMENT_FULL_NAME(26) = 'iron';
    ELEMENT_FULL_NAME(35) = 'brome';

!   ### Set the symbol name of elements ############################################################
                                                                                         !
    ELEMENT_SYMBOL_NAME(1)  = 'H';
    ELEMENT_SYMBOL_NAME(2)  = 'He';
    ELEMENT_SYMBOL_NAME(3)  = 'Li';
    ELEMENT_SYMBOL_NAME(4)  = 'Be';
    ELEMENT_SYMBOL_NAME(5)  = 'B';
    ELEMENT_SYMBOL_NAME(6)  = 'C';
    ELEMENT_SYMBOL_NAME(7)  = 'N';
    ELEMENT_SYMBOL_NAME(8)  = 'O';
    ELEMENT_SYMBOL_NAME(9)  = 'F';
    ELEMENT_SYMBOL_NAME(11) = 'Na';
    ELEMENT_SYMBOL_NAME(14) = 'Si';
    ELEMENT_SYMBOL_NAME(15) = 'P';
    ELEMENT_SYMBOL_NAME(16) = 'S';
    ELEMENT_SYMBOL_NAME(17) = 'Cl';
    ELEMENT_SYMBOL_NAME(20) = 'Ca';
    ELEMENT_SYMBOL_NAME(26) = 'Fe';
    ELEMENT_SYMBOL_NAME(35) = 'Br';

!   ### Set the atomic radius of elements (based on pair potentials of the pcff force field) #######
                                                                                         !
    ELEMENT_RADIUS(1)  = 1.098d0;                                                        ! IN [A]
                                                                                         !
    ELEMENT_RADIUS(2)  = 2.9d0;                                                          ! IN [A]
                                                                                         !
    ELEMENT_RADIUS(3)  = 3.2494d0;                                                       ! IN [A]
                                                                                         !
    ELEMENT_RADIUS(6)  = 4.01d0;                                                         ! IN [A]
                                                                                         !
    ELEMENT_RADIUS(7)  = 4.07d0;                                                         ! IN [A]
                                                                                         !
    ELEMENT_RADIUS(8)  = 3.533d0;                                                        ! IN [A]
                                                                                         !
    ELEMENT_RADIUS(9)  = 3.2d0;                                                          !
                                                                                         !
    ELEMENT_RADIUS(11) = 3.9624d0;                                                       !
                                                                                         !
    ELEMENT_RADIUS(14) = 4.45d0;                                                         !
                                                                                         !
    ELEMENT_RADIUS(15) = 4.2950;                                                         !
                                                                                         !
    ELEMENT_RADIUS(16) = 4.027d0;                                                        !
                                                                                         !
    ELEMENT_RADIUS(17) = 3.915d0;                                                        !

!   ELEMENT_RADIUS(20) = 
                                                                                         !
    ELEMENT_RADIUS(26) = 2.6595d0;                                                       !
                                                                                         !
    ELEMENT_RADIUS(35) = 4.215d0;                                                        !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine SET_ATOM_PROPERTIES
