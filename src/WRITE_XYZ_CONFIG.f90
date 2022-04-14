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

subroutine WRITE_XYZ_CONFIG(icanal,CHEXT,        &
                                   NATOM,        &
                                   CELL_AXIS,    &
                                   CELL_ANGDEG,  &
                                   CONFIG_NAT,   &
                                   CONFIG_RIJ)

!   ************************************************************************************************
!   **                                 WRITE XYZ FILES                                            **
!   ************************************************************************************************
!   **                                                                                            **
!   ** INPUTS:                                                                                    **
!   ** ------                                                                                     **
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                          **
!   ** CHEXT  : NAME OF THE LAMMPS CONFIGURATION                                                  **
!   **                                                                                            **
!   ** NATOM       : NUMBER OF ATOMS                                                              **
!   ** CONFIG_NAT  : NATURE OF ATOMS IN THE SIMULATION BOX                                        **
!   ** CONFIG_RIJ  : POSITION OF ATOMS IN THE SIMULATION BOX                                      **
!   ** CELL_AXIS   : CELL DIMENSIONS                                                              **
!   ** CELL_ANGDEG : CELL ANGLES                                                                  **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=150), intent(in) :: CHEXT;

    integer (kind=4), intent(in) :: NATOM;

    real (kind=8), dimension(1:3), intent(in) :: CELL_AXIS, CELL_ANGDEG;

    character (len=20), dimension(1:100000), intent(in) :: CONFIG_NAT;

    real (kind=8), dimension(1:3,1:100000), intent(in) :: CONFIG_RIJ;

!   ************************************************************************************************

    integer (kind=4) :: i;

!   ************************************************************************************************

    open(103,file=TRIM(CHEXT)//'.xyz');
    write(103,'(i8)') NATOM;
    write(103,'(6f15.6)') CELL_AXIS(1:3), CELL_ANGDEG(1:3);

    do i = 1, NATOM
        write(103,'(a3,3f15.6)') CONFIG_NAT(i), CONFIG_RIJ(1:3,i);
    end do

    close(103);

end subroutine WRITE_XYZ_CONFIG

