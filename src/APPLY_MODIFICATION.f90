!   //////////////////////////////////////////////////////////////////////////////////////////////////
!   //                                                                                              //
!   //  Copyright (2022) Patrick A. Bonnaud                                                         //
!   //                                                                                              //
!   //  This file is part of BUILDLIBRARY (Build Molecular Models for a Library of Molecules).      //
!   //                                                                                              //
!   //  BUILDLIBRARY is free software; you can redistribute it and/or modify it under the terms of  //
!   //  the GNU General Public License as published by the Free Software Foundation; either         //
!   //  version 2 of the License, or (at your option) any later version.                            //
!   //                                                                                              //
!   //  BUILDLIBRARY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;   //
!   //  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   //
!   //  See the GNU General Public License for more details.                                        //
!   //                                                                                              //
!   //  You should have received a copy of the GNU General Public License along with this program.  //
!   //  If not, see <http://www.gnu.org/licenses/>.                                                 //
!   //                                                                                              //
!   //////////////////////////////////////////////////////////////////////////////////////////////////

subroutine APPLY_MODIFICATION()

!   ************************************************************************************************
!   **                         APPLY MODIFICATION TO THE SIMULATION BOX                           **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: i;

    integer (kind=4) :: iespece, insp;

    real (kind=8) :: RGIJ, RGIJ_NEW;

    real (kind=8), dimension(1:3) :: DRIJ, DRGIJ, VECTU, COMA_NEW;

!   ************************************************************************************************

!   DRGIJ(1:3) = COMA_COORD(1:3,1) - COMA_COORD(1:3,ICOMA);

!   call APPLY_PBC(DRGIJ(1:3),DRGIJ(1:3));

!   RGIJ = DSQRT( DOT_PRODUCT(DRGIJ(1:3),DRGIJ(1:3)) );

!   write(99,*) DRGIJ(1:3);
!   write(99,*) 'RGIJ [A] : ', RGIJ;
!   write(99,*);

!   VECTU(1:3) = DRGIJ(1:3) / RGIJ;

!   write(99,*) VECTU(1:3);
!   write(99,*) DSQRT( DOT_PRODUCT(VECTU(1:3),VECTU(1:3)));
!   write(99,*);

!   RGIJ_NEW = RGIJ + DRI_TRANS;

!   write(99,*) 'RGIJ_NEW : ', RGIJ_NEW; 

!   DRGIJ(1:3) = RGIJ_NEW * VECTU(1:3);

!   write(99,*) DRGIJ(1:3);
!   write(99,*);

!   COMA_NEW(1:3) = COMA_COORD(1:3,ICOMA) + DRGIJ(1:3);

!   call APPLY_PBC(COMA_NEW(1:3),COMA_NEW(1:3));

!   write(99,*) COMA_NEW(1:3);
!   write(99,*) COMA_COORD(1:3,1);

!   do i = 1, NLIST_SPECIES_MODIF
!       do insp = 1, NATSP(LIST_SPECIES_MODIF(i))

!           DRIJ(1:3) = DATA(LIST_SPECIES_MODIF(i))%RI(1:3,insp) - COMA_COORD(1:3,1);
!           call APPLY_PBC(DRIJ(1:3),DRIJ(1:3));

!           DATA(LIST_SPECIES_MODIF(i))%RI(1:3,insp) = COMA_NEW(1:3) + DRIJ(1:3);
!           call APPLY_PBC(DATA(LIST_SPECIES_MODIF(i))%RI(1:3,insp), &
!                          DATA(LIST_SPECIES_MODIF(i))%RI(1:3,insp));
!           read(11,*) DATA(iespece)%NAT(insp), DATA(iespece)%RI(1:3,insp);
!       end do
!   end do

!   COMA_COORD_FINAL(1:3,1) = COMA_NEW(1:3);    
!   COMA_COORD_FINAL(1:3,ICOMA) = COMA_COORD(1:3,ICOMA);    

end subroutine APPLY_MODIFICATION



