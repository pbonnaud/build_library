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

subroutine APPLY_BUILD_CELLBOX(CELL_AXIS,CELL_ANGDEG)

!   ************************************************************************************************
!   **                               CREATE THE CELL BOX FOR SIMULATIONS                          **
!   ************************************************************************************************

    use module_data_in;

    use module_physical_constants;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3), intent(out) :: CELL_AXIS, CELL_ANGDEG;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    real (kind=8) :: XMAX, XMIN, YMAX, YMIN, ZMAX, ZMIN;

    real (kind=8) :: LENGTHA, LENGTHB, LENGTHC, LENGTHR;

    real (kind=8) :: VPART, VBOX;

    real (kind=8) :: LBOXA, LBOXB, LBOXC;

    real (kind=8) :: PART_DENSITY, PART_VOLUME_FRACTION;

!   ************************************************************************************************

    write(99,*) 'CREATION OF THE SIMULATION BOX :';
    write(99,*);

    XMAX = -100.0d0;
    XMIN =  100.0d0;
    YMAX = -100.0d0;
    YMIN =  100.0d0;
    ZMAX = -100.0d0;
    ZMIN =  100.0d0;

!   do i = 1, espece
!       do j = 1, NATSP_TMP(i)
!           if ( DATA_TMP(i)%NAT(j) == 'XXX' ) CYCLE;
!       end do
!   end do

    do i = 1, NSLAB_TMP
!       if ( SUB_TMP(i)%RI(1) > XMAX ) XMAX = SUB_TMP(i)%RI(1);
!       if ( SUB_TMP(i)%RI(1) < XMIN ) XMIN = SUB_TMP(i)%RI(1);
!       if ( SUB_TMP(i)%RI(2) > YMAX ) YMAX = SUB_TMP(i)%RI(2);
!       if ( SUB_TMP(i)%RI(2) < YMIN ) YMIN = SUB_TMP(i)%RI(2);
!       if ( SUB_TMP(i)%RI(3) > ZMAX ) ZMAX = SUB_TMP(i)%RI(3);
!       if ( SUB_TMP(i)%RI(3) < ZMIN ) ZMIN = SUB_TMP(i)%RI(3);
    end do

    write(99,*) 'XMAX [A] : ', XMAX;
    write(99,*) 'XMIN [A] : ', XMIN;
    write(99,*) 'YMAX [A] : ', YMAX;
    write(99,*) 'YMIN [A] : ', YMIN;
    write(99,*) 'ZMAX [A] : ', ZMAX;
    write(99,*) 'ZMIN [A] : ', ZMIN;
    write(99,*);

    LENGTHA = XMAX - XMIN;
    LENGTHB = YMAX - YMIN;
    LENGTHC = ZMAX - ZMIN;

    write(99,*) 'LENGTHA [A] : ', LENGTHA;
    write(99,*) 'LENGTHB [A] : ', LENGTHB;
    write(99,*) 'LENGTHC [A] : ', LENGTHC;
    write(99,*);

    LENGTHR = ( LENGTHA + LENGTHB ) * 0.5d0;

    write(99,*) 'LENGTHR [A] : ', LENGTHR;

    VPART = PI * LENGTHR**2 * LENGTHC * 0.25d0;

    VPART = VPART * 0.001d0;

    write(99,*) 'VPART [nm3] : ', VPART;
    write(99,*);

    LBOXA = LENGTHA + 60.0d0;
    LBOXB = LENGTHB + 60.0d0;
    LBOXC = LENGTHC + 40.0d0;

    write(99,*) 'LBOXA [A] : ', LBOXA;
    write(99,*) 'LBOXB [A] : ', LBOXB;
    write(99,*) 'LBOXC [A] : ', LBOXC;
    write(99,*);

    VBOX = LBOXA * LBOXB * LBOXC;
    VBOX = VBOX * 0.001d0;
   
    write(99,*) 'VBOX [nm3] : ', VBOX;
    write(99,*);

    PART_VOLUME_FRACTION = VPART / VBOX;
    PART_DENSITY = 1.0E27 / VBOX / Na;

    write(99,'(a28,f15.6)') 'VOLUME FRACTION       [1] : ', PART_VOLUME_FRACTION;
    write(99,'(a28,f15.6)') 'PARTICLE DENSITY [mmol/L] : ', PART_DENSITY;
    write(99,*);

    CELL_AXIS(1) = LBOXA;
    CELL_AXIS(2) = LBOXB;
    CELL_AXIS(3) = LBOXC;

    CELL_ANGDEG(1:3) = 90.0d0;

    call SET_PASSAGE();

!   stop;

end subroutine APPLY_BUILD_CELLBOX
