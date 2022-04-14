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

subroutine CHECK_SILICACHAINS(PASSA,PASSB)

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

    integer (kind=4) :: f, g, h, i, j, k;

    integer (kind=4) :: ipore, iatfnd;

    real (kind=8) :: grnd;

!   ************************************************************************************************

    integer (kind=4) :: INNCOUNT, NTOT_NCT, ISUPERPOSE;

    real (kind=8) :: RIJ, RIJ_WOPDB, HALF_LBOX, DIFF_RIJ;

    real (kind=8), dimension(1:3) :: DRIJ, DRIJNO, RINEW;

    integer (kind=4), allocatable, dimension(:,:) :: LABEL_FINAL; 

    integer (kind=4), allocatable, dimension(:,:,:) :: NEIGHBOR_NBRE; 
    integer (kind=4), allocatable, dimension(:,:,:,:) :: NEIGHBOR_LIST; 

    integer (kind=4) :: NADDED_FINAL, RDCHX_OS;

    integer (kind=4), dimension(1:2) :: NBRE_ADDED, NBRE_OS;

    integer (kind=4), dimension(1:2,1:100000) :: LIST_OS;

    real (kind=8), allocatable, dimension(:,:,:) :: COORD_ADDED;

    character (len=3), allocatable, dimension(:,:) :: NAT_ADDED; 

!   ************************************************************************************************

    write(*,*);

    HALF_LBOX = 100.0d0;

    do ipore = 1, 2 
        DRIJNO(1:3)   = 0.0d0;
        DRIJNO(ipore) = 1.0d0;
        DRIJ(1:3) = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));
        if ( 0.5d0 * DRIJ(ipore) < HALF_LBOX ) then
            HALF_LBOX = 0.5d0 * DRIJ(ipore);
!           write(*,*) 'XXX', HALF_LBOX;
        end if
    end do
    write(*,*) 'HALF_LBOX : ', HALF_LBOX; 
    write(*,*);

    NTOT_NCT = 0;

    allocate(NEIGHBOR_LIST(1:10000,1:2,1:2,1:10));
    allocate(NEIGHBOR_NBRE(1:10000,1:2,1:2));

    NEIGHBOR_NBRE(1:10000,1:2,1:2) = 0;
    NEIGHBOR_LIST(1:10000,1:2,1:2,1:10) = 0;

    allocate(LABEL_FINAL(1:NWORK_EXTD,1:2));
    LABEL_FINAL(1:NWORK_EXTD,1:2) = 0;

    NBRE_ADDED(1:2) = 0;

    allocate(COORD_ADDED(1:100000,1:3,1:2),NAT_ADDED(1:100000,1:2));

    COORD_ADDED(1:100000,1:3,1:2) = 0.0d0;
    NAT_ADDED(1:100000,1:2) = 'XXX';

    do h = 1, 2
        do i = 1, NBFINAL(h)
            if ( NATFINAL(i,h) == 'Si' ) then
                INNCOUNT = 0;
                do g = 1, 2
                    do j = 1, NBFINAL(g)
                        if ( NATFINAL(j,g) == 'O' ) then
                            DRIJ(1:3) = RFINAL(j,1:3,g) - RFINAL(i,1:3,h);
                            RIJ_WOPDB = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));                        
                            DRIJNO(1:3) = MATMUL(PASSB(1:3,1:3),DRIJ(1:3)); 
                            DRIJNO(1:3) = DRIJNO(1:3) - ANINT( DRIJNO(1:3) );
                            DRIJ(1:3) = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));
                            RIJ = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));
                            if ( RIJ < 2.5d0 ) then
                                INNCOUNT = INNCOUNT + 1;
                                NEIGHBOR_NBRE(i,h,g) = INNCOUNT;
                                NEIGHBOR_LIST(i,h,g,INNCOUNT) = j;
                                LABEL_FINAL(j,g) = LABEL_FINAL(j,g) + 1;
                                if ( LABEL_FINAL(j,g) == 1 ) then
                                    RFINAL(j,1:3,g) = RFINAL(i,1:3,h) + DRIJ(1:3);
                                else if ( LABEL_FINAL(j,g) == 2 ) then
                                    DIFF_RIJ = RIJ_WOPDB - RIJ;
                                    if ( DIFF_RIJ > 1E-4 ) then
!                                       write(*,*) 'XXXXXXX', RIJ, RIJ_WOPDB;
                                        NBRE_ADDED(g) = NBRE_ADDED(g) + 1;
                                        NAT_ADDED(NBRE_ADDED(g),g) = NATFINAL(j,g);
                                        COORD_ADDED(NBRE_ADDED(g),1:3,g) = RFINAL(i,1:3,h) + DRIJ(1:3);
                                    end if
                                end if
!                               write(*,'(4i6,f15.6)') h, i, g, j, RIJ;
                            end if    
                        end if
                    end do
                    if ( INNCOUNT > 0 ) write(*,'(8i6)') h, i, g, INNCOUNT, NEIGHBOR_LIST(i,h,g,1:INNCOUNT);
                end do

                if ( INNCOUNT < 4 ) then
                    write(*,*) 'INCOMPLETE TETRAHEDRA';
                    NTOT_NCT = NTOT_NCT + 1;
                else if ( INNCOUNT > 4 ) then
                    write(*,*) 'TOO MANY NN STOP';
                    stop;
                end if 
            end if
        end do
    end do

    write(*,*);
    write(*,*) 'TOTAL NUMBER OF INCOMPLETE TETRAHEDRA : ', NTOT_NCT;
    write(*,*);
    write(*,*) 'LIST OF ATOMS TO BE ADDED : ', NBRE_ADDED(1:2);
    write(*,*);

    do g = 1, 2
        NADDED_FINAL = NBFINAL(g);
        do i = 1, NBRE_ADDED(g)
!           write(*,'(2i6,2x,a3,2x,3f15.6)') g, i, NAT_ADDED(i,g), COORD_ADDED(i,1:3,g);
            NWORK = NWORK + 1;
            NADDED_FINAL = NADDED_FINAL + 1;
            NATFINAL(NADDED_FINAL,g) = NAT_ADDED(i,g);
            RFINAL(NADDED_FINAL,1:3,g) = COORD_ADDED(i,1:3,g);            
        end do
        NBFINAL(g) = NADDED_FINAL; 
!       write(*,*) NBFINAL(g), NWORK; 
    end do

    write(*,*) 'NWORK : ', NWORK;
    write(*,*);
    write(*,*) 'MAKE LIST OF Os ATOMS : ';
    write(*,*);

    NBRE_OS(1:2) = 0;
    LIST_OS(1:2,1:100000) = 0;

    do h = 1, 2
        do i = 1, NBFINAL(h)
            if ( NATFINAL(i,h) == 'Os' ) then
                NBRE_OS(h) = NBRE_OS(h) + 1;
                LIST_OS(h,NBRE_OS(h)) = i;
            end if
        end do
        write(*,*) 'NOS : ', NBRE_OS(h);
    end do

!   CHOSE RANDONLY Os ATOMS TO BE REMOVED    
    do h = 1, 2
        do i = 1, NBRE_ADDED(h)
            iatfnd = 0;
            do while ( iatfnd == 0 )
                RDCHX_OS = INT( NBRE_OS(h) * grnd() ) + 1;
                if ( NATFINAL(LIST_OS(h,RDCHX_OS),h) /= 'XXX' ) then
                    NATFINAL(LIST_OS(h,RDCHX_OS),h) = 'XXX';
                    do j = LIST_OS(h,RDCHX_OS), NBFINAL(h)-1
                        NATFINAL(j,h) = NATFINAL(j+1,h);
                        RFINAL(j,1:3,h) = RFINAL(j+1,1:3,h);
                    end do
                    NWORK = NWORK - 1;
                    NBFINAL(h) = NBFINAL(h) - 1
                    iatfnd = 1;
                end if
            end do
!            write(*,*) i, RDCHX_OS; 
        end do
    end do

    write(*,*) 'NWORK : ', NWORK;
    write(*,*);


    deallocate(NEIGHBOR_LIST,NEIGHBOR_NBRE,COORD_ADDED,NAT_ADDED);


end subroutine CHECK_SILICACHAINS

