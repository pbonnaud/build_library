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

subroutine CHECK_DISTANCES(PASSA,PASSB)

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;


    integer (kind=4) :: f, g, h, i, j, k, m;

    real (kind=8) :: grnd;

!   ************************************************************************************************

    integer (kind=4) :: ISTART, JSTART, ITOOCLOSE, NIMOL;

    real (kind=8) :: RIJ;

    real (kind=8), dimension(1:3) :: DRIJ, DRIJNO;

    real (kind=8), dimension(1:3,1:100) :: DRIJ_MOL; 

    integer (kind=4), dimension(1:2) :: HALFNBRE;

    integer (kind=4), dimension(1:100) :: IMOL;

!   ************************************************************************************************

    write(99,*);
    write(99,*) 'NBDOPLE  : ', NBDOPLE(1:2);
    HALFNBRE(1:2) = NBDOPLE(1:2) / 2;
    write(99,*) 'HALFNBRE : ', HALFNBRE(1:2); 

    do h = 1, 2           
        do i = 1, NBFINAL(h) !HALFNBRE(h) !   NBDOPLE(h)
            do g = 1, 2
                do j = 1, NBFINAL(g)    !JSTART, NBDOPLE(g)
                    DRIJ(1:3)   = RDOPLE(j,1:3,g) - RFINAL(i,1:3,h);
                    DRIJNO(1:3) = MATMUL(PASSB(1:3,1:3),DRIJ(1:3));
                    DRIJNO(1:3) = DRIJNO(1:3) - ANINT(DRIJNO(1:3));
                    DRIJ(1:3)   = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));
                    RIJ = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));
                    if ( RIJ < 0.2d0 ) then
                        write(99,*);
                        write(99,'(a70)') '!!! ATOMS ARE TOO CLOSE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!';
                        write(99,*) 'RIJ : ', RIJ;
                        write(99,'(i4,1x,i4,1x,a3,3f15.6)') h, i, NATFINAL(i,h), RFINAL(i,1:3,h);
                        write(99,'(i4,1x,i4,1x,a3,3f15.6)') g, j, NATDOPLE(j,g), RDOPLE(j,1:3,g);
                        write(99,*);
                        if ( NATDOPLE(j,g) == 'Ow' .AND. NATDOPLE(j+1,g) == 'Hw' .AND. NATDOPLE(j+2,g) == 'Hw' ) then
                            IMOL(1) = j;
                            IMOL(2) = j+1;
                            IMOL(3) = j+2;
                            NIMOL = 3;
                        else if ( NATDOPLE(j-1,g) == 'Ow' .AND. NATDOPLE(j,g) == 'Hw' .AND. NATDOPLE(j+1,g) == 'Hw' ) then
                            IMOL(1) = j-1;
                            IMOL(2) = j;
                            IMOL(3) = j+1;
                            NIMOL = 3;
                        else if ( NATDOPLE(j-2,g) == 'Ow' .AND. NATDOPLE(j-1,g) == 'Hw' .AND. NATDOPLE(j,g) == 'Hw' ) then
                            IMOL(1) = j-2;
                            IMOL(2) = j-1;
                            IMOL(3) = j;
                            NIMOL = 3;
                        else if ( NATDOPLE(j,g) == 'Cw' ) then
                            IMOL(1) = j;
                            NIMOL = 1;
                        else if ( NATDOPLE(j,g) == 'O' .OR. NATDOPLE(j,g) == 'Si' ) then

                        else
                            write(99,*) 'CYCLE';
                            CYCLE;
                        end if

                        do m = 1, NIMOL
                            write(99,*) NATDOPLE(IMOL(m),g), RDOPLE(IMOL(m),1:3,g);
                        end do
!                       write(99,*) NATDOPLE(IMOL(2),g), RDOPLE(IMOL(2),1:3,g);
!                       write(99,*) NATDOPLE(IMOL(3),g), RDOPLE(IMOL(3),1:3,g);

                        if (NIMOL > 1 ) then
                            do m = 2, NIMOL
                                DRIJ_MOL(1:3,m) = RDOPLE(IMOL(m),1:3,g) - RDOPLE(IMOL(1),1:3,g);
                            end do
!                           DRIJ_MOL(1:3,2) = RDOPLE(IMOL(3),1:3,g) - RDOPLE(IMOL(1),1:3,g);
                        end if

                        ITOOCLOSE = 0;

                        do while ( ITOOCLOSE == 0 )
                            DRIJNO(1) = ( grnd() - 0.5d0 );
                            DRIJNO(2) = ( grnd() - 0.5d0 );
                            DRIJNO(3) = ( grnd() - 0.5d0 );

                            RDOPLE(IMOL(1),1:3,g) = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));
                            if ( NIMOL > 1 ) then
                                do m = 2, NIMOL 
                                    RDOPLE(IMOL(m),1:3,g) = RDOPLE(IMOL(1),1:3,g) + DRIJ_MOL(1:3,m);
                                end do
                            end if

!                           RDOPLE(IMOL(3),1:3,g) = RDOPLE(IMOL(1),1:3,g) + DRIJ_MOL(1:3,2);

                            ITOOCLOSE = 1;
                            do f = 1, 2
                                do k = 1, NBFINAL(f)
                        !           do m = i, i+2
                                    do m = 1, NIMOL
                                        DRIJ(1:3)   = RFINAL(k,1:3,f) - RDOPLE(IMOL(m),1:3,h);
                                        DRIJNO(1:3) = MATMUL(PASSB(1:3,1:3),DRIJ(1:3));
                                        DRIJNO(1:3) = DRIJNO(1:3) - ANINT(DRIJNO(1:3));
                                        DRIJ(1:3)   = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));
                                        RIJ = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));
                                        if ( RIJ < 0.9999d0 ) then
                                            ITOOCLOSE = 0;
                                            write(99,*) 'RIJ 1 : ', RIJ;
                                        end if
                        !               if ( ) then
                                    end do
                                end do
                            end do
                        end do
                    end if
                end do
            end do
        end do
    end do

!   write(99,*) 'OK DISTANCES';

!   stop;

end subroutine CHECK_DISTANCES
