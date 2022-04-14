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

!subroutine CIbary(ici,a,b,c,zsign,init)

!   use data_in

!   ************************************************************************************************

!   implicit none

!   ************************************************************************************************

!   integer :: ici, init;

!   real (kind=8), intent(in) :: a, b, c, zsign;

!   ************************************************************************************************

!   integer :: e, f, h, i, j;

!   ### Parametres vecteur de translation ##########################################################
!   real (kind=8) :: x1, y1, z1;
!   real (kind=8) :: x2, y2, z2;
!   real (kind=8) :: x3, y3, z3;
!   real (kind=8) :: x4, y4, z4;
!   real (kind=8) :: x5, y5, z5;

!   real (kind=8) :: xtno, ytno, ztno;

!   real (kind=8) :: trx, try, trz;
!   real (kind=8) :: trxno, tryno, trzno;

!   real (kind=8) :: dx1, dy1, dz1;
!   real (kind=8) :: dx2, dy2, dz2;
!   real (kind=8) :: dx3, dy3, dz3;
!   real (kind=8) :: dx5, dy5, dz5;

!   real (kind=8) :: cost, sint;

!   real (kind=8) :: zosef;

!   real (kind=8) :: rij;

!   real (kind=8) :: grnd;

!   ### Critere de choix du site du CI #############################################################
!   integer :: presSi;

!   ************************************************************************************************

!   presSi = 2;

!   do while ( presSi > 0 )


!   if ( ici == init ) then

!       h = 0;

!       do while ( h == 0 )

!          j = int( grnd() * NMAX ) + 1;
!          if ( slab0(j)%nat == 'H' .and. slab0(j)%wt ==.true. .and. &
!                sign(1.0d0,slab0(j)%z) == zsign ) h = 1;

!       end do

!       write(99,'(2i4,2a3,i4)') ici, j, slab0(j)%nat,conio(1)%nat,s;
!       write(99,*) slab0(j)%x, slab0(j)%y, slab0(j)%z;
!       write(99,*)

!       ### Origine du vecteur de translation ######################################################
!       xt = slab0(j)%x;
!       yt = slab0(j)%y;
!       zt = slab0(j)%z;

!   else if ( ici > init ) then

!       write(99,*)
!       write(99,*) ' xt  : ', xt;
!       write(99,*) ' yt  : ', yt;
!       write(99,*) ' zt  : ', zt;
!       write(99,*) 
!       write(99,*) ' dxt : ', dxt;
!       write(99,*) ' dyt : ', dyt;
!       write(99,*)

!       ### determination du prochain point origine ####################################################
!       call passage(1,xt,yt,zt,a,b,c,x1,y1,z1);
!       call passage(1,dxt,dyt,zosef,a,b,c,dx1,dy1,dz1);

!       x1 = x1 + 2.0d0 * dx1;
!       y1 = y1 + 2.0d0 * dy1;

!       x1 = x1 - anint( x1 );
!       y1 = y1 - anint( y1 );

!       call passage(0,x1,y1,z1,a,b,c,xt,yt,zt);

!       write(99,*) ' xt trans : ', xt;
!       write(99,*) ' yt trans : ', yt;
!       write(99,*) ' zt trans : ', zt;
!       write(99,*)
!       write(99,*) ' N        : ', N;
!       write(99,*)

!       j = 0;

!       f = 0;

!       do while ( j == 0 .and. f < 2 )

!           do i = 1, NMAX

!               if ( slab0(i)%nat == 'H' .and. slab0(i)%wt == .true. ) then

!                   call passage(1,slab0(i)%x,slab0(i)%y,slab0(i)%z,a,b,c,x2,y2,z2);

!                   dx2 = x2 - x1;
!                   dy2 = y2 - y1;
!                   dz2 = z2 - z1;

!                   dx2 = dx2 - anint( dx2 ); 
!                   dy2 = dy2 - anint( dy2 );
!                   dz2 = dz2 - anint( dz2 );

!                   call passage(0,dx2,dy2,dz2,a,b,c,dx2,dy2,dz2);

!                   rij = dsqrt( dx2 * dx2 + dy2 * dy2 + dz2 * dz2 );

!                   if ( rij < 0.1d0 ) then
!                       j = i;
!                       write(99,*) ' j    : ', j;
!                       write(99,'(a3,3f15.5)') slab0(j)%nat, slab0(j)%x, &
!                                               slab0(j)%y, slab0(j)%z;
!                       slab0(j)%nat = conio(1)%nat;
!                       slab0(j)%q   = conio(1)%q;

!                   end if

!               end if

!           end do

!           if ( j == 0 ) then

!               if ( ici <= partconio ) then
!                   cost = abs( dxt / rt );
!                   sint = abs( dyt / rt );
!               else if ( ici > partconio ) then
!                   cost = dxt / rt;
!                   sint = abs( dyt / rt );
!               end if

!               trx = 5.03d0 * sint;
!               try = 5.03d0 * cost;
!               trz = 0.0d0;

!               call passage(1,trx,try,trz,a,b,c,trxno,tryno,trzno);

!               x1 = x1 + trxno;
!               y1 = y1 + tryno;

!               x1 = x1 - anint( x1 );
!               y1 = y1 - anint( y1 );

!               call passage(0,x1,y1,zosef,a,b,c,xt,yt,zosef);

!               write(99,*) ' Changement de ligne :'
!               write(99,*)
!               write(99,*) ' cos   : ', cost;
!               write(99,*) ' sin   : ', sint;
!               write(99,*)
!               write(99,*) ' trx   : ', trx;
!               write(99,*) ' try   : ', try;
!               write(99,*) ' trz   : ', trz;
!               write(99,*)
!               write(99,*) ' xtr   : ', xt;
!               write(99,*) ' ytr   : ', yt;
!               write(99,*) ' ztr   : ', zt;
!               write(99,*)

!               f = f + 1;

!           end if

!       end do

!       if ( f == 2 ) then
!           write(99,*) 'Pas de sol'
!           stop
!       end if

!   end if

!   call passage(1,slab0(j)%x,slab0(j)%y,slab0(j)%z,a,b,c,x3,y3,z3);

!   e = 0;
!   h = 0;

!   do while ( h == 0 .and. e < NMAX )

!       e = e + 1;

!       if ( e /= j .and. slab0(e)%wt == .true. .and. slab0(e)%nat == 'H' .and.  &
!                         sign(1.0d0,slab0(e)%z) == zsign ) then

!           call passage(1,slab0(e)%x,slab0(e)%y,slab0(e)%z,a,b,c,x4,y4,z4);

!           dx2 = x4 - x3; 
!           dy2 = y4 - y3;
!           dz2 = z4 - z3;

!           dx2 = dx2 - anint( dx2 );
!           dy2 = dy2 - anint( dy2 );
!           dz2 = dz2 - anint( dz2 );

!           call passage(0,dx2,dy2,dz2,a,b,c,dx3,dy3,dz3);

!           rij = dsqrt( dx3 * dx3 + dy3 * dy3 + dz3 * dz3 );

!           if ( ici == init ) then

!               if ( rij > 2.4d0 .and. rij < 2.6d0 ) then

!                   write(99,*) ' --> rij : ', rij; 
!                   write(99,*) ' x       : ', slab0(e)%x;
!                   write(99,*) ' y       : ', slab0(e)%y;
!                   write(99,*) ' z       : ', slab0(e)%z;
!                   write(99,*)

!                   h = 1;

!               end if

!           else

!               if ( rij > 2.4d0 .and. rij < 2.6d0 .and. sign(1.0d0,dx3) == sign(1.0d0,dxt) .and. &
!                    sign(1.0d0,dy3) == sign(1.0d0,dyt) ) then

!                   write(99,*) ' --> rij : ', rij; 
!                   write(99,*) ' x       : ', slab0(e)%x;
!                   write(99,*) ' y       : ', slab0(e)%y;
!                   write(99,*) ' z       : ', slab0(e)%z;
!                   write(99,*)

!                   h = 1;

!               end if

!           end if


!       end if

!       if ( e == NMAX ) stop 'Pas de sol.'

!   end do

!   write(99,*) e;
!   write(99,*) slab0(e)%x, slab0(e)%y, slab0(e)%z;
!   write(99,*)

!   ### Placement du CI au barycentre des charges ##################################################
!   x3 = x3 + 0.5d0 * dx2;
!   y3 = y3 + 0.5d0 * dy2;
!   z3 = z3 + 0.5d0 * dz2;

!   x3 = x3 - anint( x3 );
!   y3 = y3 - anint( y3 );
!   z3 = z3 - anint( z3 );

!   presSi = 0;

!   if ( ici == init ) then

!       do i = 1, NMAX

!           if ( slab0(i)%nat == 'Si' ) then

!               call passage(1,slab0(i)%x,slab0(i)%y,slab0(i)%z,a,b,c,x5,y5,z5);

!               dx5 = x5 - x3;
!               dy5 = y5 - y3;
!               dz5 = z5 - z3;

!               dx5 = dx5 - anint( dx5 );
!               dy5 = dy5 - anint( dy5 );
!               dz5 = dz5 - anint( dz5 );

!               call passage(0,dx5,dy5,dz5,a,b,c,dx5,dy5,dz5);

!               rij = dsqrt( dx5 * dx5 + dy5 * dy5 + dz5 * dz5 );

!               if ( rij < 2.5d0 ) presSi = presSi + 1;

!           end if

!       end do

!   end if     

!   if ( presSi == 0 ) then

!       ### Remplacment du premier hydrogene #######################################################
!       slab0(j)%nat = conio(1)%nat;
!       slab0(j)%q   = conio(1)%q;

!       ### Elimination du second hydrogene ########################################################
!       slab0(e)%wt = .false.;
!       s = s - 1;

!       ### Composantes du vecteur de translation ######################################################
!       dxt = dx3;
!       dyt = dy3;
!       dzt = dz3;

!       rt  = dsqrt( dxt * dxt + dyt * dyt + dzt * dzt );

!       call passage(0,x3,y3,z3,a,b,c,slab0(j)%x,slab0(j)%y,slab0(j)%z);

!       write(99,*) slab0(j)%x, slab0(j)%y, slab0(j)%z;
!       write(99,*)
!       write(99,*) ' rt   : ', rt;
!       write(99,*)
!       write(99,*) ' *********************************************'
!       write(99,*)

!   end if

!   end do

!end subroutine CIbary

subroutine VECT_PRODUCT(VECTU,VECTV,VECTW)

!   ************************************************************************************************
!   **                                COMPUTE VECTOR PRODUCT                                      **
!   ************************************************************************************************

    implicit none

!   ************************************************************************************************

    real (kind=8), dimension(1:3), intent(in) :: VECTU, VECTV;

    real (kind=8), dimension(1:3), intent(out) :: VECTW;

!   ***********************************************************************************************

    VECTW(1) = VECTU(2) * VECTV(3) - VECTU(3) * VECTV(2);
    VECTW(2) = VECTU(3) * VECTV(1) - VECTU(1) * VECTV(3);
    VECTW(3) = VECTU(1) * VECTV(2) - VECTU(2) * VECTV(1);

end subroutine VECT_PRODUCT




