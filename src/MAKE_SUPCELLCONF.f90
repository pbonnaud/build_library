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

subroutine MAKE_SUPCELLCONF()

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: g, h, i, j, k;

    integer (kind=4) :: iatsc, imax;

    integer (kind=4) :: icheck_silica_chains;

    integer (kind=4) :: ADD_NLABEL;

    integer (kind=4), dimension(1:10) :: ADD_LABEL, ADD_IONS;

    real (kind=8) :: DTRBIN, TRBIN;

    real (kind=8) :: CSTE01;

    real (kind=8) :: THETA, COS_THETA, SIN_THETA, TRANS_INIT;

    real (kind=8) :: NORM_VECTU, NORM_VECTV, NORM_VECTW;
    real (kind=8), dimension(1:3) :: VECTU, VECTV, VECTW;

    real (kind=8), dimension(1:3) :: CELL_AXIS, CELL_ANGDEG;

    real (kind=8), dimension(1:3,1:3) :: PASSC;

    real (kind=8) :: NORM_VECTU_PRIME, NORM_VECTV_PRIME, NORM_VECTW_PRIME;
    real (kind=8), dimension(1:3) :: VECTU_PRIME, VECTV_PRIME, VECTW_PRIME;

    real (kind=8), dimension(1:3) :: TRANSR, RIJ, RIJ_PRIME, RINO, RION, RCOM, RCOM_IDEAL;
    real (kind=8), dimension(1:2,1:3) :: COMRON, COMRNO;

    real (kind=8), dimension(1:2,1:2,1:3) :: PART_LIMITS;

    real (kind=8), dimension(1:3,1:3) :: PASSA, PASSB;

    character (len=3) :: CHBIN;

    character (len=50) :: CHCONF_END_NAME;

!   ### PARAMETERS THAT WE DON'T CARE ABOUT ########################################################
    integer (kind=4) :: IOSEF1;

    real (kind=8) :: ROSEF1;

!   ************************************************************************************************

    write(99,'(a70)') '**********************************************************************';
    write(99,'(a70)') '**                   CREATION OF THE SUPER CELL                     **';
    write(99,'(a70)') '**********************************************************************';
    write(99,*);

    NMULTI = MULTI(1) * MULTI(2) * MULTI(3);

    write(99,'(a9,3i4)') 'MULTI  : ', MULTI(1:3);
    write(99,'(a9,i4)')  'NMULTI : ', NMULTI;
    write(99,*);

    NWORK = NMULTI * NTOT_ATOMS;
    NBATLIST(1:NBNATIN) = NBATLIST(1:NBNATIN) * NMULTI;

    N(1:NLABEL) = N(1:NLABEL) * NMULTI;

    write(99,*) 'NEW NBRES OF ATOMS : ';
    IOSEF1 = 0;

    do i = 1, NBNATIN
        IOSEF1 = IOSEF1 + NBATLIST(i);
        write(99,*) NATLIST(i), NBATLIST(i), IOSEF1;
    end do

    write(99,*);
    write(99,*) 'NWORK : ', NWORK;
    write(99,*);

    NTOT_ATOMS_FINAL = NWORK;

    CELL_AXIS(1:3) = CELL_AXIS(1:3) * MULTI(1:3);

    call SET_PASSAGE(99,CELL_AXIS(1:3),CELL_ANGDEG(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));

    NWORK_EXTD = NWORK;
    if ( HPORE(1) > 0.0d0 .OR. HPORE(2) > 0.0d0 .OR. HPORE(3) > 0.0d0 ) NWORK_EXTD = NWORK_EXTD + 1000; 

    write(99,*) 'NWORK_EXTD : ', NWORK_EXTD;

    call APPLY_REPLICA_CREATION();

    if ( IPLATELET == 1 ) then
        call APPLY_TRANSLATION();
        call SET_CENTROSYM();
        call APPLY_PLATELET();
        call APPLY_REBUILD_WATER();
        call APPLY_REBUILD_WATER_Hw();
        call APPLY_REBUILD_SILICA_CHAINS();
        call CHECK_ELECTRONEUTRALITY(ADD_NLABEL,ADD_LABEL(1:10),ADD_IONS(1:10));
        call APPLY_BUILD_CELLBOX(CELL_AXIS(1:3),CELL_ANGDEG(1:3));
        call APPLY_ION_INSERTION(ADD_NLABEL,       &
                                 ADD_LABEL(1:10),  &
                                 ADD_IONS(1:10),   &
                                 CELL_AXIS(1:3),   &
                                 CELL_ANGDEG(1:3), &
                                 PASSA(1:3,1:3),   &
                                 PASSB(1:3,1:3));

        call WRITE_FINAL_CONFIG(CELL_AXIS(1:3),   &
                                CELL_ANGDEG(1:3), &
                                PASSA(1:3,1:3),   &
                                PASSB(1:3,1:3));

        ICOMA = 1; 
        COMA_COORD_FINAL(1:3,ICOMA) = 0.0d0;

        call WRITE_FIELD();
        call WRITE_DATA();
    end if

!   icheck_silica_chains = 0;

!   do i = 1, NLABEL
!       if ( ATOM_LABEL(i) == 'Si' ) icheck_silica_chains = icheck_silica_chains + 1; 
!       if ( ATOM_LABEL(i) == 'O' ) icheck_silica_chains = icheck_silica_chains + 1; 
!   end do

!   write(99,*) 'icheck_silica_chains : ', icheck_silica_chains;
!   write(99,*);

!   if ( icheck_silica_chains == 2 ) call CHECK_SILICACHAINS();

!   CELL_AXIS(1:3) = CELL_AXIS(1:3) + HPORE(1:3); ! ADD PORE SPACE
!   call SET_PASSAGE();                           ! SET THE NEW CELL PARAMETERS

!   CELL_AXIS(1:2) = 50.0d0 + 40.0d0;
!   CELL_AXIS(3) = 120.0d0;

!   CELL_ANGDEG(1:3) = 90.0d0;

!   call SET_PASSAGE();


!   call SET_PBC()  ! PUT BACK PERIODIC BOUNDARY CONDITIONS

!   CHCONF_END_NAME = '001';
!   if ( idople == 2 ) CHCONF_END_NAME = trim(CHCONF_END_NAME)//'-SINGLE';
!   CHCONF_END_NAME = trim(CHCONF_END_NAME)//'.xyz'; 
    
!   write(99,'(a18,a30)') 'CHCONF_END_NAME : ', trim(CHCONF_END_NAME);

!   if ( icheck_silica_chains == 2 ) call SET_CONFIG(CHCONF_END_NAME);

!   COMRON(1:2,1:3) = 0.0d0;

!   if ( icheck_silica_chains == 2 ) call WRITE_FIELD(0,CHBIN,COMRON(1:2,1:3)) ! WRITE FIELD FILE FOR MC SIMULATIONS

!   if ( idople == 2 ) then   ! DUPLICATE AND TRANSLATE THE SECOND PARTICLE

!       NBREDOPLE(1:100) = 0;

!       CENTER OF MASS OF THE TWO PARTICLES ARE IDENTICAL 
!       COMRON(1:2,1:3) = 0.0d0;

!       write(99,*);

!       allocate(NATDOPLE(1:NWORK,1:2),RDOPLE(1:NWORK,1:3,1:2));

!       NWORK = NWORK * 2;

!       do i = 1, 2  ! DUPLICATE THE PARTICLE IN A NEW TABLE      
!           NATDOPLE(1:NBFINAL(i),i) = NATFINAL(1:NBFINAL(i),i);  
!           RDOPLE(1:NBFINAL(i),1:3,i) = RFINAL(1:NBFINAL(i),1:3,i);
!       end do

!       APPLY ROTATION TO THE DUPLICATED MOLECULE
!       if ( IROTA == 1 ) then
!           do h = 1, 2
!               do i = 1, NBFINAL(h)
!                   call SET_ROTATION(RDOPLE(i,1:3,h),COMRON(2,1:3),RDOPLE(i,1:3,h));
!               end do
!           end do
!       end if

!       APPLY TRANSLATION TO THE DUPLICATED MOLECULE
!       if ( ITRAN == 1 ) then
!           do h = 1, 2
!               do i = 1, NBFINAL(h)
!                   call SET_TRANSLATION(RDOPLE(i,1:3,h),RDOPLE(i,1:3,h));
!               end do
!           end do
!           call SET_TRANSLATION(COMRON(2,1:3),COMRON(2,1:3));
!       end if 

!       write(99,'(a13,3f12.5)') 'COM 1 [A] : ', COMRON(1,1:3);
!       write(99,'(a13,3f12.5)') 'COM 2 [A] : ', COMRON(2,1:3);

!       FIND CARTESIAN COORDINATE LIMITS OF BOTH PARTICLE
!       PART_LIMITS(1:2,1:2,1:3) = 0.0d0;

!       do h = 1, 2
!           do i = 1, NBFINAL(h)
!               if ( NATFINAL(i,h) == 'O' .OR. NATFINAL(i,h) == 'Si' ) then
!                   do j = 1, 3
!                       if ( RFINAL(i,j,h) < PART_LIMITS(1,1,j) ) PART_LIMITS(1,1,j) = RFINAL(i,j,h);
!                       if ( RFINAL(i,j,h) > PART_LIMITS(1,2,j) ) PART_LIMITS(1,2,j) = RFINAL(i,j,h);
!                       if ( RDOPLE(i,j,h) < PART_LIMITS(2,1,j) ) PART_LIMITS(2,1,j) = RDOPLE(i,j,h);
!                       if ( RDOPLE(i,j,h) > PART_LIMITS(2,2,j) ) PART_LIMITS(2,2,j) = RDOPLE(i,j,h);
!                   end do
!               end if
!           end do
!       end do

!       do i = 1, 2
!           write(99,'(a20,i4)') 'SPACE LIMITS PART : ', i;
!           do j = 1, 3
!               write(99,*) j, PART_LIMITS(i,1:2,j)
!           end do
!           write(99,*);
!       end do

!       EXTEND SIMULATION CELL
!       CELL_AXIS(3) = CELL_AXIS(3) * 2;        
!       call SET_PASSAGE();
!       TRANSR(1:3) = 0.0d0;

!       TRANSLATE THE FIRST PARTICLE
!       TRANSR(1:3) = (/ 0.0d0, 0.0d0, -0.5d0/);
!       COMRNO(1,1:3) = MATMUL(PASSB(1:3,1:3),COMRON(1,1:3));
!       COMRNO(1,1:3) = COMRNO(1,1:3) + TRANSR(1:3);
!       COMRON(1,1:3) = MATMUL(PASSA(1:3,1:3),COMRNO(1,1:3));
!       write(99,'(a13,3f12.5)') 'COM 1 [A] : ', COMRON(1,1:3);

!       do h = 1, 2
!           do i = 1, NBFINAL(h)
!               RINO(1:3)  = MATMUL(PASSB(1:3,1:3),RFINAL(i,1:3,h));
!               RINO(1:3)  = RINO(1:3) + TRANSR(1:3);
!               RINO(1:3)  = RINO(1:3) - ANINT( RINO(1:3) );
!               RFINAL(i,1:3,h) = MATMUL(PASSA(1:3,1:3),RINO(1:3));  
!           end do
!       end do

!       SET INITIAL POSITION OF COM OF PARTICLE 2
!       COMRNO(2,1:3) = MATMUL(PASSB(1:3,1:3),COMRON(2,1:3));

!       CSTE01 = ( PART_LIMITS(1,2,3) + PART_LIMITS(2,2,3) ) / CELL_AXIS(3);

!       COMRNO(2,3) = - 0.5d0 + CSTE01;
!       COMRON(2,1:3) = MATMUL(PASSA(1:3,1:3),COMRNO(2,1:3)); 
!       write(99,'(a13,3f12.5)') 'COM 2 [A] : ', COMRON(2,1:3);
!       write(99,*);

!       VECTU_PRIME(1:3) = COMRON(2,1:3) - COMRON(1,1:3);
!       NORM_VECTU_PRIME = DSQRT( DOT_PRODUCT(VECTU_PRIME(1:3),VECTU_PRIME(1:3)) );
!       write(99,*) 'NORM_VECTU_PRIME (BRUT) : ', NORM_VECTU_PRIME;
        
!       TRANS_INIT = NORM_VECTU_PRIME; ! SAVE VALUE FOR LATER TRANSLATION OF PARTICLE 2;

!       RINO(1:2) = 0.0d0;
!       RINO(3) = -0.5d0 + ( TRANS_INIT - 2.0d0 ) / CELL_AXIS(3);

!       RCOM_IDEAL(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!       write(99,'(a13,3f12.5)') 'RCOM_IDEAL : ',  RCOM_IDEAL(1:3);

!       VECTU_PRIME(1:3) = RCOM_IDEAL(1:3) - COMRON(1,1:3);
!       NORM_VECTU_PRIME = DSQRT( DOT_PRODUCT(VECTU_PRIME(1:3),VECTU_PRIME(1:3)) );
!       write(99,*) 'NORM_VECTU_PRIME (BRUT) : ', NORM_VECTU_PRIME;

!       VECTU_PRIME(1:3) = VECTU_PRIME(1:3) / NORM_VECTU_PRIME;
!       NORM_VECTU_PRIME = DSQRT( DOT_PRODUCT(VECTU_PRIME(1:3),VECTU_PRIME(1:3)) );
!       write(99,*) 'NORM_VECTU_PRIME (NORM) : ', NORM_VECTU_PRIME;

!       ROTATION : USE OF THE INVARIANT ROTATION METHOD
!       VECTU(1:3) = COMRON(2,1:3) - COMRON(1,1:3);
!       NORM_VECTU = DSQRT( DOT_PRODUCT(VECTU(1:3),VECTU(1:3)) );
!       write(99,*) 'NORM_VECTU (BRUT) : ', NORM_VECTU;

!       VECTU(1:3) = VECTU(1:3) / NORM_VECTU;
!       NORM_VECTU = DSQRT( DOT_PRODUCT(VECTU(1:3),VECTU(1:3)) );
!       write(99,*) 'NORM_VECTU (NORM) : ', NORM_VECTU;

!       COS_THETA = DOT_PRODUCT(VECTU(1:3),VECTU_PRIME(1:3));
!       THETA = ACOS( COS_THETA );

!       write(99,'(a41,3f15.6)') 'COS(THETA) / THETA (RAD) / THETA (DEG) : ', COS_THETA, THETA, THETA * 360.0d0 / twopi;

!       call VECT_PRODUCT(VECTU(1:3),VECTU_PRIME(1:3),VECTW(1:3));
!       NORM_VECTW = DSQRT( DOT_PRODUCT(VECTW(1:3),VECTW(1:3)) );

!       write(99,*) 'NORM_VECTW (BRUT) : ', NORM_VECTW;
!       VECTW(1:3) = VECTW(1:3) / NORM_VECTW;
!       NORM_VECTW = DSQRT( DOT_PRODUCT(VECTW(1:3),VECTW(1:3)) );
!       write(99,*) 'NORM_VECTW (NORM) : ', NORM_VECTW;

!       VECTW_PRIME(1:3) = VECTW(1:3);

!       call VECT_PRODUCT(VECTW(1:3),VECTU(1:3),VECTV(1:3));
!       NORM_VECTV = DSQRT( DOT_PRODUCT(VECTV(1:3),VECTV(1:3)) );
!       write(99,*) 'NORM_VECTV (NORM) : ', NORM_VECTV;

!       call VECT_PRODUCT(VECTW_PRIME(1:3),VECTU_PRIME(1:3),VECTV_PRIME(1:3));
!       NORM_VECTV_PRIME = DSQRT( DOT_PRODUCT(VECTV_PRIME(1:3),VECTV_PRIME(1:3)) );
!       write(99,*) 'NORM_VECTV_PRIME (NORM) : ', NORM_VECTV_PRIME;

!       if ( DOT_PRODUCT(VECTU_PRIME(1:3),VECTV(1:3)) > 0.0d0 ) THETA = - THETA;

!       if ( THETA /= 0.0d0 ) then
!           SIN_THETA = DSIN( THETA );

!           PASSC(1:3,1) = VECTU(1:3); 
!           PASSC(1:3,2) = VECTV(1:3); 
!           PASSC(1:3,3) = VECTW(1:3); 

!           do h = 1, 2 ! ROTATE THE FIRST PARTICLE
!               do i = 1, NBFINAL(h)
!                   RION(1:3) = RFINAL(i,1:3,h) - COMRON(1,1:3); 
!                   RINO(1:3) = MATMUL(PASSB(1:3,1:3),RION(1:3));
!                   RINO(1:3) = RINO(1:3) - ANINT( RINO(1:3) );
!                   RION(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));

!                   RIJ(1) = DOT_PRODUCT(RION(1:3),VECTU(1:3));
!                   RIJ(2) = DOT_PRODUCT(RION(1:3),VECTV(1:3));
!                   RIJ(3) = DOT_PRODUCT(RION(1:3),VECTW(1:3));

!                   RIJ_PRIME(1) =   RIJ(1) * COS_THETA + RIJ(2) * SIN_THETA;
!                   RIJ_PRIME(2) = - RIJ(1) * SIN_THETA + RIJ(2) * COS_THETA;
!                   RIJ_PRIME(3) =   RIJ(3);

!                   PASSAGE IN THE GENERAL AXIS
!                   RION(1:3) = MATMUL(PASSC(1:3,1:3),RIJ_PRIME(1:3));

!                   RFINAL(i,1:3,h) = COMRON(1,1:3) + RION(1:3);

!                   RINO(1:3) = MATMUL(PASSB(1:3,1:3),RFINAL(i,1:3,h));
!                   RINO(1:3) = RINO(1:3) - ANINT( RINO(1:3) );
!                   RFINAL(i,1:3,h) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!               end do
!           end do
!       end if

!       if ( THETA /= 0.0d0 ) then
!           COMRON(2,1:3) = RCOM_IDEAL(1:3) ! UPDATE COM OF THE SECOND PARTICLE (AFTER ROTATION)
!           COMRNO(2,1:3) = MATMUL(PASSB(1:3,1:3),COMRON(2,1:3));
!       end if

!       DEFINE TRANSLATION STEP
!       DTRBIN = ABS( COMRNO(2,3) ) / REAL( nbins );

!       write(99,*) 'DTRBIN : ', DTRBIN, DTRBIN * CELL_AXIS(3);

!       TRBIN = TRANS_INIT / CELL_AXIS(3) !+  DTRBIN;

!       TRANSR(1:3) = COMRNO(2,1:3);
!       TRANSR(3) = TRANSR(3) !-  DTRBIN;

!       write(99,*) 'TRANS(3) : ', TRANSR(3), TRANSR(3) * CELL_AXIS(3);

!       TRANSLATE THE SECOND PARTICLE        
!       do h = 1, 2
!           do i = 1, NBFINAL(h)
!               RINO(1:3)  = MATMUL(PASSB(1:3,1:3),RDOPLE(i,1:3,h));
!               RINO(1:3)  = RINO(1:3) + TRANSR(1:3);
!               RINO(1:3)  = RINO(1:3) - ANINT( RINO(1:3) );
!               RDOPLE(i,1:3,h) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!           end do
!       end do

!       if ( THETA /= 0.0d0 ) then
!           do h = 1, 2 ! ROTATE ATOMS OF THE SECOND PARTICLE
!               do i = 1, NBFINAL(h)
!                   RION(1:3) = RDOPLE(i,1:3,h) - COMRON(1,1:3);
!                   RINO(1:3) = MATMUL(PASSB(1:3,1:3),RION(1:3));
!                   RINO(1:3) = RINO(1:3) - ANINT( RINO(1:3) );
!                   RION(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));

!                   RIJ(1) = DOT_PRODUCT(RION(1:3),VECTU(1:3));
!                   RIJ(2) = DOT_PRODUCT(RION(1:3),VECTV(1:3));
!                   RIJ(3) = DOT_PRODUCT(RION(1:3),VECTW(1:3));

!                   RIJ_PRIME(1) =   RIJ(1) * COS_THETA + RIJ(2) * SIN_THETA;
!                   RIJ_PRIME(2) = - RIJ(1) * SIN_THETA + RIJ(2) * COS_THETA;
!                   RIJ_PRIME(3) =   RIJ(3);

!                   PASSAGE IN THE GENERAL AXIS
!                   RION(1:3) = MATMUL(PASSC(1:3,1:3),RIJ_PRIME(1:3));

!                   RDOPLE(i,1:3,h) = COMRON(1,1:3) + RION(1:3);

!                   RINO(1:3) = MATMUL(PASSB(1:3,1:3),RDOPLE(i,1:3,h));
!                   RINO(1:3) = RINO(1:3) - ANINT( RINO(1:3) );
!                   RDOPLE(i,1:3,h) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!               end do
!           end do
!       end if

!       do j = 1, nbins
!           write(99,*) '*********************************************************************';
!           write(99,*) 'TRBIN : ', TRBIN, j;
!           write(99,'(a13,3f12.5)') 'COM 2 [A] : ', COMRON(2,1:3);

!           do h = 1, 2 ! TRANSLATION OF PARTICLE 2
!               do g = 1, NBFINAL(h)
!                   RINO(1:3)  = MATMUL(PASSB(1:3,1:3),RDOPLE(g,1:3,h));
!                   RINO(3)    = RINO(3) + DTRBIN;
!                   RINO(1:3)  = RINO(1:3) - ANINT( RINO(1:3) );
!                   RDOPLE(g,1:3,h) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
!               end do
!           end do

!           call CHECK_DISTANCES();

!           write(99,*);
!           write(CHBIN,'(i3.3)') j; 

!           open(104,file='001-'//trim(CHBIN)//'.xyz');

!           if ( iboxc == 1 ) then
!               write(104,*) NWORK+8;
!           else
!               write(104,*) NWORK;
!           end if

!           write(104,'(6f15.6)') CELL_AXIS(1:3), CELL_ANGDEG(1:3);

!           do h = 1, 2           ! WRITE COORDINATES OF H2O MOLECULES
!               do i = 1, NBFINAL(h)
!                   if ( NATFINAL(i,h) == 'Ow' .OR. NATFINAL(i,h) == 'Hw' ) then
!                       write(104,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(1) = NBREDOPLE(1) + 1;
!                   end if
!               end do
!           end do

!           do h = 1, 2           ! WRITE COORDINATES OF H2O MOLECULES
!               do i = 1, NBFINAL(h)
!                   if ( NATDOPLE(i,h) == 'Ow' .OR. NATDOPLE(i,h) == 'Hw' ) then
!                       write(104,'(a3,3f15.6)') NATDOPLE(i,h), RDOPLE(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(1) = NBREDOPLE(1) + 1;
!                   end if
!               end do
!           end do

!           do h = 1, 2           ! WRITE COORDINATES OF CW IONS
!               do i = 1, NBFINAL(h)
!                   if ( NATFINAL(i,h) == 'Cw' ) then
!                       write(104,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(2) = NBREDOPLE(2) + 1;
!                   end if
!               end do
!           end do

!           do h = 1, 2           ! WRITE COORDINATES OF CW IONS
!               do i = 1, NBFINAL(h)
!                   if ( NATDOPLE(i,h) == 'Cw' ) then
!                       write(104,'(a3,3f15.6)') NATDOPLE(i,h), RDOPLE(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(2) = NBREDOPLE(2) + 1;
!                   end if
!               end do
!           end do

!           WRITE THE CSH PARTICLE THAT IS ALLOWED TO MOVE  (MIDDLE OF THE BOX)
!           do h = 1, 2
!               do i = 1, NBFINAL(h)
!                   if ( NATDOPLE(i,h) == 'Ca' ) then
!                       write(104,'(a3,3f15.6)') NATDOPLE(i,h), RDOPLE(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(3) = NBREDOPLE(3) + 1;
!                   end if
!               end do
!           end do
!           do h = 1, 2
!               do i = 1, NBFINAL(h)
!                   if ( NATDOPLE(i,h) == 'Os' ) then
!                       write(104,'(a3,3f15.6)') NATDOPLE(i,h), RDOPLE(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(4) = NBREDOPLE(4) + 1;
!                   end if
!               end do
!           end do
!           do h = 1, 2
!               do i = 1, NBFINAL(h)
!                   if ( NATDOPLE(i,h) == 'Si' ) then
!                       write(104,'(a3,3f15.6)') NATDOPLE(i,h), RDOPLE(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(5) = NBREDOPLE(5) + 1;
!                   end if
!               end do
!           end do
!           do h = 1, 2
!               do i = 1, NBFINAL(h)
!                   if ( NATDOPLE(i,h) == 'O' ) then
!                       write(104,'(a3,3f15.6)') NATDOPLE(i,h), RDOPLE(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(6) = NBREDOPLE(6) + 1;
!                   end if
!               end do
!           end do

!           WRITE THE CSH PARTICLE THAT IS FIXED (TOP & DOWN OF THE BOX)
!           do h = 1, 2
!               do i = 1, NBFINAL(h)
!                   if ( NATFINAL(i,h) == 'Ca' ) then
!                       write(104,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(7) = NBREDOPLE(7) + 1;
!                   end if
!               end do
!           end do
!           do h = 1, 2
!               do i = 1, NBFINAL(h)
!                   if ( NATFINAL(i,h) == 'Os' ) then 
!                       write(104,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(8) = NBREDOPLE(8) + 1;
!                   end if
!               end do
!           end do
!           do h = 1, 2
!               do i = 1, NBFINAL(h)
!                   if ( NATFINAL(i,h) == 'Si' ) then
!                       write(104,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(9) = NBREDOPLE(9) + 1;
!                   end if
!               end do
!           end do
!           do h = 1, 2
!               do i = 1, NBFINAL(h)
!                   if ( NATFINAL(i,h) == 'O' ) then
!                       write(104,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
!                       if ( j == nbins ) NBREDOPLE(10) = NBREDOPLE(10) + 1;
!                   end if
!               end do
!           end do

!           if ( iboxc == 1 ) then
!               write(104,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/-0.5, -0.5, -0.5/));
!               write(104,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/-0.5, -0.5,  0.5/));
!               write(104,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/-0.5,  0.5, -0.5/));
!               write(104,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/-0.5,  0.5,  0.5/));
!               write(104,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/ 0.5, -0.5, -0.5/));
!               write(104,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/ 0.5, -0.5,  0.5/));
!               write(104,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/ 0.5,  0.5, -0.5/));
!               write(104,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/ 0.5,  0.5,  0.5/));
!           end if

!           close(104);

!           call WRITE_FIELD(1,CHBIN,COMRON(1:2,1:3)) ! WRITE FIELD FILE FOR MC SIMULATIONS

!           if ( j == nbins ) then
!               write(99,*) 'NH2O : ', NBREDOPLE(1);
!               write(99,*) 'Cw   : ', NBREDOPLE(2);
!               write(99,*) 'PART1 Ca : ', NBREDOPLE(3);
!               write(99,*) 'PART1 Os : ', NBREDOPLE(4);
!               write(99,*) 'PART1 Si : ', NBREDOPLE(5);
!               write(99,*) 'PART1 O  : ', NBREDOPLE(6);
!               write(99,*);
!               write(99,*) 'PART2 Ca : ', NBREDOPLE(7);
!               write(99,*) 'PART2 Os : ', NBREDOPLE(8);
!               write(99,*) 'PART2 Si : ', NBREDOPLE(9);
!               write(99,*) 'PART2 O  : ', NBREDOPLE(10);
!               write(99,*);

!               call system('cp 001-'//trim(CHBIN)//'.xyz   001.xyz'); 
!               call system('cp 001-'//trim(CHBIN)//'_FIELD.dat   001_FIELD.dat'); 
!           end if

!           UPDATE TRANSLATION
!           TRBIN = TRBIN + DTRBIN;
!           TRANSR(3) = TRANSR(3) + DTRBIN;

!           UPDATE CENTER OF MASS
!           COMRNO(2,1:3) = MATMUL(PASSB(1:3,1:3),COMRON(2,1:3));           
!           COMRNO(2,3) = COMRNO(2,3) + DTRBIN;
!           COMRON(2,1:3) = MATMUL(PASSA(1:3,1:3),COMRNO(2,1:3));
!       end do

!       deallocate(NATDOPLE,RDOPLE);

!   else
!       ### WRITE FIELD FILE FOR MC SIMULATIONS ####################################################
!   end if

!   if ( icheck_silica_chains == 2 ) deallocate(NATFINAL,RFINAL);

end subroutine MAKE_SUPCELLCONF
