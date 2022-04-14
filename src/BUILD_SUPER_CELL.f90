subroutine BUILD_SUPER_CELL(icanal,PATH_DIRECTORY,CH_WORKING_FILE)

!   ************************************************************************************************
!   **                     APPLY REPLICA CREATION OF THE ORIGINAL UNIT CELL                       **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=200), intent(in) :: PATH_DIRECTORY;

    character (len=5), intent(in) :: CH_WORKING_FILE;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m;

    integer (kind=4) :: iespece, insp, islab, IATOM;

    real (kind=8), dimension(1:3) :: TRANSR, RINO;

    integer (kind=4) :: INIT_NATOM, NATOM;

    real (kind=8), dimension(1:3,1:3) :: PASSA, PASSB;

    real (kind=8), dimension(1:3) :: INIT_CELL_AXIS, INIT_CELL_ANGDEG, CELL_AXIS, CELL_ANGDEG;

    real (kind=8), dimension(1:3,1:1000) :: INIT_RI;

    character (len=3), dimension(1:1000) :: INIT_NAT;

    real (kind=8), dimension(1:3,1:10000) :: CONFIG_RI;

    character (len=3), dimension(1:10000) :: CONFIG_NAT;

!   ************************************************************************************************

    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '| BUILD SUPER CELL '//REPEAT(' ',38)//'|';
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    open(10,file=TRIM(PATH_DIRECTORY)//TRIM(CH_WORKING_FILE)//'_Primitive_Cell.xyz');
    read(10,*) INIT_NATOM;
    read(10,*) INIT_CELL_AXIS(1:3), INIT_CELL_ANGDEG(1:3);
    do i = 1, INIT_NATOM;
        read(10,*) INIT_NAT(i), INIT_RI(1:3,i);
    end do
    close(10);

    write(icanal,*) INIT_CELL_AXIS(1:3);
    write(icanal,*) INIT_CELL_ANGDEG(1:3);

    CELL_AXIS(1:3)   = INIT_CELL_AXIS(1:3) * MULTI(1:3);
    CELL_ANGDEG(1:3) = INIT_CELL_ANGDEG(1:3);

    call SET_PASSAGE(99,CELL_AXIS(1:3),CELL_ANGDEG(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));

    NATOM = INIT_NATOM * MULTI(1) * MULTI(2) * MULTI(3);

    write(icanal,*) '| NATOM : ', NATOM;

    IATOM = 0;

    TRANSR(1) = - ( REAL(MULTI(1)) - 1.0d0 ) / ( 2.0d0 * REAL(MULTI(1)) );

    do k = 1, MULTI(1)
        TRANSR(2) = - ( REAL(MULTI(2)) - 1.0d0 ) / ( 2.0d0 * REAL(MULTI(2)) );
        do j = 1, MULTI(2)
           TRANSR(3) = - ( REAL(MULTI(3)) - 1.0d0 ) / ( 2.0d0 * REAL(MULTI(3)) );
            do i = 1, MULTI(3)

                do m = 1, INIT_NATOM;
                    IATOM = IATOM + 1;
                    CONFIG_NAT(IATOM)    = INIT_NAT(m); 
                    CONFIG_RI(1:3,IATOM) = INIT_RI(1:3,m);

                    RINO(1:3) = MATMUL(PASSB(1:3,1:3),CONFIG_RI(1:3,IATOM));
                    RINO(1:3) = RINO(1:3) + TRANSR(1:3);
                    CONFIG_RI(1:3,IATOM) = MATMUL(PASSA(1:3,1:3),RINO(1:3));
                end do

                if ( MULTI(3) /= 1 ) TRANSR(3) = TRANSR(3) + 1.0d0 / REAL( MULTI(3) );
            end do
            TRANSR(2) = TRANSR(2) + 1.0d0 / REAL( MULTI(2) );
        end do
        TRANSR(1) = TRANSR(1) + 1.0d0 / REAL( MULTI(1) );
    end do

    write(icanal,*) 'XXX 1', IATOM;

    open(103,file=TRIM(PATH_DIRECTORY)//TRIM(CH_WORKING_FILE)//'_Super_Cell.xyz');
    write(103,*) NATOM;
    write(103,'(6f15.6)') CELL_AXIS(1:3), CELL_ANGDEG(1:3);
    do i = 1, NATOM;
        write(103,'(a3,3f15.6)') CONFIG_NAT(i), CONFIG_RI(1:3,i);
    end do
    close(103);

    write(icanal,*);
    write(icanal,*) 'END REPLICA CREATION';
    write(icanal,*);

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine BUILD_SUPER_CELL
