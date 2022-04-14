subroutine BUILD_PRIMITIVE_CELL(icanal, &
                                PATH_DIRECTORY, &
                                CH_WORKING_FILE, &
                                PRIMITIVE_CELL_AXIS, &
                                PRIMITIVE_CELL_ANGDEG, &
                                PASSA, &
                                PASSB, &
                                PRIMITIVE_NATOM, &
                                PRIMITIVE_NAT, &
                                PRIMITIVE_RI, &
                                NBRE_REPEAT_PATTERN, &
                                REPEAT_PATTERN_SI)

    use module_data_in;

!   ***********************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: PRIMITIVE_NATOM, NBRE_REPEAT_PATTERN;

    real (kind=8), dimension(1:3,1:20), intent(in) :: PRIMITIVE_RI, REPEAT_PATTERN_SI;

    character (len=3), dimension(1:20), intent(in) :: PRIMITIVE_NAT;

    character (len=200), intent(in) :: PATH_DIRECTORY;

    character (len=5), intent(in) :: CH_WORKING_FILE;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

    real (kind=8), dimension(1:3), intent(in) :: PRIMITIVE_CELL_AXIS, PRIMITIVE_CELL_ANGDEG;


!   ***********************************************************************************************

    integer (kind=4) :: i, j;

    real (kind=8), dimension(1:3) :: RINO, RI;

!   ***********************************************************************************************

    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '| BUILD THE PRIMITIVE CELL '//REPEAT(' ',42)//'|';
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    open(102,file=TRIM(PATH_DIRECTORY)//TRIM(CH_WORKING_FILE)//'_Primitive_Cell.xyz');
    write(102,*) NBRE_REPEAT_PATTERN * PRIMITIVE_NATOM;
    write(102,'(6f15.6)') PRIMITIVE_CELL_AXIS(1:3), PRIMITIVE_CELL_ANGDEG(1:3);
    do i = 1, NBRE_REPEAT_PATTERN    
        do j = 1, PRIMITIVE_NATOM
            RINO(1:3) = PRIMITIVE_RI(1:3,j) + REPEAT_PATTERN_SI(1:3,i);
            RI(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));

            call APPLY_PBC(RI(1:3),RI(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));

            write(102,'(a3,3f15.6)') PRIMITIVE_NAT(j), RI(1:3);
        end do
    end do
    close(102);

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine BUILD_PRIMITIVE_CELL
