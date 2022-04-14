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

subroutine WRITE_FINAL_CONFIG(CELL_AXIS,CELL_ANGDEG,PASSA,PASSB)!PATH_DIRECTORY,CH_WORKING_FILE);

!   ************************************************************************************************
!   **                           WRITE THE MOLECULAR CONFIGURATION                                **
!   ************************************************************************************************

    use module_data_in;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3), intent(in) :: CELL_AXIS, CELL_ANGDEG;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   character (len=200), intent(in) :: PATH_DIRECTORY;

!   character (len=5), intent(in) :: CH_WORKING_FILE;

!   ************************************************************************************************

    integer (kind=4) :: i;

    integer (kind=4) :: iespece, insp, jespece, jnsp, itmp;

    real (kind=8) :: RIJ;

    real (kind=8), dimension(1:3) :: DRIJ, DRIJNO;

    character (len=200) :: CH_DATA_FILE_NAME;

    logical :: PROBE1;

!   ************************************************************************************************

    write(99,*) 'WRITE FINAL CONFIGURATION : ';
    write(99,*);

!   CH_DATA_FILE_NAME = trim(PATH_DIRECTORY)//trim(CH_WORKING_FILE)//'.xyz';

!   write(99,*) CH_DATA_FILE_NAME;

!   inquire(FILE=trim(CH_DATA_FILE_NAME),EXIST=PROBE1);

!   if (PROBE1 == .TRUE. ) then
!       write(99,*) trim(CH_DATA_FILE_NAME)//' ALREADY EXIST !!!';
!       write(99,*) 'END OF PROGRAM';
!       close(99);
!       stop;
!   end if

!   deallocate(DATA_FINAL,SUB_FINAL);

!   allocate(DATA_FINAL(1:espece));

!   NATSP_FINAL(1:espece) = NATSP_TMP(1:espece);
    NSLAB_FINAL = NSLAB_TMP;

!   do i = 1, espece
!       allocate(DATA_FINAL(i)%NAT(1:NATSP_TMP(i)));
!       allocate(DATA_FINAL(i)%RI(1:3,1:NATSP_TMP(i)));
!   end do

!   allocate(SUB_FINAL(1:NSLAB_FINAL));

    NATSP_FINAL(1) = 0;

!   do iespece = 1, espece
    do insp = 1, NATSP_TMP(1)
!       if ( DATA_TMP(1)%NAT(insp) /= 'Ow' ) CYCLE;

        NATSP_FINAL(1) = NATSP_FINAL(1) + 1;

!       DATA_FINAL(1)%NAT(NATSP_FINAL(1))    = DATA_TMP(1)%NAT(insp);
!       DATA_FINAL(1)%RI(1:3,NATSP_FINAL(1)) = DATA_TMP(1)%RI(1:3,insp);

!       DATA_TMP(1)%NAT(insp) = 'XXX';

        itmp = NATSP_FINAL(1);

        do jnsp = 1, NATSP_TMP(1)
!           if ( DATA_TMP(1)%NAT(jnsp) /= 'Hw' ) CYCLE;
!           DRIJ(1:3)   = DATA_TMP(1)%RI(1:3,jnsp) - DATA_FINAL(1)%RI(1:3,itmp);
!           DRIJNO(1:3) = MATMUL(PASSB(1:3,1:3),DRIJ(1:3));
!           DRIJNO(1:3) = DRIJNO(1:3) - ANINT( DRIJNO(1:3) );
!           DRIJ(1:3)   = MATMUL(PASSA(1:3,1:3),DRIJNO(1:3));
!           RIJ = DSQRT(DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)));

            if ( RIJ < 1.1d0 ) then
                NATSP_FINAL(1) = NATSP_FINAL(1) + 1;

!               DATA_FINAL(1)%NAT(NATSP_FINAL(1))    = DATA_TMP(1)%NAT(jnsp);
!               DATA_FINAL(1)%RI(1:3,NATSP_FINAL(1)) = DATA_TMP(1)%RI(1:3,jnsp);

!               DATA_TMP(1)%NAT(jnsp) = 'XXX';
            end if

        end do
    end do
!   end do

    do i = 2, espece
        NATSP_FINAL(i) = NATSP_TMP(i);
!       DATA_FINAL(i)%NAT(1:NATSP_FINAL(i))    = DATA_TMP(i)%NAT(1:NATSP_TMP(i));
!       DATA_FINAL(i)%RI(1:3,1:NATSP_FINAL(i)) = DATA_TMP(i)%RI(1:3,1:NATSP_TMP(i));
    end do

    do i = 1, NSLAB_FINAL
!       SUB_FINAL(i)%NAT     = SUB_TMP(i)%NAT;
!       SUB_FINAL(i)%RI(1:3) = SUB_TMP(i)%RI(1:3);
    end do

    NATOM_FINAL = NSLAB_FINAL;

    do i = 1, espece
        NATOM_FINAL = NATOM_FINAL + NATSP_FINAL(i);
    end do

!   open(110,file=trim(CH_DATA_FILE_NAME),status='new');
    open(110,file='001.xyz');
    write(110,*) NATOM_FINAL;
    write(110,'(6f20.13)') CELL_AXIS(1:3), CELL_ANGDEG(1:3);

!   do iespece = 1, espece
!       do insp = 1, NATSP_FINAL(iespece)
!           write(110,'(a3,3f17.8)') DATA_FINAL(iespece)%NAT(insp), DATA_FINAL(iespece)%RI(1:3,insp);
!       end do
!   end do

!   if ( chslab == 1 ) then
!       do i = 1, NSLAB_FINAL    ! READ SUBSTRATE PARAMETERS
!           write(110,'(a3,3f17.8)') SUB_FINAL(i)%NAT, SUB_FINAL(i)%RI(1:3);
!       end do
!   end if

    close(110);

end subroutine WRITE_FINAL_CONFIG
