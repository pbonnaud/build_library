subroutine BUILD_TEMPLATED_CONFIG(icanal) 

!   ************************************************************************************************
!   **                                 READ LAMMPS FILES                                          **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal        : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                   **
!   **                                                                                            **
!   ************************************************************************************************

    use module_physical_constants;

    use module_size_arrays;

    use module_config;

    use module_library;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: IMAX_ATOMID;

!   integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4;

    integer (kind=4) :: ILENGTH_TITLE, ILENGTH_CHEXT;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Build templated config';                                                  !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of variables for the library ################################################
                                                                                         ! 
    NATOM_LIBRARY     = 0;                                                               !
                                                                                         !
    NBOND_LIBRARY     = 0;                                                               !
                                                                                         !
    NANGLE_LIBRARY    = 0;                                                               !
                                                                                         !
    NDIHEDRAL_LIBRARY = 0;                                                               !
                                                                                         !
    NIMPROPER_LIBRARY = 0;                                                               !
                                                                                         !
!   ### Set the number of atoms and the number of atom types for the current molecule ##############
                                                                                         !
    NTYPE_ATOM_LIBRARY = NTYPE_ATOM;                                                     !
                                                                                         !
    IMAX_ATOMID = 0;                                                                     !
                                                                                         !
    do i = 1, NATOM;                                                                     !
                                                                                         !
        if ( CONFIG_MOLECULEID(i) == 1 ) NATOM_LIBRARY = NATOM_LIBRARY + 1;              !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a18,i8,a44)') '| NATOM_LIBRARY : ', &                                 !
                                  NATOM_LIBRARY,       &                                 !
                                  REPEAT(' ',43)//'|';                                   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays with the size of the total number of atoms to generate library files #######
                                                                                         !
    allocate(CONFIG_ATOMID_LIBRARY(1:NATOM_LIBRARY));                                    !
                                                                                         !
    allocate(CONFIG_ATOM_TYPE_LIBRARY(1:NATOM_LIBRARY));                                 !
                                                                                         !
    allocate(CONFIG_QI_LIBRARY(1:NATOM_LIBRARY));                                        !
                                                                                         !
    allocate(CONFIG_RI_LIBRARY(1:3,1:NATOM_LIBRARY));                                    !
                                                                                         !
    allocate(CONFIG_NAT_LIBRARY(1:NATOM_LIBRARY));                                       !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of arrays with the size of the total number of atoms ########################
                                                                                         !
    CONFIG_ATOMID_LIBRARY(1:NATOM_LIBRARY) = 0;                                          !
                                                                                         !
    CONFIG_ATOM_TYPE_LIBRARY(1:NATOM_LIBRARY) = 0;                                       !
                                                                                         !
    CONFIG_QI_LIBRARY(1:NATOM_LIBRARY) = 0.0d0;                                          ! 
                                                                                         !
    CONFIG_RI_LIBRARY(1:3,1:NATOM_LIBRARY) = 0.0d0;                                      !
                                                                                         !
    CONFIG_NAT_LIBRARY(1:NATOM_LIBRARY) = 'XXX';                                         !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Setting atom properties ####################################################################
                                                                                         !
    NATOM_LIBRARY = 0;                                                                   !

    do i = 1, NATOM;

        if ( CONFIG_MOLECULEID(i) == 1 ) then;

            NATOM_LIBRARY = NATOM_LIBRARY + 1;

            CONFIG_NAT_LIBRARY(NATOM_LIBRARY) = CONFIG_NAT(i)

            CONFIG_RI_LIBRARY(1:3,NATOM_LIBRARY) = CONFIG_RI(1:3,i);

            CONFIG_QI_LIBRARY(NATOM_LIBRARY) = CONFIG_QI(i);

            CONFIG_ATOMID_LIBRARY(NATOM_LIBRARY) = CONFIG_ATOMID(i);
        
            CONFIG_ATOM_TYPE_LIBRARY(NATOM_LIBRARY) = CONFIG_ATOM_TYPE(i);

            if ( CONFIG_ATOMID_LIBRARY(NATOM_LIBRARY) > IMAX_ATOMID ) then;

                IMAX_ATOMID = CONFIG_ATOMID_LIBRARY(NATOM_LIBRARY);

            end if

        end if

    end do
                                                                                         !
    write(icanal,'(a16,i8,a46)') '| IMAX_ATOMID : ', IMAX_ATOMID, REPEAT(' ',45)//'|';   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the number of bonds in the molecule ####################################################

    do i = 1, NBOND;

        if ( BOND_ATOMID(1,i) > IMAX_ATOMID ) CYCLE;

        if ( BOND_ATOMID(2,i) > IMAX_ATOMID ) CYCLE;

        NBOND_LIBRARY = NBOND_LIBRARY + 1;

    end do

    write(icanal,'(a18,i8,a44)') '| NBOND_LIBRARY : ', &
                                 NBOND_LIBRARY,        &
                                 REPEAT(' ',43)//'|';

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays containing bond properties of the library molecule to build ################
                                                                                         !
    allocate(BOND_TYPE_LIBRARY(1:NBOND_LIBRARY));

    allocate(BOND_ATOMID_LIBRARY(1:2,1:NBOND_LIBRARY));

!   stop;

!   ### Initialization of arrays containing bond properties of the library molecule to build #######

    BOND_TYPE_LIBRARY(1:NBOND_LIBRARY) = 0;

    BOND_ATOMID_LIBRARY(1:2,1:NBOND_LIBRARY) = 0;

!   stop;

!   ### Set bond properties of the molecule to build ###############################################

    NBOND_LIBRARY = 0;

    do i = 1, NBOND;

        if ( BOND_ATOMID(1,i) > IMAX_ATOMID ) CYCLE;

        if ( BOND_ATOMID(2,i) > IMAX_ATOMID ) CYCLE;
  
        NBOND_LIBRARY = NBOND_LIBRARY + 1;

        BOND_TYPE_LIBRARY(NBOND_LIBRARY) = BOND_TYPE(i);

        BOND_ATOMID_LIBRARY(1:2,NBOND_LIBRARY) = BOND_ATOMID(1:2,i);

    end do

!   stop;

!   ### Set the number of angles in the molecule ###################################################

    do i = 1, NANGLE;

        if ( ANGLE_ATOMID(1,i) > IMAX_ATOMID ) CYCLE;

        if ( ANGLE_ATOMID(2,i) > IMAX_ATOMID ) CYCLE;

        if ( ANGLE_ATOMID(3,i) > IMAX_ATOMID ) CYCLE;

        NANGLE_LIBRARY = NANGLE_LIBRARY + 1;

    end do

    write(icanal,'(a19,i8,a43)') '| NANGLE_LIBRARY : ', &
                                 NANGLE_LIBRARY,        & 
                                 REPEAT(' ',42)//'|';

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!   stop;

!   ### Allocate arrays containing angle properties ################################################

    allocate(ANGLE_TYPE_LIBRARY(1:NANGLE_LIBRARY));

    allocate(ANGLE_ATOMID_LIBRARY(1:3,1:NANGLE_LIBRARY));

!   stop;

!   ### Initialization of arrays containing angle properties #######################################

    ANGLE_TYPE_LIBRARY(1:NANGLE_LIBRARY) = 0;

    ANGLE_ATOMID_LIBRARY(1:3,1:NANGLE_LIBRARY) = 0;

!   stop;

!   ### Set angle properties of the library molecule to build ######################################

    NANGLE_LIBRARY = 0;

    do i = 1, NANGLE;

        if ( ANGLE_ATOMID(1,i) > IMAX_ATOMID ) CYCLE;

        if ( ANGLE_ATOMID(2,i) > IMAX_ATOMID ) CYCLE;

        if ( ANGLE_ATOMID(3,i) > IMAX_ATOMID ) CYCLE;

        NANGLE_LIBRARY = NANGLE_LIBRARY + 1;

        ANGLE_TYPE_LIBRARY(NANGLE_LIBRARY) = ANGLE_TYPE(i);

        ANGLE_ATOMID_LIBRARY(1:3,NANGLE_LIBRARY) = ANGLE_ATOMID(1:3,i);

    end do

!   stop;

!   ### Set the number of dihedrals in the library molecule to build ###############################
                                                                                         !
    do i = 1, NDIHEDRAL;

        if ( DIHEDRAL_ATOMID(1,i) > IMAX_ATOMID ) CYCLE;

        if ( DIHEDRAL_ATOMID(2,i) > IMAX_ATOMID ) CYCLE;

        if ( DIHEDRAL_ATOMID(3,i) > IMAX_ATOMID ) CYCLE;

        if ( DIHEDRAL_ATOMID(4,i) > IMAX_ATOMID ) CYCLE;

        NDIHEDRAL_LIBRARY = NDIHEDRAL_LIBRARY + 1;

    end do

    write(icanal,'(a22,i8,a40)') '| NDIHEDRAL_LIBRARY : ', &
                                 NDIHEDRAL_LIBRARY,        &
                                 REPEAT(' ',39)//'|'; 

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
 
!   stop;

!   ### Allocate arrays for dihedral properties of the molecule to build ###########################

    allocate(DIHEDRAL_TYPE_LIBRARY(1:NDIHEDRAL_LIBRARY));

    allocate(DIHEDRAL_ATOMID_LIBRARY(1:4,1:NDIHEDRAL_LIBRARY));

!   stop;

!   ### Initialization of arrays for dihedral properties of the molecule to build ##################
 
    DIHEDRAL_TYPE_LIBRARY(1:NDIHEDRAL_LIBRARY) = 0;

    DIHEDRAL_ATOMID_LIBRARY(1:4,1:NDIHEDRAL_LIBRARY) = 0;

!   stop;

!   ### Set dihedral properties of the library molecule to build ###################################

    NDIHEDRAL_LIBRARY = 0;

    do i = 1, NDIHEDRAL;

        if ( DIHEDRAL_ATOMID(1,i) > IMAX_ATOMID ) CYCLE;

        if ( DIHEDRAL_ATOMID(2,i) > IMAX_ATOMID ) CYCLE;

        if ( DIHEDRAL_ATOMID(3,i) > IMAX_ATOMID ) CYCLE;

        if ( DIHEDRAL_ATOMID(4,i) > IMAX_ATOMID ) CYCLE;

        NDIHEDRAL_LIBRARY = NDIHEDRAL_LIBRARY + 1;

        DIHEDRAL_TYPE_LIBRARY(NDIHEDRAL_LIBRARY) = DIHEDRAL_TYPE(i);

        DIHEDRAL_ATOMID_LIBRARY(1:4,NDIHEDRAL_LIBRARY) = DIHEDRAL_ATOMID(1:4,i);

    end do

!   stop;

    NIMPROPER_LIBRARY = 0;

    if ( NIMPROPER > 0 ) then;

        do i = 1, NIMPROPER;

            if ( IMPROPER_ATOMID(1,i) > IMAX_ATOMID ) CYCLE;

            if ( IMPROPER_ATOMID(2,i) > IMAX_ATOMID ) CYCLE;

            if ( IMPROPER_ATOMID(3,i) > IMAX_ATOMID ) CYCLE;

            if ( IMPROPER_ATOMID(4,i) > IMAX_ATOMID ) CYCLE;

            NIMPROPER_LIBRARY = NIMPROPER_LIBRARY + 1;

        end do

        write(icanal,'(a22,i8,a40)') '| NIMPROPER_LIBRARY : ', &
                                     NIMPROPER_LIBRARY,        &
                                     REPEAT(' ',39)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        allocate(IMPROPER_TYPE_LIBRARY(1:NIMPROPER_LIBRARY));

        allocate(IMPROPER_ATOMID_LIBRARY(1:4,1:NIMPROPER_LIBRARY));

        IMPROPER_TYPE_LIBRARY(1:NIMPROPER_LIBRARY) = 0;

        IMPROPER_ATOMID_LIBRARY(1:4,1:NIMPROPER_LIBRARY) = 0;

        NIMPROPER_LIBRARY = 0;

        NTYPE_IMPROPER_LIBRARY = 0;

        do i = 1, NIMPROPER;

            if ( IMPROPER_ATOMID(1,i) > IMAX_ATOMID ) CYCLE;

            if ( IMPROPER_ATOMID(2,i) > IMAX_ATOMID ) CYCLE;

            if ( IMPROPER_ATOMID(3,i) > IMAX_ATOMID ) CYCLE;

            if ( IMPROPER_ATOMID(4,i) > IMAX_ATOMID ) CYCLE;

            NIMPROPER_LIBRARY = NIMPROPER_LIBRARY + 1;

            IMPROPER_TYPE_LIBRARY(NIMPROPER_LIBRARY) = IMPROPER_TYPE(i);

            IMPROPER_ATOMID_LIBRARY(1:4,NIMPROPER_LIBRARY) = IMPROPER_ATOMID(1:4,i);

            if ( IMPROPER_TYPE_LIBRARY(NIMPROPER_LIBRARY) > NTYPE_IMPROPER_LIBRARY ) &
                    NTYPE_IMPROPER_LIBRARY = IMPROPER_TYPE_LIBRARY(NIMPROPER_LIBRARY);
        end do

        write(icanal,'(a27,i8,a35)') '| NTYPE_IMPROPER_LIBRARY : ', &
                                     NTYPE_IMPROPER_LIBRARY,        &
                                     REPEAT(' ',34)//'|';

    end if

!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine BUILD_TEMPLATED_CONFIG
