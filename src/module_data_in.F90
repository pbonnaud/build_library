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

module module_data_in

    implicit none;

!   ************************************************************************************************
                                                                                         !
!   real (kind=8), parameter :: PI    = 3.14159265358979323846d0;                        !
                                                                                         !
!   real (kind=8), parameter :: TWOPI = 2.0d0 * PI;                                      !
                                                                                         !
!   real (kind=8), parameter :: Na    = 6.02214179E23;                                   ! IN [at/mol]
                                                                                         !
!   ************************************************************************************************
                                                                                         !
!   ### Set arrays for properties of atoms in the Mendeleiev table #################################
                                                                                         !
    real (kind=8), dimension(1:200) :: MOLAR_MASS_ELEMENT, &                             !
                                       ELEMENT_RADIUS;                                   !
                                                                                         !
    character (len=150), dimension(1:200) :: ELEMENT_FULL_NAME,   &                      !
                                             ELEMENT_SYMBOL_NAME;                        !
                                                                                         !
!   ### INPUT FILE PARAMETERS ######################################################################

    integer (kind=4) :: iadjust_charges;

    integer (kind=4) :: imodif, ipsub, ifsub, idople, iboxc, IROTA, ITRAN;

    integer (kind=4) :: ch_amorph;

    real (kind=8) ::  lmax;

    real (kind=8) :: TOTAL_NET_CHARGE;

    integer (kind=4) :: nbins;

    integer (kind=4) :: IKEEPGATH, ICENTROSYM, IOXYTETRA, IPLATELET;

    integer (kind=4) :: IUMBRELLA, NLUMBRELLA;

    real (kind=8) :: DELTAZ_KEEPGATH;

    character (len=3) :: CHFSUB;

    integer (kind=4) :: NLIST_SPECIES_MODIF;

    integer (kind=4), dimension(1:10) :: LIST_SPECIES_MODIF;

    real (kind=4) :: DRI_TRANS;

!   ### READ DATA PARAMETERS #######################################################################

    integer (kind=4) :: chslab;

    integer (kind=4) :: espece, NSLAB, NSLAB_FINAL, NSLAB_PLATE;

    integer (kind=4), dimension(1:10) :: nsite, N, NATSP, FIXSITE, NSITE_TOT;

    integer (kind=4), dimension(1:10) :: NATSP_FINAL, NATSP_PLATE;

    integer (kind=4) :: NLABEL;

    real (kind=8), dimension(1:10) :: ATOM_SIG, ATOM_EPS, ATOM_CHAR, ATOM_MASS;

!   ### READ FIELD PARAMETERS ######################################################################

    integer (kind=4) :: ICOMA;

    real (kind=8), dimension(1:3,1:10) :: COMA_COORD, COMA_COORD_FINAL;

!   ### READ CONFIGURATION PARAMETERS ##############################################################

    integer (kind=4) :: NTOT_ATOMS, NTOT_ATOMS_FINAL;

    integer (kind=4) :: NATOM_FINAL, NATOM_TMP, NSLAB_TMP;

    integer (kind=4), dimension(1:10) :: NATSP_TMP;

!   ### STARTING CELL PARAMETERS ###################################################################

    integer (kind=4) :: NMAX;

    integer (kind=4) :: NVDWIWS;

    integer (kind=4) :: NMOTIF, NMAILLE, NBCAT;

    integer (kind=4) :: NMULTI, NWORK, NWORK_EXTD;

    integer (kind=4) :: NBNATIN;

    integer (kind=4), dimension(1:2) :: NBINIT, NBFINAL, NBDOPLE; 

    integer (kind=4), dimension(1:20) :: NBATLIST;

    real (kind=8), dimension(1:20) :: QLIST;

    character (len=3), dimension(1:20) :: NATLIST; 

    integer (kind=4), dimension(1:3) :: MULTI;

    integer (kind=4), dimension(1:100) :: NMOL;

    integer (kind=4) :: nisp, nsusp;

    character (len=3), allocatable, dimension(:,:) :: NATINIT, NATFINAL, NATDOPLE;

    integer (kind=4), allocatable, dimension(:,:) :: ADSTATUS;

    real (kind=8), allocatable, dimension(:,:,:) :: RINIT, RFINAL, RDOPLE;

    real (kind=8), dimension(1:3) :: DOPLE_ANGDEG, DOPLE_ANGRAD, DOPLE_TRANS;

!   ### FINAL CELL PARAMETERS ######################################################################

    integer (kind=4) :: NESPECE, NVDWI, NSTOT_WALL;

    integer (kind=4) :: ILJONES, IPNTRAZ;

    integer (kind=4) :: FLAG_BONDS, FLAG_ANGLES;

    real (kind=8), dimension(1:3,1:10) :: RCDM_MOL, ORIENT_MOL;

    integer (kind=4), dimension(1:10) :: ICDM_MOL;

    real (kind=8), dimension(1:20,1:20) :: MASSMO, BONDS_KR, ANGLES_KR, BONDS_R0, ANGLES_THETA0;

    integer (kind=4), dimension(1:20,1:20) :: MULAT;

    real (kind=8), dimension(1:20) :: MTOT;

    real (kind=8), dimension(1:10,1:10) :: SIGATO, EPSATO, CHARATO;

    real (kind=8), dimension(1:10,1:10) :: AIJATO, BIJATO, C06ATO, C08ATO, C10ATO;

!   ### Parametres de la maille ####################################################################

    real (kind=8) :: VOL;

    real (kind=8), dimension(1:3) :: HPORE;

end module module_data_in
