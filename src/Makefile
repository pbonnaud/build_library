# Master makefile for BUILDLIBRARY                                      
# Author: P. Bonnaud SEPTEMBER 2017                                     
#                                                                       
#=======================================================================
# Define default settings                                               
#=======================================================================
                                                                        
.SUFFIXES: .F90                                                         
.SUFFIXES: .f90                                                         
                                                                        
#=======================================================================
# Define object files                                                   
                                                                        
OBJS = APPLY_BUILD_CELLBOX.o APPLY_ION_INSERTION.o APPLY_MODIFICATION.o \
APPLY_PBC.o APPLY_PLATELET.o APPLY_REBUILD_SILICA_CHAINS.o \
APPLY_REBUILD_WATER.o APPLY_REBUILD_WATER_Hw.o APPLY_REPLICA_CREATION.o \
APPLY_TRANSLATION.o BUILDLIBRARY.o BUILD_PARAMETER_AMBER_TO_LAMMPS.o \
BUILD_PRIMITIVE_CELL.o BUILD_SUPER_CELL.o BUILD_TEMPLATED_CONFIG.o \
CHECK_DISTANCES.o CHECK_ELECTRONEUTRALITY.o CHECK_FILE_LIBRARY.o \
CHECK_OVERLAPPING_ATOMS.o CHECK_SILICACHAINS.o MAKE_PROGRAM_TITLE.o \
MAKE_ROUTINE_TITLE.o MAKE_SUPCELLCONF.o MANA_ERROR.o \
MOVE_TO_LIBRARY_DIRECTORY.o RANDOM_INSERTION_MOLECULE.o READ_AMBER_INPCRD.o \
READ_AMBER_PRMTOP.o READ_CONFIG.o READ_DATA.o \
READ_FIELD.o READ_FILE_LIBRARY.o READ_INPUT.o \
READ_INPUT_LAMMPS.o READ_INTERATOMIC_POTENTIALS_TEMPLATE.o READ_LAMMPS_CONFIG.o \
READ_LAMMPS_TEMPLATE.o SET_ATOM_PROPERTIES.o SET_CENTROSYM.o \
SET_CONFIG.o SET_PASSAGE.o SET_PBC.o \
SET_ROTATION_MATRIX.o SET_TRANSLATION.o TOOLS.o \
WRITE_CONFIG.o WRITE_DATA.o WRITE_ERROR_MESSAGE.o \
WRITE_FIELD.o WRITE_FINAL_CONFIG.o WRITE_INTERATOMIC_POTENTIALS_TEMPLATE.o \
WRITE_LAMMPS_CONFIG.o WRITE_LAMMPS_INFO.o WRITE_LAMMPS_TEMPLATE.o \
WRITE_XYZ_CONFIG.o module_config.o module_data_in.o \
module_frc_arrays.o module_library.o module_osef.o \
module_physical_constants.o module_size_arrays.o
                                                                      
#=====================================================================
                                                                      
MODULE = module_config.mod module_data_in.mod module_frc_arrays.mod\
module_library.mod module_osef.mod module_physical_constants.mod\
module_size_arrays.mod
                                                                      
#=====================================================================
# Define targets                                                      
                                                                      
all: BUILDLIBRARY_AMD.exe                                             
                                                                      
BUILDLIBRARY_AMD.exe: $(MODULE) $(OBJS)                               
	ifort -O3 -o BUILDLIBRARY_AMD.exe  $(OBJS)                          
                                                                      
module_config.mod: $(MODULE) module_config.F90
	ifort -c module_config.F90
                                                                      
module_data_in.mod: $(MODULE) module_data_in.F90
	ifort -c module_data_in.F90
                                                                      
module_frc_arrays.mod: $(MODULE) module_frc_arrays.F90
	ifort -c module_frc_arrays.F90
                                                                      
module_library.mod: $(MODULE) module_library.F90
	ifort -c module_library.F90
                                                                      
module_osef.mod: $(MODULE) module_osef.F90
	ifort -c module_osef.F90
                                                                      
module_physical_constants.mod: $(MODULE) module_physical_constants.F90
	ifort -c module_physical_constants.F90
                                                                      
module_size_arrays.mod: $(MODULE) module_size_arrays.F90
	ifort -c module_size_arrays.F90
                                                                      
                                                                      
.f90.o:                                                               
	ifort -c $<                                                       
                                                                      
.F90.o:                                                               
	ifort -c $<                                                         
                                                                      
clean:                                                                
	rm -f *.o *.mod                                                     
