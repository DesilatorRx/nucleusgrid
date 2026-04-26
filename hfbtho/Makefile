
SHELL := /bin/bash
COMPILE_DIR = $(PWD)

# Defines version number of BLAS/LAPACK library package
#  MKL (Quartz) .........: mkl-2017.1
#  MKL (Cab) ............: mkl-11.3.2
#  ESSL (Vulcan) ........: 5.1
VERSION_LIBRARY = mkl-2017.1

COMPILER    = GFORTRAN
FORTRAN_MPI = gfortran

# Compilation options
#   - DEBUG .............: Activates debugging options for each compiler
#   - PEDANTIC ..........: Optional pedantic flag, set to -pedantic to use
DEBUG          = FALSE
PEDANTIC       =
# Preprocessor options
#   - SWITCH_ESSL .......: Use IBM ESSL library for LAPACK and BLAS (different routine calls)
#   - USE_OPENMP ........: Use OpenMP multithreading
SWITCH_ESSL    = 0
USE_OPENMP     = 1

# Physics options
#    - USE_MPI ..........: 1 MPI parallelism for DRIP_lINES, DO_MASSTABLE or DO_PES modes
#                          0 Inactive
#    - DRIP_LINES .......: 1 Calculates a mass table from dripline to dripline
#                          0 Inactive
#    - DO_MASSTABLE .....: 1 Calculates a section of the mass table
#                          0 Inactive
#    - READ_FUNCTIONAL ..: 1 Reads the parameters of the functional from a file
#                          0 The energy functional is defined in the code based on the value of the
#                            input keyword 'functional'
#    - GOGNY_SYMMETRIES..: 1 For production runs; Assumes several symmetries in the finite range matrix
#                            elements.
#                          0 For debugging the finite range matrix elements; Makes no assumption of any
#                            symmetry in the finite range matrix elements.
#                            Requires compiler compatible with 2008 standard of up to 15 dimensional arrays.
#    - GOGNY_HYPER.......: 1 Uses hypergeometric function to calculate finite range matrix elements
#                            Accurate for big and small basis size
#                          0 Uses the direct Gogny transformation to calculate the finite range matrix
#                            elements; Accurate for small basis size ONLY (N<20)
#    - USE_QRPA .........: 1 Records HFB solutions in a format readable by the pnFAM-QRPA code of Mustonen & Shafer
#                          0 Inactive
USE_MPI          = 0
DRIP_LINES       = 0
DO_MASSTABLE     = 0
READ_FUNCTIONAL  = 0
GOGNY_SYMMETRIES = 1
GOGNY_HYPER      = 1
USE_QRPA         = 0

# Names and paths to BLAS and LAPACK Libraries
#  * BG/Q ARCHITECTURES (ESSL)
#      ALCF (OpenMP) ....: -L/soft/libraries/alcf/current/xl/LAPACK/lib -llapack \
#                          -L/soft/libraries/alcf/current/xl/BLAS/lib -lblas
#      ALCF (ESSL) ......: -L/soft/libraries/alcf/current/xl/LAPACK/lib -llapack \
#                          -L/soft/libraries/essl/current/lib64 -lesslsmpbg \
#                          -L/soft/compilers/ibmcmp-may2014/xlf/bg/14.1/bglib64 \
#                          -lxlf90_r -lxlfmath -lxlopt -lxl -Wl,-E
#      Vulcan (ESSL) ....: -L/usr/local/tools/essl/$(VERSION_LIBRARY)/lib -lesslsmpbg
#  * LINUX CLUSTER(MKL OR ACML)
#      LC (Cab,OpenMP) ..: -L/usr/local/tools/$(VERSION_LIBRARY)/lib -I/usr/local/tools/$(VERSION_LIBRARY)/include \
#                          -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#      LC (Quartz,OpenMP): -L/usr/tce/packages/mkl/$(VERSION_LIBRARY)/lib -I/usr/tce/packages/mkl/$(VERSION_LIBRARY)/include \
#                          -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#      LC (Quartz,serial): -L/usr/tce/packages/mkl/$(VERSION_LIBRARY)/lib -I/usr/tce/packages/mkl/$(VERSION_LIBRARY)/include \
#                          -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#      NERSC ............: -mkl (after: module swap PrgEnv-pgi PrgEnv-intel)
#      OLCF (OpenMP).....: -L$(ACML_DIR)/ifort64_mp/lib -lacml_mp -liomp5 -lifcoremt_pic -limf -lirc -lsvml
#      OLCF .............: -L$(ACML_DIR)/ifort64/lib -lacml -lifcoremt_pic -limf -lirc -lsvml
#  * LINUX DESKTOP
#      local ............: -L$(HOME)/local -llapack -lblas
LINEAR_ALGEBRA = -L$(HOME)/local -llapack -lblas

#======================================================================#
#  Nothing beyond this line should be changed, in principle            #
#======================================================================#

# Check if compile directory exists, returns error code and exit otherwise
EXISTS = $(firstword $(wildcard $(COMPILE_DIR)) )
$(info $$COMPILE_DIR is [${COMPILE_DIR}] )
ifeq ($(EXISTS), )
      $(error ${COMPILE_DIR} directory does not exist)
endif

# Consistency checks between MPI option and dripline or mass table modes
ifeq ($(DRIP_LINES),1)
     ifeq ($(DO_MASSTABLE),1)
           $(error DRIP_LINES and DO_MASSTABLE modes are incompatible)
     endif
     ifneq ($(USE_MPI),1)
            $(error DRIP_LINES mode requires USE_MPI=1)
     endif
endif
ifeq ($(USE_MPI),1)
      ifeq ($(DRIP_LINES),0)
            ifeq ($(DO_MASSTABLE),0)
                  $(error USE_MPI=1 not allowed for regular HFBTHO calculation)
            endif
      endif
endif

# Names
HFBTHO_EXE     = hfbtho_main
HFBTHO_SOURCE  = hfbtho_main.f90
HFBTHO_OBJ     = hfbtho_main.o

# Defining compiler options for: IFORT FORTRAN COMPILER (ifort)
ifeq ($(COMPILER),IFORT)

      FORMAT_F90   = -free -extend_source
      STATIC       =
      PREPROCESSOR = -fpp -DUSE_OPENMP=$(USE_OPENMP) \
                          -DUSE_MPI=$(USE_MPI)  \
                          -DUSE_PETSC=$(USE_PETSC) \
                          -DDO_MASSTABLE=$(DO_MASSTABLE) \
                          -DDRIP_LINES=$(DRIP_LINES) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL) \
                          -DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) \
                          -DGOGNY_HYPER=$(GOGNY_HYPER) \
                          -DUSE_QRPA=$(USE_QRPA) \
                          -DSWITCH_ESSL=$(SWITCH_ESSL) $(FFLAGS)

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -O3
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -check all -g -traceback
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS = $(OPTIONS_FC) -qopenmp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

# Defining compiler options for: GNU FORTRAN COMPILER (gfortran)
ifeq ($(COMPILER),GFORTRAN)

      FORMAT_F90   = -ffree-form -ffree-line-length-none
      STATIC       =
      PREPROCESSOR = -cpp -DUSE_OPENMP=$(USE_OPENMP) \
                          -DUSE_MPI=$(USE_MPI)  \
                          -DUSE_PETSC=$(USE_PETSC) \
                          -DDO_MASSTABLE=$(DO_MASSTABLE) \
                          -DDRIP_LINES=$(DRIP_LINES) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL) \
                          -DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) \
                          -DGOGNY_HYPER=$(GOGNY_HYPER) \
                          -DUSE_QRPA=$(USE_QRPA) \
                          -DSWITCH_ESSL=$(SWITCH_ESSL) $(FFLAGS)

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -O2
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC)  -g -O0 -Wall \
                         -Warray-bounds -Wunderflow -Warray-temporaries \
                         -Wcharacter-truncation -Wtabs -Wintrinsic-shadow -Walign-commons -frange-check \
                         -fbounds-check -Wconversion -Wuninitialized $(PEDANTIC) \
                         -finit-real=nan \
                         -ftrapv
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS = $(OPTIONS_FC) -fopenmp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

# Defining compiler options for: CRAY COMPILER
ifeq ($(COMPILER),CRAY)

      FORMAT_F90   = -f free
      STATIC       =
      PREPROCESSOR = -e Z -DUSE_OPENMP=$(USE_OPENMP) \
                          -DUSE_MPI=$(USE_MPI)  \
                          -DUSE_PETSC=$(USE_PETSC) \
                          -DDO_MASSTABLE=$(DO_MASSTABLE) \
                          -DDRIP_LINES=$(DRIP_LINES) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL) \
                          -DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) \
                          -DGOGNY_HYPER=$(GOGNY_HYPER) \
                          -DUSE_QRPA=$(USE_QRPA) \
                          -DSWITCH_ESSL=$(SWITCH_ESSL) $(FFLAGS)

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -O3
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -e c -e D
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS = $(OPTIONS_FC)
      else
            OPTIONS = $(OPTIONS_FC) -h noomp
      endif

endif

# Defining compiler options for: CRAY COMPILER
ifeq ($(COMPILER),IBM)

      FORMAT_F90   = -qstrict -qfree=f90 -qsuffix=cpp=f90
      STATIC       =
      PREPROCESSOR = '-WF,-DUSE_OPENMP=$(USE_OPENMP)' \
                     '-WF,-DUSE_MPI=$(USE_MPI)' \
                     '-WF,-DUSE_PETSC=$(USE_PETSC)'\
                     '-WF,-DDO_MASSTABLE=$(DO_MASSTABLE) '\
                     '-WF,-DDRIP_LINES=$(DRIP_LINES) '\
                     '-WF,-DREAD_FUNCTIONAL=$(READ_FUNCTIONAL) '\
                     '-WF,-DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) '\
                     '-WF,-DGOGNY_HYPER=$(GOGNY_HYPER) '\
                     '-WF,-DUSE_QRPA=$(USE_QRPA)' \
                     '-WF,-DSWITCH_ESSL=$(SWITCH_ESSL) ' $(FFLAGS)

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -qhot -O2
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -g -C -qflttrap
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS = $(OPTIONS_FC) -qsmp=omp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

#=========================#
# Beginning of the action #
#=========================#

HFBTHO_DIR = $(COMPILE_DIR)
DOC_DIR    = $(HFBTHO_DIR)/doc

# Object files
HFBTHO_VERSION_OBJ      = hfbtho_version.o
HFBTHO_UTILITIES_OBJ    = hfbtho_utilities.o
HFBTHO_UNEDF_OBJ        = hfbtho_unedf.o
HFBTHO_VARIABLES_OBJ    = hfbtho_variables.o
HFBTHO_ELLIPINT_OBJ     = hfbtho_elliptic_integrals.o
HFBTHO_BESSEL_OBJ       = hfbtho_bessel.o
HFBTHO_GAUSS_OBJ        = hfbtho_gauss.o
HFBTHO_LINALGEB_OBJ     = hfbtho_linear_algebra.o
HFBTHO_THO_OBJ          = hfbtho_tho.o
HFBTHO_MULTIPOLE_OBJ    = hfbtho_multipole_moments.o
HFBTHO_COLLECTIVE_OBJ   = hfbtho_collective.o
HFBTHO_FISSION_OBJ      = hfbtho_fission.o
HFBTHO_SOLVER_OBJ       = hfbtho_solver.o
HFBTHO_STORAGE_OBJ      = hfbtho_storage.o
HFBTHO_MPIMT_OBJ        = hfbtho_large_scale.o
HFBTHO_FUNCTIONAL_OBJ   = hfbtho_read_functional.o
HFBTHO_GOGNY_OBJ        = hfbtho_gogny.o
HFBTHO_LOCALIZATION_OBJ = hfbtho_localization.o
HFBTHO_IO_OBJ           = hfbtho_io.o


# Options are: anl, $(HFBTHO_EXE), doc (using Doxygen)
all: $(HFBTHO_EXE) doc

# Doxygen documentation
doc: $(HFBTHO_EXE)
	( cd $(HFBTHO_DIR); $(MAKE) -C $(DOC_DIR); cd .. )

# HFBTHO solver
$(HFBTHO_VERSION_OBJ) : hfbtho_version.f90
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $<
$(HFBTHO_UTILITIES_OBJ) : hfbtho_utilities.f90
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $<
$(HFBTHO_UNEDF_OBJ) : hfbtho_unedf.f90 $(HFBTHO_UTILITIES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_VARIABLES_OBJ) : hfbtho_variables.f90 $(HFBTHO_VERSION_OBJ) $(HFBTHO_UNEDF_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_IO_OBJ) : hfbtho_io.f90 $(HFBTHO_VERSION_OBJ) $(HFBTHO_UNEDF_OBJ) $(HFBTHO_COLLECTIVE_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_ELLIPINT_OBJ) : hfbtho_elliptic_integrals.f90 $(HFBTHO_UTILITIES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_BESSEL_OBJ) : hfbtho_bessel.f90 $(HFBTHO_UTILITIES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_GAUSS_OBJ) : hfbtho_gauss.f90 $(HFBTHO_UTILITIES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_GOGNY_OBJ) : hfbtho_gogny.f90 $(HFBTHO_UTILITIES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_MPIMT_OBJ) : hfbtho_large_scale.f90 $(HFBTHO_UTILITIES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_FUNCTIONAL_OBJ) : hfbtho_read_functional.f90 $(HFBTHO_UTILITIES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_LINALGEB_OBJ) : hfbtho_linear_algebra.f90 $(HFBTHO_UTILITIES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $^
$(HFBTHO_THO_OBJ) : hfbtho_tho.f90 $(HFBTHO_LINALGEB_OBJ) $(HFBTHO_VARIABLES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_MULTIPOLE_OBJ) : hfbtho_multipole_moments.f90 $(HFBTHO_UTILITIES_OBJ) $(HFBTHO_VARIABLES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_FISSION_OBJ) : hfbtho_fission.f90 $(HFBTHO_UTILITIES_OBJ) $(HFBTHO_VARIABLES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_COLLECTIVE_OBJ) : hfbtho_collective.f90 $(HFBTHO_VARIABLES_OBJ) $(LAPACK_OBJ) $(BLAS_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_SOLVER_OBJ) : hfbtho_solver.f90 $(HFBTHO_THO_OBJ) $(HFBTHO_GAUSS_OBJ) $(HFBTHO_BESSEL_OBJ) \
                                         $(HFBTHO_ELLIPINT_OBJ) $(HFBTHO_VARIABLES_OBJ) \
                                         $(HFBTHO_MULTIPOLE_OBJ) $(HFBTHO_FISSION_OBJ) $(HFBTHO_COLLECTIVE_OBJ) \
                                         $(LAPACK_OBJ) $(BLAS_OBJ) $(HFBTHO_GOGNY_OBJ) $(HFBTHO_IO_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^

# Extensions: Interface with pnFAM-QRPA and localization functions
$(HFBTHO_STORAGE_OBJ) : hfbtho_storage.f90 $(HFBTHO_VARIABLES_OBJ) $(HFBTHO_SOLVER_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^
$(HFBTHO_LOCALIZATION_OBJ) : hfbtho_localization.f90 $(HFBTHO_VARIABLES_OBJ) $(HFBTHO_UTILITIES_OBJ) \
                                                     $(HFBTHO_SOLVER_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^

# HFBTHO (1/3): object file
$(HFBTHO_OBJ) : $(HFBTHO_SOURCE) $(HFBTHO_SOLVER_OBJ) $(HFBTHO_MPIMT_OBJ) \
                $(HFBTHO_FUNCTIONAL_OBJ) $(HFBTHO_LOCALIZATION_OBJ) $(HFBTHO_STORAGE_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -c $^

# HFBTHO (2/3): all object files from Fortran modules
GLOBAL_OBJECTS = $(HFBTHO_VERSION_OBJ) $(HFBTHO_UTILITIES_OBJ) $(HFBTHO_UNEDF_OBJ) \
                 $(HFBTHO_VARIABLES_OBJ) $(HFBTHO_ELLIPINT_OBJ) $(HFBTHO_BESSEL_OBJ) \
                 $(HFBTHO_GAUSS_OBJ) $(HFBTHO_LINALGEB_OBJ) $(HFBTHO_THO_OBJ) \
                 $(HFBTHO_MULTIPOLE_OBJ) $(HFBTHO_COLLECTIVE_OBJ) $(HFBTHO_FISSION_OBJ) \
                 $(HFBTHO_SOLVER_OBJ) $(HFBTHO_STORAGE_OBJ) $(HFBTHO_MPIMT_OBJ) \
                 $(HFBTHO_FUNCTIONAL_OBJ) $(HFBTHO_GOGNY_OBJ) $(HFBTHO_LOCALIZATION_OBJ) \
                 $(HFBTHO_IO_OBJ)

# HFBTHO (3/3): executable
$(HFBTHO_EXE) : $(HFBTHO_OBJ) $(GLOBAL_OBJECTS)
	$(FORTRAN_MPI) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LINEAR_ALGEBRA)

# Cleaning
clean ::
	$(MAKE) clean -C $(DOC_DIR)
	-rm -f *.o *.oo *.ipo *.mod
