HEALPIX = /Users/guillermo/Documents/Healpix_3.82
SHARPLDIR = /Users/guillermo/Documents/Healpix_3.82/lib
MKLROOT = /opt/intel/oneapi/mkl/2022.1.0

F90_BINDIR = $(HEALPIX)/bin
F90_INCDIR = $(HEALPIX)/include
F90_LIBDIR = $(HEALPIX)/lib
F90_BUILDDIR = $(HEALPIX)/build
FITSDIR = /opt/local/lib
LIBFITS = cfitsio

F90_FC = ifort
F90_FFLAGS = -O3 -r8 -I$(F90_INCDIR) -cm -w -sox -qopt-report=0 -qopenmp -fPIC -heap-arrays 256
F90_LDFLAGS = -L$(F90_LIBDIR) -L$(FITSDIR) -L$(SHARPLDIR) -lhealpix -lhpxgif -lsharp -l$(LIBFITS) -Wl,-rpath,$(FITSDIR) -Wl,-rpath,$(SHARPLDIR) -Wl,-rpath,$(F90_LIBDIR) -lcurl -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/../../compiler/2022.1.0/lib -mkl

# Program name
PROGRAM = hilc

# Source file
SRCS = hilc.f90

$(PROGRAM): $(SRCS)
	$(F90_FC) $(F90_FFLAGS) $^ -o $@ $(F90_LDFLAGS)
