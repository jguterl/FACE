# makefile for FACE20
# compiler
#FC     = ifort
.SUFFIXES:
FC     = gfortran
fopenmp=-fopenmp 
FMKL=
# PROGRAM NAME: FACE
PROGRAM = FACE20
# build directory
BUILDDIR=build
# source directory
SRCDIR=src
# directory where binary are created
BINDIR=bin
OBJDIR=$(BUILDDIR)
MODDIR=modules

# compile flags
FFLAGS =-fbounds-check -fcheck=all -Wall -g #-qopenmp -mkl=parallel#-Wall -Werror -Wextra -fno-align-commons 
DBGFLAGS =  -fbacktrace -g
RLSFLAGS = -O3 


# list source files
SOURCES_f90 = $(wildcard $(SRCDIR)/*.f90)
SOURCES_f = $(wildcard $(SRCDIR)/*.f)

# list objects
DBGOBJECTS_f  := $(SOURCES_f:$(SRCDIR)/%.f=$(OBJDIR)/%.o)
DBGOBJECTS_f90  := $(SOURCES_f90:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
RLSOBJECTS_f  := $(SOURCES_f:$(SRCDIR)/%.f=$(OBJDIR)/%.o)
RLSOBJECTS_f90  := $(SOURCES_f90:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)



# program name
EXE_RELEASE=$(PROGRAM)_rls.exe
EXE_DEBUG=$(PROGRAM)_dbg.exe

all: prep release

# debug rules
$(DBGOBJECTS_f90): $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) -cpp ${fopenmp} -g $(FFLAGS) $(DBGFLAGS)  -c $< -o $@ 


$(DBGOBJECTS_f): $(OBJDIR)/%.o : $(SRCDIR)/%.f
	$(FC) -cpp ${fopenmp} -g $(FFLAGS)  $(DBGFLAGS) -c $< -o $@ 

debug: prep debug_exe

debug_exe: $(DBGOBJECTS_f) $(DBGOBJECTS_f90)
	$(FC)  -cpp -g ${fopenmp} $(FFLAGS) $(DBGFLAGS) -o $(BINDIR)/$(EXE_DEBUG) $^


# release rules

$(RLSOBJECTS_f90): $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC)  -cpp ${fopenmp}  $(FFLAGS) $(RLSFLAGS) ${FMKL} -c $< -o $@ -J$(MODDIR)


$(RLSOBJECTS_f): $(OBJDIR)/%.o : $(SRCDIR)/%.f
	$(FC)  -cpp ${fopenmp}  -cpp $(FFLAGS) $(RLSFLAGS) ${FMKL} -c $< -o $@ -J$(MODDIR)

release: prep release_exe

release_exe: $(RLSOBJECTS_f) $(RLSOBJECTS_f90)
	$(FC) -cpp ${fopenmp} $(FFLAGS) $(RLSFLAGS) ${FMKL} -o $(BINDIR)/$(EXE_RELEASE) $^

# default
default: release

# cleaning rules
clean:
	@rm -f $(BUILDDIR)/*.o $(MODDIR)/*.mod $(BINDIR)/$(EXE_RELEASE) $(BINDIR)/$(EXE_DEBUG) 

prep:
	@mkdir -p $(BUILDDIR) ${BINDIR} ${MODDIR}



.PHONY : all clean release release_exe debug debug_exe prep


include $(SRCDIR)/dependency.mk
