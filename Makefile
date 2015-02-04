########################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
########################################################################

########################################################################
MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
USER        := $(shell whoami)
MYHOST      := $(shell hostname -s)

########################################################################
EXECFILE    := g2psim
RECEXECFILE := g2prec
LIBFILE     := libG2PSim.so
LIBNAME     := g2psim
USERDICT    := $(LIBNAME)_dict
VERSION     := 1.9.0

########################################################################
SRCDIR      := src
INCDIR      := include
OBJDIR      := obj.$(ARCH)

########################################################################
MODELDIR    := HRSTrans:G2PPhys
MODELLIST   := HRSTrans G2PPhys
INCDIR      := $(INCDIR):$(MODELDIR)

########################################################################
# Compiler
AR          := ar
CC          := gcc
CXX         := g++
FF          := gfortran
LD          := g++

########################################################################
# Flags
ifeq ($(ARCH),i686)
    MODE    := -m32
else
    MODE    := -m64
endif
INCDIRS     := $(patsubst %,-I%,$(subst :, ,$(INCDIR)))
CFLAGS      := -Wall -fPIC -O3 -g $(MODE)
CXXFLAGS    := -Wall -fPIC -O3 -g $(MODE)
FFLAGS      := -Wall -fPIC -O3 -g $(MODE)
ifeq ($(MYOS),Darwin)
#in Darwin, do not use -fno-leading-underscore
    FFLAGS  += -fno-second-underscore -fno-automatic -fbounds-check \
               -funroll-all-loops -fdollar-ok -ffixed-line-length-none \
               -fno-range-check
else
    FFLAGS  += -fno-leading-underscore -fno-second-underscore \
               -fno-automatic -fbounds-check -funroll-all-loops \
               -fdollar-ok -ffixed-line-length-none -fno-range-check
endif
GPPFLAGS    := -MM
LDFLAGS     := -O3 -g $(MODE)

########################################################################
# Generate obj file list
FSOURCES    := $(wildcard $(SRCDIR)/*.[Ff])
CSOURCES    := $(wildcard $(SRCDIR)/*.[Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][XxPp][XxPp])
SOURCES     := $(FSOURCES) $(CSOURCES)
# header files
HEADERS     := $(foreach n,$(subst :, ,$(INCDIR)),$(wildcard $(n)/*.hh))
# add .o to all the source files
OBJS        := $(addsuffix .o, $(basename $(SOURCES)))
OBJS        := $(patsubst  $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))
DEPS        := $(subst .o,.d,$(OBJS))

########################################################################
# Libs
SYSLIBS     := -lstdc++ -lgfortran

ifdef LIBCONFIG
    ifneq ($(LIBCONFIG),/usr)
        INCDIRS += -I$(LIBCONFIG)/include
    endif
    SYSLIBS += -L$(LIBCONFIG)/lib -lconfig
else
$(error $$LIBCONFIG environment variable not defined)
endif

ifdef ANALYZER
ANADIRS     := $(wildcard $(addprefix $(ANALYZER)/, src hana_decode hana_scaler))
ANACXXFLAGS := $(addprefix -I, $(ANADIRS))
ANALIBS     := -L$(ANALYZER) -lHallA -ldc -lscaler
else
ANACXXFLAGS :=
ANALIBS     :=
endif

CFLAGS      += $(INCDIRS)
CXXFLAGS    += $(INCDIRS)
FFLAGS      += $(INCDIRS)

OTHERLIBS   := -LHRSTrans -lHRSTrans -LG2PPhys -lG2PPhys

########################################################################
# ROOT configure
ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs)
ROOTGLIBS   := $(shell root-config --glibs) -lMinuit

CXXFLAGS    += $(ROOTCFLAGS)
LIBS        := $(SYSLIBS) $(ROOTLIBS)
GLIBS       := $(SYSLIBS) $(ROOTGLIBS)

########################################################################
# You can specify the .SUFFIXES
.SUFFIXES: .c .C .cc .CC .cpp .cxx .f .F
.PHONY: all clean test
VPATH       := $(SRCDIR)

########################################################################
all: lib script exe

########################################################################
# Make the $(TARGET).d file and include it.
$(OBJDIR)/%.d: %.c
	@echo Making dependency for file $< ......
	@set -e; \
	$(CC) $(GPPFLAGS) $(CFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.C
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cc
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.CC
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cpp
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cxx
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

#$(OBJDIR)/%.d: %.f
#	@echo Making dependency for file $< ......
#	@set -e; \
#	$(FF) -cpp $(GPPFLAGS) $(FFLAGS) $< | \
#	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
#	[ -s $@ ] || rm -f $@

#$(OBJDIR)/%.d: %.F
#	@echo Making dependency for file $< ......
#	@set -e; \
#	$(FF) -cpp $(GPPFLAGS) $(FFLAGS) $< | \
#	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
#	[ -s $@ ] || rm -f $@

ifneq ($(DEPS),)
-include $(DEPS)
endif

########################################################################
exe: dir $(OBJS) $(OBJDIR)/Main.o $(OBJDIR)/$(USERDICT).o
	@$(LD) $(LDFLAGS) -o $(EXECFILE) $(OBJDIR)/Main.o $(OBJS)\
           $(OBJDIR)/$(USERDICT).o $(LIBS) $(OTHERLIBS)
	@echo "Linking $(EXECFILE) ...... done!"

$(OBJDIR)/Main.o: Main.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

########################################################################
g2prec: dir $(OBJS) $(OBJDIR)/Rec.o $(OBJDIR)/$(USERDICT).o
	@$(LD) $(LDFLAGS) -o $(RECEXECFILE) $(OBJDIR)/Rec.o $(OBJS)\
           $(OBJDIR)/$(USERDICT).o $(LIBS) $(OTHERLIBS) $(ANALIBS)
	@echo "Linking $(RECEXECFILE) ...... done!"

$(OBJDIR)/Rec.o: Rec.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS) $(ANACXXFLAGS)

########################################################################
lib: dir $(OBJS) $(OBJDIR)/$(USERDICT).o
	@for model in $(MODELLIST); do \
		if [ -d $$model ]; then \
			make -s -C $$model; \
		fi; \
	done;
	@$(LD) -shared $(LDFLAGS) -o $(LIBFILE).$(VERSION) \
           $(OBJS) $(OBJDIR)/$(USERDICT).o $(LIBS) $(OTHERLIBS)
	@ln -sf $(LIBFILE).$(VERSION) $(LIBFILE)
	@echo "Linking $(LIBFILE) ...... done!"

script:
	@if [ ! -e "Run.C" ] ; then cp -pr "templates/Run.C" "Run.C" ; \
	echo "Generate Run.C ...... done!"; fi
	@if [ ! -e "sim.cfg" ] ; then cp -pr "templates/sim.cfg" "sim.cfg" ; \
	echo "Generate sim.cfg ...... done!"; fi

$(USERDICT).cxx: $(HEADERS) $(LIBNAME)_linkdef.h
	@echo "Generating dictionary $(USERDICT).cxx ......"
	@$(ROOTSYS)/bin/rootcint -f $@ -c $(INCDIRS) $^

$(OBJDIR)/$(USERDICT).o: $(USERDICT).cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

########################################################################
$(OBJDIR)/%.o: %.c
	@echo Compiling $< ......
	@$(CC) -c $< -o $@  $(CFLAGS)

$(OBJDIR)/%.o: %.C
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.CC
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cpp
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.f
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

$(OBJDIR)/%.o: %.F
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

dir:
	@if [ ! -d $(OBJDIR) ] ; then mkdir -p $(OBJDIR) ;fi

########################################################################
clean: dir
	@rm -f $(OBJDIR)/*
	@rm -f $(USERDICT).cxx $(USERDICT).h
	@rm -f $(EXECFILE) $(RECEXECFILE)
	@rm -f $(LIBFILE) $(LIBFILE).$(VERSION)
	@rm -f *~ *# */*~ */*#

distclean: clean
	@for model in $(MODELLIST); do \
		if [[ -d $$model ]]; then \
			make distclean -s -C $$model; \
		fi; \
	done;

test:
	@echo \\MYOS\:$(MYOS) \\ARCH\:$(ARCH)
	@echo \\CFLAGS\:$(CFLAGS)
	@echo \\CXXFLAGS\:$(CXXFLAGS)
	@echo \\FFLAGS\:$(FFLAGS)
	@echo \\LDFLAGS\:$(LDFLAGS)
	@echo \\SYSLIBS\:$(SYSLIBS)
	@echo \\fsources\: $(FSOURCES)
	@echo \\sources\: $(SOURCES)
	@echo \\headers\: $(HEADERS)
	@echo \\objs\: $(OBJS)
	@echo \\dependencies: \$(DEPS)

help: test

env: test
