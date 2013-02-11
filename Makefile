########################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
########################################################################

########################################################################
MYOS        := $(shell uname)
ARCH        := $(shell arch)
USER        := $(shell whoami)
MYHOST      := $(shell hostname -s)

########################################################################
VERSION     := 1.3.1
EXECFILE    := G2PSim
LIBFILE     := libG2PSim.so
LIBNAME     := G2PSim
USERDICT    := $(LIBNAME)_Dict

########################################################################
SRCDIR      := src
INCDIR      := include
OBJDIR      := obj.$(ARCH)

########################################################################
# Compiler
AR          := ar
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
CFLAGS      := -Wall -fPIC -O3 -g $(MODE) -I$(INCDIR)
CXXFLAGS    := -Wall -fPIC -O3 -g $(MODE) -I$(INCDIR) 
FFLAGS      := -Wall -fPIC -O3 -g $(MODE) -I$(INCDIR)
ifeq ($(MYOS),Darwin) 
#in Darwin, do not use -fno-leading-underscore
    FFLAGS  += -fno-second-underscore -fno-automatic -fbounds-check \
               -fno-range-check -funroll-all-loops -fdollar-ok \
               -ffixed-line-length-none
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
HEADERS     := $(wildcard $(INCDIR)/*.hh)
HEADERS     += $(wildcard HRSTransport/*.hh)
HEADERS     += $(wildcard G2PXSection/*.hh)
# add .o to all the source files
OBJS        := $(addsuffix .o, $(basename $(SOURCES)))
OBJS        := $(patsubst  $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))
DEPS        := $(subst .o,.d,$(OBJS))

########################################################################
# Libs
SYSLIBS     := -lstdc++
ifeq ($(MYOS),Darwin) # Assume using homebrew
    SYSLIBS += -L/usr/local/Cellar/gfortran/4.7.2/gfortran/lib -lgfortran
else
    SYSLIBS += -lgfortran
endif
OTHERLIBS   := -LHRSTransport -lHRSTransport -LG2PXSection -lG2PXS

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
VPATH       := $(SRCDIR):$(INCDIR) 

########################################################################
all: exe

########################################################################
# Make the $(TARGET).d file and include it.
$(OBJDIR)/%.d: %.c 
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
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

$(OBJDIR)/%.d: %.f
	@echo Making dependency for file $< ......
	@set -e; \
	$(FF) -cpp $(GPPFLAGS) $(FFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.F
	@echo Making dependency for file $< ......
	@set -e; \
	$(FF) -cpp $(GPPFLAGS) $(FFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

ifneq ($(DEPS),)
-include $(DEPS)
endif

########################################################################
exe: lib script $(OBJDIR)/Main.o
	@$(LD) $(LDFLAGS) -o $(EXECFILE) $(OBJDIR)/Main.o ./$(LIBFILE) $(LIBS)
	@echo "Linking $(EXECFILE) ... done!"

$(OBJDIR)/Main.o: Main.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

script:
	@if [ ! -e "Run.C" ] ; then cp -pr "scripts/RunTmp.C" "Run.C" ; \
	echo "Generate Run.C ... done!"; fi

########################################################################
lib: $(OBJDIR) $(OBJS) $(OBJDIR)/$(USERDICT).o
	@make -s -C HRSTransport
	@make -s -C G2PXSection
	@$(LD) -shared $(LDFLAGS) -o $(LIBFILE).$(VERSION) \
           $(OBJS) $(OBJDIR)/$(USERDICT).o $(LIBS) $(OTHERLIBS)
	@ln -sf $(LIBFILE).$(VERSION) $(LIBFILE)
	@echo "Linking $(LIBFILE) ... done!"

$(USERDICT).cxx: $(HEADERS) $(LIBNAME)_LinkDef.h
		@echo "Generating dictionary $(USERDICT)..."
		@$(ROOTSYS)/bin/rootcint -f $@ -c $(CXXFLAGS) $^

$(OBJDIR)/$(USERDICT).o: $(USERDICT).cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

########################################################################
$(OBJDIR)/%.o: %.c
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

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

$(OBJDIR):
	@if [ ! -d $(OBJDIR) ] ; then mkdir -p $(OBJDIR) ;fi

########################################################################
clean: $(OBJDIR)
	@rm -f $(OBJDIR)/*
	@rm -f $(USERDICT).cxx $(USERDICT).h
	@rm -f $(EXECFILE) $(LIBFILE) $(LIBFILE).*
	@rm -f *~ *# */*~ */*#

distclean: clean
	@make distclean -s -C HRSTransport
	@make distclean -s -C G2PXSection

test:	
	@echo \\MYOS\:$(MYOS) \\ARCH\:$(ARCH)
	@echo \\CFLAGS\:$(CFLAGS)	
	@echo \\CXXFLAGS\:$(CXXFLAGS)        
	@echo \\FFLAGS\:$(FFLAGS)
	@echo \\LDFLAGS\:$(LDFLAGS)
	@echo \\SYSLIBS\:$(SYSLIBS)
	@echo \\fsources\: $(FSOURCES)	
	@echo \\sources\: $(SOURCES)	
	@echo \\objs\: $(OBJS)	
	@echo \\dependencies: \$(DEPS)

help: test

env: test
