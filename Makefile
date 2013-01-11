##################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
##################################################################

##################################################################
MYOS        := $(shell uname)
ARCH        := $(shell arch)
USER        := $(shell whoami)
MYHOST      := $(shell hostname -s)

##################################################################
VERSION     := 1.0.0
EXECFILE    := g2pSim

##################################################################
SRCDIR      := src
INCDIR      := include
OBJDIR      := obj.$(ARCH)
OTHERINC    := 

###################################################################
# Compiler
ifeq ($(MYOS),Darwin)
    AR      := libtool -static -o
else
    AR      := ar r
endif
CXX         := g++
FF          := gfortran
LD          := g++

###################################################################
# Flags
ifeq ($(ARCH),i686)
    MODE    := -m32
else
    MODE    := -m64
endif
CFLAGS      := -Wall -Wno-unused-variable -O3 -g \
               $(MODE) -I$(INCDIR) $(OTHERINC) 
CXXFLAGS    := -Wall -Wno-unused-variable -O3 -g -Wno-deprecated \
               $(MODE) -I$(INCDIR) $(OTHERINC) 
FFLAGS      := -Wall -Wno-unused-variable -O3 -g \
               $(MODE) -I$(INCDIR) $(OTHERINC)
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
GPPFLAGS    := -M
LDFLAGS     := -O3 -g $(MODE)

###################################################################
# Generate obj file list
FSOURCES    := $(wildcard $(SRCDIR)/*.[Ff])
CSOURCES    := $(wildcard $(SRCDIR)/*.[Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][XxPp][XxPp])
SOURCES     := $(FSOURCES) $(CSOURCES)
#add .o to all the source files
#######OBJS  = $(OBJS:.f=.o)#######
OBJS        := $(addsuffix .o, $(basename $(SOURCES)))
OBJS        := $(patsubst  $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))
#do not make dependence for fortran code
COBJS       := $(addsuffix .o, $(basename $(CSOURCES)))
COBJS       := $(patsubst  $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(COBJS))
DEPS        := $(subst .o,.d,$(COBJS))

###################################################################
# Libs
SYSLIBS     := -lstdc++
ifeq ($(MYOS),Darwin) # Assume using homebrew
    SYSLIBS += -L/usr/local/Cellar/gfortran/4.7.2/gfortran/lib -lgfortran
else
    SYSLIBS += -lgfortran
endif
OTHERLIBS   := -LHRSTransport/obj.${ARCH} -lHRSTransport \
               -LCrossSection/obj.${ARCH} -lCrossSection

###################################################################
# ROOT configure
ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs)
ROOTGLIBS   := $(shell root-config --glibs) -lMinuit

CXXFLAGS    += $(ROOTCFLAGS) 
LIBS        := $(SYSLIBS) $(ROOTLIBS) $(OTHERLIBS)
GLIBS       := $(SYSLIBS) $(ROOTGLIBS) $(OTHERLIBS)


###################################################################
# You can specify the .SUFFIXES
.SUFFIXES: .C .c .CC .cc .CPP .cpp .CXX .cxx .F .f
.PHONY: all clean test
VPATH    := $(SRCDIR):$(INCDIR):$(OTHERINC) 

###################################################################
all: exe

###################################################################
# Make the $(TARGET).d file and include it.
$(OBJDIR)/%.d: %.c 
	@echo Making dependency for file $< ......
	@set -e;\
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) -w $< |\
	sed 's!$*\.o!$(OBJDIR)/& $@!' >$@;\
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.C 
	@echo Making dependency for file $< ......
	@set -e;\
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) -w $< |\
	sed 's!$*\.o!$(OBJDIR)/& $@!' >$@;\
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cc
	@echo Making dependency for file $< ......
	@set -e;\
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) -w $< |\
	sed 's!$*\.o!$(OBJDIR)/& $@!' >$@;\
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cpp 
	@echo Making dependency for file $< ......
	@set -e;\
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) -w $< |\
	sed 's!$*\.o!$(OBJDIR)/& $@!' >$@;\
	[ -s $@ ] || rm -f $@
$(OBJDIR)/%.d: %.cxx
	@echo Making dependency for file $< ......
	@set -e;\
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) -w $< |\
	sed 's!$*\.o!$(OBJDIR)/& $@!' >$@;\
	[ -s $@ ] || rm -f $@

#there is a problem to build dependency for fortran code...
$(OBJDIR)/%.d: %.f 
	@echo Making dependency for file $< ......
	@set -e;\
	$(FF) $(GPPFLAGS) $(FFLAGS) -w  $< |\
	sed 's!$*\.o!$(OBJDIR)/& $@!' >$@;\
	[ -s $@ ] || rm -f $@ 
	@-rm -f $*.o 

$(OBJDIR)/%.d: %.F 
	@echo Making dependency for file $< ......
	@set -e;\
	$(FF) $(GPPFLAGS) $(FFLAGS) -w $< |\
	sed 's!$*\.o!$(OBJDIR)/& $@!' >$@;\
	[ -s $@ ] || rm -f $@ 
	@-rm -f $*.o 

ifneq ($(DEPS),)
-include $(DEPS)
endif

##########################################################
exe: $(OBJDIR) $(OBJS)
	@make -C HRSTransport lib
	@make -C CrossSection lib
	@$(LD) $(LDFLAGS) -o $(EXECFILE) $(OBJS) $(LIBS)
	@echo "Linking $(EXECFILE) ... done!"

##########################################################
$(OBJDIR)/%.o: %.F
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

$(OBJDIR)/%.o: %.f
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

$(OBJDIR)/%.o: %.c
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.C
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cpp
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR):
	@if [ ! -d $(OBJDIR) ] ; then mkdir -p $(OBJDIR) ;fi

##########################################################
clean:
	@rm -f $(OBJS) $(DEPS)
	@rm -f *~ *# $(EXECFILE) */*~ */*#

distclean: clean
	@cd HRSTransport; make clean; cd ..
	@cd CrossSection; make clean; cd ..

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

