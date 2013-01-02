###################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $CERN_ROOT and $ROOTSYS defined
###################################################################
MYOS     := $(shell uname)
ARCH     := $(shell arch)
USER     := $(shell whoami)
MYHOST   := $(shell hostname -s)

##################################################################
VERSION   = 1.1
TARGET    = SNAKE
EXECFILE  = TestSNAKE

SRCDIR   := src
INCDIR   := include 

OTHERINC := 
OTHERLIBS := -lgfortran -LHRSTransport/obj.${ARCH} -lHRSTransport

###################################################################
TARGETLIB = $(TARGET)_v$(VERSION)
OBJDIR   := obj

###################################################################
#in the darwin machines, ar cmd has something wrong!!!
#I switch to libtool 
#USAGE: libtool -static -o out.lib ObjList.o
AR       := ar r
ifeq ($(MYOS),Darwin)
AR       := libtool -static -o
endif

###################################################################
ifeq ($(ARCH),i686)
MODE     := -m32
else
MODE     := -m64
endif
GPPFLAGS := -M
CXX      := g++
FF       := gfortran
CXXFLAGS := -Wall -O3 -g -Wno-deprecated $(MODE) -I$(INCDIR) $(OTHERINC) 
CFLAGS   := -Wall -O3 -g $(MODE) -I$(INCDIR) $(OTHERINC) 
FFLAGS   := -Wall -O3 -g $(MODE) -I$(INCDIR) $(OTHERINC) 
LD       := g++
LDFLAGS  := -O3 -g $(MODE)  
SYSLIBS  := -lstdc++

###################################################################
#you can specify the .SUFFIXES
.SUFFIXES: .C .c .CC .cc .CPP .cpp .CXX .cxx .F .f
.PHONY: all clean test
VPATH := $(SRCDIR):$(INCDIR):$(OTHERINC) 

###################################################################
#generate obj file list
FSOURCES  := $(wildcard $(SRCDIR)/*.[Ff])
CSOURCES  := $(wildcard $(SRCDIR)/*.[Cc])
CSOURCES  += $(wildcard $(SRCDIR)/*.[Cc][Cc])
CSOURCES  += $(wildcard $(SRCDIR)/*.[Cc][XxPp][XxPp])
SOURCES   := $(FSOURCES) $(CSOURCES)
#add .o to all the source files
#an example of replace/patrern replace statement
#######OBJS  = $(OBJS:.f=.o)#######
OBJS  := $(addsuffix .o, $(basename $(SOURCES)))
OBJS  := $(patsubst  $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))
#do not make dependence for fortran code
COBJS := $(addsuffix .o, $(basename $(CSOURCES)))
COBJS := $(patsubst  $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(COBJS))
DEPS  := $(subst .o,.d,$(COBJS))

###################################################################
ifneq ($(FSOURCES),)
ifeq ($(FF),gfortran)
SYSLIBS   += -lgfortran
endif
endif

###################################################################
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) -lMinuit $(SYSLIBS) 

ifneq ($(FSOURCES),)
ifdef CERN_ROOT
FFLAGS       += -I$(CERN_ROOT)/include 
GLIBS        += -L$(CERN_ROOT)/lib -lpdflib804 -lmathlib -lphtools \
 -lgeant321 -lpawlib -lgraflib -lgrafX11 -lpacklib -lkernlib 
endif
endif

CXXFLAGS     += $(ROOTCFLAGS) 
LDFLAGS      += $(GLIBS) $(OTHERLIBS) 

###################################################################
all: exe

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
lib: $(OBJDIR) $(OBJS)
	@make -C HRSTransport lib
	@rm -f $(OBJDIR)/lib$(TARGETLIB).a
	@$(AR) $(OBJDIR)/lib$(TARGETLIB).a $(OBJS)
	@if [ -e $(OBJDIR)/lib$(TARGETLIB).a ] ; then echo "Finish creating target lib"  $(OBJDIR)/lib$(TARGETLIB).a ;fi
#	@libtool -static  -o $(OBJDIR)/lib$(TARGETLIB).a $(OBJS)
#	@ar r $(OBJDIR)/lib$(TARGETLIB).a $(OBJS)

exe: lib
#	@if [ ! -d Graph ]; then mkdir -p Graph; fi
	$(LD) -o $(EXECFILE) -L$(OBJDIR) -l$(TARGETLIB) $(LDFLAGS)
	@echo "Linking $(EXECFILE) ...done!"

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

##########################################################
$(OBJDIR):
	@if [ ! -d $(OBJDIR) ] ; then mkdir -p $(OBJDIR) ;fi

clean:
	@rm -f $(OBJS) $(DEPS) $(OBJDIR)/lib$(TARGETLIB).a
	@rm -f *~ *# $(EXECFILE) */*~ */*#

test:	
	@echo \\MYOS\:$(MYOS) \\ARCH\:$(ARCH)
	@echo \\CFLAGS\:$(CFLAGS)	
	@echo \\CXXFLAGS\:$(CXXFLAGS)        
	@echo \\FFLAGS\:$(FFLAGS)
	@echo \\LDFLAGS\:$(LDFLAGS)
	@echo \\SYSLIBS\:$(SYSLIBS)
	@echo \\Fsources\: $(FSOURCES)	
	@echo \\sources\: $(SOURCES)	
	@echo \\objs\: $(OBJS)	
	@echo \\dependencies: \$(DEPS)

help: test

env: test

#@echo \\sources1\: $(SOURCES1)
#@echo \\obj1\: $(OBJS1)
##########################################################
