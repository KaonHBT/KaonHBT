
#OPT           =-02
OPT           =
 
CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -Wall -fPIC
LD            = gfortran
LDFLAGS       =
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS) -lEG 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)


#-DCERNLIB_TYPE -DCERNLIB_DOUBLE -DCERNLIB_NOQUAD -DCERNLIB_LINUX -Dsl64_gcc482 -m32 -fd-lines-as-code -std=legacy -fno-second-underscore -w -fno-automatic -Wall -W -Wsurprising -fPIC -g -ffixed-line-length-132 -I.sl64_gcc482/include -I. -I/afs/rhic.bnl.gov/star/packages/SL15e/.sl64_gcc482/include -I/afs/rhic.bnl.gov/star/packages/SL15e/asps/Simulation/geant321/include -I/cern/pro/include -I/afs/rhic.bnl.gov/star/packages/SL15e -c .sl64_gcc482/obj/StRoot/LendFSI/cernlib.F -o .sl64_gcc482/obj/StRoot/LendFSI/

# Fortran:
F77           = gfortran
#Per - puvodne F77FLAGS      = -m32 -fd-lines-as-code  -std=legacy -fno-second-underscore -w -fno-automatic -Wall -W -Wsurprising -fPIC
F77FLAGS      =  -m32 -fd-lines-as-code -std=legacy -fno-second-underscore -w -fno-automatic -Wall -W -Wsurprising -fPIC -g -ffixed-line-length-132 
F77OPT        = $(OPT)
CXXOUT        = -o # keep whitespace after "-o"

#linker
#Petr - poznamka - pokud pouziju g++ a ne gfortran bude problem
LD           =gfortran
# Petr - puvodne LD_FLAGS     =-g -m32 -shared -Wl,-Bdynamic 
LD_FLAGS     =-g -m32 -shared -Wl,-Bdynamic -lgfortran  
#----------------------------

SHAREDSO       =fsiw110707_test.so 
 
 
#------------------------------------------------------------------------------
all:     $(SHAREDSO) 

clean:
	@rm -f $(SHAREDSO) fsiw110707_test.o cernlib.o

fsiw110707_test.so: fsiw110707_test.o cernlib.o
	$(LD) $(LD_FLAGS) -o $@  $^

%.o : %.cxx
	@echo cxx: $@
	$(CXX) -O2 $(CXXFLAGS) -c $<

%.o: %.f
	$(F77) $(F77OPT) $(F77FLAGS) $(CXXOUT)$@ -c $<


