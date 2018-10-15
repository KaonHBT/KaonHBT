# _________________________________________________________
# --- FILE
#
# Directories
ifdef NODEBUG
DEBUG = 0
else
DEBUG= 1
endif
ifeq ($(DEBUG),0)
DBin =  BIN
DDObj = OBJ
DObj = OBJ
OptFlag = -O4
else
DBin =  bin
DDObj = obj
DObj = obj
OptFlag = -g
endif
DSource = Code

# Source 
# main
CFBWP = $(DSource)/FitBWPaper.cxx
CFBW = $(DSource)/FitBW.cxx
CCBW = $(DSource)/CalcBWPaper.cxx
CHBTBW = $(DSource)/FitHbtRPBW.cxx
CXiBW = $(DSource)/FitBWXi.cxx
# common
CBWF = $(DSource)/Common/BlastWaveFitter.cxx
HBWF = $(DSource)/Common/BlastWaveFitter.h
CBW = $(DSource)/Common/BlastWave.cxx
HBW = $(DSource)/Common/BlastWave.h
CDFF = $(DSource)/Common/DataForFit.cxx
HDFF = $(DSource)/Common/DataForFit.h
# Multicent
CMulti = $(DSource)/MultiCent/MultiCentFitBW.cxx
CIO = $(DSource)/MultiCent/IOHandler.cxx
HIO = $(DSource)/MultiCent/IOHandler.h
CPHSpec = $(DSource)/MultiCent/PHENIXAuAu200SpectraHandler.cxx 
HPHSpec = $(DSource)/MultiCent/PHENIXAuAu200SpectraHandler.h
CSTSpec = $(DSource)/MultiCent/STARdEdx200SpectraHandler.cxx 
HSTSpec = $(DSource)/MultiCent/STARdEdx200SpectraHandler.h
CSTV2 = $(DSource)/MultiCent/STAR200V2Handler.cxx
HSTV2 = $(DSource)/MultiCent/STAR200V2Handler.h
CSTV2V4 = $(DSource)/MultiCent/STAR200V2V4Handler.cxx
HSTV2V4 = $(DSource)/MultiCent/STAR200V2V4Handler.h
CPKV2V4 = $(DSource)/MultiCent/PKolbSpectraV2V4Handler.cxx
HPKV2V4 = $(DSource)/MultiCent/PKolbSpectraV2V4Handler.h
CSTV2Spec = $(DSource)/MultiCent/STAR200V2SpectraHandler.cxx
HSTV2Spec = $(DSource)/MultiCent/STAR200V2SpectraHandler.h
CSTXi = $(DSource)/MultiCent/STAR200XiHandler.cxx
HSTXi = $(DSource)/MultiCent/STAR200XiHandler.h
CPHV2Spec = $(DSource)/MultiCent/PHENIX200SpectraV2Handler.cxx 
HPHV2Spec = $(DSource)/MultiCent/PHENIX200SpectraV2Handler.h
CSTOm = $(DSource)/MultiCent/STAR200OmSpectraHandler.cxx
HSTOm = $(DSource)/MultiCent/STAR200OmSpectraHandler.h
CSTPhi = $(DSource)/MultiCent/STAR200PhiSpectraHandler.cxx
HSTPhi = $(DSource)/MultiCent/STAR200PhiSpectraHandler.h
CSTasHBT = $(DSource)/MultiCent/STAR200asHBTHandler.cxx
HSTasHBT = $(DSource)/MultiCent/STAR200asHBTHandler.h
CNASpec = $(DSource)/MultiCent/NA49SpectraHandler.cxx 
HNASpec = $(DSource)/MultiCent/NA49SpectraHandler.h
CSTasHBTSpec = $(DSource)/MultiCent/STAR200asHBTSpectraHandler.cxx
HSTasHBTSpec = $(DSource)/MultiCent/STAR200asHBTSpectraHandler.h
CSTasHBTV2Spec = $(DSource)/MultiCent/STAR200asHBTV2SpectraHandler.cxx
HSTasHBTV2Spec = $(DSource)/MultiCent/STAR200asHBTV2SpectraHandler.h
CPHV2HBTSpec = $(DSource)/MultiCent/PHENIX200SpectraV2HBTHandler.cxx 
HPHV2HBTSpec = $(DSource)/MultiCent/PHENIX200SpectraV2HBTHandler.h
CSTHBTSpec = $(DSource)/MultiCent/STAR200HBTSpectraHandler.cxx
HSTHBTSpec = $(DSource)/MultiCent/STAR200HBTSpectraHandler.h
CPHHBTSpec = $(DSource)/MultiCent/PHENIX200SpectraHBTHandler.cxx 
HPHHBTSpec = $(DSource)/MultiCent/PHENIX200SpectraHBTHandler.h
CSTHBT = $(DSource)/MultiCent/STAR200HBTHandler.cxx
HSTHBT = $(DSource)/MultiCent/STAR200HBTHandler.h
HDT = $(DSource)/MultiCent/DataTestHandler.h
CDT = $(DSource)/MultiCent/DataTestHandler.cxx
HHbt = $(DSource)/MultiCent/HbtPRCDataHandler.h
CHbt = $(DSource)/MultiCent/HbtPRCDataHandler.cxx
CSTPXHBTSpec = $(DSource)/MultiCent/STARKplusKplusHBTSpectraHandler.cxx
HSTPXHBTSpec = $(DSource)/MultiCent/STARKplusKplusHBTSpectraHandler.h

# Object

OFBWP = $(DObj)/FitBWPaper.o
OCBW = $(DObj)/CalcBWPaper.o
OFBW = $(DObj)/FitBW.o
OHBTBW = $(DObj)/FitHbtRPBW.o
OBWF = $(DObj)/BlastWaveFitter.o
OBW = $(DObj)/BlastWave.o
ODFF = $(DObj)/DataForFit.o
OXiBW = $(DObj)/FitXiBW.o
OMulti = $(DObj)/MultiCentFitBW.o
OIO = $(DObj)/IOHandler.o
OPHSpec = $(DObj)/PHENIXAuAu200SpectraHandler.o
OSTSpec = $(DObj)/STARdEdx200SpectraHandler.o
OSTV2 = $(DObj)/STAR200V2Handler.o
OSTV2V4 = $(DObj)/STAR200V2V4Handler.o
OPKV2V4 = $(DObj)/PKolbSpectraV2V4Handler.o
OSTV2Spec = $(DObj)/STAR200V2SpectraHandler.o
OSTXi = $(DObj)/STAR200XiHandler.o
OPHV2Spec = $(DObj)/PHENIX200SpectraV2Handler.o
OPHV2HBTSpec = $(DObj)/PHENIX200SpectraV2HBTHandler.o
OPHHBTSpec = $(DObj)/PHENIX200SpectraHBTHandler.o
OSTOm = $(DObj)/STAR200OmSpectraHandler.o
OSTPhi = $(DObj)/STAR200PhiSpectraHandler.o
OSTasHBT = $(DObj)/STAR200asHBTHandler.o
ONASpec = $(DObj)/NA49SpectraHandler.o
OSTasHBTSpec = $(DObj)/STAR200asHBTSpectraHandler.o
OSTasHBTV2Spec = $(DObj)/STAR200asHBTV2SpectraHandler.o
OSTHBTSpec = $(DObj)/STAR200HBTSpectraHandler.o
OSTPXHBTSpec = $(DObj)/STARKplusKplusHBTSpectraHandler.o
OSTHBT = $(DObj)/STAR200HBTHandler.o
ODT = $(DObj)/DataTestHandler.o
OHbt = $(DObj)/HbtPRCDataHandler.o

# Bin

EFBWP = $(DBin)/FitBWPaper.exe
ECBW = $(DBin)/CalcBWPaper.exe
EFBW = $(DBin)/FitBW.exe
EHBTBW = $(DBin)/FitHbtRPBW.exe
EXiBW = $(DBin)/FitXiBW.exe
EMulti = $(DBin)/MultiCentFitBW.exe
#
# _________________________________________________________
# --- OPTION
#

#COPT = -c $(OptFlag) -I $(DSource)/Common -I$(DSource)/MultiCent -I ${ROOTSYS}/include 
COPT = -c $(OptFlag) -I $(DSource)/Common -I$(DSource)/MultiCent $(shell root-config --cflags)

#LLIB = -LC:\cygwin\root/lib -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lm
LLIB = $(shell root-config --libs) -lMinuit
LOPT = $(OptFlag) $(shell root-config --ldflags)

#
# _________________________________________________________
# --- Ruling
#

paper : $(EFBWP)

fit : $(EFBW)

calc: $(ECBW)

hbt : $(EHBTBW)

xi : $(EXiBW)

multi : $(EMulti)

clean : 
	rm $(OFBW) $(OBWF) $(OBW) $(ODFF)

$(EFBWP) : $(DBin) $(DObj) $(OFBWP) $(OBWF) $(OBW) $(ODFF) 
	g++ $(LOPT) -o $(EFBWP) $(OFBWP) $(OBWF) $(OBW) $(ODFF) $(LLIB)

$(EFBW) : $(DBin) $(DObj) $(OFBW) $(OBWF) $(OBW) $(ODFF) 
	g++ $(LOPT) -o $(EFBW) $(OFBW) $(OBWF) $(OBW) $(ODFF) $(LLIB)

$(ECBW) : $(DBin) $(DObj) $(OCBW) $(OBWF) $(OBW) $(ODFF) 
	g++ $(LOPT) -o $(ECBW) $(OCBW) $(OBWF) $(OBW) $(ODFF) $(LLIB)	

$(EHBTBW) : $(DBin) $(DObj) $(OHBTBW) $(OBWF) $(OBW) $(ODFF) 
	g++ $(LOPT) -o $(EHBTBW) $(OHBTBW) $(OBWF) $(OBW) $(ODFF) $(LLIB)

$(EXiBW) : $(DBin) $(DObj) $(OXiBW) $(OBWF) $(OBW) $(ODFF) 
	g++ $(LOPT) -o $(EXiBW) $(OXiBW) $(OBWF) $(OBW) $(ODFF) $(LLIB)	

$(EMulti) : $(DBin) $(DObj) $(OMulti) $(OBWF) $(OBW) $(ODFF) $(OIO) $(OPHSpec) $(OSTSpec) $(OSTV2) $(OSTV2Spec) $(OSTV2V4) $(OPKV2V4) $(OSTXi) $(OPHV2Spec) $(OPHV2HBTSpec) $(OPHHBTSpec) $(OSTOm) $(OSTasHBT) $(OSTasHBTSpec) $(OSTPhi) $(ONASpec) $(OSTasHBTV2Spec) $(OSTHBTSpec) $(OSTHBT) $(ODT) $(OHbt) $(OSTPXHBTSpec)
	g++ $(LOPT) -o $(EMulti) $(OMulti) $(OBWF) $(OBW) $(ODFF) $(OIO) $(OPHSpec) $(OSTSpec) $(OSTV2) $(OSTV2Spec) $(OSTV2V4) $(OPKV2V4) $(OSTXi) $(OPHV2Spec) $(OPHV2HBTSpec) $(OPHHBTSpec) $(OSTOm) $(OSTPhi) $(OSTasHBT) $(ONASpec) $(OSTasHBTSpec) $(OSTasHBTV2Spec) $(OSTHBTSpec)  $(OSTPXHBTSpec) $(OSTHBT)  $(ODT) $(OHbt) $(LLIB) 	

$(OFBWP) : $(CFBWP) $(HFBW) $(HBWF)
	g++ $(COPT) -o $(OFBWP) $(CFBWP)

$(OFBW) : $(CFBW) $(HFBW) $(HBWF)
	g++ $(COPT) -o $(OFBW) $(CFBW)

$(OCBW) : $(CCBW) $(HFBW) $(HBWF)
	g++ $(COPT) -o $(OCBW) $(CCBW)	

$(OHBTBW) : $(CHBTBW) $(HFBW) $(HBWF)
	g++ $(COPT) -o $(OHBTBW) $(CHBTBW)		

$(OBWF) : $(CBWF) $(HBWF) $(HBW) $(HDFF)
	g++ $(COPT) -o $(OBWF) $(CBWF)

$(OBW) : $(CBW) $(HBW)
	g++ $(COPT) -o $(OBW) $(CBW)

$(ODFF) : $(CDFF) $(HDFF) $(HBW)
	g++ $(COPT) -o $(ODFF) $(CDFF)

$(DBin) : 
	tcsh -c 'if(! -d $(DBin)) mkdir $(DBin)'

$(DObj) : 
	tcsh -c 'if(! -d $(DDObj)) mkdir $(DDObj)'
	tcsh -c 'if(! -d $(DObj)) mkdir $(DObj)'

$(OXiBW) : $(CXiBW) $(HFBW) $(HBWF)
	g++ $(COPT) -o $(OXiBW) $(CXiBW)

$(OMulti) : $(CMulti)  $(HFBW) $(HBWF) $(HIO) $(HPhSpec) $(HSTSpec) $(HSTV2) $(HSTV2V4) $(HPKV2V4) $(HSTV2Spec) $(HStXi) $(HPHV2Spec) $(HPHV2HBTSpec) $(HPHHBTSpec) $(HSTOm) $(HSTasHBT) $(HSTasHBTSpec) $(HSTHBTSpec) $(HSTPhi) $(HNASpec) $(HSTasHBTV2Spec)  $(HSTHBT) $(HDT) $(HHbt) $(HSTPXHBTSpec)
	g++ $(COPT) -o $(OMulti) $(CMulti)

$(OIO) : $(CIO) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OIO) $(CIO)

$(OPHSpec) : $(CPHSpec) $(HPHSpec) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OPHSpec) $(CPHSpec)

$(OSTSpec) : $(CSTSpec) $(HSTSpec) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTSpec) $(CSTSpec)

$(OSTV2) : $(CSTV2) $(HSTV2) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTV2) $(CSTV2)

$(OSTV2V4) : $(CSTV2V4) $(HSTV2V4) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTV2V4) $(CSTV2V4)

$(OPKV2V4) : $(CPKV2V4) $(HPKV2V4) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OPKV2V4) $(CPKV2V4)

$(OSTV2Spec) : $(CSTV2Spec) $(HSTV2Spec) $(HIO) $(HSTV2) $(HSTSpec) $(HBWF)
	g++ $(COPT) -o $(OSTV2Spec) $(CSTV2Spec)

$(OSTXi) : $(CSTXi) $(HSTXi) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTXi) $(CSTXi)

$(OPHV2Spec) : $(CPHV2Spec) $(HPHV2Spec) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OPHV2Spec) $(CPHV2Spec)

$(OPHV2HBTSpec) : $(CPHV2HBTSpec) $(HPHV2HBTSpec) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OPHV2HBTSpec) $(CPHV2HBTSpec)

$(OSTOm) : $(CSTOm) $(HSTOm) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTOm) $(CSTOm)

$(OSTPhi) : $(CSTPhi) $(HSTPhi) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTPhi) $(CSTPhi)

$(OSTasHBT) : $(CSTasHBT) $(HSTasHBT) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTasHBT) $(CSTasHBT)

$(ONASpec) : $(CNASpec) $(HNASpec) $(HIO) $(HBWF)
	g++ $(COPT) -o $(ONASpec) $(CNASpec)

$(OSTasHBTSpec) : $(CSTasHBTSpec) $(HSTasHBTSpec) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTasHBTSpec) $(CSTasHBTSpec)

$(OSTasHBTV2Spec) : $(CSTasHBTV2Spec) $(HSTasHBTV2Spec) $(HSTasHBTSpec) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTasHBTV2Spec) $(CSTasHBTV2Spec)

$(OSTHBTSpec) : $(CSTHBTSpec) $(HSTHBTSpec) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTHBTSpec) $(CSTHBTSpec)

$(OSTPXHBTSpec) : $(CSTPXHBTSpec) $(HSTPXHBTSpec) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTPXHBTSpec) $(CSTPXHBTSpec)

$(OPHHBTSpec) : $(CPHHBTSpec) $(HPHHBTSpec) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OPHHBTSpec) $(CPHHBTSpec)

$(OSTHBT) : $(CSTHBT) $(HSTHBT) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OSTHBT) $(CSTHBT)

$(ODT) : $(CDT) $(HDT) $(HIO) $(HBWF)
	g++ $(COPT) -o $(ODT) $(CDT)

$(OHbt) : $(CHbt) $(HHbt) $(HIO) $(HBWF)
	g++ $(COPT) -o $(OHbt) $(CHbt)









