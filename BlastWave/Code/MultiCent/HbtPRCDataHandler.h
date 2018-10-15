#ifndef HbtPRCDataHandler_h
#define HbtPRCDataHandler_h
//#include "STAR200asHBTSpectraHandler.h"
#include "IOHandler.h"
class TFile;

class HbtPRCDataHandler : public IOHandler{//STAR200asHBTSpectraHandler{
 public:
  HbtPRCDataHandler(BlastWaveFitter* aBWFitter, 
		    const char* aInputFileName=0, 
		    int aAsHbt=0, // 0=no, 1=Fourier, 2=Broken 
		    int aV2=0,
		    int aSpectraFit=1);
  ~HbtPRCDataHandler(){}
 protected:
  virtual void initLoad();
  virtual void load();
  virtual void fixParameters();
  virtual void loadSpectra(int iPart);
  virtual void close();
  TFile* mSpectraInputFile;
  TFile* mV2InputFile;
  int mAsHbt;
  int mSpectraFit;
  double* mFixedT;
  double* mFixedRho0;
};

#endif
