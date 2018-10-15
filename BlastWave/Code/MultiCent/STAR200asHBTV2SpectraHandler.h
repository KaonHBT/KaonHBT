#ifndef STAR200asHBTV2SpectraHandler_h
#define STAR200asHBTV2SpectraHandler_h
//#include "STAR200asHBTSpectraHandler.h"
#include "IOHandler.h"
class TFile;

class STAR200asHBTV2SpectraHandler : public IOHandler{//STAR200asHBTSpectraHandler{
 public:
  STAR200asHBTV2SpectraHandler(BlastWaveFitter* aBWFitter, 
			       const char* aInputFileName=0,
			       int aFourier=0, // Fit Fourier coef
			       int aV2=1 // fit v2 as well
			       );
  ~STAR200asHBTV2SpectraHandler(){}
 protected:
  virtual void initLoad();
  virtual void load();
  virtual void fixParameters();
  virtual void loadSpectra(int iPart);
  virtual void close();
  TFile* mSpectraInputFile;
  TFile* mV2InputFile;
  int mFourier;
};

#endif
