#ifndef STAR200asHBTSpectraHandler_h
#define STAR200asHBTSpectraHandler_h
#include "IOHandler.h"
class TFile;

class STAR200asHBTSpectraHandler : public IOHandler{
 public:
  STAR200asHBTSpectraHandler(BlastWaveFitter* aBWFitter, 
		   const char* aInputFileName=0);
  ~STAR200asHBTSpectraHandler(){}
 protected:
  virtual void initLoad();
  virtual void load();
  virtual void fixParameters();
  virtual void loadSpectra(int iPart);
  virtual void close();
  TFile* mSpectraInputFile;
};

#endif
