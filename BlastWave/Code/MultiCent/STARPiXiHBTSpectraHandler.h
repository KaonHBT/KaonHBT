#ifndef STARPiXiHBTSpectraHandler_h
#define STARPiXiHBTSpectraHandler_h
//#include "STARPiXiasHBTSpectraHandler.h"
#include "IOHandler.h"
class TFile;

class STARPiXiHBTSpectraHandler : public IOHandler{//STARPiXiasHBTSpectraHandler{
 public:
  STARPiXiHBTSpectraHandler(BlastWaveFitter* aBWFitter, 
			       const char* aInputFileName=0);
  ~STARPiXiHBTSpectraHandler(){}
 protected:
  virtual void initLoad();
  virtual void load();
  //virtual void fixParameters();
  virtual void loadSpectra(int iPart);
  virtual void close();
  TFile* mSpectraInputFile;
  //virtual void fixParameters();
};

#endif
