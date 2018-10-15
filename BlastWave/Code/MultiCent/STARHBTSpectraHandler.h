#ifndef STAR200HBTSpectraHandler_h
#define STAR200HBTSpectraHandler_h
//#include "STAR200asHBTSpectraHandler.h"
#include "IOHandler.h"
class TFile;

class STAR200HBTSpectraHandler : public IOHandler{//STAR200asHBTSpectraHandler{
 public:
  STAR200HBTSpectraHandler(BlastWaveFitter* aBWFitter, 
			       const char* aInputFileName=0);
  ~STAR200HBTSpectraHandler(){}
 protected:
  virtual void initLoad();
  virtual void load();
  //virtual void fixParameters();
  virtual void loadSpectra(int iPart);
  virtual void close();
  TFile* mSpectraInputFile;
};

#endif
