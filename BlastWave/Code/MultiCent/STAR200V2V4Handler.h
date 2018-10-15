#ifndef STAR200V2V4Handler_h
#define STAR200V2V4Handler_h
#include "IOHandler.h"
class TFile;

class STAR200V2V4Handler : public IOHandler{
 public:
  STAR200V2V4Handler(BlastWaveFitter* aBWFitter,
		     int aV4=1,
		     const char* aInputFileName=0);
  ~STAR200V2V4Handler(){}
 protected:
  TFile* mSpectraInputFile;
  virtual void fixParameters();
  virtual void initLoad();
  virtual void load();
  int mV4;
};

#endif
