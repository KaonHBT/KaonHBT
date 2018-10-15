#ifndef STAR200V2SpectraHandler_h
#define STAR200V2SpectraHandler_h
#include "STAR200V2Handler.h"
class STARdEdx200SpectraHandler;

class STAR200V2SpectraHandler : public STAR200V2Handler{
 public:
  STAR200V2SpectraHandler(BlastWaveFitter* aBWFitter, 
			      const char* aInputFileName=0);
  ~STAR200V2SpectraHandler(){}
  void writeSpectra();
 private:
  STARdEdx200SpectraHandler* mSpectraHandler;
  void initLoad();
  void load();
};

#endif
