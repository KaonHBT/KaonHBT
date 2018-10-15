#ifndef STARdEdx200SpectraHandler_h
#define STARdEdx200SpectraHandler_h
#include "IOHandler.h"

class STARdEdx200SpectraHandler : public IOHandler{
 public:
  STARdEdx200SpectraHandler(BlastWaveFitter* aBWFitter, 
			      const char* aInputFileName=0);
  ~STARdEdx200SpectraHandler(){}

  virtual void initLoad();
  virtual void load();
};

#endif
