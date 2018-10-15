#ifndef STAR200V2Handler_h
#define STAR200V2Handler_h
#include "IOHandler.h"

class STAR200V2Handler : public IOHandler{
 public:
  STAR200V2Handler(BlastWaveFitter* aBWFitter, 
			      const char* aInputFileName=0);
  ~STAR200V2Handler(){}
 protected:
  virtual void initLoad();
  virtual void load();
};

#endif
