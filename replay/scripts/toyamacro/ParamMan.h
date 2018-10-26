/*
  ParamMan.h
*/

#ifndef ParamMan_hh
#define ParamMan_hh 1

#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include "define.h"
//

class ParamMan
{
private:
//   std::string ParamFileName;
  const char* ParamFileName;
  char* OutputFileName;

public:
//   explicit ParamMan( const std::string & filename );
  ParamMan();
  explicit ParamMan( const char* filename);
  ~ParamMan(){}

private:
  struct S2param{
    int cid;
    int lr;
    int tb;
    double tdcOffset[nS2*4];//LR,TB
    double tdcGain[nS2*4];
  };

  struct S0param{
    int cid;
    int lr;
    int tb;
    double tdcOffset[nS0*4];
    double tdcGain[nS0*4];
  };

  struct RFparam{
    int cid;
    int lr;
    int tb;
    double tdcOffset[nRF*2];//LR
    double tdcGain[nRF*2];
  };


  S2param F1S2;
  S0param F1S0;
  RFparam RF;
  S2param FbS2;
  S0param FbS0;

public:
  void SetFileName( const char* filename )
  { ParamFileName=filename; } 
  void SetTdcOffset( int cid, int seg, int lr, int tb, double tdcOffset );
  void SetTdcGain(   int cid, int seg, int lr, int tb, double tdcGain   );
  void SetTimeTune(  int cid, int seg, int lr, int tb, double tdcGain   );
  void SetTimeTune(  int cid, int seg, int lr, double tdcGain   );
  
  bool SetVal( void );
  double GetTdcOffset( int cid, int seg, int lr, int tb );
  double GetTdcGain(   int cid, int seg, int lr, int tb );
  double time(         int cid, int seg, int lr, int tb, double tdc);
  void WriteToFile( const char* OutputFileName );
};


#endif
