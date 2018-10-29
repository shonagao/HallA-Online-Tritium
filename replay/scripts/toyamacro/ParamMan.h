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
using namespace std;
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
    double tdcOffset[nS2];//LR,TB
    double tdcGain[nS2];
    double timeWalk[nS2];
  };

  struct S0param{
    int cid;
    int lr;
    int tb;
    double tdcOffset[nS0];
    double tdcGain[nS0];
    double timeWalk[nS0];
  };

  struct RFparam{
    int cid;
    int lr;
    int tb;
    double tdcOffset[nRF];//LR
    double tdcGain[nRF];
  };


  S2param L_F1S2T, L_F1S2B, R_F1S2T, R_F1S2B;
  S0param L_F1S0T, L_F1S0B, R_F1S0T, R_F1S0B;
  RFparam L_RF   ,          R_RF;
  S2param L_FbS2T, L_FbS2B, R_FbS2T, R_FbS2B;
  S0param L_FbS0T, L_FbS0B, R_FbS0T, R_FbS0B;

public:
  void SetFileName( const char* filename )
  { ParamFileName=filename; } 
  void SetAParam( string &sline, string &name, double *par, int &npar );
  void SetTdcOffset( int cid, int seg, int lr, int tb, double tdcOffset );
  void SetTdcGain(   int cid, int seg, int lr, int tb, double tdcGain   );
  void SetTimeTune(  int cid, int seg, int lr, int tb, double tdcGain   );
  void SetTimeTune(  int cid, int seg, int lr, double tdcGain   );
  
  bool SetVal( void );
  double GetTdcOffset( int cid, int seg, int lr, int tb );
  double GetTdcGain(   int cid, int seg, int lr, int tb );
  double time(         int cid, int seg, int lr, int tb, double tdc, double adc=-1.);
  double TimewalkCorrection(double par, double adc, double peak=2000.);
  void WriteToFile( const char* OutputFileName );
};


#endif
