/*
  ParamMan.cc
*/


#include "ParamMan.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <map>

ParamMan::ParamMan( const char* filename )
  : ParamFileName(filename)
{
  L_F1S2T.cid=CID_F1S2;L_F1S2B.cid=CID_F1S2;
  R_F1S2T.cid=CID_F1S2;R_F1S2B.cid=CID_F1S2;
  L_F1S2T.lr=0;L_F1S2B.lr=0;
  R_F1S2T.lr=1;R_F1S2B.lr=1;
  L_F1S2T.tb=0;L_F1S2B.tb=1;
  R_F1S2T.tb=0;R_F1S2B.tb=1;
  for(int i=0; i<nS2; i++){
    L_F1S2T.tdcOffset[i]=0.00; L_F1S2T.timeWalk[i]=0.00;   L_F1S2T.tdcGain[i]=0.05623;
    L_F1S2B.tdcOffset[i]=0.00; L_F1S2B.timeWalk[i]=0.00;   L_F1S2B.tdcGain[i]=0.05623;
    R_F1S2T.tdcOffset[i]=0.00; R_F1S2T.timeWalk[i]=0.00;   R_F1S2T.tdcGain[i]=0.05623;
    R_F1S2B.tdcOffset[i]=0.00; R_F1S2B.timeWalk[i]=0.00;   R_F1S2B.tdcGain[i]=0.05623;
  }

  L_F1S0T.cid=CID_F1S0;L_F1S0B.cid=CID_F1S0;
  R_F1S0T.cid=CID_F1S0;R_F1S0B.cid=CID_F1S0;
  L_F1S0T.lr=0;L_F1S0B.lr=0;
  R_F1S0T.lr=1;R_F1S0B.lr=1;
  L_F1S0T.tb=0;L_F1S0B.tb=1;
  R_F1S0T.tb=0;R_F1S0B.tb=1;
  for(int i=0; i<nS0; i++){
    L_F1S0T.tdcOffset[i]=0.00; L_F1S0T.timeWalk[i]=0.00;   L_F1S0T.tdcGain[i]=0.05623;
    L_F1S0B.tdcOffset[i]=0.00; L_F1S0B.timeWalk[i]=0.00;   L_F1S0B.tdcGain[i]=0.05623;
    R_F1S0T.tdcOffset[i]=0.00; R_F1S0T.timeWalk[i]=0.00;   R_F1S0T.tdcGain[i]=0.05623;
    R_F1S0B.tdcOffset[i]=0.00; R_F1S0B.timeWalk[i]=0.00;   R_F1S0B.tdcGain[i]=0.05623;
  }

  L_RF.cid=CID_RF;
  R_RF.cid=CID_RF;
  L_RF.lr=0;
  R_RF.lr=1;
  for(int i=0; i<nRF; i++){
    L_RF.tdcOffset[i]=0.00;   L_RF.tdcGain[i]=0.05623;
    R_RF.tdcOffset[i]=0.00;   R_RF.tdcGain[i]=0.05623;
  }

  L_FbS2T.cid=CID_FbS2;L_FbS2B.cid=CID_FbS2;
  R_FbS2T.cid=CID_FbS2;R_FbS2B.cid=CID_FbS2;
  L_FbS2T.lr=0;L_FbS2B.lr=0;
  R_FbS2T.lr=1;R_FbS2B.lr=1;
  L_FbS2T.tb=0;L_FbS2B.tb=1;
  R_FbS2T.tb=0;R_FbS2B.tb=1;
  for(int i=0; i<nS2; i++){
    L_FbS2T.tdcOffset[i]=0.00; L_FbS2T.timeWalk[i]=0.00;  L_FbS2T.tdcGain[i]=0.5;
    L_FbS2B.tdcOffset[i]=0.00; L_FbS2B.timeWalk[i]=0.00;  L_FbS2B.tdcGain[i]=0.5;
    R_FbS2T.tdcOffset[i]=0.00; R_FbS2T.timeWalk[i]=0.00;  R_FbS2T.tdcGain[i]=0.5;
    R_FbS2B.tdcOffset[i]=0.00; R_FbS2B.timeWalk[i]=0.00;  R_FbS2B.tdcGain[i]=0.5;
  }

  L_FbS0T.cid=CID_FbS0;L_FbS0B.cid=CID_FbS0;
  R_FbS0T.cid=CID_FbS0;R_FbS0B.cid=CID_FbS0;
  L_FbS0T.lr=0;L_FbS0B.lr=0;
  R_FbS0T.lr=1;R_FbS0B.lr=1;
  L_FbS0T.tb=0;L_FbS0B.tb=1;
  R_FbS0T.tb=0;R_FbS0B.tb=1;
  for(int i=0; i<nS0; i++){
    L_FbS0T.tdcOffset[i]=0.00; L_FbS0T.timeWalk[i]=0.00;   L_FbS0T.tdcGain[i]=0.5;
    L_FbS0B.tdcOffset[i]=0.00; L_FbS0B.timeWalk[i]=0.00;   L_FbS0B.tdcGain[i]=0.5;
    R_FbS0T.tdcOffset[i]=0.00; R_FbS0T.timeWalk[i]=0.00;   R_FbS0T.tdcGain[i]=0.5;
    R_FbS0B.tdcOffset[i]=0.00; R_FbS0B.timeWalk[i]=0.00;   R_FbS0B.tdcGain[i]=0.5;
  }

}
////////////////////////////////
//ParamMan::~ParamMan()
//{
//}
//____________________________________________________________________________________________
bool ParamMan::SetVal( void )
{
  static const std::string funcname = "ParamMan::SetVal";
  
  ifstream ifs(ParamFileName);
  if ( ifs.fail() ) {
    cout << "file open fail : " << ParamFileName << endl;
    return false;
  }

  string line;
  string name;
  int npar;
  while(1){
    getline(ifs,line);
    if(line[0]=='-') {cout<<line<<endl;continue;}
    if(line[0]=='#') continue;
    if( ifs.eof() ) break;
    name = "R.s2.R.off =" ;             npar=16; SetAParam(line,name,R_FbS2T.tdcOffset,npar);
    name = "R.s2.L.off =" ;             npar=16; SetAParam(line,name,R_FbS2B.tdcOffset,npar);
    name = "R.s0.R.off =" ;             npar=1;  SetAParam(line,name,R_FbS0T.tdcOffset,npar);
    name = "R.s0.L.off =" ;             npar=1;  SetAParam(line,name,R_FbS0B.tdcOffset,npar);
    name = "R.s2.R.off_F1 =" ;          npar=16; SetAParam(line,name,R_F1S2T.tdcOffset,npar);//tmp
    name = "R.s2.L.off_F1 =" ;          npar=16; SetAParam(line,name,R_F1S2B.tdcOffset,npar);//tmp
    name = "R.s0.R.off_F1 =" ;          npar=1;  SetAParam(line,name,R_F1S0T.tdcOffset,npar);//tmp
    name = "R.s0.L.off_F1 =" ;          npar=1;  SetAParam(line,name,R_F1S0B.tdcOffset,npar);//tmp
    name = "R.RF.off =" ;               npar=1;  SetAParam(line,name,R_RF.tdcOffset,npar);//tmp
    name = "L.s2.R.off =" ;             npar=16; SetAParam(line,name,L_FbS2T.tdcOffset,npar);
    name = "L.s2.L.off =" ;             npar=16; SetAParam(line,name,L_FbS2B.tdcOffset,npar);
    name = "L.s0.R.off =" ;             npar=1;  SetAParam(line,name,L_FbS0T.tdcOffset,npar);
    name = "L.s0.L.off =" ;             npar=1;  SetAParam(line,name,L_FbS0B.tdcOffset,npar);
    name = "L.s2.R.off_F1 =" ;          npar=16; SetAParam(line,name,L_F1S2T.tdcOffset,npar);//tmp
    name = "L.s2.L.off_F1 =" ;          npar=16; SetAParam(line,name,L_F1S2B.tdcOffset,npar);//tmp
    name = "L.s0.R.off_F1 =" ;          npar=1;  SetAParam(line,name,L_F1S0T.tdcOffset,npar);//tmp
    name = "L.s0.L.off_F1 =" ;          npar=1;  SetAParam(line,name,L_F1S0B.tdcOffset,npar);//tmp
    name = "L.RF.off =" ;               npar=1;  SetAParam(line,name,L_RF.tdcOffset,npar);//tmp
    
    //time walk correction parameters style is different from replay DB. Sorry.
    name = "R.s2.R.timewalk_params =" ;    npar=16; SetAParam(line,name,R_FbS2T.timeWalk ,npar);
    name = "R.s2.L.timewalk_params =" ;    npar=16; SetAParam(line,name,R_FbS2B.timeWalk ,npar);
    name = "R.s0.R.timewalk_params =" ;    npar=1;  SetAParam(line,name,R_FbS0T.timeWalk ,npar);
    name = "R.s0.L.timewalk_params =" ;    npar=1;  SetAParam(line,name,R_FbS0B.timeWalk ,npar);
    name = "R.s2.R.timewalk_params_F1 =" ; npar=16; SetAParam(line,name,R_F1S2T.timeWalk,npar);//tmp
    name = "R.s2.L.timewalk_params_F1 =" ; npar=16; SetAParam(line,name,R_F1S2B.timeWalk,npar);//tmp
    name = "R.s0.R.timewalk_params_F1 =" ; npar=1;  SetAParam(line,name,R_F1S0T.timeWalk,npar);//tmp
    name = "R.s0.L.timewalk_params_F1 =" ; npar=1;  SetAParam(line,name,R_F1S0B.timeWalk,npar);//tmp
    name = "L.s2.R.timewalk_params =" ;    npar=16; SetAParam(line,name,L_FbS2T.timeWalk,npar);
    name = "L.s2.L.timewalk_params =" ;    npar=16; SetAParam(line,name,L_FbS2B.timeWalk,npar);
    name = "L.s0.R.timewalk_params =" ;    npar=1;  SetAParam(line,name,L_FbS0T.timeWalk,npar);
    name = "L.s0.L.timewalk_params =" ;    npar=1;  SetAParam(line,name,L_FbS0B.timeWalk,npar);
    name = "L.s2.R.timewalk_params_F1 =" ; npar=16; SetAParam(line,name,L_F1S2T.timeWalk,npar);//tmp
    name = "L.s2.L.timewalk_params_F1 =" ; npar=16; SetAParam(line,name,L_F1S2B.timeWalk,npar);//tmp
    name = "L.s0.R.timewalk_params_F1 =" ; npar=1;  SetAParam(line,name,L_F1S0T.timeWalk,npar);//tmp
    name = "L.s0.L.timewalk_params_F1 =" ; npar=1;  SetAParam(line,name,L_F1S0B.timeWalk,npar);//tmp
  }
  
//   std::cout << "[" << funcname << "]: Initialization finished" << std::endl;
  return true;
}

///////////////////////////////////
void ParamMan::SetAParam( string &line, string &name, double *par, int &npar )
{
  if ( line.compare(0,name.size(),name) != 0 ) {
   return;
  }
    cout<<line<<endl;
  istringstream sline(line);
  sline >> name;
  sline >> name;
  for ( int i=0; i<npar; i++ ) {
    sline >> par[i];
    cout<<i<<" "<<par[i]<<endl;
  }
}
//____________________________________________________________________________________________
double ParamMan::time( int cid, int seg, int lr, int tb, double tdc , double adc)
{
  static const std::string funcname = "ParamMan::time";

  if     (cid==L_F1S2T.cid && lr==L_F1S2T.lr && tb==L_F1S2T.tb) return L_F1S2T.tdcGain[seg]*(tdc-L_F1S2T.tdcOffset[seg]) - TimewalkCorrection(L_F1S2T.timeWalk[seg], adc);
  else if(cid==L_F1S2B.cid && lr==L_F1S2B.lr && tb==L_F1S2B.tb) return L_F1S2B.tdcGain[seg]*(tdc-L_F1S2B.tdcOffset[seg]) - TimewalkCorrection(L_F1S2B.timeWalk[seg], adc);
  else if(cid==R_F1S2T.cid && lr==R_F1S2T.lr && tb==R_F1S2T.tb) return R_F1S2T.tdcGain[seg]*(tdc-R_F1S2T.tdcOffset[seg]) - TimewalkCorrection(R_F1S2T.timeWalk[seg], adc);
  else if(cid==R_F1S2B.cid && lr==R_F1S2B.lr && tb==R_F1S2B.tb) return R_F1S2B.tdcGain[seg]*(tdc-R_F1S2B.tdcOffset[seg]) - TimewalkCorrection(R_F1S2B.timeWalk[seg], adc);
  else if(cid==L_F1S0T.cid && lr==L_F1S0T.lr && tb==L_F1S0T.tb) return L_F1S0T.tdcGain[seg]*(tdc-L_F1S0T.tdcOffset[seg]) - TimewalkCorrection(L_F1S0T.timeWalk[seg], adc);
  else if(cid==L_F1S0B.cid && lr==L_F1S0B.lr && tb==L_F1S0B.tb) return L_F1S0B.tdcGain[seg]*(tdc-L_F1S0B.tdcOffset[seg]) - TimewalkCorrection(L_F1S0B.timeWalk[seg], adc);
  else if(cid==R_F1S0T.cid && lr==R_F1S0T.lr && tb==R_F1S0T.tb) return R_F1S0T.tdcGain[seg]*(tdc-R_F1S0T.tdcOffset[seg]) - TimewalkCorrection(R_F1S0T.timeWalk[seg], adc);
  else if(cid==R_F1S0B.cid && lr==R_F1S0B.lr && tb==R_F1S0B.tb) return R_F1S0B.tdcGain[seg]*(tdc-R_F1S0B.tdcOffset[seg]) - TimewalkCorrection(R_F1S0B.timeWalk[seg], adc);
  else if(cid==L_RF.cid    && lr==L_RF.lr)                      return L_RF.tdcGain[seg]   *(tdc-L_RF.tdcOffset[seg]);
  else if(cid==R_RF.cid    && lr==R_RF.lr)                      return R_RF.tdcGain[seg]   *(tdc-R_RF.tdcOffset[seg]);
  else if(cid==L_FbS2T.cid && lr==L_FbS2T.lr && tb==L_FbS2T.tb) return L_FbS2T.tdcGain[seg]*(tdc-L_FbS2T.tdcOffset[seg]) - TimewalkCorrection(L_FbS2T.timeWalk[seg] , adc);
  else if(cid==L_FbS2B.cid && lr==L_FbS2B.lr && tb==L_FbS2B.tb) return L_FbS2B.tdcGain[seg]*(tdc-L_FbS2B.tdcOffset[seg]) - TimewalkCorrection(L_FbS2B.timeWalk[seg] , adc);
  else if(cid==R_FbS2T.cid && lr==R_FbS2T.lr && tb==R_FbS2T.tb) return R_FbS2T.tdcGain[seg]*(tdc-R_FbS2T.tdcOffset[seg]) - TimewalkCorrection(R_FbS2T.timeWalk[seg] , adc);
  else if(cid==R_FbS2B.cid && lr==R_FbS2B.lr && tb==R_FbS2B.tb) return R_FbS2B.tdcGain[seg]*(tdc-R_FbS2B.tdcOffset[seg]) - TimewalkCorrection(R_FbS2B.timeWalk[seg] , adc);
  else if(cid==L_FbS0T.cid && lr==L_FbS0T.lr && tb==L_FbS0T.tb) return L_FbS0T.tdcGain[seg]*(tdc-L_FbS0T.tdcOffset[seg]) - TimewalkCorrection(L_FbS0T.timeWalk[seg] , adc);
  else if(cid==L_FbS0B.cid && lr==L_FbS0B.lr && tb==L_FbS0B.tb) return L_FbS0B.tdcGain[seg]*(tdc-L_FbS0B.tdcOffset[seg]) - TimewalkCorrection(L_FbS0B.timeWalk[seg] , adc);
  else if(cid==R_FbS0T.cid && lr==R_FbS0T.lr && tb==R_FbS0T.tb) return R_FbS0T.tdcGain[seg]*(tdc-R_FbS0T.tdcOffset[seg]) - TimewalkCorrection(R_FbS0T.timeWalk[seg] , adc);
  else if(cid==R_FbS0B.cid && lr==R_FbS0B.lr && tb==R_FbS0B.tb) return R_FbS0B.tdcGain[seg]*(tdc-R_FbS0B.tdcOffset[seg]) - TimewalkCorrection(R_FbS0B.timeWalk[seg] , adc);
  else   cerr << "[" << funcname << "]: unknown id" << endl;

  return -1.;
}

///////////////////////////////////
void ParamMan::SetTdcOffset( int cid, int seg, int lr, int tb,
				 double tdcOffset )
{
  static const std::string funcname = "ParamMan::SettdcOffset";
  
  if     (cid==L_F1S2T.cid && lr==L_F1S2T.lr && tb==L_F1S2T.tb) L_F1S2T.tdcOffset[seg] = tdcOffset;
  else if(cid==L_F1S2B.cid && lr==L_F1S2B.lr && tb==L_F1S2B.tb) L_F1S2B.tdcOffset[seg] = tdcOffset;
  else if(cid==R_F1S2T.cid && lr==R_F1S2T.lr && tb==R_F1S2T.tb) R_F1S2T.tdcOffset[seg] = tdcOffset;
  else if(cid==R_F1S2B.cid && lr==R_F1S2B.lr && tb==R_F1S2B.tb) R_F1S2B.tdcOffset[seg] = tdcOffset;
  else if(cid==L_F1S0T.cid && lr==L_F1S0T.lr && tb==L_F1S0T.tb) L_F1S0T.tdcOffset[seg] = tdcOffset;
  else if(cid==L_F1S0B.cid && lr==L_F1S0B.lr && tb==L_F1S0B.tb) L_F1S0B.tdcOffset[seg] = tdcOffset;
  else if(cid==R_F1S0T.cid && lr==R_F1S0T.lr && tb==R_F1S0T.tb) R_F1S0T.tdcOffset[seg] = tdcOffset;
  else if(cid==R_F1S0B.cid && lr==R_F1S0B.lr && tb==R_F1S0B.tb) R_F1S0B.tdcOffset[seg] = tdcOffset;
  else if(cid==L_RF.cid    && lr==L_RF.lr)                      L_RF.tdcOffset[seg]    = tdcOffset;
  else if(cid==R_RF.cid    && lr==R_RF.lr)                      R_RF.tdcOffset[seg]    = tdcOffset;
  else if(cid==L_FbS2T.cid && lr==L_FbS2T.lr && tb==L_FbS2T.tb) L_FbS2T.tdcOffset[seg] = tdcOffset;
  else if(cid==L_FbS2B.cid && lr==L_FbS2B.lr && tb==L_FbS2B.tb) L_FbS2B.tdcOffset[seg] = tdcOffset;
  else if(cid==R_FbS2T.cid && lr==R_FbS2T.lr && tb==R_FbS2T.tb) R_FbS2T.tdcOffset[seg] = tdcOffset;
  else if(cid==R_FbS2B.cid && lr==R_FbS2B.lr && tb==R_FbS2B.tb) R_FbS2B.tdcOffset[seg] = tdcOffset;
  else if(cid==L_FbS0T.cid && lr==L_FbS0T.lr && tb==L_FbS0T.tb) L_FbS0T.tdcOffset[seg] = tdcOffset;
  else if(cid==L_FbS0B.cid && lr==L_FbS0B.lr && tb==L_FbS0B.tb) L_FbS0B.tdcOffset[seg] = tdcOffset;
  else if(cid==R_FbS0T.cid && lr==R_FbS0T.lr && tb==R_FbS0T.tb) R_FbS0T.tdcOffset[seg] = tdcOffset;
  else if(cid==R_FbS0B.cid && lr==R_FbS0B.lr && tb==R_FbS0B.tb) R_FbS0B.tdcOffset[seg] = tdcOffset;
  else   cerr << "[" << funcname << "]: unknown id" << endl;
}
///////////////////////////////////
void ParamMan::SetTdcGain( int cid, int seg, int lr, int tb,
			       double tdcGain )
{
  static const std::string funcname = "ParamMan::SettdcGain";
  
  if     (cid==L_F1S2T.cid && lr==L_F1S2T.lr && tb==L_F1S2T.tb) L_F1S2T.tdcGain[seg] = tdcGain;
  else if(cid==L_F1S2B.cid && lr==L_F1S2B.lr && tb==L_F1S2B.tb) L_F1S2B.tdcGain[seg] = tdcGain;
  else if(cid==R_F1S2T.cid && lr==R_F1S2T.lr && tb==R_F1S2T.tb) R_F1S2T.tdcGain[seg] = tdcGain;
  else if(cid==R_F1S2B.cid && lr==R_F1S2B.lr && tb==R_F1S2B.tb) R_F1S2B.tdcGain[seg] = tdcGain;
  else if(cid==L_F1S0T.cid && lr==L_F1S0T.lr && tb==L_F1S0T.tb) L_F1S0T.tdcGain[seg] = tdcGain;
  else if(cid==L_F1S0B.cid && lr==L_F1S0B.lr && tb==L_F1S0B.tb) L_F1S0B.tdcGain[seg] = tdcGain;
  else if(cid==R_F1S0T.cid && lr==R_F1S0T.lr && tb==R_F1S0T.tb) R_F1S0T.tdcGain[seg] = tdcGain;
  else if(cid==R_F1S0B.cid && lr==R_F1S0B.lr && tb==R_F1S0B.tb) R_F1S0B.tdcGain[seg] = tdcGain;
  else if(cid==L_RF.cid    && lr==L_RF.lr)                      L_RF.tdcGain[seg]    = tdcGain;
  else if(cid==R_RF.cid    && lr==R_RF.lr)                      R_RF.tdcGain[seg]    = tdcGain;
  else if(cid==L_FbS2T.cid && lr==L_FbS2T.lr && tb==L_FbS2T.tb) L_FbS2T.tdcGain[seg] = tdcGain;
  else if(cid==L_FbS2B.cid && lr==L_FbS2B.lr && tb==L_FbS2B.tb) L_FbS2B.tdcGain[seg] = tdcGain;
  else if(cid==R_FbS2T.cid && lr==R_FbS2T.lr && tb==R_FbS2T.tb) R_FbS2T.tdcGain[seg] = tdcGain;
  else if(cid==R_FbS2B.cid && lr==R_FbS2B.lr && tb==R_FbS2B.tb) R_FbS2B.tdcGain[seg] = tdcGain;
  else if(cid==L_FbS0T.cid && lr==L_FbS0T.lr && tb==L_FbS0T.tb) L_FbS0T.tdcGain[seg] = tdcGain;
  else if(cid==L_FbS0B.cid && lr==L_FbS0B.lr && tb==L_FbS0B.tb) L_FbS0B.tdcGain[seg] = tdcGain;
  else if(cid==R_FbS0T.cid && lr==R_FbS0T.lr && tb==R_FbS0T.tb) R_FbS0T.tdcGain[seg] = tdcGain;
  else if(cid==R_FbS0B.cid && lr==R_FbS0B.lr && tb==R_FbS0B.tb) R_FbS0B.tdcGain[seg] = tdcGain;
  else   cerr << "[" << funcname << "]: unknown id" << endl;
}
///////////////////////////////////
void ParamMan::SetTimeTune( int cid, int seg, int lr, int tb,
				double time )
{
  static const std::string funcname = "ParamMan::SetNpeTune";
  
  if     (cid==L_F1S2T.cid && lr==L_F1S2T.lr && tb==L_F1S2T.tb) L_F1S2T.tdcOffset[seg] += time/L_F1S2T.tdcGain[seg];
  else if(cid==L_F1S2B.cid && lr==L_F1S2B.lr && tb==L_F1S2B.tb) L_F1S2B.tdcOffset[seg] += time/L_F1S2B.tdcGain[seg];
  else if(cid==R_F1S2T.cid && lr==R_F1S2T.lr && tb==R_F1S2T.tb) R_F1S2T.tdcOffset[seg] += time/R_F1S2T.tdcGain[seg];
  else if(cid==R_F1S2B.cid && lr==R_F1S2B.lr && tb==R_F1S2B.tb) R_F1S2B.tdcOffset[seg] += time/R_F1S2B.tdcGain[seg];
  else if(cid==L_F1S0T.cid && lr==L_F1S0T.lr && tb==L_F1S0T.tb) L_F1S0T.tdcOffset[seg] += time/L_F1S0T.tdcGain[seg];
  else if(cid==L_F1S0B.cid && lr==L_F1S0B.lr && tb==L_F1S0B.tb) L_F1S0B.tdcOffset[seg] += time/L_F1S0B.tdcGain[seg];
  else if(cid==R_F1S0T.cid && lr==R_F1S0T.lr && tb==R_F1S0T.tb) R_F1S0T.tdcOffset[seg] += time/R_F1S0T.tdcGain[seg];
  else if(cid==R_F1S0B.cid && lr==R_F1S0B.lr && tb==R_F1S0B.tb) R_F1S0B.tdcOffset[seg] += time/R_F1S0B.tdcGain[seg];
  else if(cid==L_RF.cid    && lr==L_RF.lr)                      L_RF.tdcOffset[seg]    += time/L_RF.tdcGain[seg];
  else if(cid==R_RF.cid    && lr==R_RF.lr)                      R_RF.tdcOffset[seg]    += time/R_RF.tdcGain[seg];
  else if(cid==L_FbS2T.cid && lr==L_FbS2T.lr && tb==L_FbS2T.tb) L_FbS2T.tdcOffset[seg] += time/L_FbS2T.tdcGain[seg];
  else if(cid==L_FbS2B.cid && lr==L_FbS2B.lr && tb==L_FbS2B.tb) L_FbS2B.tdcOffset[seg] += time/L_FbS2B.tdcGain[seg];
  else if(cid==R_FbS2T.cid && lr==R_FbS2T.lr && tb==R_FbS2T.tb) R_FbS2T.tdcOffset[seg] += time/R_FbS2T.tdcGain[seg];
  else if(cid==R_FbS2B.cid && lr==R_FbS2B.lr && tb==R_FbS2B.tb) R_FbS2B.tdcOffset[seg] += time/R_FbS2B.tdcGain[seg];
  else if(cid==L_FbS0T.cid && lr==L_FbS0T.lr && tb==L_FbS0T.tb) L_FbS0T.tdcOffset[seg] += time/L_FbS0T.tdcGain[seg];
  else if(cid==L_FbS0B.cid && lr==L_FbS0B.lr && tb==L_FbS0B.tb) L_FbS0B.tdcOffset[seg] += time/L_FbS0B.tdcGain[seg];
  else if(cid==R_FbS0T.cid && lr==R_FbS0T.lr && tb==R_FbS0T.tb) R_FbS0T.tdcOffset[seg] += time/R_FbS0T.tdcGain[seg];
  else if(cid==R_FbS0B.cid && lr==R_FbS0B.lr && tb==R_FbS0B.tb) R_FbS0B.tdcOffset[seg] += time/R_FbS0B.tdcGain[seg];
  else   cerr << "[" << funcname << "]: unknown id" << endl;
}

///////////////////////////////////
double ParamMan::GetTdcOffset( int cid, int seg, int lr, int tb )
{
  static const std::string funcname = "ParamMan::GetTdcOffset";

  if     (cid==L_F1S2T.cid && lr==L_F1S2T.lr && tb==L_F1S2T.tb)return L_F1S2T.tdcOffset[seg];
  else if(cid==L_F1S2B.cid && lr==L_F1S2B.lr && tb==L_F1S2B.tb)return L_F1S2B.tdcOffset[seg];
  else if(cid==R_F1S2T.cid && lr==R_F1S2T.lr && tb==R_F1S2T.tb)return R_F1S2T.tdcOffset[seg];
  else if(cid==R_F1S2B.cid && lr==R_F1S2B.lr && tb==R_F1S2B.tb)return R_F1S2B.tdcOffset[seg];
  else if(cid==L_F1S0T.cid && lr==L_F1S0T.lr && tb==L_F1S0T.tb)return L_F1S0T.tdcOffset[seg];
  else if(cid==L_F1S0B.cid && lr==L_F1S0B.lr && tb==L_F1S0B.tb)return L_F1S0B.tdcOffset[seg];
  else if(cid==R_F1S0T.cid && lr==R_F1S0T.lr && tb==R_F1S0T.tb)return R_F1S0T.tdcOffset[seg];
  else if(cid==R_F1S0B.cid && lr==R_F1S0B.lr && tb==R_F1S0B.tb)return R_F1S0B.tdcOffset[seg];
  else if(cid==L_RF.cid    && lr==L_RF.lr)                     return L_RF.tdcOffset[seg]   ;
  else if(cid==R_RF.cid    && lr==R_RF.lr)                     return R_RF.tdcOffset[seg]   ;
  else if(cid==L_FbS2T.cid && lr==L_FbS2T.lr && tb==L_FbS2T.tb)return L_FbS2T.tdcOffset[seg];
  else if(cid==L_FbS2B.cid && lr==L_FbS2B.lr && tb==L_FbS2B.tb)return L_FbS2B.tdcOffset[seg];
  else if(cid==R_FbS2T.cid && lr==R_FbS2T.lr && tb==R_FbS2T.tb)return R_FbS2T.tdcOffset[seg];
  else if(cid==R_FbS2B.cid && lr==R_FbS2B.lr && tb==R_FbS2B.tb)return R_FbS2B.tdcOffset[seg];
  else if(cid==L_FbS0T.cid && lr==L_FbS0T.lr && tb==L_FbS0T.tb)return L_FbS0T.tdcOffset[seg];
  else if(cid==L_FbS0B.cid && lr==L_FbS0B.lr && tb==L_FbS0B.tb)return L_FbS0B.tdcOffset[seg];
  else if(cid==R_FbS0T.cid && lr==R_FbS0T.lr && tb==R_FbS0T.tb)return R_FbS0T.tdcOffset[seg];
  else if(cid==R_FbS0B.cid && lr==R_FbS0B.lr && tb==R_FbS0B.tb)return R_FbS0B.tdcOffset[seg];
  else   cerr << "[" << funcname << "]: unknown id" << endl;

  return -1.;
}

///////////////////////////////////
double ParamMan::GetTdcGain( int cid, int seg, int lr, int tb )
{
  static const std::string funcname = "ParamMan::GetTdcGain";

  if     (cid==L_F1S2T.cid && lr==L_F1S2T.lr && tb==L_F1S2T.tb)return L_F1S2T.tdcGain[seg];
  else if(cid==L_F1S2B.cid && lr==L_F1S2B.lr && tb==L_F1S2B.tb)return L_F1S2B.tdcGain[seg];
  else if(cid==R_F1S2T.cid && lr==R_F1S2T.lr && tb==R_F1S2T.tb)return R_F1S2T.tdcGain[seg];
  else if(cid==R_F1S2B.cid && lr==R_F1S2B.lr && tb==R_F1S2B.tb)return R_F1S2B.tdcGain[seg];
  else if(cid==L_F1S0T.cid && lr==L_F1S0T.lr && tb==L_F1S0T.tb)return L_F1S0T.tdcGain[seg];
  else if(cid==L_F1S0B.cid && lr==L_F1S0B.lr && tb==L_F1S0B.tb)return L_F1S0B.tdcGain[seg];
  else if(cid==R_F1S0T.cid && lr==R_F1S0T.lr && tb==R_F1S0T.tb)return R_F1S0T.tdcGain[seg];
  else if(cid==R_F1S0B.cid && lr==R_F1S0B.lr && tb==R_F1S0B.tb)return R_F1S0B.tdcGain[seg];
  else if(cid==L_RF.cid    && lr==L_RF.lr)                     return L_RF.tdcGain[seg]   ;
  else if(cid==R_RF.cid    && lr==R_RF.lr)                     return R_RF.tdcGain[seg]   ;
  else if(cid==L_FbS2T.cid && lr==L_FbS2T.lr && tb==L_FbS2T.tb)return L_FbS2T.tdcGain[seg];
  else if(cid==L_FbS2B.cid && lr==L_FbS2B.lr && tb==L_FbS2B.tb)return L_FbS2B.tdcGain[seg];
  else if(cid==R_FbS2T.cid && lr==R_FbS2T.lr && tb==R_FbS2T.tb)return R_FbS2T.tdcGain[seg];
  else if(cid==R_FbS2B.cid && lr==R_FbS2B.lr && tb==R_FbS2B.tb)return R_FbS2B.tdcGain[seg];
  else if(cid==L_FbS0T.cid && lr==L_FbS0T.lr && tb==L_FbS0T.tb)return L_FbS0T.tdcGain[seg];
  else if(cid==L_FbS0B.cid && lr==L_FbS0B.lr && tb==L_FbS0B.tb)return L_FbS0B.tdcGain[seg];
  else if(cid==R_FbS0T.cid && lr==R_FbS0T.lr && tb==R_FbS0T.tb)return R_FbS0T.tdcGain[seg];
  else if(cid==R_FbS0B.cid && lr==R_FbS0B.lr && tb==R_FbS0B.tb)return R_FbS0B.tdcGain[seg];
  else   cerr << "[" << funcname << "]: unknown id" << endl;

  return -1.;
}
///////////////////////////////////
double ParamMan::TimewalkCorrection( double par ,double adc,double peak)
{
  double tw=0.,tw_ref=0.;

  //if(par<=0.)return 0.;
  if(adc<=0.)return 0.;
  tw = par/sqrt(adc);  
  tw_ref = par/sqrt(2000.);  
  return tw-tw_ref;
}

///////////////////////////////////
void ParamMan::WriteToFile(const char* OutputFileName)   //wrinting param file
{
  ofstream fout;
  if( fout.is_open() ) fout.close();
  fout.open(OutputFileName, ios::out|ios::trunc);
  fout.setf(ios_base::fixed);
  //fout.open(name.str().c_str(), std::ios::out|std::ios::trunc);
  //fout.setf(std::ios_base::fixed);

  fout<<"R.s2.R.off = "<<flush;
  for(int i=0; i<nS2;i++){
    fout << std::setprecision(2)
         <<R_FbS2T.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"R.s2.L.off = "<<flush;
  for(int i=0; i<nS2;i++){
    fout << std::setprecision(2)
         <<R_FbS2B.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"R.s0.R.off = "<<flush;
  for(int i=0; i<nS0;i++){
    fout << std::setprecision(2)
         <<R_FbS0T.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"R.s0.L.off = "<<flush;
  for(int i=0; i<nS0;i++){
    fout << std::setprecision(2)
         <<R_FbS0B.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"R.s2.R.off_F1 = "<<flush;
  for(int i=0; i<nS2;i++){
    fout << std::setprecision(2)
         <<R_F1S2T.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"R.s2.L.off_F1 = "<<flush;
  for(int i=0; i<nS2;i++){
    fout << std::setprecision(2)
         <<R_F1S2B.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"R.s0.R.off_F1 = "<<flush;
  for(int i=0; i<nS0;i++){
    fout << std::setprecision(2)
         <<R_F1S0T.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"R.s0.L.off_F1 = "<<flush;
  for(int i=0; i<nS0;i++){
    fout << std::setprecision(2)
         <<R_F1S0B.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;


  //Left
  fout<<"L.s2.R.off = "<<flush;
  for(int i=0; i<nS2;i++){
    fout << std::setprecision(2)
         <<L_FbS2T.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"L.s2.L.off = "<<flush;
  for(int i=0; i<nS2;i++){
    fout << std::setprecision(2)
         <<L_FbS2B.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"L.s0.R.off = "<<flush;
  for(int i=0; i<nS0;i++){
    fout << std::setprecision(2)
         <<L_FbS0T.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"L.s0.L.off = "<<flush;
  for(int i=0; i<nS0;i++){
    fout << std::setprecision(2)
         <<L_FbS0B.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"L.s2.R.off_F1 = "<<flush;
  for(int i=0; i<nS2;i++){
    fout << std::setprecision(2)
         <<L_F1S2T.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"L.s2.L.off_F1 = "<<flush;
  for(int i=0; i<nS2;i++){
    fout << std::setprecision(2)
         <<L_F1S2B.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"L.s0.R.off_F1 = "<<flush;
  for(int i=0; i<nS0;i++){
    fout << std::setprecision(2)
         <<L_F1S0T.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;

  fout<<"L.s0.L.off_F1 = "<<flush;
  for(int i=0; i<nS0;i++){
    fout << std::setprecision(2)
         <<L_F1S0B.tdcOffset[i]<<" "<<flush;
  }
  fout<<endl;


  if(fout.is_open()) fout.close();
  cout << OutputFileName << " was written"<<endl;
}
