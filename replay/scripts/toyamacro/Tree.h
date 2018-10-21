#ifndef Tree_h
#define Tree_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"

class Tree
{
public:
  Tree();
  ~Tree();

//private:
public:

  int ENum;

  ////////
  //Left//
  ////////
  //S2L
  double L_s2_dedx[16]      , L_s2_dtime[16] , L_s2_la[16]      , L_s2_la_c[16]     , L_s2_la_p[16]       , L_s2_lbadped[16]   , L_s2_lnhits[16]    ;
  double L_s2_loverflow[16] , L_s2_lpeak[16] , L_s2_lt[16]      , L_s2_lt_c[16]     , L_s2_lt_fadc[16]    , L_s2_ltc_fadc[16]  , L_s2_lunderflow[16];
  double L_s2_ra[16]        , L_s2_ra_c[16]  , L_s2_ra_p[16]    , L_s2_rbadped[16]  , L_s2_rnhits[16]     , L_s2_roverflow[16] , L_s2_rpeak[16]     ;
  double L_s2_rt[16]        , L_s2_rt_c[16]  , L_s2_rt_fadc[16] , L_s2_rtc_fadc[16] , L_s2_runderflow[16] ;
  double L_s2_t_pads;
  double L_s2_time[16], L_s2_troff[16], L_s2_y_adc[16], L_s2_y_t[16];
  double L_s2_trdx      ;//not filled
  double L_s2_trpath    ;//not filled
  double L_s2_trpad     ;//not filled
  double L_s2_trx       ;//not filled
  double L_s2_try       ;//not filled

  /////////
  //Right//
  /////////
//private:
public:
  TChain *tree;

public:
  void chain_tree(string ifname);
  ////////
  //Left//
  ////////
  void readtreeS0L();
  void readtreeS2L();
  /////////
  //Right//
  /////////
  void readtreeS0R();
  void readtreeS2R();

};

#endif
