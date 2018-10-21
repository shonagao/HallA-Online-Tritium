#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "Tree.h"

Tree::Tree()
{
  tree = new TChain("T");
  //  cout<<"start tree class"<<endl;
}

Tree::~Tree()
{
  //cout<<"end tree class"<<endl;
}
void Tree::chain_tree(string ifname)
{
  tree->Add(Form("%s",ifname.c_str()));
  cout<<ifname<<endl;

  ENum = tree->GetEntries();
  tree->SetBranchStatus("*" ,0);
}
//////////////////////////////////////////////////
void Tree::readtreeS0L()
{
}
//////////////////////////////////////////////////
void Tree::readtreeS2L()
{
  tree->SetBranchStatus("L.s2.dedx"       ,1);tree->SetBranchAddress("L.s2.dedx"      , L_s2_dedx      );
  tree->SetBranchStatus("L.s2.dtime"      ,1);tree->SetBranchAddress("L.s2.dtime"     , L_s2_dtime     );
  tree->SetBranchStatus("L.s2.la"         ,1);tree->SetBranchAddress("L.s2.la"        , L_s2_la        );
  tree->SetBranchStatus("L.s2.la_c"       ,1);tree->SetBranchAddress("L.s2.la_c"      , L_s2_la_c      );
  tree->SetBranchStatus("L.s2.la_p"       ,1);tree->SetBranchAddress("L.s2.la_p"      , L_s2_la_p      );
  tree->SetBranchStatus("L.s2.lbadped"    ,1);tree->SetBranchAddress("L.s2.lbadped"   , L_s2_lbadped   );
  tree->SetBranchStatus("L.s2.lnhits"     ,1);tree->SetBranchAddress("L.s2.lnhits"    , L_s2_lnhits    );
  tree->SetBranchStatus("L.s2.loverflow"  ,1);tree->SetBranchAddress("L.s2.loverflow" , L_s2_loverflow );
  tree->SetBranchStatus("L.s2.lpeak"      ,1);tree->SetBranchAddress("L.s2.lpeak"     , L_s2_lpeak     );
  tree->SetBranchStatus("L.s2.lt"         ,1);tree->SetBranchAddress("L.s2.lt"        , L_s2_lt        );
  tree->SetBranchStatus("L.s2.lt_c"       ,1);tree->SetBranchAddress("L.s2.lt_c"      , L_s2_lt_c      );
  tree->SetBranchStatus("L.s2.lt_fadc"    ,1);tree->SetBranchAddress("L.s2.lt_fadc"   , L_s2_lt_fadc   );
  tree->SetBranchStatus("L.s2.ltc_fadc"   ,1);tree->SetBranchAddress("L.s2.ltc_fadc"  , L_s2_ltc_fadc  );
  tree->SetBranchStatus("L.s2.lunderflow" ,1);tree->SetBranchAddress("L.s2.lunderflow", L_s2_lunderflow);
  tree->SetBranchStatus("L.s2.ra"         ,1);tree->SetBranchAddress("L.s2.ra"        , L_s2_ra        );
  tree->SetBranchStatus("L.s2.ra_c"       ,1);tree->SetBranchAddress("L.s2.ra_c"      , L_s2_ra_c      );
  tree->SetBranchStatus("L.s2.ra_p"       ,1);tree->SetBranchAddress("L.s2.ra_p"      , L_s2_ra_p      );
  tree->SetBranchStatus("L.s2.rbadped"    ,1);tree->SetBranchAddress("L.s2.rbadped"   , L_s2_rbadped   );
  tree->SetBranchStatus("L.s2.rnhits"     ,1);tree->SetBranchAddress("L.s2.rnhits"    , L_s2_rnhits    );
  tree->SetBranchStatus("L.s2.roverflow"  ,1);tree->SetBranchAddress("L.s2.roverflow" , L_s2_roverflow );
  tree->SetBranchStatus("L.s2.rpeak"      ,1);tree->SetBranchAddress("L.s2.rpeak"     , L_s2_rpeak     );
  tree->SetBranchStatus("L.s2.rt"         ,1);tree->SetBranchAddress("L.s2.rt"        , L_s2_rt        );
  tree->SetBranchStatus("L.s2.rt_c"       ,1);tree->SetBranchAddress("L.s2.rt_c"      , L_s2_rt_c      );
  tree->SetBranchStatus("L.s2.rt_fadc"    ,1);tree->SetBranchAddress("L.s2.rt_fadc"   , L_s2_rt_fadc   );
  tree->SetBranchStatus("L.s2.rtc_fadc"   ,1);tree->SetBranchAddress("L.s2.rtc_fadc"  , L_s2_rtc_fadc  );
  tree->SetBranchStatus("L.s2.runderflow" ,1);tree->SetBranchAddress("L.s2.runderflow", L_s2_runderflow);
  tree->SetBranchStatus("L.s2.t_pads"     ,1);tree->SetBranchAddress("L.s2.t_pads"    ,&L_s2_t_pads    );
  tree->SetBranchStatus("L.s2.time"       ,1);tree->SetBranchAddress("L.s2.time"      , L_s2_time      );
  tree->SetBranchStatus("L.s2.trdx"       ,1);tree->SetBranchAddress("L.s2.trdx"      ,&L_s2_trdx      );
  tree->SetBranchStatus("L.s2.troff"      ,1);tree->SetBranchAddress("L.s2.troff"     , L_s2_troff     );
  tree->SetBranchStatus("L.s2.trpad"      ,1);tree->SetBranchAddress("L.s2.trpad"     ,&L_s2_trpad     );
  tree->SetBranchStatus("L.s2.trpath"     ,1);tree->SetBranchAddress("L.s2.trpath"    ,&L_s2_trpath    );
  tree->SetBranchStatus("L.s2.trx"        ,1);tree->SetBranchAddress("L.s2.trx"       ,&L_s2_trx       );
  tree->SetBranchStatus("L.s2.try"        ,1);tree->SetBranchAddress("L.s2.try"       ,&L_s2_try       );
  tree->SetBranchStatus("L.s2.y_adc"      ,1);tree->SetBranchAddress("L.s2.y_adc"     ,L_s2_y_adc     );
  tree->SetBranchStatus("L.s2.y_t"        ,1);tree->SetBranchAddress("L.s2.y_t"       ,L_s2_y_t       );
}                                 
//////////////////////////////////////////////////
void Tree::readtreeS2R()
{
}                                 
//////////////////////////////////////////////////
void Tree::readtreeS0R()
{
}                                 
//////////////////////////////////////////////////
