#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "Tree.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
Tree::Tree()
{
  tree = new TChain("T");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
Tree::~Tree()
{
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Tree::merge(string ifname){
  tree->Add(ifname.c_str());
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Tree::readtree_COMN()
{
  tree->SetBranchStatus("rbax"   ,1);  tree->SetBranchAddress("rbax"   ,&rbax   ); 
  tree->SetBranchStatus("rbay"   ,1);  tree->SetBranchAddress("rbay"   ,&rbay   ); 
  tree->SetBranchStatus("rbbx"   ,1);  tree->SetBranchAddress("rbbx"   ,&rbbx   ); 
  tree->SetBranchStatus("rbby"   ,1);  tree->SetBranchAddress("rbby"   ,&rbby   ); 
  tree->SetBranchStatus("rbx"    ,1);  tree->SetBranchAddress("rbx"    ,&rbx    ); 
  tree->SetBranchStatus("rby"    ,1);  tree->SetBranchAddress("rby"    ,&rby    ); 
  tree->SetBranchStatus("bpmaws" ,1);  tree->SetBranchAddress("bpmaws" ,&bpmaws ); 
  tree->SetBranchStatus("bpmbws" ,1);  tree->SetBranchAddress("bpmbws" ,&bpmbws ); 
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Tree::readtree_LHRS()
{
  tree->SetBranchStatus("Ndata.DL.T1"     ,1);  tree->SetBranchAddress("Ndata.DL.T1"     ,&L_T1            );
  tree->SetBranchStatus("Ndata.DL.T2"     ,1);  tree->SetBranchAddress("Ndata.DL.T2"     ,&L_T2            );
  tree->SetBranchStatus("Ndata.DL.T3"     ,1);  tree->SetBranchAddress("Ndata.DL.T3"     ,&L_T3            );
  tree->SetBranchStatus("Ndata.DL.T4"     ,1);  tree->SetBranchAddress("Ndata.DL.T4"     ,&L_T4            );
  tree->SetBranchStatus("Ndata.DL.T5"     ,1);  tree->SetBranchAddress("Ndata.DL.T5"     ,&L_T5            );
  tree->SetBranchStatus("Ndata.DL.T6"     ,1);  tree->SetBranchAddress("Ndata.DL.T6"     ,&L_T6            );
  tree->SetBranchStatus("Ndata.DL.T7"     ,1);  tree->SetBranchAddress("Ndata.DL.T7"     ,&L_T7            );
  tree->SetBranchStatus("Ndata.DL.T8"     ,1);  tree->SetBranchAddress("Ndata.DL.T8"     ,&L_T8            );

  tree->SetBranchStatus("L.s0.la"         ,1);  tree->SetBranchAddress("L.s0.la"         , L_s0_la         );
  tree->SetBranchStatus("L.s0.la_c"       ,1);  tree->SetBranchAddress("L.s0.la_c"       , L_s0_la_c       );
  tree->SetBranchStatus("L.s0.la_p"       ,1);  tree->SetBranchAddress("L.s0.la_p"       , L_s0_la_p       );
  tree->SetBranchStatus("L.s0.lt"         ,1);  tree->SetBranchAddress("L.s0.lt"         , L_s0_lt         );
  tree->SetBranchStatus("L.s0.lt_c"       ,1);  tree->SetBranchAddress("L.s0.lt_c"       , L_s0_lt_c       );
  tree->SetBranchStatus("L.s0.ra"         ,1);  tree->SetBranchAddress("L.s0.ra"         , L_s0_ra         );
  tree->SetBranchStatus("L.s0.ra_c"       ,1);  tree->SetBranchAddress("L.s0.ra_c"       , L_s0_ra_c       );
  tree->SetBranchStatus("L.s0.ra_p"       ,1);  tree->SetBranchAddress("L.s0.ra_p"       , L_s0_ra_p       );
  tree->SetBranchStatus("L.s0.rt"         ,1);  tree->SetBranchAddress("L.s0.rt"         , L_s0_rt         );
  tree->SetBranchStatus("L.s0.rt_c"       ,1);  tree->SetBranchAddress("L.s0.rt_c"       , L_s0_rt_c       );
  tree->SetBranchStatus("L.s0.t_pads"     ,1);  tree->SetBranchAddress("L.s0.t_pads"     , L_s0_t_pads     );
  tree->SetBranchStatus("L.s0.time"       ,1);  tree->SetBranchAddress("L.s0.time"       , L_s0_time       );
  tree->SetBranchStatus("L.s0.dedx"       ,1);  tree->SetBranchAddress("L.s0.dedx"       , L_s0_dedx       );
  tree->SetBranchStatus("L.s0.trdy"       ,1);  tree->SetBranchAddress("L.s0.trdy"       , L_s0_trdy       );
  tree->SetBranchStatus("L.s0.troff"      ,1);  tree->SetBranchAddress("L.s0.troff"      , L_s0_troff      );
  tree->SetBranchStatus("L.s0.trpad"      ,1);  tree->SetBranchAddress("L.s0.trpad"      , L_s0_trpad      );
  tree->SetBranchStatus("L.s0.trpath"     ,1);  tree->SetBranchAddress("L.s0.trpath"     , L_s0_trpath     );
  tree->SetBranchStatus("L.s0.trx"        ,1);  tree->SetBranchAddress("L.s0.trx"        , L_s0_trx        );
  tree->SetBranchStatus("L.s0.try"        ,1);  tree->SetBranchAddress("L.s0.try"        , L_s0_try        );
  tree->SetBranchStatus("L.s0.lnhits"     ,1);  tree->SetBranchAddress("L.s0.lnhits"     , L_s0_lnhits     );
  tree->SetBranchStatus("L.s0.lpeak"      ,1);  tree->SetBranchAddress("L.s0.lpeak"      , L_s0_lpeak      );
  tree->SetBranchStatus("L.s0.lt_fadc"    ,1);  tree->SetBranchAddress("L.s0.lt_fadc"    , L_s0_lt_fadc    );
  tree->SetBranchStatus("L.s0.ltc_fadc"   ,1);  tree->SetBranchAddress("L.s0.ltc_fadc"   , L_s0_ltc_fadc   );
  tree->SetBranchStatus("L.s0.lbadped"    ,0);  tree->SetBranchAddress("L.s0.lbadped"    , L_s0_lbadped    );
  tree->SetBranchStatus("L.s0.lunderflow" ,0);  tree->SetBranchAddress("L.s0.lunderflow" , L_s0_lunderflow );
  tree->SetBranchStatus("L.s0.loverflow"  ,0);  tree->SetBranchAddress("L.s0.loverflow"  , L_s0_loverflow  );
  tree->SetBranchStatus("L.s0.rnhits"     ,1);  tree->SetBranchAddress("L.s0.rnhits"     , L_s0_rnhits     );
  tree->SetBranchStatus("L.s0.rpeak"      ,1);  tree->SetBranchAddress("L.s0.rpeak"      , L_s0_rpeak      );
  tree->SetBranchStatus("L.s0.rt_fadc"    ,1);  tree->SetBranchAddress("L.s0.rt_fadc"    , L_s0_rt_fadc    );
  tree->SetBranchStatus("L.s0.rtc_fadc"   ,1);  tree->SetBranchAddress("L.s0.rtc_fadc"   , L_s0_rtc_fadc   );
  tree->SetBranchStatus("L.s0.rbadped"    ,0);  tree->SetBranchAddress("L.s0.rbadped"    , L_s0_rbadped    );
  tree->SetBranchStatus("L.s0.runderflow" ,0);  tree->SetBranchAddress("L.s0.runderflow" , L_s0_runderflow );
  tree->SetBranchStatus("L.s0.roverflow"  ,0);  tree->SetBranchAddress("L.s0.roverflow"  , L_s0_roverflow  );

  tree->SetBranchStatus("L.s2.la"         ,1);  tree->SetBranchAddress("L.s2.la"         , L_s2_la         );
  tree->SetBranchStatus("L.s2.la_c"       ,1);  tree->SetBranchAddress("L.s2.la_c"       , L_s2_la_c       );
  tree->SetBranchStatus("L.s2.la_p"       ,1);  tree->SetBranchAddress("L.s2.la_p"       , L_s2_la_p       );
  tree->SetBranchStatus("L.s2.lt"         ,1);  tree->SetBranchAddress("L.s2.lt"         , L_s2_lt         );
  tree->SetBranchStatus("L.s2.lt_c"       ,1);  tree->SetBranchAddress("L.s2.lt_c"       , L_s2_lt_c       );
  tree->SetBranchStatus("L.s2.ra"         ,1);  tree->SetBranchAddress("L.s2.ra"         , L_s2_ra         );
  tree->SetBranchStatus("L.s2.ra_c"       ,1);  tree->SetBranchAddress("L.s2.ra_c"       , L_s2_ra_c       );
  tree->SetBranchStatus("L.s2.ra_p"       ,1);  tree->SetBranchAddress("L.s2.ra_p"       , L_s2_ra_p       );
  tree->SetBranchStatus("L.s2.rt"         ,1);  tree->SetBranchAddress("L.s2.rt"         , L_s2_rt         );
  tree->SetBranchStatus("L.s2.rt_c"       ,1);  tree->SetBranchAddress("L.s2.rt_c"       , L_s2_rt_c       );
  tree->SetBranchStatus("L.s2.t_pads"     ,1);  tree->SetBranchAddress("L.s2.t_pads"     , L_s2_t_pads     );
  tree->SetBranchStatus("L.s2.time"       ,1);  tree->SetBranchAddress("L.s2.time"       , L_s2_time       );
  tree->SetBranchStatus("L.s2.dedx"       ,1);  tree->SetBranchAddress("L.s2.dedx"       , L_s2_dedx       );
  tree->SetBranchStatus("L.s2.trdx"       ,1);  tree->SetBranchAddress("L.s2.trdx"       , L_s2_trdx       );
  tree->SetBranchStatus("L.s2.troff"      ,1);  tree->SetBranchAddress("L.s2.troff"      , L_s2_troff      );
  tree->SetBranchStatus("L.s2.trpad"      ,1);  tree->SetBranchAddress("L.s2.trpad"      , L_s2_trpad      );
  tree->SetBranchStatus("L.s2.trpath"     ,1);  tree->SetBranchAddress("L.s2.trpath"     , L_s2_trpath     );
  tree->SetBranchStatus("L.s2.trx"        ,1);  tree->SetBranchAddress("L.s2.trx"        , L_s2_trx        );
  tree->SetBranchStatus("L.s2.try"        ,1);  tree->SetBranchAddress("L.s2.try"        , L_s2_try        );
  tree->SetBranchStatus("L.s2.lnhits"     ,1);  tree->SetBranchAddress("L.s2.lnhits"     , L_s2_lnhits     );
  tree->SetBranchStatus("L.s2.lpeak"      ,1);  tree->SetBranchAddress("L.s2.lpeak"      , L_s2_lpeak      );
  tree->SetBranchStatus("L.s2.lt_fadc"    ,1);  tree->SetBranchAddress("L.s2.lt_fadc"    , L_s2_lt_fadc    );
  tree->SetBranchStatus("L.s2.ltc_fadc"   ,1);  tree->SetBranchAddress("L.s2.ltc_fadc"   , L_s2_ltc_fadc   );
  tree->SetBranchStatus("L.s2.lbadped"    ,0);  tree->SetBranchAddress("L.s2.lbadped"    , L_s2_lbadped    );
  tree->SetBranchStatus("L.s2.lunderflow" ,0);  tree->SetBranchAddress("L.s2.lunderflow" , L_s2_lunderflow );
  tree->SetBranchStatus("L.s2.loverflow"  ,0);  tree->SetBranchAddress("L.s2.loverflow"  , L_s2_loverflow  );
  tree->SetBranchStatus("L.s2.rnhits"     ,1);  tree->SetBranchAddress("L.s2.rnhits"     , L_s2_rnhits     );
  tree->SetBranchStatus("L.s2.rpeak"      ,1);  tree->SetBranchAddress("L.s2.rpeak"      , L_s2_rpeak      );
  tree->SetBranchStatus("L.s2.rt_fadc"    ,1);  tree->SetBranchAddress("L.s2.rt_fadc"    , L_s2_rt_fadc    );
  tree->SetBranchStatus("L.s2.rtc_fadc"   ,1);  tree->SetBranchAddress("L.s2.rtc_fadc"   , L_s2_rtc_fadc   );
  tree->SetBranchStatus("L.s2.rbadped"    ,0);  tree->SetBranchAddress("L.s2.rbadped"    , L_s2_rbadped    );
  tree->SetBranchStatus("L.s2.runderflow" ,0);  tree->SetBranchAddress("L.s2.runderflow" , L_s2_runderflow );
  tree->SetBranchStatus("L.s2.roverflow"  ,0);  tree->SetBranchAddress("L.s2.roverflow"  , L_s2_roverflow  );

  tree->SetBranchStatus("L.tr.n"          ,1);  tree->SetBranchAddress("L.tr.n"          ,&L_tr_n          );
  tree->SetBranchStatus("L.tr.flag"       ,1);  tree->SetBranchAddress("L.tr.flag"       , L_tr_flag       );
  tree->SetBranchStatus("L.tr.ndof"       ,1);  tree->SetBranchAddress("L.tr.ndof"       , L_tr_ndof       );
  tree->SetBranchStatus("L.tr.chi2"       ,1);  tree->SetBranchAddress("L.tr.chi2"       , L_tr_chi2       );
  tree->SetBranchStatus("L.tr.beta"       ,1);  tree->SetBranchAddress("L.tr.beta"       , L_tr_beta       );
  tree->SetBranchStatus("L.tr.d_x"        ,1);  tree->SetBranchAddress("L.tr.d_x"        , L_tr_d_x        );
  tree->SetBranchStatus("L.tr.d_y"        ,1);  tree->SetBranchAddress("L.tr.d_y"        , L_tr_d_y        );
  tree->SetBranchStatus("L.tr.d_th"       ,1);  tree->SetBranchAddress("L.tr.d_th"       , L_tr_d_th       );
  tree->SetBranchStatus("L.tr.d_ph"       ,1);  tree->SetBranchAddress("L.tr.d_ph"       , L_tr_d_ph       );
  tree->SetBranchStatus("L.tr.r_x"        ,1);  tree->SetBranchAddress("L.tr.r_x"        , L_tr_r_x        );
  tree->SetBranchStatus("L.tr.r_y"        ,1);  tree->SetBranchAddress("L.tr.r_y"        , L_tr_r_y        );
  tree->SetBranchStatus("L.tr.r_th"       ,1);  tree->SetBranchAddress("L.tr.r_th"       , L_tr_r_th       );
  tree->SetBranchStatus("L.tr.r_ph"       ,1);  tree->SetBranchAddress("L.tr.r_ph"       , L_tr_r_ph       );
  tree->SetBranchStatus("L.tr.x"          ,1);  tree->SetBranchAddress("L.tr.x"          , L_tr_x          );
  tree->SetBranchStatus("L.tr.y"          ,1);  tree->SetBranchAddress("L.tr.y"          , L_tr_y          );
  tree->SetBranchStatus("L.tr.th"         ,1);  tree->SetBranchAddress("L.tr.th"         , L_tr_th         );
  tree->SetBranchStatus("L.tr.ph"         ,1);  tree->SetBranchAddress("L.tr.ph"         , L_tr_ph         );
  tree->SetBranchStatus("L.tr.time"       ,1);  tree->SetBranchAddress("L.tr.time"       , L_tr_time       );
  tree->SetBranchStatus("L.tr.p"          ,1);  tree->SetBranchAddress("L.tr.p"          , L_tr_p          );
  tree->SetBranchStatus("L.tr.pathl"      ,1);  tree->SetBranchAddress("L.tr.pathl"      , L_tr_pathl      );
  tree->SetBranchStatus("L.tr.px"         ,1);  tree->SetBranchAddress("L.tr.px"         , L_tr_px         );
  tree->SetBranchStatus("L.tr.py"         ,1);  tree->SetBranchAddress("L.tr.py"         , L_tr_py         );
  tree->SetBranchStatus("L.tr.pz"         ,1);  tree->SetBranchAddress("L.tr.pz"         , L_tr_pz         );
  tree->SetBranchStatus("L.tr.tg_dp"      ,1);  tree->SetBranchAddress("L.tr.tg_dp"      , L_tr_tg_dp      );
  tree->SetBranchStatus("L.tr.tg_y"       ,1);  tree->SetBranchAddress("L.tr.tg_y"       , L_tr_tg_y       );
  tree->SetBranchStatus("L.tr.tg_th"      ,1);  tree->SetBranchAddress("L.tr.tg_th"      , L_tr_tg_th      );
  tree->SetBranchStatus("L.tr.tg_ph"      ,1);  tree->SetBranchAddress("L.tr.tg_ph"      , L_tr_tg_ph      );
  tree->SetBranchStatus("L.tr.vx"         ,1);  tree->SetBranchAddress("L.tr.vx"         , L_tr_vx         );
  tree->SetBranchStatus("L.tr.vy"         ,1);  tree->SetBranchAddress("L.tr.vy"         , L_tr_vy         );
  tree->SetBranchStatus("L.tr.vz"         ,1);  tree->SetBranchAddress("L.tr.vz"         , L_tr_vz         );

  tree->SetBranchStatus("DL.ltRFtime"     ,1);  tree->SetBranchAddress("DL.ltRFtime"     , L_rtRFtime      );
  tree->SetBranchStatus("LTDC.F1FirstHit" ,1);  tree->SetBranchAddress("LTDC.F1FirstHit" , L_F1Fhit        );

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Tree::readtree_RHRS()
{
  tree->SetBranchStatus("Ndata.DR.T1"     ,1);  tree->SetBranchAddress("Ndata.DR.T1"     ,&R_T1            );
  tree->SetBranchStatus("Ndata.DR.T2"     ,1);  tree->SetBranchAddress("Ndata.DR.T2"     ,&R_T2            );
  tree->SetBranchStatus("Ndata.DR.T3"     ,1);  tree->SetBranchAddress("Ndata.DR.T3"     ,&R_T3            );
  tree->SetBranchStatus("Ndata.DR.T4"     ,1);  tree->SetBranchAddress("Ndata.DR.T4"     ,&R_T4            );
  tree->SetBranchStatus("Ndata.DR.T5"     ,1);  tree->SetBranchAddress("Ndata.DR.T5"     ,&R_T5            );
  tree->SetBranchStatus("Ndata.DR.T6"     ,1);  tree->SetBranchAddress("Ndata.DR.T6"     ,&R_T6            );
  tree->SetBranchStatus("Ndata.DR.T7"     ,1);  tree->SetBranchAddress("Ndata.DR.T7"     ,&R_T7            );
  tree->SetBranchStatus("Ndata.DR.T8"     ,1);  tree->SetBranchAddress("Ndata.DR.T8"     ,&R_T8            );

  tree->SetBranchStatus("R.s0.la"         ,1);  tree->SetBranchAddress("R.s0.la"         , R_s0_la         );
  tree->SetBranchStatus("R.s0.la_c"       ,1);  tree->SetBranchAddress("R.s0.la_c"       , R_s0_la_c       );
  tree->SetBranchStatus("R.s0.la_p"       ,1);  tree->SetBranchAddress("R.s0.la_p"       , R_s0_la_p       );
  tree->SetBranchStatus("R.s0.lt"         ,1);  tree->SetBranchAddress("R.s0.lt"         , R_s0_lt         );
  tree->SetBranchStatus("R.s0.lt_c"       ,1);  tree->SetBranchAddress("R.s0.lt_c"       , R_s0_lt_c       );
  tree->SetBranchStatus("R.s0.ra"         ,1);  tree->SetBranchAddress("R.s0.ra"         , R_s0_ra         );
  tree->SetBranchStatus("R.s0.ra_c"       ,1);  tree->SetBranchAddress("R.s0.ra_c"       , R_s0_ra_c       );
  tree->SetBranchStatus("R.s0.ra_p"       ,1);  tree->SetBranchAddress("R.s0.ra_p"       , R_s0_ra_p       );
  tree->SetBranchStatus("R.s0.rt"         ,1);  tree->SetBranchAddress("R.s0.rt"         , R_s0_rt         );
  tree->SetBranchStatus("R.s0.rt_c"       ,1);  tree->SetBranchAddress("R.s0.rt_c"       , R_s0_rt_c       );
  tree->SetBranchStatus("R.s0.t_pads"     ,1);  tree->SetBranchAddress("R.s0.t_pads"     , R_s0_t_pads     );
  tree->SetBranchStatus("R.s0.time"       ,1);  tree->SetBranchAddress("R.s0.time"       , R_s0_time       );
  tree->SetBranchStatus("R.s0.dedx"       ,1);  tree->SetBranchAddress("R.s0.dedx"       , R_s0_dedx       );
  tree->SetBranchStatus("R.s0.trdy"       ,1);  tree->SetBranchAddress("R.s0.trdy"       , R_s0_trdy       );
  tree->SetBranchStatus("R.s0.troff"      ,1);  tree->SetBranchAddress("R.s0.troff"      , R_s0_troff      );
  tree->SetBranchStatus("R.s0.trpad"      ,1);  tree->SetBranchAddress("R.s0.trpad"      , R_s0_trpad      );
  tree->SetBranchStatus("R.s0.trpath"     ,1);  tree->SetBranchAddress("R.s0.trpath"     , R_s0_trpath     );
  tree->SetBranchStatus("R.s0.trx"        ,1);  tree->SetBranchAddress("R.s0.trx"        , R_s0_trx        );
  tree->SetBranchStatus("R.s0.try"        ,1);  tree->SetBranchAddress("R.s0.try"        , R_s0_try        );
  tree->SetBranchStatus("R.s0.lnhits"     ,1);  tree->SetBranchAddress("R.s0.lnhits"     , R_s0_lnhits     );
  tree->SetBranchStatus("R.s0.lpeak"      ,1);  tree->SetBranchAddress("R.s0.lpeak"      , R_s0_lpeak      );
  tree->SetBranchStatus("R.s0.lt_fadc"    ,1);  tree->SetBranchAddress("R.s0.lt_fadc"    , R_s0_lt_fadc    );
  tree->SetBranchStatus("R.s0.ltc_fadc"   ,1);  tree->SetBranchAddress("R.s0.ltc_fadc"   , R_s0_ltc_fadc   );
  tree->SetBranchStatus("R.s0.lbadped"    ,0);  tree->SetBranchAddress("R.s0.lbadped"    , R_s0_lbadped    );
  tree->SetBranchStatus("R.s0.lunderflow" ,0);  tree->SetBranchAddress("R.s0.lunderflow" , R_s0_lunderflow );
  tree->SetBranchStatus("R.s0.loverflow"  ,0);  tree->SetBranchAddress("R.s0.loverflow"  , R_s0_loverflow  );
  tree->SetBranchStatus("R.s0.rnhits"     ,1);  tree->SetBranchAddress("R.s0.rnhits"     , R_s0_rnhits     );
  tree->SetBranchStatus("R.s0.rpeak"      ,1);  tree->SetBranchAddress("R.s0.rpeak"      , R_s0_rpeak      );
  tree->SetBranchStatus("R.s0.rt_fadc"    ,1);  tree->SetBranchAddress("R.s0.rt_fadc"    , R_s0_rt_fadc    );
  tree->SetBranchStatus("R.s0.rtc_fadc"   ,1);  tree->SetBranchAddress("R.s0.rtc_fadc"   , R_s0_rtc_fadc   );
  tree->SetBranchStatus("R.s0.rbadped"    ,0);  tree->SetBranchAddress("R.s0.rbadped"    , R_s0_rbadped    );
  tree->SetBranchStatus("R.s0.runderflow" ,0);  tree->SetBranchAddress("R.s0.runderflow" , R_s0_runderflow );
  tree->SetBranchStatus("R.s0.roverflow"  ,0);  tree->SetBranchAddress("R.s0.roverflow"  , R_s0_roverflow  );

  tree->SetBranchStatus("R.s2.la"         ,1);  tree->SetBranchAddress("R.s2.la"         , R_s2_la         );
  tree->SetBranchStatus("R.s2.la_c"       ,1);  tree->SetBranchAddress("R.s2.la_c"       , R_s2_la_c       );
  tree->SetBranchStatus("R.s2.la_p"       ,1);  tree->SetBranchAddress("R.s2.la_p"       , R_s2_la_p       );
  tree->SetBranchStatus("R.s2.lt"         ,1);  tree->SetBranchAddress("R.s2.lt"         , R_s2_lt         );
  tree->SetBranchStatus("R.s2.lt_c"       ,1);  tree->SetBranchAddress("R.s2.lt_c"       , R_s2_lt_c       );
  tree->SetBranchStatus("R.s2.ra"         ,1);  tree->SetBranchAddress("R.s2.ra"         , R_s2_ra         );
  tree->SetBranchStatus("R.s2.ra_c"       ,1);  tree->SetBranchAddress("R.s2.ra_c"       , R_s2_ra_c       );
  tree->SetBranchStatus("R.s2.ra_p"       ,1);  tree->SetBranchAddress("R.s2.ra_p"       , R_s2_ra_p       );
  tree->SetBranchStatus("R.s2.rt"         ,1);  tree->SetBranchAddress("R.s2.rt"         , R_s2_rt         );
  tree->SetBranchStatus("R.s2.rt_c"       ,1);  tree->SetBranchAddress("R.s2.rt_c"       , R_s2_rt_c       );
  tree->SetBranchStatus("R.s2.t_pads"     ,1);  tree->SetBranchAddress("R.s2.t_pads"     , R_s2_t_pads     );
  tree->SetBranchStatus("R.s2.time"       ,1);  tree->SetBranchAddress("R.s2.time"       , R_s2_time       );
  tree->SetBranchStatus("R.s2.dedx"       ,1);  tree->SetBranchAddress("R.s2.dedx"       , R_s2_dedx       );
  tree->SetBranchStatus("R.s2.trdx"       ,1);  tree->SetBranchAddress("R.s2.trdx"       , R_s2_trdx       );
  tree->SetBranchStatus("R.s2.troff"      ,1);  tree->SetBranchAddress("R.s2.troff"      , R_s2_troff      );
  tree->SetBranchStatus("R.s2.trpad"      ,1);  tree->SetBranchAddress("R.s2.trpad"      , R_s2_trpad      );
  tree->SetBranchStatus("R.s2.trpath"     ,1);  tree->SetBranchAddress("R.s2.trpath"     , R_s2_trpath     );
  tree->SetBranchStatus("R.s2.trx"        ,1);  tree->SetBranchAddress("R.s2.trx"        , R_s2_trx        );
  tree->SetBranchStatus("R.s2.try"        ,1);  tree->SetBranchAddress("R.s2.try"        , R_s2_try        );
  tree->SetBranchStatus("R.s2.lnhits"     ,1);  tree->SetBranchAddress("R.s2.lnhits"     , R_s2_lnhits     );
  tree->SetBranchStatus("R.s2.lpeak"      ,1);  tree->SetBranchAddress("R.s2.lpeak"      , R_s2_lpeak      );
  tree->SetBranchStatus("R.s2.lt_fadc"    ,1);  tree->SetBranchAddress("R.s2.lt_fadc"    , R_s2_lt_fadc    );
  tree->SetBranchStatus("R.s2.ltc_fadc"   ,1);  tree->SetBranchAddress("R.s2.ltc_fadc"   , R_s2_ltc_fadc   );
  tree->SetBranchStatus("R.s2.lbadped"    ,0);  tree->SetBranchAddress("R.s2.lbadped"    , R_s2_lbadped    );
  tree->SetBranchStatus("R.s2.lunderflow" ,0);  tree->SetBranchAddress("R.s2.lunderflow" , R_s2_lunderflow );
  tree->SetBranchStatus("R.s2.loverflow"  ,0);  tree->SetBranchAddress("R.s2.loverflow"  , R_s2_loverflow  );
  tree->SetBranchStatus("R.s2.rnhits"     ,1);  tree->SetBranchAddress("R.s2.rnhits"     , R_s2_rnhits     );
  tree->SetBranchStatus("R.s2.rpeak"      ,1);  tree->SetBranchAddress("R.s2.rpeak"      , R_s2_rpeak      );
  tree->SetBranchStatus("R.s2.rt_fadc"    ,1);  tree->SetBranchAddress("R.s2.rt_fadc"    , R_s2_rt_fadc    );
  tree->SetBranchStatus("R.s2.rtc_fadc"   ,1);  tree->SetBranchAddress("R.s2.rtc_fadc"   , R_s2_rtc_fadc   );
  tree->SetBranchStatus("R.s2.rbadped"    ,0);  tree->SetBranchAddress("R.s2.rbadped"    , R_s2_rbadped    );
  tree->SetBranchStatus("R.s2.runderflow" ,0);  tree->SetBranchAddress("R.s2.runderflow" , R_s2_runderflow );
  tree->SetBranchStatus("R.s2.roverflow"  ,0);  tree->SetBranchAddress("R.s2.roverflow"  , R_s2_roverflow  );

  tree->SetBranchStatus("R.a1.a"          ,1);  tree->SetBranchAddress("R.a1.a"          , R_a1_a          ); 
  tree->SetBranchStatus("R.a1.a_c"        ,1);  tree->SetBranchAddress("R.a1.a_c"        , R_a1_a_c        ); 
  tree->SetBranchStatus("R.a1.a_p"        ,1);  tree->SetBranchAddress("R.a1.a_p"        , R_a1_a_p        ); 
  tree->SetBranchStatus("R.a1.t"          ,1);  tree->SetBranchAddress("R.a1.t"          , R_a1_t          ); 
  tree->SetBranchStatus("R.a1.t_c"        ,1);  tree->SetBranchAddress("R.a1.t_c"        , R_a1_t_c        ); 
  tree->SetBranchStatus("R.a1.asum_p"     ,1);  tree->SetBranchAddress("R.a1.asum_p"     ,&R_a1_asum_p     ); 
  tree->SetBranchStatus("R.a1.asum_c"     ,1);  tree->SetBranchAddress("R.a1.asum_c"     ,&R_a1_asum_c     ); 
  tree->SetBranchStatus("R.a1.trpath"     ,1);  tree->SetBranchAddress("R.a1.trpath"     , R_a1_trpath     ); 
  tree->SetBranchStatus("R.a1.trx"        ,1);  tree->SetBranchAddress("R.a1.trx"        , R_a1_trx        ); 
  tree->SetBranchStatus("R.a1.try"        ,1);  tree->SetBranchAddress("R.a1.try"        , R_a1_try        ); 
  tree->SetBranchStatus("R.a1.nhits"      ,1);  tree->SetBranchAddress("R.a1.nhits"      , R_a1_nhits      ); 
  tree->SetBranchStatus("R.a1.peak"       ,1);  tree->SetBranchAddress("R.a1.peak"       , R_a1_peak       ); 
  tree->SetBranchStatus("R.a1.t_fadc"     ,1);  tree->SetBranchAddress("R.a1.t_fadc"     , R_a1_t_fadc     ); 
  tree->SetBranchStatus("R.a1.tc_fadc"    ,1);  tree->SetBranchAddress("R.a1.tc_fadc"    , R_a1_tc_fadc    ); 
  tree->SetBranchStatus("R.a1.nbadped"    ,0);  tree->SetBranchAddress("R.a1.nbadped"    , R_a1_nbadped    ); 
  tree->SetBranchStatus("R.a1.noverflow"  ,0);  tree->SetBranchAddress("R.a1.noverflow"  , R_a1_noverflow  ); 
  tree->SetBranchStatus("R.a1.nunderflow" ,0);  tree->SetBranchAddress("R.a1.nunderflow" , R_a1_nunderflow ); 

  tree->SetBranchStatus("R.a2.a"          ,1);  tree->SetBranchAddress("R.a2.a"          , R_a2_a          ); 
  tree->SetBranchStatus("R.a2.a_c"        ,1);  tree->SetBranchAddress("R.a2.a_c"        , R_a2_a_c        ); 
  tree->SetBranchStatus("R.a2.a_p"        ,1);  tree->SetBranchAddress("R.a2.a_p"        , R_a2_a_p        ); 
  tree->SetBranchStatus("R.a2.t"          ,1);  tree->SetBranchAddress("R.a2.t"          , R_a2_t          ); 
  tree->SetBranchStatus("R.a2.t_c"        ,1);  tree->SetBranchAddress("R.a2.t_c"        , R_a2_t_c        ); 
  tree->SetBranchStatus("R.a2.asum_p"     ,1);  tree->SetBranchAddress("R.a2.asum_p"     ,&R_a2_asum_p     ); 
  tree->SetBranchStatus("R.a2.asum_c"     ,1);  tree->SetBranchAddress("R.a2.asum_c"     ,&R_a2_asum_c     ); 
  tree->SetBranchStatus("R.a2.trpath"     ,1);  tree->SetBranchAddress("R.a2.trpath"     , R_a2_trpath     ); 
  tree->SetBranchStatus("R.a2.trx"        ,1);  tree->SetBranchAddress("R.a2.trx"        , R_a2_trx        ); 
  tree->SetBranchStatus("R.a2.try"        ,1);  tree->SetBranchAddress("R.a2.try"        , R_a2_try        ); 
  tree->SetBranchStatus("R.a2.nhits"      ,1);  tree->SetBranchAddress("R.a2.nhits"      , R_a2_nhits      ); 
  tree->SetBranchStatus("R.a2.peak"       ,1);  tree->SetBranchAddress("R.a2.peak"       , R_a2_peak       ); 
  tree->SetBranchStatus("R.a2.t_fadc"     ,1);  tree->SetBranchAddress("R.a2.t_fadc"     , R_a2_t_fadc     ); 
  tree->SetBranchStatus("R.a2.tc_fadc"    ,1);  tree->SetBranchAddress("R.a2.tc_fadc"    , R_a2_tc_fadc    ); 
  tree->SetBranchStatus("R.a2.nbadped"    ,0);  tree->SetBranchAddress("R.a2.nbadped"    , R_a2_nbadped    ); 
  tree->SetBranchStatus("R.a2.noverflow"  ,0);  tree->SetBranchAddress("R.a2.noverflow"  , R_a2_noverflow  ); 
  tree->SetBranchStatus("R.a2.nunderflow" ,0);  tree->SetBranchAddress("R.a2.nunderflow" , R_a2_nunderflow ); 

  tree->SetBranchStatus("R.cer.a"         ,1);  tree->SetBranchAddress("R.cer.a"         , R_cer_a         ); 
  tree->SetBranchStatus("R.cer.a_c"       ,1);  tree->SetBranchAddress("R.cer.a_c"       , R_cer_a_c       ); 
  tree->SetBranchStatus("R.cer.a_p"       ,1);  tree->SetBranchAddress("R.cer.a_p"       , R_cer_a_p       ); 
  tree->SetBranchStatus("R.cer.t"         ,1);  tree->SetBranchAddress("R.cer.t"         , R_cer_t         ); 
  tree->SetBranchStatus("R.cer.t_c"       ,1);  tree->SetBranchAddress("R.cer.t_c"       , R_cer_t_c       ); 
  tree->SetBranchStatus("R.cer.asum_p"    ,1);  tree->SetBranchAddress("R.cer.asum_p"    ,&R_cer_asum_p    ); 
  tree->SetBranchStatus("R.cer.asum_c"    ,1);  tree->SetBranchAddress("R.cer.asum_c"    ,&R_cer_asum_c    ); 
  tree->SetBranchStatus("R.cer.trpath"    ,1);  tree->SetBranchAddress("R.cer.trpath"    , R_cer_trpath    ); 
  tree->SetBranchStatus("R.cer.trx"       ,1);  tree->SetBranchAddress("R.cer.trx"       , R_cer_trx       ); 
  tree->SetBranchStatus("R.cer.try"       ,1);  tree->SetBranchAddress("R.cer.try"       , R_cer_try       ); 
  tree->SetBranchStatus("R.cer.nhits"     ,1);  tree->SetBranchAddress("R.cer.nhits"     , R_cer_nhits     ); 
  tree->SetBranchStatus("R.cer.peak"      ,1);  tree->SetBranchAddress("R.cer.peak"      , R_cer_peak      ); 
  tree->SetBranchStatus("R.cer.t_fadc"    ,1);  tree->SetBranchAddress("R.cer.t_fadc"    , R_cer_t_fadc    ); 
  tree->SetBranchStatus("R.cer.tc_fadc"   ,1);  tree->SetBranchAddress("R.cer.tc_fadc"   , R_cer_tc_fadc   ); 
  tree->SetBranchStatus("R.cer.nbadped"   ,0);  tree->SetBranchAddress("R.cer.nbadped"   , R_cer_nbadped   ); 
  tree->SetBranchStatus("R.cer.noverflow" ,0);  tree->SetBranchAddress("R.cer.noverflow" , R_cer_noverflow ); 
  tree->SetBranchStatus("R.cer.nunderflow",0);  tree->SetBranchAddress("R.cer.nunderflow", R_cer_nunderflow); 

  tree->SetBranchStatus("R.ps.a"          ,1);  tree->SetBranchAddress("R.ps.a"          , R_ps_a          );
  tree->SetBranchStatus("R.ps.a_c"        ,1);  tree->SetBranchAddress("R.ps.a_c"        , R_ps_a_c        );
  tree->SetBranchStatus("R.ps.a_p"        ,1);  tree->SetBranchAddress("R.ps.a_p"        , R_ps_a_p        );
  tree->SetBranchStatus("R.ps.t"          ,1);  tree->SetBranchAddress("R.ps.t"          , R_ps_t          );
  tree->SetBranchStatus("R.ps.t_c"        ,1);  tree->SetBranchAddress("R.ps.t_c"        , R_ps_t_c        );
  tree->SetBranchStatus("R.ps.asum_p"     ,1);  tree->SetBranchAddress("R.ps.asum_p"     ,&R_ps_asum_p     ); 
  tree->SetBranchStatus("R.ps.asum_c"     ,1);  tree->SetBranchAddress("R.ps.asum_c"     ,&R_ps_asum_c     ); 

  tree->SetBranchStatus("R.sh.a"          ,1);  tree->SetBranchAddress("R.sh.a"          , R_sh_a          );
  tree->SetBranchStatus("R.sh.a_c"        ,1);  tree->SetBranchAddress("R.sh.a_c"        , R_sh_a_c        );
  tree->SetBranchStatus("R.sh.a_p"        ,1);  tree->SetBranchAddress("R.sh.a_p"        , R_sh_a_p        );
  tree->SetBranchStatus("R.sh.t"          ,1);  tree->SetBranchAddress("R.sh.t"          , R_sh_t          );
  tree->SetBranchStatus("R.sh.t_c"        ,1);  tree->SetBranchAddress("R.sh.t_c"        , R_sh_t_c        );
  tree->SetBranchStatus("R.sh.asum_p"     ,1);  tree->SetBranchAddress("R.sh.asum_p"     ,&R_sh_asum_p     ); 
  tree->SetBranchStatus("R.sh.asum_c"     ,1);  tree->SetBranchAddress("R.sh.asum_c"     ,&R_sh_asum_c     ); 

  tree->SetBranchStatus("R.tr.n"          ,1);  tree->SetBranchAddress("R.tr.n"          ,&R_tr_n          );
  tree->SetBranchStatus("R.tr.flag"       ,1);  tree->SetBranchAddress("R.tr.flag"       , R_tr_flag       );
  tree->SetBranchStatus("R.tr.ndof"       ,1);  tree->SetBranchAddress("R.tr.ndof"       , R_tr_ndof       );
  tree->SetBranchStatus("R.tr.chi2"       ,1);  tree->SetBranchAddress("R.tr.chi2"       , R_tr_chi2       );
  tree->SetBranchStatus("R.tr.beta"       ,1);  tree->SetBranchAddress("R.tr.beta"       , R_tr_beta       );
  tree->SetBranchStatus("R.tr.d_x"        ,1);  tree->SetBranchAddress("R.tr.d_x"        , R_tr_d_x        );
  tree->SetBranchStatus("R.tr.d_y"        ,1);  tree->SetBranchAddress("R.tr.d_y"        , R_tr_d_y        );
  tree->SetBranchStatus("R.tr.d_th"       ,1);  tree->SetBranchAddress("R.tr.d_th"       , R_tr_d_th       );
  tree->SetBranchStatus("R.tr.d_ph"       ,1);  tree->SetBranchAddress("R.tr.d_ph"       , R_tr_d_ph       );
  tree->SetBranchStatus("R.tr.r_x"        ,1);  tree->SetBranchAddress("R.tr.r_x"        , R_tr_r_x        );
  tree->SetBranchStatus("R.tr.r_y"        ,1);  tree->SetBranchAddress("R.tr.r_y"        , R_tr_r_y        );
  tree->SetBranchStatus("R.tr.r_th"       ,1);  tree->SetBranchAddress("R.tr.r_th"       , R_tr_r_th       );
  tree->SetBranchStatus("R.tr.r_ph"       ,1);  tree->SetBranchAddress("R.tr.r_ph"       , R_tr_r_ph       );
  tree->SetBranchStatus("R.tr.x"          ,1);  tree->SetBranchAddress("R.tr.x"          , R_tr_x          );
  tree->SetBranchStatus("R.tr.y"          ,1);  tree->SetBranchAddress("R.tr.y"          , R_tr_y          );
  tree->SetBranchStatus("R.tr.th"         ,1);  tree->SetBranchAddress("R.tr.th"         , R_tr_th         );
  tree->SetBranchStatus("R.tr.ph"         ,1);  tree->SetBranchAddress("R.tr.ph"         , R_tr_ph         );
  tree->SetBranchStatus("R.tr.time"       ,1);  tree->SetBranchAddress("R.tr.time"       , R_tr_time       );
  tree->SetBranchStatus("R.tr.p"          ,1);  tree->SetBranchAddress("R.tr.p"          , R_tr_p          );
  tree->SetBranchStatus("R.tr.pathl"      ,1);  tree->SetBranchAddress("R.tr.pathl"      , R_tr_pathl      );
  tree->SetBranchStatus("R.tr.px"         ,1);  tree->SetBranchAddress("R.tr.px"         , R_tr_px         );
  tree->SetBranchStatus("R.tr.py"         ,1);  tree->SetBranchAddress("R.tr.py"         , R_tr_py         );
  tree->SetBranchStatus("R.tr.pz"         ,1);  tree->SetBranchAddress("R.tr.pz"         , R_tr_pz         );
  tree->SetBranchStatus("R.tr.tg_dp"      ,1);  tree->SetBranchAddress("R.tr.tg_dp"      , R_tr_tg_dp      );
  tree->SetBranchStatus("R.tr.tg_y"       ,1);  tree->SetBranchAddress("R.tr.tg_y"       , R_tr_tg_y       );
  tree->SetBranchStatus("R.tr.tg_th"      ,1);  tree->SetBranchAddress("R.tr.tg_th"      , R_tr_tg_th      );
  tree->SetBranchStatus("R.tr.tg_ph"      ,1);  tree->SetBranchAddress("R.tr.tg_ph"      , R_tr_tg_ph      );
  tree->SetBranchStatus("R.tr.vx"         ,1);  tree->SetBranchAddress("R.tr.vx"         , R_tr_vx         );
  tree->SetBranchStatus("R.tr.vy"         ,1);  tree->SetBranchAddress("R.tr.vy"         , R_tr_vy         );
  tree->SetBranchStatus("R.tr.vz"         ,1);  tree->SetBranchAddress("R.tr.vz"         , R_tr_vz         );

  tree->SetBranchStatus("DR.rtRFtime"     ,1);  tree->SetBranchAddress("DR.rtRFtime"     , R_rtRFtime      );
  tree->SetBranchStatus("RTDC.F1FirstHit" ,1);  tree->SetBranchAddress("RTDC.F1FirstHit" , R_F1Fhit        );

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Tree::readtree_COIN(){

}

