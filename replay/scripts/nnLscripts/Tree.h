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

#define MAX 10  // Maximum No. of Tracks
#define RS0 1   // No. of Segments of R-S0
#define RS2 16  // No. of Segments of R-S2
#define RA1 24  // No. of Segments of R-AC1
#define RA2 26  // No. of Segments of R-AC2
#define RCR 10  // No. of Segments of R-GC
#define RPS 48  // No. of Segments of R-Pre-Shower
#define RSH 75  // No. of Segments of R-Shower
#define LS0 1   // No. of Segments of L-S0
#define LS2 16  // No. of Segments of L-S2

class Tree
{
public:
  Tree();
  ~Tree();

public:
////////////
// Common //
////////////
// Beam Position
  double rbax, rbay, rbbx, rbby;  // BPM
  double rbx, rby;  // from raster
  double bpmaws, bpmbws;  // Sum of BPM raw cur

///////////////////
// HRS Right Arm //
///////////////////
//// Event Flag ////
  double DR_T1, DR_T2, DR_T3, DR_T4, DR_T5, DR_T6, DR_T7, DR_T8;

//// S0 ////
  double R_s0_da[RS0][RS0], R_s0_da_c[RS0], R_s0_da_p[RS0], R_s0_dedx[RS0];
  double R_s0_dt[RS0], R_s0_dt_c[RS0], R_s0_dtime[RS0];
  double R_s0_la[RS0], R_s0_la_c[RS0], R_s0_la_p[RS0];
  double R_s0_lbadped[RS0], R_s0_lnhits[RS0], R_s0_loverflow[RS0], R_s0_lpeak[RS0];
  double R_s0_lt[RS0], R_s0_lt_c[RS0], R_s0_lt_fadc[RS0], R_s0_ltc_fadc[RS0], R_s0_lunderflow[RS0];
  double R_s0_ra[RS0], R_s0_ra_c[RS0], R_s0_ra_p[RS0];
  double R_s0_rbadped[RS0], R_s0_rnhits[RS0], R_s0_roverflow[RS0], R_s0_rpeak[RS0];
  double R_s0_rt[RS0], R_s0_rt_c[RS0], R_s0_rt_fadc[RS0], R_s0_rtc_fadc[RS0], R_s0_runderflow[RS0];
  double R_s0_t_pads[MAX];
  double R_s0_time[RS0];
  double R_s0_trdy[MAX];
  double R_s0_troff[RS0];
  double R_s0_trpad[MAX], R_s0_trpath[MAX], R_s0_trx[MAX], R_s0_try[MAX];
  double R_s0_ua[RS0], R_s0_ua_c[RS0], R_s0_ua_p[RS0];
  double R_s0_ut[RS0], R_s0_ut_c[RS0], R_s0_x_adc[RS0], R_s0_x_t[RS0];

//// S2 ////
  double R_s2_dedx[RS2], R_s2_dtime[RS2];
  double R_s2_la[RS2], R_s2_la_c[RS2], R_s2_la_p[RS2];
  double R_s2_lbadped[RS2], R_s2_lnhits[RS2], R_s2_loverflow[RS2], R_s2_lpeak[RS2];
  double R_s2_lt[RS2], R_s2_lt_c[RS2], R_s2_lt_fadc[RS2], R_s2_ltc_fadc[RS2], R_s2_lunderflow[RS2];
  double R_s2_ra[RS2], R_s2_ra_c[RS2], R_s2_ra_p[RS2];
  double R_s2_rbadped[RS2], R_s2_rnhits[RS2], R_s2_roverflow[RS2], R_s2_rpeak[RS2];
  double R_s2_rt[RS2], R_s2_rt_c[RS2], R_s2_rt_fadc[RS2], R_s2_rtc_fadc[RS2], R_s2_runderflow[RS2];
  double R_s2_t_pads[MAX];
  double R_s2_time[RS2];
  double R_s2_trdx[MAX];
  double R_s2_troff[RS2];
  double R_s2_trpad[MAX], R_s2_trpath[MAX], R_s2_trx[MAX], R_s2_try[MAX];
  double R_s2_y_adc[RS2], R_s2_y_t[RS2];

//// AC ////
  double R_a1_a[RA1], R_a1_a_c[RA1], R_a1_a_p[RA1];
  double R_a1_nbadped[RA1], R_a1_nhits[RA1], R_a1_noverflow[RA1], R_a1_nunderflow[RA1], R_a1_peak[RA1];
  double R_a1_t[RA1], R_a1_t_c[RA1], R_a1_t_fadc[RA1], R_a1_tc_fadc[RA1];
  double R_a1_trpath[MAX], R_a1_trx[MAX], R_a1_try[MAX];
  double R_a2_a[RA2], R_a2_a_c[RA2], R_a2_a_p[RA2];
  double R_a2_nbadped[RA2], R_a2_nhits[RA2], R_a2_noverflow[RA2], R_a2_nunderflow[RA2], R_a2_peak[RA2];
  double R_a2_t[RA2], R_a2_t_c[RA2], R_a2_t_fadc[RA2], R_a2_tc_fadc[RA2];
  double R_a2_trpath[MAX], R_a2_trx[MAX], R_a2_try[MAX];

//// Pre-shower ////
  double R_ps_a[RPS], R_ps_a_c[RPS], R_ps_a_p[RPS];
  double R_ps_eblk[6], R_ps_nblk[6];
  double R_ps_trpath[MAX], R_ps_trx[MAX], R_ps_try[MAX];

//// Gas Cherenkov ////
  double R_cer_a[RCR], R_cer_a_c[RCR], R_cer_a_p[RCR];
  double R_cer_nbadped[RCR], R_cer_nhits[RCR], R_cer_noverflow[RCR], R_cer_nunderflow[RCR], R_cer_peak[RCR];
  double R_cer_t[RCR], R_cer_t_c[RCR], R_cer_t_fadc[RCR], R_cer_tc_fadc[RCR];
  double R_cer_trpath[MAX], R_cer_trx[MAX], R_cer_try[MAX];

//// Shower ////
  double R_sh_a[RSH], R_sh_a_c[RSH], R_sh_a_p[RSH];
  double R_sh_eblk[9], R_sh_nblk[9];
  double R_sh_trpath[MAX], R_sh_trx[MAX], R_sh_try[MAX];

//// Tracking ////
  double R_tr_beta[MAX];
  double R_tr_chi2[MAX];
  double R_tr_d_ph[MAX], R_tr_d_th[MAX], R_tr_d_x[MAX], R_tr_d_y[MAX];
  double R_tr_dbeta[MAX], R_tr_dtime[MAX];
  double R_tr_flag[MAX], R_tr_ndof[MAX];
  double R_tr_p[MAX], R_tr_pathl[MAX], R_tr_ph[MAX], R_tr_px[MAX], R_tr_py[MAX], R_tr_pz[MAX];
  double R_tr_r_ph[MAX], R_tr_r_th[MAX], R_tr_r_x[MAX], R_tr_r_y[MAX];
  double R_tr_tg_dp[MAX], R_tr_tg_ph[MAX], R_tr_tg_th[MAX], R_tr_tg_y[MAX];
  double R_tr_th[MAX], R_tr_time[MAX];
  double R_tr_vx[MAX], R_tr_vy[MAX], R_tr_vz[MAX];
  double R_tr_x[MAX], R_tr_y[MAX];

///////////////////
// HRS Light Arm //
///////////////////
//// S0 ////
  double L_s0_da[LS0], L_s0_da_c[LS0], L_s0_da_p[LS0], L_s0_dedx[LS0];
  double L_s0_dt[LS0], L_s0_dt_c[LS0], L_s0_dtime[LS0];
  double L_s0_la[LS0], L_s0_la_c[LS0], L_s0_la_p[LS0];
  double L_s0_lbadped[LS0], L_s0_lnhits[LS0], L_s0_loverflow[LS0], L_s0_lpeak[LS0];
  double L_s0_lt[LS0], L_s0_lt_c[LS0], L_s0_lt_fadc[LS0], L_s0_ltc_fadc[LS0], L_s0_lunderflow[LS0];
  double L_s0_ra[LS0], L_s0_ra_c[LS0], L_s0_ra_p[LS0];
  double L_s0_rbadped[LS0], L_s0_rnhits[LS0], L_s0_roverflow[LS0], L_s0_rpeak[LS0];
  double L_s0_rt[LS0], L_s0_rt_c[LS0], L_s0_rt_fadc[LS0], L_s0_rtc_fadc[LS0], L_s0_runderflow[LS0]; 
  double L_s0_t_pads[MAX];
  double L_s0_time[LS0];
  double L_s0_trdy[MAX];
  double L_s0_troff[LS0];
  double L_s0_trpad[MAX], L_s0_trpath[MAX], L_s0_trx[MAX], L_s0_try[MAX];
  double L_s0_ua[LS0], L_s0_ua_c[LS0], L_s0_ua_p[LS0];
  double L_s0_ut[LS0], L_s0_ut_c[LS0], L_s0_x_adc[LS0], L_s0_x_t[LS0];

//// S2 ////
  double L_s2_dedx[LS2], L_s2_dtime[LS2];
  double L_s2_la[LS2], L_s2_la_c[LS2], L_s2_la_p[LS2];
  double L_s2_lbadped[LS2], L_s2_lnhits[LS2], L_s2_loverflow[LS2], L_s2_lpeak[LS2];
  double L_s2_lt[LS2], L_s2_lt_c[LS2], L_s2_lt_fadc[LS2], L_s2_ltc_fadc[LS2], L_s2_lunderflow[LS2];
  double L_s2_ra[LS2], L_s2_ra_c[LS2], L_s2_ra_p[LS2];
  double L_s2_rbadped[LS2], L_s2_rnhits[LS2], L_s2_roverflow[LS2], L_s2_rpeak[LS2];
  double L_s2_rt[LS2], L_s2_rt_c[LS2], L_s2_rt_fadc[LS2], L_s2_rtc_fadc[LS2], L_s2_runderflow[LS2];
  double L_s2_t_pads[MAX];
  double L_s2_time[LS2];
  double L_s2_trdx[MAX];
  double L_s2_troff[LS2];
  double L_s2_trpad[MAX], L_s2_trpath[MAX], L_s2_trx[MAX], L_s2_try[MAX];
  double L_s2_y_adc[LS2], L_s2_y_t[LS2];

// Tracking
  double L_tr_beta[MAX];
  double L_tr_chi2[MAX];
  double L_tr_d_ph[MAX], L_tr_d_th[MAX], L_tr_d_x[MAX], L_tr_d_y[MAX];
  double L_tr_dbeta[MAX], L_tr_dtime[MAX];
  double L_tr_flag[MAX], L_tr_ndof[MAX];
  double L_tr_p[MAX], L_tr_pathl[MAX], L_tr_ph[MAX], L_tr_px[MAX], L_tr_py[MAX], L_tr_pz[MAX];
  double L_tr_r_ph[MAX], L_tr_r_th[MAX], L_tr_r_x[MAX], L_tr_r_y[MAX];
  double L_tr_tg_dp[MAX], L_tr_tg_ph[MAX], L_tr_tg_th[MAX], L_tr_tg_y[MAX];
  double L_tr_th[MAX], L_tr_time[MAX];
  double L_tr_vx[MAX], L_tr_vy[MAX], L_tr_vz[MAX];
  double L_tr_x[MAX], L_tr_y[MAX];

public:
  TChain *tree;

  void readtree_COMN();
  void readtree_RHRS();
  void readtree_LHRS();
  void readtree_COIN();
  void merge(string ifname);
  int GetEntries()    const { return tree->GetEntries(); }

private:

};

#endif
