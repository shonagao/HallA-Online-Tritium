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

#include "ParamMan.h"

#define MAX 10     // Maximum No. of Tracks
#define RS0 1      // No. of Segments of R-S0
#define RS2 16     // No. of Segments of R-S2
#define RA1 24     // No. of Segments of R-AC1
#define RA2 26     // No. of Segments of R-AC2
#define RCR 10     // No. of Segments of R-GC
#define RPS 48     // No. of Segments of R-Pre-Shower
#define RSH 75     // No. of Segments of R-Shower
#define RF1TDC 64  // No. of ch of R-F1TDC
#define LS0 1      // No. of Segments of L-S0
#define LS2 16     // No. of Segments of L-S2
#define LF1TDC 64  // No. of ch of L-F1TDC

class Tree
{
public:
  Tree();
  ~Tree();

//private:
public:

  int ENum;
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
  double DR_evtypebits;
  double DR_T1, DR_T2, DR_T3, DR_T4, DR_T5, DR_T6, DR_T7, DR_T8;

//// S0 ////
  double R_s0_la[RS0], R_s0_la_c[RS0], R_s0_la_p[RS0];                 // Bottom ADC, decADC, ADC(ped-supp)
  double R_s0_lt[RS0], R_s0_lt_c[RS0];                                 // Bottom TDC, decTDC
  double R_s0_ra[RS0], R_s0_ra_c[RS0], R_s0_ra_p[RS0];                 // Top ADC, decADC, ADC(ped-supp)
  double R_s0_rt[RS0], R_s0_rt_c[RS0];                                 // Top TDC, decTDC
  double R_s0_t_pads[MAX];                                             // Paddle T&B TDC coin
  double R_s0_time[RS0];                                               // Time of hit (sec)
  double R_s0_dedx[RS0];                                               // dE/dx (no angle correction)
  double R_s0_trdy[MAX], R_s0_troff[RS0];                              // Track Deviation (m), Trigger offset
  double R_s0_trpad[MAX], R_s0_trpath[MAX];                            // hit-pad associated with track, track path length
  double R_s0_trx[MAX], R_s0_try[MAX];                                 // track x-pos, track y-pos
// FADC
  double R_s0_lnhits[RS0], R_s0_lpeak[RS0];                            // Hits of Bottom, FADC peak val
  double R_s0_lt_fadc[RS0], R_s0_ltc_fadc[RS0];                        // FADC TDC, FADC decTDC
  double R_s0_lbadped[RS0], R_s0_lunderflow[RS0], R_s0_loverflow[RS0]; // FADC ped-quality, underflow, overflow bit
  double R_s0_rnhits[RS0], R_s0_rpeak[RS0];                            // Hits of Bottom, FADC peak val
  double R_s0_rt_fadc[RS0], R_s0_rtc_fadc[RS0];                        // FADC TDC, FADC decTDC
  double R_s0_rbadped[RS0], R_s0_runderflow[RS0], R_s0_roverflow[RS0]; // FADC ped-quality, underflow, overflow bit

//// S2 ////
  double R_s2_la[RS2], R_s2_la_c[RS2], R_s2_la_p[RS2];                 // Bottom ADC, decADC, ADC(ped-supp)
  double R_s2_lt[RS2], R_s2_lt_c[RS2];                                 // Bottom TDC, decTDC
  double R_s2_ra[RS2], R_s2_ra_c[RS2], R_s2_ra_p[RS2];                 // Top ADC, decADC, ADC(ped-supp)
  double R_s2_rt[RS2], R_s2_rt_c[RS2];                                 // Top TDC, decTDC
  double R_s2_t_pads[RS2];                                             // Paddle T&B TDC coin
  double R_s2_time[RS2];                                               // Time of hit (sec)
  double R_s2_dedx[RS2];                                               // dE/dx (no angle correction)
  double R_s2_trdx[MAX], R_s2_troff[RS2];                              // Track Deviation (m), Trigger offset
  double R_s2_trpad[MAX], R_s2_trpath[MAX];                            // hit-pad associated with track, track path length
  double R_s2_trx[MAX], R_s2_try[MAX];                                 // track x-pos, track y-pos
// FADC
  double R_s2_lnhits[RS2], R_s2_lpeak[RS2];                            // Hits of Bottom, FADC peak val
  double R_s2_lt_fadc[RS2], R_s2_ltc_fadc[RS2];                        // FADC TDC, FADC decTDC
  double R_s2_lbadped[RS2], R_s2_lunderflow[RS2], R_s2_loverflow[RS2]; // FADC ped-quality, underflow, overflow bit
  double R_s2_rnhits[RS2], R_s2_rpeak[RS2];                            // Hits of Bottom, FADC peak val
  double R_s2_rt_fadc[RS2], R_s2_rtc_fadc[RS2];                        // FADC TDC, FADC decTDC
  double R_s2_rbadped[RS2], R_s2_runderflow[RS2], R_s2_roverflow[RS2]; // FADC ped-quality, underflow, overflow bit

//// AC ////
  double R_a1_a[RA1], R_a1_a_c[RA1], R_a1_a_p[RA1];                     // ADC, decADC, ADC(ped-supp)
  double R_a1_t[RA1], R_a1_t_c[RA1];                                    // TDC, decTDC
  double R_a1_asum_p, R_a1_asum_c;                                      // sumADC(ped-supp), sum-decADC
  double R_a1_trpath[MAX], R_a1_trx[MAX], R_a1_try[MAX];                // track path length, track x-pos, track y-pos
  double R_a1_nhits[RA1], R_a1_peak[RA1];                               // Hits, FADC peak val
  double R_a1_t_fadc[RA1], R_a1_tc_fadc[RA1];                           // FADC TDC, FADC decTDC
  double R_a1_nbadped[RA1], R_a1_nunderflow[RA1], R_a1_noverflow[RA1];  // ped-quality, underflow, overflow bit

  double R_a2_a[RA2], R_a2_a_c[RA2], R_a2_a_p[RA2];                     // ADC, decADC, ADC(ped-supp)
  double R_a2_t[RA2], R_a2_t_c[RA2];                                    // TDC, decTDC
  double R_a2_asum_p, R_a2_asum_c;                                      // sumADC(ped-supp), sum-decADC
  double R_a2_trpath[MAX], R_a2_trx[MAX], R_a2_try[MAX];                // track path length, track x-pos, track y-pos
  double R_a2_nhits[RA2], R_a2_peak[RA2];                               // Hits, FADC peak val
  double R_a2_t_fadc[RA2], R_a2_tc_fadc[RA2];                           // FADC TDC, FADC decTDC
  double R_a2_nbadped[RA2], R_a2_nunderflow[RA2], R_a2_noverflow[RA2];  // ped-quality, underflow, overflow bit

//// Gas Cherenkov ////
  double R_cer_a[RCR], R_cer_a_c[RCR], R_cer_a_p[RCR];                    // ADC, decADC, ADC(ped-supp)
  double R_cer_t[RCR], R_cer_t_c[RCR];                                    // TDC, decTDC
  double R_cer_asum_p, R_cer_asum_c;                                      // sumADC(ped-supp), sum-decADC
  double R_cer_trpath[MAX], R_cer_trx[MAX], R_cer_try[MAX];               // track path length, track x-pos, track y-pos
  double R_cer_nhits[RCR], R_cer_peak[RCR];                               // Hits, FADC peak val
  double R_cer_t_fadc[RCR], R_cer_tc_fadc[RCR];                           // FADC TDC, FADC decTDC
  double R_cer_nbadped[RCR], R_cer_nunderflow[RCR], R_cer_noverflow[RCR]; // ped-quality, underflow, overflow bit

//// Pre-shower ////
  double R_ps_a[RCR], R_ps_a_c[RCR], R_ps_a_p[RCR];                     // ADC, decADC, ADC(ped-supp)
  double R_ps_t[RCR], R_ps_t_c[RCR];                                    // TDC, decTDC
  double R_ps_asum_p, R_ps_asum_c;                                      // sumADC(ped-supp), sum-decADC

//// Shower ////
  double R_sh_a[RCR], R_sh_a_c[RCR], R_sh_a_p[RCR];                     // ADC, decADC, ADC(ped-supp)
  double R_sh_t[RCR], R_sh_t_c[RCR];                                    // TDC, decTDC
  double R_sh_asum_p, R_sh_asum_c;                                      // sumADC(ped-supp), sum-decADC

//// Tracking ////
  double R_tr_n;                                                                 // No. of Tracks
  double R_tr_flag[MAX], R_tr_ndof[MAX];                                         // track status, track NDoF
  double R_tr_chi2[MAX];                                                         // track chi2
  double R_tr_beta[MAX];                                                         // beta of track
  double R_tr_d_x[MAX], R_tr_d_y[MAX], R_tr_d_th[MAX], R_tr_d_ph[MAX];           // x, y, theta, phi at Detector
  double R_tr_r_x[MAX], R_tr_r_y[MAX], R_tr_r_th[MAX], R_tr_r_ph[MAX];           // x, y, theta, phi at Rot-Coordinate
  double R_tr_x[MAX], R_tr_y[MAX], R_tr_th[MAX], R_tr_ph[MAX];                   // x, y, theta, phi
  double R_tr_time[MAX];                                                         // theta, time@RF
  double R_tr_p[MAX], R_tr_pathl[MAX], R_tr_px[MAX], R_tr_py[MAX], R_tr_pz[MAX]; // mom(unit 10GeV?), length(TtoP), momx, momy, momz
  double R_tr_tg_dp[MAX], R_tr_tg_y[MAX], R_tr_tg_th[MAX], R_tr_tg_ph[MAX];      // dp, y, theta, phi at target
  double R_tr_vx[MAX], R_tr_vy[MAX], R_tr_vz[MAX]; // vertex X, Y, Z

  double R_rtRFtime[6];

///////////////////
// HRS Light Arm //
///////////////////
//// S0 ////
  double L_s0_la[LS0], L_s0_la_c[LS0], L_s0_la_p[LS0];                 // Bottom ADC, decADC, ADC(ped-supp)
  double L_s0_lt[LS0], L_s0_lt_c[LS0];                                 // Bottom TDC, decTDC
  double L_s0_ra[LS0], L_s0_ra_c[LS0], L_s0_ra_p[LS0];                 // Top ADC, decADC, ADC(ped-supp)
  double L_s0_rt[LS0], L_s0_rt_c[LS0];                                 // Top TDC, decTDC
  double L_s0_t_pads[MAX];                                             // Paddle T&B TDC coin
  double L_s0_time[LS0];                                               // Time of hit (sec)
  double L_s0_dedx[LS0];                                               // dE/dx (no angle correction)
  double L_s0_trdy[MAX], L_s0_troff[LS0];                              // Track Deviation (m), Trigger offset
  double L_s0_trpad[MAX], L_s0_trpath[MAX];                            // hit-pad associated with track, track path length
  double L_s0_trx[MAX], L_s0_try[MAX];                                 // track x-pos, track y-pos
// FADC
  double L_s0_lnhits[LS0], L_s0_lpeak[LS0];                            // Hits of Bottom, FADC peak val
  double L_s0_lt_fadc[LS0], L_s0_ltc_fadc[LS0];                        // FADC TDC, FADC decTDC
  double L_s0_lbadped[LS0], L_s0_lunderflow[LS0], L_s0_loverflow[LS0]; // FADC ped-quality, underflow, overflow bit
  double L_s0_rnhits[LS0], L_s0_rpeak[LS0];                            // Hits of Bottom, FADC peak val
  double L_s0_rt_fadc[LS0], L_s0_rtc_fadc[LS0];                        // FADC TDC, FADC decTDC
  double L_s0_rbadped[LS0], L_s0_runderflow[LS0], L_s0_roverflow[LS0]; // FADC ped-quality, underflow, overflow bit

//// S2 ////
  double L_s2_la[LS2], L_s2_la_c[LS2], L_s2_la_p[LS2];                 // Bottom ADC, decADC, ADC(ped-supp)
  double L_s2_lt[LS2], L_s2_lt_c[LS2];                                 // Bottom TDC, decTDC
  double L_s2_ra[LS2], L_s2_ra_c[LS2], L_s2_ra_p[LS2];                 // Top ADC, decADC, ADC(ped-supp)
  double L_s2_rt[LS2], L_s2_rt_c[LS2];                                 // Top TDC, decTDC
  double L_s2_t_pads[LS2];                                             // Paddle T&B TDC coin
  double L_s2_time[LS2];                                               // Time of hit (sec)
  double L_s2_dedx[LS2];                                               // dE/dx (no angle correction)
  double L_s2_trdx[MAX], L_s2_troff[LS2];                               // Track Deviation (m), Trigger offset
  double L_s2_trpad[MAX], L_s2_trpath[MAX];                            // hit-pad associated with track, track path length
  double L_s2_trx[MAX], L_s2_try[MAX];                                 // track x-pos, track y-pos
// FADC
  double L_s2_lnhits[LS2], L_s2_lpeak[LS2];                            // Hits of Bottom, FADC peak val
  double L_s2_lt_fadc[LS2], L_s2_ltc_fadc[LS2];                        // FADC TDC, FADC decTDC
  double L_s2_lbadped[LS2], L_s2_lunderflow[LS2], L_s2_loverflow[LS2]; // FADC ped-quality, underflow, overflow bit
  double L_s2_rnhits[LS2], L_s2_rpeak[LS2];                            // Hits of Bottom, FADC peak val
  double L_s2_rt_fadc[LS2], L_s2_rtc_fadc[LS2];                        // FADC TDC, FADC decTDC
  double L_s2_rbadped[LS2], L_s2_runderflow[LS2], L_s2_roverflow[LS2]; // FADC ped-quality, underflow, overflow bit

//// Tracking ////
  double L_tr_n;                                                                 // No. of Tracks
  double L_tr_flag[MAX], L_tr_ndof[MAX];                                         // track status, track NDoF
  double L_tr_chi2[MAX];                                                         // track chi2
  double L_tr_beta[MAX];                                                         // beta of track
  double L_tr_d_x[MAX], L_tr_d_y[MAX], L_tr_d_th[MAX], L_tr_d_ph[MAX];           // x, y, theta, phi at Detector
  double L_tr_r_x[MAX], L_tr_r_y[MAX], L_tr_r_th[MAX], L_tr_r_ph[MAX];           // x, y, theta, phi at Rot-Coordinate
  double L_tr_x[MAX], L_tr_y[MAX], L_tr_th[MAX], L_tr_ph[MAX];                   // x, y, theta, phi
  double L_tr_time[MAX];                                                         // theta, time@RF
  double L_tr_p[MAX], L_tr_pathl[MAX], L_tr_px[MAX], L_tr_py[MAX], L_tr_pz[MAX]; // mom, length(TtoP), momx, momy, momz
  double L_tr_tg_dp[MAX], L_tr_tg_y[MAX], L_tr_tg_th[MAX], L_tr_tg_ph[MAX];      // dp, y, theta, phi at target
  double L_tr_vx[MAX], L_tr_vy[MAX], L_tr_vz[MAX]; // vertex X, Y, Z

  double L_rtRFtime[6];


  /////////
  //F1TDC//
  /////////
  double RTDC_F1FirstHit[RF1TDC];
  double RS2T_F1TDC[RS2], RS2B_F1TDC[RS2];
  double RS2_F1time[RS2], RS2T_F1time[RS2], RS2B_F1time[RS2];
  double RS0T_F1TDC[RS2], RS0B_F1TDC[RS2];
  double RS0_F1time[RS0], RS0T_F1time[RS0], RS0B_F1time[RS0];
  double RRF_F1time;
  double LTDC_F1FirstHit[LF1TDC];
  double LS2_F1time[LS2], LS2T_F1time[LS2], LS2B_F1time[LS2];
  double LS0_F1time[LS0], LS0T_F1time[LS0], LS0B_F1time[LS0];
  double LRF_F1time;

public:
  TChain *tree;

public:
  void add_tree(string ifname){tree->Add(Form("%s",ifname.c_str()));};
  void pack_tree();
  void chain_tree(string ifname);
  void readtreeHRSLR();
  ////////
  //Left//
  ////////
  void readtreeHRSL();
  void readtreeS0L();
  void readtreeS2L();
  void readtreeTrackL();
  void readtreeF1TDCL();
  void convertF1TDCL(ParamMan *param);
  /////////
  //Right//
  /////////
  void readtreeHRSR();
  void readtreeS0R();
  void readtreeS2R();
  void readtreeA1R();
  void readtreeA2R();
  void readtreeTrackR();
  void readtreeF1TDCR();
  void convertF1TDCR(ParamMan *param);
  double GetBeta_S0S2wF1TDCR(int itrack);
};

#endif
