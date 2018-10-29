#ifndef ana_h
#define ana_h 1

using namespace std;

#include "Setting.h"
#include "Tree.h"

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
class ana
{
  public:
    ana();
    ~ana();

  private:
    Setting *set;
    Tree *tr;
    TFile *ofp;

    int GetMaxEvent()     { return ENumMax; }

  public:
    void ReadTree(string name);
    void ReadParam(string name);
    void Roop();
    void Draw();

    bool Close();

    void MakeHist();

    void SetMaxEvent( int N )  { ENumMax = N; }

  private:
    int ENumMax;
    int ENum;

    bool GetEntry(int n);

    double CalcF1TDC(double tdc, double toff)   { return (tdc-toff)*0.05623;  }
    double CalcFADCTDC(double tdc, double toff) { return (tdc-toff);  }

    int LineID;
    int TextID;
    void DrawLine(double xmin, double ymin, double xmax, double ymax, int color=2, int style=1);
    void DrawText(double x, double y, string str, double size=0.05, int color=1);

  private:
// Lines, Textx
    TLine  *line[2000];
    TLatex *text[2000];

    TH2D *h_rbay_rbax, *h_rbby_rbbx;
    TH2D *h_rby_rbx;

//// LHRS ////
    TH1D *h_L_trig;

    TH1D *h_L_tr_n, *h_L_tr_ch2;
    TH1D *h_L_p, *h_L_pathl, *h_L_px, *h_L_py, *h_L_pz;
    TH1D *h_L_tgy, *h_L_tgth, *h_L_tgph;
    TH1D *h_L_vx, *h_L_vy, *h_L_vz;
    TH2D *h_L_y_x, *h_L_th_x, *h_L_ph_y;
    TH2D *h_L_tgph_tgth;

    TH1D *h_L_beta, *h_L_m2;
    TH2D *h_L_beta_p , *h_L_beta_m2;
    TH2D *h_L_dedx_p, *h_L_dedx_m2;
    TH1D *h_L_s0_dedx;
    TH2D *h_L_s0_dedx_x, *h_L_s0_beta_x;
    TH1D *h_L_s2_pad;
    TH1D *h_L_s2_dedx;
    TH2D *h_L_s2_dedx_x, *h_L_s2_beta_x;
    TH2D *h_L_s2_dedx_pad, *h_L_s2_beta_pad;

    TH1D *h_L_tgt;
    TH2D *h_L_s2pad_tgt;
    TH2D *h_L_p_tgt, *h_L_pathl_tgt, *h_L_tgy_tgt, *h_L_tgth_tgt, *h_L_tgph_tgt;
    TH2D *h_L_x_tgt, *h_L_y_tgt;

//// RHRS ////
    TH1D *h_R_trig;

    TH1D *h_R_tr_n, *h_R_tr_ch2;
    TH1D *h_R_p, *h_R_pathl, *h_R_px, *h_R_py, *h_R_pz;
    TH1D *h_R_tgy, *h_R_tgth, *h_R_tgph;
    TH1D *h_R_vx, *h_R_vy, *h_R_vz;
    TH2D *h_R_y_x, *h_R_th_x, *h_R_ph_y;
    TH2D *h_R_tgph_tgth;

    TH1D *h_R_beta, *h_R_m2;
    TH2D *h_R_beta_p , *h_R_beta_m2;
    TH2D *h_R_dedx_p, *h_R_dedx_m2;
    TH1D *h_R_s0_dedx;
    TH2D *h_R_s0_dedx_x, *h_R_s0_beta_x;
    TH1D *h_R_s2_pad;
    TH1D *h_R_s2_dedx;
    TH2D *h_R_s2_dedx_x, *h_R_s2_beta_x;
    TH2D *h_R_s2_dedx_pad, *h_R_s2_beta_pad;
    TH1D *h_R_a1_sum, *h_R_a2_sum;
    TH2D *h_R_a1_sum_x, *h_R_a2_sum_x;
    TH2D *h_R_a1_sum_p, *h_R_a2_sum_p;
    TH2D *h_R_a1_sum_m2, *h_R_a2_sum_m2;

    TH1D *h_R_tgt;
    TH2D *h_R_s2pad_tgt;
    TH2D *h_R_p_tgt, *h_R_pathl_tgt, *h_R_tgy_tgt, *h_R_tgth_tgt, *h_R_tgph_tgt;
    TH2D *h_R_x_tgt, *h_R_y_tgt;

//// Coin ////
    TH1D *h_ct;
    TH2D *h_Ls2x_ct;
    TH2D *h_Rs2x_ct;
    TH2D *h_a1sum_ct, *h_a2sum_ct;
    TH1D *h_mm, *h_mmall, *h_mmfoil;
    TH2D *h_Lp_mm, *h_Ll_mm, *h_Ltgy_mm, *h_Ltgth_mm, *h_Ltgph_mm;
    TH2D *h_Lvx_mm, *h_Lvy_mm, *h_Lvz_mm;
    TH2D *h_Lx_mm, *h_Ly_mm, *h_Lth_mm, *h_Lph_mm;
    TH2D *h_Rp_mm, *h_Rl_mm, *h_Rtgy_mm, *h_Rtgth_mm, *h_Rtgph_mm;
    TH2D *h_Rvx_mm, *h_Rvy_mm, *h_Rvz_mm;
    TH2D *h_Rx_mm, *h_Ry_mm, *h_Rth_mm, *h_Rph_mm;

  private:
    double L_s0l_toff    , L_s0r_toff;
    double L_s2l_toff[16], L_s2r_toff[16];
    double R_s0l_toff    , R_s0r_toff;
    double R_s2l_toff[16], R_s2r_toff[16];

    double L_s0l_t    , L_s0r_t    , L_s0_t;
    double L_s2l_t[16], L_s2r_t[16], L_s2_t[16];
    double R_s0l_t    , R_s0r_t    , R_s0_t;
    double R_s2l_t[16], R_s2r_t[16], R_s2_t[16];

};

#endif
