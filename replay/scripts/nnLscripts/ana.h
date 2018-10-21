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

    int GetMaxEvent()     { return ENumMax; }

  public:
    void ReadTree(string name);
    void Roop();
    void Draw();

    bool Close();

    void MakeHist();

    void SetMaxEvent( int N )  { ENumMax = N; }

  private:
    int ENumMax;
    int ENum;

    bool GetEntry(int n);

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

//// RHRS ////
    TH1D *h_R_T1, *h_R_T2, *h_R_T3, *h_R_T4, *h_R_T5, *h_R_T6, *h_R_T7, *h_R_T8;

    TH1D *h_R_tr_n, *h_R_tr_chi2;
    TH1D *h_R_p, *h_R_pathl, *h_R_px, *h_R_py, *h_R_pz;
    TH1D *h_R_tgy, *h_R_tgth, *h_R_tgph, *h_R_tgz;
    TH2D *h_R_y_x, *h_R_th_x, *h_R_ph_y;
    TH2D *h_R_tgph_tgth;

    TH1D *h_R_dedx, *h_R_beta;
    TH2D *h_R_s2_dedx_x, *h_R_s2_beta_x;
    TH2D *h_R_s2_dedx_beta;
    TH2D *h_R_a1_asum_x, *h_R_a2_asum_x;

    TH1D *h_R_mass2;
    TH2D *h_R_a1_asum_mass2, *h_R_a2_asum_mass2;

    TH1D *h_R_tgtime;

//// LHRS ////
    TH1D *h_L_T1, *h_L_T2, *h_L_T3, *h_L_T4, *h_L_T5, *h_L_T6, *h_L_T7, *h_L_T8;

    TH1D *h_L_tr_n, *h_L_tr_chi2;
    TH1D *h_L_p, *h_L_pathl, *h_L_px, *h_L_py, *h_L_pz;
    TH1D *h_L_tgy, *h_L_tgth, *h_L_tgph, *h_L_tgz;
    TH2D *h_L_y_x, *h_L_th_x, *h_L_ph_y;
    TH2D *h_L_tgph_tgth;

    TH1D *h_L_dedx, *h_L_beta;
    TH2D *h_L_s2_dedx_x, *h_L_s2_beta_x;
    TH2D *h_L_s2_dedx_beta;

    TH1D *h_L_mass2;

    TH1D *h_L_tgtime;

//// Coin ////
    TH1D *h_ctime;
    TH1D *h_missmass;

};

#endif
