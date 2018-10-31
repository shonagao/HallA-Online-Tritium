using namespace std;

#include "ana.h"

bool batch = false;
bool RHRS = false;
bool LHRS = false;

const double c = 0.299792458; // speed of light (m/ns)
const double MK = 0.493677;   // Kaon Mass (GeV)
const double Me = 0.000511;   // Electron Mass (GeV)
const double Mp = 0.938272;
const double ML = 1.115683;
const double TDCtoT = 0.05623;  //  (ns/ch)
const double PI = atan(1)*4.;

bool pdf_out = false;
bool root_out = false;
string ofname("output.pdf");
string ofroot("output.root");

#define F1TDC
//#define FADC

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
ana::ana()
{
  set = new Setting();
  set->Initialize();

  if( root_out ) ofp = new TFile(Form("%s",ofroot.c_str()),"recreate");

  tr = new Tree();

  tr->tree->SetBranchStatus("*",0);
  tr->readtree_COMN();
  if(RHRS)
    tr->readtree_RHRS();
  if(LHRS)
    tr->readtree_LHRS();

}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
ana::~ana(){
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::ReadTree(string name){
  tr->merge(name);
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
bool ana::GetEntry(int n){
  tr->tree->GetEntry(n);
  return true;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::Roop(){
  ENum = tr->tree->GetEntries();
  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();

  bool L_Tr = false; // LHRS Tracking Chi2 cut
  bool L_FP = false; // LHRS FP plane cut
  bool R_Tr = false; // RHRS Tracking Chi2 cut
  bool R_FP = false; // RHRS FP plane cut

  for(int n=0;n<ENum;n++){
    GetEntry(n);

    h_rbay_rbax->Fill( tr->rbax, tr->rbay );
    h_rbby_rbbx->Fill( tr->rbbx, tr->rbby );
    h_rby_rbx  ->Fill( tr->rbx , tr->rby );
//////////////
//// LHRS ////
//////////////
    if(LHRS){
      int NLtr = (int)tr->L_tr_n;  if(NLtr>MAX) NLtr = MAX;
      if(tr->L_T1>0){ h_L_trig->Fill(1); }
      if(tr->L_T2>0){ h_L_trig->Fill(2); }
      if(tr->L_T3>0){ h_L_trig->Fill(3); }

#ifdef F1TDC
      L_s0l_t = CalcF1TDC( tr->L_F1Fhit[27] - tr->L_F1Fhit[30], L_s0l_toff );
      L_s0r_t = CalcF1TDC( tr->L_F1Fhit[28] - tr->L_F1Fhit[30], L_s0r_toff );
      if( tr->L_F1Fhit[27]>0 && tr->L_F1Fhit[28]>0 ){
        L_s0_t = (L_s0l_t + L_s0r_t) / 2.;
      } else{ L_s0_t = -99; }

      for(int i=0;i<16;i++){
        L_s2l_t[i] = CalcF1TDC( tr->L_F1Fhit[i]    - tr->L_F1Fhit[30], L_s2l_toff[i] );
        L_s2r_t[i] = CalcF1TDC( tr->L_F1Fhit[i+48] - tr->L_F1Fhit[40], L_s2r_toff[i] );
        if( tr->L_F1Fhit[i]>0 && tr->L_F1Fhit[i+48]>0 ){
          L_s2_t[i] = (L_s2l_t[i] + L_s2r_t[i]) / 2.;
        } else{ L_s2_t[i] = -9999; }
      }
#endif
#ifdef FADC
      L_s0l_t = CalcFADCTDC( tr->L_s0_ltc_fadc[0], L_s0l_toff );
      L_s0r_t = CalcFADCTDC( tr->L_s0_rtc_fadc[0], L_s0r_toff );
      if( tr->L_s0_ltc_fadc[0]>0 && tr->L_s0_rtc_fadc[0]>0 ){
        L_s0_t = (L_s0l_t + L_s0r_t) / 2.;
      } else{ L_s0_t = -99; }

      for(int i=0;i<16;i++){
        L_s2l_t[i] = CalcFADCTDC( tr->L_s2_ltc_fadc[i], L_s2l_toff[i] );
        L_s2r_t[i] = CalcFADCTDC( tr->L_s2_rtc_fadc[i], L_s2r_toff[i] );
        if( tr->L_s2_ltc_fadc[i]>0 && tr->L_s2_rtc_fadc[i]>0 ){
          L_s2_t[i] = (L_s2l_t[i] + L_s2r_t[i]) / 2.;
        } else{ L_s2_t[i] = -9999; }
      }
#endif

      h_L_tr_n->Fill( tr->L_tr_n );
      for(int t=0;t<NLtr;t++){
        L_Tr = L_FP = false;
        // Cuts
        if( tr->L_tr_chi2[t]<0.01 ) L_Tr = true;
        if( tr->L_tr_th[t]<0.17*tr->L_tr_x[t]+0.025
         && tr->L_tr_th[t]>0.17*tr->L_tr_x[t]-0.035
         && tr->L_tr_th[t]<0.4*tr->L_tr_x[t]+0.13 ) L_FP = true;
      
        int s2pad = (int)tr->L_s2_trpad[t];
        double p    = tr->L_tr_p[t];
        double path = tr->L_s2_trpath[t] - tr->L_s0_trpath[t];
        double beta = -99, m2 = -99;
        if( L_s2_t[s2pad]>0 && L_s0_t>0 && s2pad>=0 ){
          beta = path / ( L_s2_t[s2pad] - L_s0_t ) / c;
          m2 = ( 1./beta/beta - 1. ) * p * p;
        }
        double betae = p / sqrt(Me*Me + p*p);

        h_L_tr_ch2   ->Fill( tr->L_tr_chi2[t] );
      
        if( L_Tr && L_FP && s2pad>=0 ){
          h_L_p        ->Fill( tr->L_tr_p[t] );
          h_L_pathl    ->Fill( tr->L_tr_pathl[t] );
          h_L_px       ->Fill( tr->L_tr_px[t] );
          h_L_py       ->Fill( tr->L_tr_py[t] );
          h_L_pz       ->Fill( tr->L_tr_pz[t] );
          h_L_tgy      ->Fill( tr->L_tr_tg_y[t] );
          h_L_tgth     ->Fill( tr->L_tr_tg_th[t] );
          h_L_tgph     ->Fill( tr->L_tr_tg_ph[t] );
          h_L_vx       ->Fill( tr->L_tr_vx[t] );
          h_L_vy       ->Fill( tr->L_tr_vy[t] );
          h_L_vz       ->Fill( tr->L_tr_vz[t] );
          h_L_y_x      ->Fill( tr->L_tr_x[t]    , tr->L_tr_y[t] );
          h_L_th_x     ->Fill( tr->L_tr_x[t]    , tr->L_tr_th[t] );
          h_L_ph_y     ->Fill( tr->L_tr_y[t]    , tr->L_tr_ph[t] );
          h_L_tgph_tgth->Fill( tr->L_tr_tg_th[t], tr->L_tr_tg_ph[t] );
      
          h_L_beta->Fill( beta );
          h_L_m2  ->Fill( m2 );
          h_L_beta_p ->Fill(  p, beta );
          h_L_beta_m2->Fill( m2, beta );
          h_L_dedx_p     ->Fill(  p, (tr->L_s0_dedx[0] + tr->L_s2_dedx[s2pad])/2. );
          h_L_dedx_m2    ->Fill( m2, (tr->L_s0_dedx[0] + tr->L_s2_dedx[s2pad])/2. );
          h_L_s0_dedx    ->Fill( tr->L_s0_dedx[0] );
          h_L_s0_dedx_x  ->Fill( tr->L_s0_trx[t], tr->L_s0_dedx[0] );
          h_L_s0_beta_x  ->Fill( tr->L_s0_trx[t], beta );
          h_L_s2_pad     ->Fill( s2pad );
          h_L_s2_dedx    ->Fill( tr->L_s2_dedx[s2pad] );
          h_L_s2_dedx_x  ->Fill( tr->L_s2_trx[t], tr->L_s2_dedx[s2pad] );
          h_L_s2_beta_x  ->Fill( tr->L_s2_trx[t], beta );
          h_L_s2_dedx_pad->Fill( s2pad, tr->L_s2_dedx[s2pad] );
          h_L_s2_beta_pad->Fill( s2pad, beta );

          double rftime = (tr->L_F1Fhit[47] - tr->L_F1Fhit[40]) * TDCtoT;
          double tgt = (L_s2_t[s2pad] - rftime) -  (tr->L_tr_pathl[t] + tr->L_s2_trpath[t])/betae/c;
          h_L_tgt      ->Fill( tgt );
          h_L_s2pad_tgt->Fill( tgt, s2pad );
          h_L_p_tgt    ->Fill( tgt, p );
          h_L_pathl_tgt->Fill( tgt, tr->L_tr_pathl[t] );
          h_L_tgy_tgt  ->Fill( tgt, tr->L_tr_tg_y[t] );
          h_L_tgth_tgt ->Fill( tgt, tr->L_tr_tg_th[t] );
          h_L_tgph_tgt ->Fill( tgt, tr->L_tr_tg_ph[t] );
          h_L_x_tgt    ->Fill( tgt, tr->L_tr_x[t] );
          h_L_y_tgt    ->Fill( tgt, tr->L_tr_y[t] );
        } // if L_Tr && L_FP
      } // for NLtr
    } // if LHRS

//////////////
//// RHRS ////
//////////////
    if(RHRS){
      int NRtr = (int)tr->R_tr_n;  if(NRtr>MAX) NRtr = MAX;
      if(tr->R_T4>0){ h_R_trig->Fill(4); }
      if(tr->R_T5>0){ h_R_trig->Fill(5); }
      if(tr->R_T6>0){ h_R_trig->Fill(6); }

#ifdef F1TDC
      R_s0l_t = CalcF1TDC( tr->R_F1Fhit[27] - tr->R_F1Fhit[30], R_s0l_toff );
      R_s0r_t = CalcF1TDC( tr->R_F1Fhit[28] - tr->R_F1Fhit[30], R_s0r_toff );
      if( tr->R_F1Fhit[27]>0 && tr->R_F1Fhit[28]>0 ){
        R_s0_t = (R_s0l_t + R_s0r_t) / 2.;
      } else{ R_s0_t = -99; }

      for(int i=0;i<16;i++){
        R_s2l_t[i] = CalcF1TDC( tr->R_F1Fhit[i]    - tr->R_F1Fhit[30], R_s2l_toff[i] );
        R_s2r_t[i] = CalcF1TDC( tr->R_F1Fhit[i+48] - tr->R_F1Fhit[40], R_s2r_toff[i] );
        if( tr->R_F1Fhit[i]>0 && tr->R_F1Fhit[i+48]>0 ){
          R_s2_t[i] = (R_s2l_t[i] + R_s2r_t[i]) / 2.;
        } else{ R_s2_t[i] = -9999; }
      }
#endif
#ifdef FADC
      R_s0l_t = CalcFADCTDC( tr->R_s0_ltc_fadc[0], R_s0l_toff );
      R_s0r_t = CalcFADCTDC( tr->R_s0_rtc_fadc[0], R_s0r_toff );
      if( tr->R_s0_ltc_fadc[0]>0 && tr->R_s0_rtc_fadc[0]>0 ){
        R_s0_t = (R_s0l_t + R_s0r_t) / 2.;
      } else{ R_s0_t = -99; }

      for(int i=0;i<16;i++){
        R_s2l_t[i] = CalcFADCTDC( tr->R_s2_ltc_fadc[i], R_s2l_toff[i] );
        R_s2r_t[i] = CalcFADCTDC( tr->R_s2_rtc_fadc[i], R_s2r_toff[i] );
        if( tr->R_s2_ltc_fadc[i]>0 && tr->R_s2_rtc_fadc[i]>0 ){
          R_s2_t[i] = (R_s2l_t[i] + R_s2r_t[i]) / 2.;
        } else{ R_s2_t[i] = -9999; }
      }
#endif


      h_R_tr_n->Fill( tr->R_tr_n );
      for(int t=0;t<NRtr;t++){
        R_Tr = R_FP = false;
        // Cuts
        if( tr->R_tr_chi2[t]<0.01 ) R_Tr = true;
        if( tr->R_tr_th[t]<0.17*tr->R_tr_x[t]+0.025
         && tr->R_tr_th[t]>0.17*tr->R_tr_x[t]-0.035
         && tr->R_tr_th[t]<0.4*tr->R_tr_x[t]+0.13 ) R_FP = true;
      
        int s2pad = (int)tr->R_s2_trpad[t];
        double p    = tr->R_tr_p[t];
        double path = tr->R_s2_trpath[t] - tr->R_s0_trpath[t];
        double beta = 0, m2 = 0;
        if( R_s2_t[s2pad]>0 && R_s0_t>0 && s2pad>=0 ){
          beta = path / ( R_s2_t[s2pad] - R_s0_t ) / c;
          m2 = ( 1./beta/beta - 1. ) * p * p;
        } 
        double betaK = p / sqrt(MK*MK + p*p);

        h_R_tr_ch2   ->Fill( tr->R_tr_chi2[t] );
      
        if( R_Tr && R_FP && s2pad>=0 ){
          h_R_p        ->Fill( tr->R_tr_p[t] );
          h_R_pathl    ->Fill( tr->R_tr_pathl[t] );
          h_R_px       ->Fill( tr->R_tr_px[t] );
          h_R_py       ->Fill( tr->R_tr_py[t] );
          h_R_pz       ->Fill( tr->R_tr_pz[t] );
          h_R_tgy      ->Fill( tr->R_tr_tg_y[t] );
          h_R_tgth     ->Fill( tr->R_tr_tg_th[t] );
          h_R_tgph     ->Fill( tr->R_tr_tg_ph[t] );
          h_R_vx       ->Fill( tr->R_tr_vx[t] );
          h_R_vy       ->Fill( tr->R_tr_vy[t] );
          h_R_vz       ->Fill( tr->R_tr_vz[t] );
          h_R_y_x      ->Fill( tr->R_tr_x[t]    , tr->R_tr_y[t] );
          h_R_th_x     ->Fill( tr->R_tr_x[t]    , tr->R_tr_th[t] );
          h_R_ph_y     ->Fill( tr->R_tr_y[t]    , tr->R_tr_ph[t] );
          h_R_tgph_tgth->Fill( tr->R_tr_tg_th[t], tr->R_tr_tg_ph[t] );
      
          h_R_beta       ->Fill( beta );
          h_R_m2         ->Fill( m2 );
          h_R_beta_p     ->Fill(  p, beta );
          h_R_beta_m2    ->Fill( m2, beta );
          h_R_dedx_p     ->Fill(  p, (tr->R_s0_dedx[0] + tr->R_s2_dedx[s2pad])/2. );
          h_R_dedx_m2    ->Fill( m2, (tr->R_s0_dedx[0] + tr->R_s2_dedx[s2pad])/2. );
          h_R_s0_dedx    ->Fill( tr->R_s0_dedx[0] );
          h_R_s0_dedx_x  ->Fill( tr->R_s0_trx[t], tr->R_s0_dedx[0] );
          h_R_s0_beta_x  ->Fill( tr->R_s0_trx[t], beta );
          h_R_s2_pad     ->Fill( s2pad );
          h_R_s2_dedx    ->Fill( tr->R_s2_dedx[s2pad] );
          h_R_s2_dedx_x  ->Fill( tr->R_s2_trx[t], tr->R_s2_dedx[s2pad] );
          h_R_s2_beta_x  ->Fill( tr->R_s2_trx[t], beta );
          h_R_s2_dedx_pad->Fill( s2pad, tr->R_s2_dedx[s2pad] );
          h_R_s2_beta_pad->Fill( s2pad, beta );

          h_R_a1_sum    ->Fill( tr->R_a1_asum_c );
          h_R_a2_sum    ->Fill( tr->R_a2_asum_c );
          h_R_a1_sum_x  ->Fill( tr->R_a1_trx[t], tr->R_a1_asum_c );
          h_R_a2_sum_x  ->Fill( tr->R_a2_trx[t], tr->R_a2_asum_c );
          h_R_a1_sum_p  ->Fill(               p, tr->R_a1_asum_c );
          h_R_a2_sum_p  ->Fill(               p, tr->R_a2_asum_c );
          h_R_a1_sum_m2 ->Fill(              m2, tr->R_a1_asum_c );
          h_R_a2_sum_m2 ->Fill(              m2, tr->R_a2_asum_c );

          double rftime = (tr->R_F1Fhit[15] - tr->R_F1Fhit[9]) * TDCtoT;
          double tgt = (R_s2_t[s2pad] - rftime) -  (tr->R_tr_pathl[t] + tr->R_s2_trpath[t])/betaK/c;
          h_R_tgt      ->Fill( tgt );
          h_R_s2pad_tgt->Fill( tgt, s2pad );
          h_R_p_tgt    ->Fill( tgt, p );
          h_R_pathl_tgt->Fill( tgt, tr->R_tr_pathl[t] );
          h_R_tgy_tgt  ->Fill( tgt, tr->R_tr_tg_y[t] );
          h_R_tgth_tgt ->Fill( tgt, tr->R_tr_tg_th[t] );
          h_R_tgph_tgt ->Fill( tgt, tr->R_tr_tg_ph[t] );
          h_R_x_tgt    ->Fill( tgt, tr->R_tr_x[t] );
          h_R_y_tgt    ->Fill( tgt, tr->R_tr_y[t] );
        } // if R_Tr && R_FP
      } // for NRtr
    } // if RHRS

/////////////////////
//// Coincidence ////
/////////////////////
    if(LHRS && RHRS){
      int NLtr = (int)tr->L_tr_n;  if(NLtr>MAX) NLtr = MAX;
      int NRtr = (int)tr->R_tr_n;  if(NRtr>MAX) NRtr = MAX;
      
      for(int lt=0;lt<NLtr;lt++){
        L_Tr = L_FP = false;
        if( tr->L_tr_chi2[lt]<0.01 ) L_Tr = true;
        if( tr->L_tr_th[lt]<0.17*tr->L_tr_x[lt]+0.025
         && tr->L_tr_th[lt]>0.17*tr->L_tr_x[lt]-0.035
         && tr->L_tr_th[lt]<0.40*tr->L_tr_x[lt]+0.130 ) L_FP = true;

        for(int rt=0;rt<NRtr;rt++){
          R_Tr = R_FP = false;
          if( tr->R_tr_chi2[rt]<0.01 ) R_Tr = true;
          if( tr->R_tr_th[rt]<0.17*tr->R_tr_x[rt]+0.025
           && tr->R_tr_th[rt]>0.17*tr->R_tr_x[rt]-0.035
           && tr->R_tr_th[rt]<0.40*tr->R_tr_x[rt]+0.130 ) R_FP = true;
        
          if( L_Tr && L_FP && R_Tr && R_FP ){

            int L_s2pad = (int)tr->L_s2_trpad[lt];
            double L_p     = tr->L_tr_p[lt];
            double L_E     = sqrt( Me*Me + L_p*L_p );
            double L_betae = L_p / sqrt(Me*Me + L_p*L_p);

            int R_s2pad    = (int)tr->R_s2_trpad[rt];
            double R_p     = tr->R_tr_p[rt];
            double R_E     = sqrt( MK*MK + R_p*R_p );
            double R_betaK = R_p / sqrt(MK*MK + R_p*R_p);

            double L_rftime = (tr->L_F1Fhit[47] - tr->L_F1Fhit[40]) * TDCtoT;
            double R_rftime = (tr->R_F1Fhit[15] - tr->R_F1Fhit[9] ) * TDCtoT;
            double L_tgt = (L_s2_t[L_s2pad] - L_rftime) -  (tr->L_tr_pathl[lt] + tr->L_s2_trpath[lt])/L_betae/c;
            double R_tgt = (R_s2_t[R_s2pad] - R_rftime) -  (tr->R_tr_pathl[rt] + tr->R_s2_trpath[rt])/R_betaK/c;
            double ct = L_tgt - R_tgt;
            h_ct->Fill( ct );
            h_Ls2x_ct ->Fill( tr->L_s2_trx[lt], ct );
            h_Rs2x_ct ->Fill( tr->R_s2_trx[rt], ct );
            h_a1sum_ct->Fill( tr->R_a1_asum_c, ct );
            h_a2sum_ct->Fill( tr->R_a2_asum_c, ct );

            TVector3 L_v, R_v, B_v;
            double Ee = 4.3185;
            L_v.SetMagThetaPhi( L_p, tr->L_tr_tg_th[lt], tr->L_tr_tg_ph[lt] );
            R_v.SetMagThetaPhi( R_p, tr->R_tr_tg_th[rt], tr->R_tr_tg_ph[rt] );
            B_v.SetMagThetaPhi( sqrt(Ee*Ee-Me*Me), 0, 0 );
            L_v.RotateZ( -13.2 / 180. * PI );
            R_v.RotateZ(  13.2 / 180. * PI );
            double mass = sqrt( (Ee + Mp - L_E - R_E)*(Ee + Mp - L_E - R_E)
                              - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
            double mm = mass - ML;

            h_mmall ->Fill( mm );
            if( fabs( tr->L_tr_vz[lt] + 0.115 ) < 0.015 || fabs( tr->L_tr_vz[lt] - 0.135 ) < 0.015 ){ 
              h_mmfoil->Fill( mm );
            }
            if( fabs( tr->L_tr_vz[lt] - 0.01 ) < 0.1 ){ 
              h_mm      ->Fill( mm );
              h_Lp_mm   ->Fill( mm, tr->L_tr_p[lt] );
              h_Ll_mm   ->Fill( mm, tr->L_tr_pathl[lt] );
              h_Ltgy_mm ->Fill( mm, tr->L_tr_tg_y[lt] );
              h_Ltgth_mm->Fill( mm, tr->L_tr_tg_th[lt] );
              h_Ltgph_mm->Fill( mm, tr->L_tr_tg_ph[lt] );
              h_Lvx_mm  ->Fill( mm, tr->L_tr_vx[lt] );
              h_Lvy_mm  ->Fill( mm, tr->L_tr_vy[lt] );
              h_Lvz_mm  ->Fill( mm, tr->L_tr_vz[lt] );
              h_Lx_mm   ->Fill( mm, tr->L_tr_x[lt] );
              h_Ly_mm   ->Fill( mm, tr->L_tr_y[lt] );
              h_Lth_mm  ->Fill( mm, tr->L_tr_th[lt] );
              h_Lph_mm  ->Fill( mm, tr->L_tr_ph[lt] );
            }
            if( fabs( tr->R_tr_vz[rt] - 0.01 ) < 0.1 ){ 
              h_Rp_mm   ->Fill( mm, tr->R_tr_p[rt] );
              h_Rl_mm   ->Fill( mm, tr->R_tr_pathl[rt] );
              h_Rtgy_mm ->Fill( mm, tr->R_tr_tg_y[rt] );
              h_Rtgth_mm->Fill( mm, tr->R_tr_tg_th[rt] );
              h_Rtgph_mm->Fill( mm, tr->R_tr_tg_ph[rt] );
              h_Rvx_mm  ->Fill( mm, tr->R_tr_vx[rt] );
              h_Rvy_mm  ->Fill( mm, tr->R_tr_vy[rt] );
              h_Rvz_mm  ->Fill( mm, tr->R_tr_vz[rt] );
              h_Rx_mm   ->Fill( mm, tr->R_tr_x[rt] );
              h_Ry_mm   ->Fill( mm, tr->R_tr_y[rt] );
              h_Rth_mm  ->Fill( mm, tr->R_tr_th[rt] );
              h_Rph_mm  ->Fill( mm, tr->R_tr_ph[rt] );
            }

          } // if L_Tr && L_FP && R_Tr && R_FP
        } // for NRtr
      } // for NLtr
    } // if LHRS && RHRS

    if(n % 10000 == 0){ cout<<n<<" / "<<ENum<<endl; }
  } // for ENum

  delete tr;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::Draw(){

  TCanvas *c0;
  TPaveText *p1 = new TPaveText(0.2,0.5,0.8,0.7,"NDC");
  p1->SetTextSize(0.05);  p1->SetFillColor(10);  p1->SetBorderSize(1);
  TText *t1;
  if( batch && pdf_out ){
    c0 = new TCanvas("c0","c0",800,300);
    c0->Print(Form("%s[",ofname.c_str()));
  }

  TCanvas *c1 = new TCanvas("c1","BPM",1000,800);
  c1->Divide(3,3,1E-5,1E-5);
  c1->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_rbay_rbax->Draw("colz");
  c1->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_rbay_rbax->ProjectionX()->Draw();
  c1->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_rbay_rbax->ProjectionY()->Draw();
  c1->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_rbby_rbbx->Draw("colz");
  c1->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_rbby_rbbx->ProjectionX()->Draw();
  c1->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_rbby_rbbx->ProjectionY()->Draw();
  c1->cd(7)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_rby_rbx  ->Draw("colz");
  c1->cd(8)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_rby_rbx  ->ProjectionX()->Draw();
  c1->cd(9)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_rby_rbx  ->ProjectionY()->Draw();
  if( batch && pdf_out ){
    c0->cd();
    t1 = p1->AddText("BPM"); p1->Draw();
    c0->Print(Form("%s",ofname.c_str()));
    p1->Clear(); c0->Clear();
    c1->Print(Form("%s",ofname.c_str()));
  }

// LHRS
  if( LHRS ){
    TCanvas *c2 = new TCanvas("c2","LHRS Track 1",1000,800);
    c2->Divide(3,3,1E-5,1E-5);
    c2->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_trig     ->Draw();
    c2->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_tr_n     ->Draw();
    c2->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_tr_ch2   ->Draw();
    c2->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_p        ->Draw();
    c2->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_pathl    ->Draw();
    c2->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_tgph_tgth->Draw("colz");
    c2->cd(7)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_tgy      ->Draw();
    c2->cd(8)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_tgth     ->Draw();
    c2->cd(9)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_tgph     ->Draw();
  
    TCanvas *c3 = new TCanvas("c3","LHRS Track 2",1000,800);
    c3->Divide(3,2,1E-5,1E-5);
    c3->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_vx  ->Draw();
    c3->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_vy  ->Draw();
    c3->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_vz  ->Draw();
    c3->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_y_x ->Draw("colz");
    c3->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_th_x->Draw("colz");
    c3->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_ph_y->Draw("colz");
  
    TCanvas *c4 = new TCanvas("c4","LHRS beta M2",1000,800);
    c4->Divide(3,2,1E-5,1E-5);
    c4->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_beta   ->Draw();
    c4->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_m2     ->Draw();
    c4->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_beta_p ->Draw("colz");
    c4->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_beta_m2->Draw("colz");
    c4->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_dedx_p ->Draw("colz");
    c4->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_dedx_m2->Draw("colz");
  
    TCanvas *c5 = new TCanvas("c5","LHRS Scin",1000,800);
    c5->Divide(3,3,1E-5,1E-5);
    c5->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_s0_dedx    ->Draw();
    c5->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_s0_dedx_x  ->Draw("colz");
    c5->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_s0_beta_x  ->Draw("colz");
    c5->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_s2_dedx    ->Draw();
    c5->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_s2_dedx_x  ->Draw("colz");
    c5->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_s2_beta_x  ->Draw("colz");
    c5->cd(7)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_s2_pad     ->Draw();
    c5->cd(8)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_s2_dedx_pad->Draw("colz");
    c5->cd(9)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_s2_beta_pad->Draw("colz");
  
    TCanvas *c6 = new TCanvas("c6","LHRS Time at Target",1000,800);
    c6->Divide(3,3,1E-5,1E-5);
    c6->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_L_tgt      ->Draw();
    c6->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_s2pad_tgt->Draw("colz");
    c6->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_p_tgt    ->Draw("colz");
    c6->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_pathl_tgt->Draw("colz");
    c6->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_tgy_tgt  ->Draw("colz");
    c6->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_tgth_tgt ->Draw("colz");
    c6->cd(7)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_tgph_tgt ->Draw("colz");
    c6->cd(8)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_x_tgt    ->Draw("colz");
    c6->cd(9)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_L_y_tgt    ->Draw("colz");
  
    if( batch && pdf_out ){
      c0->cd();
      t1 = p1->AddText("LHRS"); p1->Draw();
      c0->Print(Form("%s",ofname.c_str()));
      p1->Clear(); c0->Clear();
      c2->Print(Form("%s",ofname.c_str()));
      c3->Print(Form("%s",ofname.c_str()));
      c4->Print(Form("%s",ofname.c_str()));
      c5->Print(Form("%s",ofname.c_str()));
      c6->Print(Form("%s",ofname.c_str()));
    }
  }

// RHRS
  if( RHRS ){
    TCanvas *c7 = new TCanvas("c7","RHRS Track 1",1000,800);
    c7->Divide(3,3,1E-5,1E-5);
    c7->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_trig     ->Draw();
    c7->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_tr_n     ->Draw();
    c7->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_tr_ch2   ->Draw();
    c7->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_p        ->Draw();
    c7->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_pathl    ->Draw();
    c7->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_tgph_tgth->Draw("colz");
    c7->cd(7)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_tgy      ->Draw();
    c7->cd(8)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_tgth     ->Draw();
    c7->cd(9)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_tgph     ->Draw();
    
    TCanvas *c8 = new TCanvas("c8","RHRS Track 2",1000,800);
    c8->Divide(3,2,1E-5,1E-5);
    c8->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_vx  ->Draw();
    c8->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_vy  ->Draw();
    c8->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_vz  ->Draw();
    c8->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_y_x ->Draw("colz");
    c8->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_th_x->Draw("colz");
    c8->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_ph_y->Draw("colz");
    
    TCanvas *c9 = new TCanvas("c9","RHRS beta M2",1000,800);
    c9->Divide(3,2,1E-5,1E-5);
    c9->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_beta   ->Draw();
    c9->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_m2     ->Draw();
    c9->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_beta_p ->Draw("colz");
    c9->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_beta_m2->Draw("colz");
    c9->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_dedx_p ->Draw("colz");
    c9->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_dedx_m2->Draw("colz");
    
    TCanvas *c10 = new TCanvas("c10","RHRS Scin",1000,800);
    c10->Divide(3,3,1E-5,1E-5);
    c10->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_s0_dedx    ->Draw();
    c10->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_s0_dedx_x  ->Draw("colz");
    c10->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_s0_beta_x  ->Draw("colz");
    c10->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_s2_dedx    ->Draw();
    c10->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_s2_dedx_x  ->Draw("colz");
    c10->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_s2_beta_x  ->Draw("colz");
    c10->cd(7)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_s2_pad     ->Draw();
    c10->cd(8)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_s2_dedx_pad->Draw("colz");
    c10->cd(9)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_s2_beta_pad->Draw("colz");
    
    TCanvas *c11 = new TCanvas("c10","RHRS AC",1000,800);
    c11->Divide(4,2,1E-5,1E-5);
    c11->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_a1_sum   ->Draw();
    c11->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_a1_sum_x ->Draw("colz");
    c11->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_a1_sum_p ->Draw("colz");
    c11->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_a1_sum_m2->Draw("colz");
    c11->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_a2_sum   ->Draw();
    c11->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_a2_sum_x ->Draw("colz");
    c11->cd(7)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_a2_sum_p ->Draw("colz");
    c11->cd(8)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_a2_sum_m2->Draw("colz");
    
    TCanvas *c12 = new TCanvas("c12","RHRS Time at Target",1000,800);
    c12->Divide(3,3,1E-5,1E-5);
    c12->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_R_tgt      ->Draw();
    c12->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_s2pad_tgt->Draw("colz");
    c12->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_p_tgt    ->Draw("colz");
    c12->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_pathl_tgt->Draw("colz");
    c12->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_tgy_tgt  ->Draw("colz");
    c12->cd(6)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_tgth_tgt ->Draw("colz");
    c12->cd(7)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_tgph_tgt ->Draw("colz");
    c12->cd(8)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_x_tgt    ->Draw("colz");
    c12->cd(9)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_R_y_tgt    ->Draw("colz");

    if( batch && pdf_out ){
      c0->cd();
      t1 = p1->AddText("RHRS"); p1->Draw();
      c0->Print(Form("%s",ofname.c_str()));
      p1->Clear(); c0->Clear();
      c7 ->Print(Form("%s",ofname.c_str()));
      c8 ->Print(Form("%s",ofname.c_str()));
      c9 ->Print(Form("%s",ofname.c_str()));
      c10->Print(Form("%s",ofname.c_str()));
      c11->Print(Form("%s",ofname.c_str()));
      c12->Print(Form("%s",ofname.c_str()));
    }
  }

// COIN
  if( LHRS && RHRS ){
    TCanvas *c13 = new TCanvas("c13","Cointime",1000,800);
    c13->Divide(3,2,1E-5,1E-5);
    c13->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_ct      ->Draw();
    c13->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ls2x_ct ->Draw("colz");
    c13->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rs2x_ct ->Draw("colz");
    c13->cd(4)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_a1sum_ct->Draw("colz");
    c13->cd(5)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_a2sum_ct->Draw("colz");
  
    TCanvas *c14 = new TCanvas("c14","Missing Mass 1",1000,800);
    c14->Divide(3,1,1E-5,1E-5);
    c14->cd(1)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_mmall   ->Draw();
    c14->cd(2)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_mmfoil  ->Draw();
    c14->cd(3)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogy(); h_mm      ->Draw();
  
    TCanvas *c15 = new TCanvas("c15","Missing Mass ",1000,800);
    c15->Divide(4,3,1E-5,1E-5);
    c15->cd(1) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lp_mm   ->Draw("colz");
    c15->cd(2) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ll_mm   ->Draw("colz");
    c15->cd(3) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ltgy_mm ->Draw("colz");
    c15->cd(4) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ltgth_mm->Draw("colz");
    c15->cd(5) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ltgph_mm->Draw("colz");
    c15->cd(6) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lvx_mm  ->Draw("colz");
    c15->cd(7) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lvy_mm  ->Draw("colz");
    c15->cd(8) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lvz_mm  ->Draw("colz");
    c15->cd(9) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lx_mm   ->Draw("colz");
    c15->cd(10)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ly_mm   ->Draw("colz");
    c15->cd(11)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lth_mm  ->Draw("colz");
    c15->cd(12)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Lph_mm  ->Draw("colz");
  
    TCanvas *c16 = new TCanvas("c16","Missing Mass ",1000,800);
    c16->Divide(4,3,1E-5,1E-5);
    c16->cd(1) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rp_mm   ->Draw("colz");
    c16->cd(2) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rl_mm   ->Draw("colz");
    c16->cd(3) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rtgy_mm ->Draw("colz");
    c16->cd(4) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rtgth_mm->Draw("colz");
    c16->cd(5) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rtgph_mm->Draw("colz");
    c16->cd(6) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rvx_mm  ->Draw("colz");
    c16->cd(7) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rvy_mm  ->Draw("colz");
    c16->cd(8) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rvz_mm  ->Draw("colz");
    c16->cd(9) ->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rx_mm   ->Draw("colz");
    c16->cd(10)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Ry_mm   ->Draw("colz");
    c16->cd(11)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rth_mm  ->Draw("colz");
    c16->cd(12)->SetMargin(0.15,0.15,0.15,0.10); gPad->SetLogz(); h_Rph_mm  ->Draw("colz");

    if( batch && pdf_out ){
      c0->cd();
      t1 = p1->AddText("COIN"); p1->Draw();
      c0->Print(Form("%s",ofname.c_str()));
      p1->Clear(); c0->Clear();
      c13->Print(Form("%s",ofname.c_str()));
      c14->Print(Form("%s",ofname.c_str()));
      c15->Print(Form("%s",ofname.c_str()));
      c16->Print(Form("%s",ofname.c_str()));
    }
  } // if LHRS && RHRS

  if( batch && pdf_out ){
    c0->Print(Form("%s]",ofname.c_str()));
  }


}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::MakeHist(){
// Hist name is defined by "h + LorR + variable" for TH1.
//                         "h + LorR + variableY + variableX" for TH2.
/////////////
//// BPM ////
/////////////
  h_rbay_rbax = new TH2D("h_rbay_rbax","h_rbay_rbax",200,-4,4,200,-3,7);
  h_rbby_rbbx = new TH2D("h_rbby_rbbx","h_rbby_rbbx",200,-4,4,200,-3,7);
  h_rby_rbx   = new TH2D("h_rby_rbx"  ,"h_rby_rbx"  ,200,-6,4,200,-6,4);
  set->SetTH2(h_rbay_rbax,"BPM A"         ,"X","Y");
  set->SetTH2(h_rbby_rbbx,"BPM B"         ,"X","Y");
  set->SetTH2(h_rby_rbx  ,"Raster Pattern","X","Y");

//////////////
//// LHRS ////
//////////////
  h_L_trig = new TH1D("h_L_trig","h_L_trig",10,0,10);
  set->SetTH1(h_L_trig,"Trigger Flag","Trig No.","Counts");

  h_L_tr_n      = new TH1D("h_L_tr_n"     ,"h_L_tr_n"     ,30 ,    0,  30);
  h_L_tr_ch2    = new TH1D("h_L_tr_ch2"   ,"h_L_tr_ch2"   ,400,    0,0.05);
  h_L_p         = new TH1D("h_L_p"        ,"h_L_p"        ,400,    1,   3);
  h_L_pathl     = new TH1D("h_L_pathl"    ,"h_L_pathl"    ,400,   20,  30);
  h_L_px        = new TH1D("h_L_px"       ,"h_L_px"       ,400,    0,   2);
  h_L_py        = new TH1D("h_L_py"       ,"h_L_px"       ,400,   -1,   1);
  h_L_pz        = new TH1D("h_L_pz"       ,"h_L_px"       ,400,    1,   3);
  h_L_tgy       = new TH1D("h_L_tgy"      ,"h_L_tgy"      ,400, -0.5, 0.5);
  h_L_tgth      = new TH1D("h_L_tgth"     ,"h_L_tgth"     ,400, -0.2, 0.2);
  h_L_tgph      = new TH1D("h_L_tgph"     ,"h_L_tgph"     ,400, -0.2, 0.2);
  h_L_vx        = new TH1D("h_L_vx"       ,"h_L_vx"       ,400,-0.01,0.01);
  h_L_vy        = new TH1D("h_L_vy"       ,"h_L_vy"       ,400,-0.01,0.01);
  h_L_vz        = new TH1D("h_L_vz"       ,"h_L_vz"       ,400, -0.5, 0.5);
  h_L_y_x       = new TH2D("h_L_y_x"      ,"h_L_y_x"      ,200,   -1,  1 ,200,-0.2,0.2);
  h_L_th_x      = new TH2D("h_L_th_x"     ,"h_L_th_x"     ,200,   -1,  1 ,200,-0.3,0.3);
  h_L_ph_y      = new TH2D("h_L_ph_y"     ,"h_L_ph_y"     ,200, -0.2, 0.2,200,-0.2,0.2);
  h_L_tgph_tgth = new TH2D("h_L_tgph_tgth","h_L_tgph_tgth",200, -0.2, 0.2,200,-0.2,0.2);
  set->SetTH1(h_L_tr_n     ,"No. of Tracks"           ,"No. of Tracks"   ,"Counts");
  set->SetTH1(h_L_tr_ch2   ,"Tracking #chi^{2}"       ,"#chi^{2}"        ,"Counts");
  set->SetTH1(h_L_p        ,"Track Momentum"          ,"p (GeV/#it{c})"  ,"Counts");
  set->SetTH1(h_L_pathl    ,"Track Path Length"       ,"l (m)"           ,"Counts");
  set->SetTH1(h_L_px       ,"Momentum X"              ,"px (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_L_py       ,"Momentum Y"              ,"py (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_L_pz       ,"Momentum Z"              ,"pz (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_L_tgy      ,"Target Plane Y"          ,"y_{t} (m)"       ,"Counts");
  set->SetTH1(h_L_tgth     ,"Target Plane #theta"     ,"#theta_{t} (rad)","Counts");
  set->SetTH1(h_L_tgph     ,"Target Plane #phi"       ,"#phi_{t} (rad)"  ,"Counts");
  set->SetTH1(h_L_vx       ,"Vertex X"                ,"x_{v} (m)"       ,"Counts");
  set->SetTH1(h_L_vy       ,"Vertex Y"                ,"y_{v} (m)"       ,"Counts");
  set->SetTH1(h_L_vz       ,"Vertex Z"                ,"z_{v} (m)"       ,"Counts");
  set->SetTH2(h_L_y_x      ,"Focal Plane Y v.s X"     ,"X (m)"           ,"Y (m)");
  set->SetTH2(h_L_th_x     ,"Focal Plane #theta v.s X","X (m)"           ,"#theta (rad)");
  set->SetTH2(h_L_ph_y     ,"Focal Plane #phi v.s Y"  ,"Y (m)"           ,"#phi (rad)");
  set->SetTH2(h_L_tgph_tgth,"Target #phi v.s #theta"  ,"#theta_{t} (rad)","#phi_{t} (rad)");

  h_L_beta        = new TH1D("h_L_beta"       ,"h_L_beta"       ,400,   0,  2); 
  h_L_m2          = new TH1D("h_L_m2"         ,"h_L_m2"         ,400,-0.5,2.5); 
  h_L_beta_p      = new TH2D("h_L_beta_p"     ,"h_L_beta_p"     ,200,   1,  3,200,   0,  2); 
  h_L_beta_m2     = new TH2D("h_L_beta_m2"    ,"h_L_beta_m2"    ,200,   1,  3,200,-0.5,2.5); 
  h_L_dedx_p      = new TH2D("h_L_dedx_p"     ,"h_L_dedx_p"     ,200,   0, 10,200,   0,  2); 
  h_L_dedx_m2     = new TH2D("h_L_dedx_m2"    ,"h_L_dedx_m2"    ,200,   0, 10,200,-0.5,2.5); 
  h_L_s0_dedx     = new TH1D("h_L_s0_dedx"    ,"h_L_s0_dedx"    ,400,   0, 10); 
  h_L_s0_beta_x   = new TH2D("h_L_s0_beta_x"  ,"h_L_s0_beta_x"  ,200,  -1,  1,200,   0, 22); 
  h_L_s0_dedx_x   = new TH2D("h_L_s0_dedx_x"  ,"h_L_s0_dedx_x"  ,200,  -1,  1,200,   0, 10); 
  h_L_s2_pad      = new TH1D("h_L_s2_pad"     ,"h_L_s2_pad"     , 16,   0, 16); 
  h_L_s2_dedx     = new TH1D("h_L_s2_dedx"    ,"h_L_s2_dedx"    ,400,   0, 10); 
  h_L_s2_beta_x   = new TH2D("h_L_s2_beta_x"  ,"h_L_s2_beta_x"  ,200,  -1,  1,200,   0,  2); 
  h_L_s2_dedx_x   = new TH2D("h_L_s2_dedx_x"  ,"h_L_s2_dedx_x"  ,200,  -1,  1,200,   0,10 ); 
  h_L_s2_beta_pad = new TH2D("h_L_s2_beta_pad","h_L_s2_beta_pad", 16,  0, 16,200,    0,10 ); 
  h_L_s2_dedx_pad = new TH2D("h_L_s2_dedx_pad","h_L_s2_dedx_pad", 16,  0, 16,200,    0, 2 ); 
  set->SetTH1(h_L_beta       ,"Track beta"                    ,"#beta"                ,"Counts");
  set->SetTH1(h_L_m2         ,"Mass Square"                   ,"M^{2} (GeV^{2}/c^{4})","Counts");
  set->SetTH2(h_L_beta_p     ,"#beta v.s Momentum"            ,"p (GeV/c)"            ,"#beta",0.0);
  set->SetTH2(h_L_beta_m2    ,"#beta v.s Mass Square"         ,"M^{2} (GeV^{2}/c^{4})","#beta");
  set->SetTH2(h_L_dedx_p     ,"Energy Deposit v.s Momentum"   ,"p (GeV/c)"            ,"dE/dx ()");
  set->SetTH2(h_L_dedx_m2    ,"Energy Deposit v.s Mass Square","M^{2} (GeV^{2}/c^{4})","dE/dx ()");
  set->SetTH1(h_L_s0_dedx    ,"Energy Deposit (S0)"           ,"dE/dx"                ,"Counts");
  set->SetTH2(h_L_s0_beta_x  ,"#beta v.s X-pos (S0)"          ,"X (m)"                ,"#beta");
  set->SetTH2(h_L_s0_dedx_x  ,"Energy Deposit (S0) v.s X-pos" ,"X (m)"                ,"dE/dx ()");
  set->SetTH1(h_L_s2_pad     ,"Hit Paddle (S2)"               ,"Paddle No."           ,"Counts");
  set->SetTH1(h_L_s2_dedx    ,"Energy Deposit (S2)"           ,"dE/dx"                ,"Counts");
  set->SetTH2(h_L_s2_beta_x  ,"#beta v.s X-pos (S2)"          ,"X (m)"                ,"#beta");
  set->SetTH2(h_L_s2_dedx_x  ,"Energy Deposit (S2) v.s X-pos" ,"X (m)"                ,"dE/dx ()");
  set->SetTH2(h_L_s2_beta_pad,"#beta v.s Paddle (S2)"         ,"Paddle No."           ,"#beta");
  set->SetTH2(h_L_s2_dedx_pad,"Energy Deposit (S2) v.s Paddle","Paddle No."           ,"dE/dx ()");

  h_L_tgt       = new TH1D("h_L_tgt"      ,"h_L_tgt"       ,400,0,1000);
  h_L_s2pad_tgt = new TH2D("h_L_s2pad_tgt","h_L_s2pad_tgt" ,200,0,1000, 16,   0, 16);
  h_L_p_tgt     = new TH2D("h_L_p_tgt"    ,"h_L_p_tgt"     ,200,0,1000,200,   1,  3);
  h_L_pathl_tgt = new TH2D("h_L_pathl_tgt","h_L_pathl_tgt" ,200,0,1000,200,  20, 30);
  h_L_tgy_tgt   = new TH2D("h_L_tgy_tgt"  ,"h_L_tgy_tgt"   ,200,0,1000,200,-0.5,0.5);
  h_L_tgth_tgt  = new TH2D("h_L_tgth_tgt" ,"h_L_tgth_tgt"  ,200,0,1000,200,-0.2,0.2);
  h_L_tgph_tgt  = new TH2D("h_L_tgph_tgt" ,"h_L_tgph_tgt"  ,200,0,1000,200,-0.2,0.2);
  h_L_x_tgt     = new TH2D("h_L_x_tgt"    ,"h_L_x_tgt"     ,200,0,1000,200,  -1,  1);
  h_L_y_tgt     = new TH2D("h_L_y_tgt"    ,"h_L_y_tgt"     ,200,0,1000,200,-0.2,0.2);
  set->SetTH1(h_L_tgt      ,"Time at Target (S2-RF)","Time (ns)"         ,"Counts");
  set->SetTH2(h_L_s2pad_tgt,"S2 Paddle v.s Time at Target (S2-RF)"       ,"Time (ns)","Paddle No.");
  set->SetTH2(h_L_p_tgt    ,"Momentum v.s Time at Target (S2-RF)"        ,"Time (ns)","Momentum (GeV/c)");
  set->SetTH2(h_L_pathl_tgt,"Path Length v.s Time at Target (S2-RF)"     ,"Time (ns)","L (m)");
  set->SetTH2(h_L_tgy_tgt  ,"Y at Target v.s Time at Target (S2-RF)"     ,"Time (ns)","Y_{t}");
  set->SetTH2(h_L_tgth_tgt ,"#theta at Target v.s Time at Target (S2-RF)","Time (ns)","#theta_{t}");
  set->SetTH2(h_L_tgph_tgt ,"#phi at Target v.s Time at Target (S2-RF)"  ,"Time (ns)","#phi_{t}");
  set->SetTH2(h_L_x_tgt    ,"X at FP v.s Time at Target (S2-RF)"         ,"Time (ns)","X");
  set->SetTH2(h_L_y_tgt    ,"Y at FP v.s Time at Target (S2-RF)"         ,"Time (ns)","Y");

//////////////
//// RHRS ////
//////////////
  h_R_trig = new TH1D("h_R_trig","h_R_trig",10,0,10);
  set->SetTH1(h_R_trig,"Trigger Flag","Trig No.","Counts");

  h_R_tr_n      = new TH1D("h_R_tr_n"     ,"h_R_tr_n"     , 30,    0,  30); 
  h_R_tr_ch2    = new TH1D("h_R_tr_ch2"   ,"h_R_tr_ch2"   ,400,    0,0.05); 
  h_R_p         = new TH1D("h_R_p"        ,"h_R_p"        ,400,    1,   3); 
  h_R_pathl     = new TH1D("h_R_pathl"    ,"h_R_pathl"    ,400,   20,  30); 
  h_R_px        = new TH1D("h_R_px"       ,"h_R_px"       ,400,    0,   2); 
  h_R_py        = new TH1D("h_R_py"       ,"h_R_px"       ,400,   -1,   1); 
  h_R_pz        = new TH1D("h_R_pz"       ,"h_R_px"       ,400,    1,   3); 
  h_R_tgy       = new TH1D("h_R_tgy"      ,"h_R_tgy"      ,400, -0.5, 0.5); 
  h_R_tgth      = new TH1D("h_R_tgth"     ,"h_R_tgth"     ,400, -0.2, 0.2); 
  h_R_tgph      = new TH1D("h_R_tgph"     ,"h_R_tgph"     ,400, -0.2, 0.2); 
  h_R_vx        = new TH1D("h_R_vx"       ,"h_R_vx"       ,400,-0.01,0.01); 
  h_R_vy        = new TH1D("h_R_vy"       ,"h_R_vy"       ,400,-0.01,0.01); 
  h_R_vz        = new TH1D("h_R_vz"       ,"h_R_vz"       ,400, 0.5, 0.5); 
  h_R_y_x       = new TH2D("h_R_y_x"      ,"h_R_y_x"      ,200,   -1,  1,200,-0.2,0.2);
  h_R_th_x      = new TH2D("h_R_th_x"     ,"h_R_th_x"     ,200,   -1,  1,200,-0.3,0.3);
  h_R_ph_y      = new TH2D("h_R_ph_y"     ,"h_R_ph_y"     ,200,-0.2 ,0.2,200,-0.2,0.2);
  h_R_tgph_tgth = new TH2D("h_R_tgph_tgth","h_R_tgph_tgth",200, -0.2, 0.2,200,-0.2,0.2);
  set->SetTH1(h_R_tr_n  ,"No. of Tracks"           ,"No. of Tracks"   ,"Counts");
  set->SetTH1(h_R_tr_ch2,"Tracking #chi^{2}"       ,"#chi^{2}"        ,"Counts");
  set->SetTH1(h_R_p     ,"Track Momentum"          ,"p (GeV/#it{c})"  ,"Counts");
  set->SetTH1(h_R_pathl ,"Track Path Length"       ,"l (m)"           ,"Counts");
  set->SetTH1(h_R_px    ,"Momentum X"              ,"px (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_R_py    ,"Momentum Y"              ,"py (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_R_pz    ,"Momentum Z"              ,"pz (GeV/#it{c})" ,"Counts");
  set->SetTH1(h_R_tgy   ,"Target Plane Y"          ,"y_{t} (m)"       ,"Counts");
  set->SetTH1(h_R_tgth  ,"Target Plane #theta"     ,"#theta_{t} (rad)","Counts");
  set->SetTH1(h_R_tgph  ,"Target Plane #phi"       ,"#phi_{t} (rad)"  ,"Counts");
  set->SetTH1(h_R_vx    ,"Vertex X"                ,"x_{v} (m)"       ,"Counts");
  set->SetTH1(h_R_vy    ,"Vertex Y"                ,"y_{v} (m)"       ,"Counts");
  set->SetTH1(h_R_vz    ,"Vertex Z"                ,"z_{v} (m)"       ,"Counts");
  set->SetTH2(h_R_y_x   ,"Focal Plane Y v.s X"     ,"X (m)"           ,"Y (m)");
  set->SetTH2(h_R_th_x  ,"Focal Plane #theta v.s X","X (m)"           ,"#theta (rad)");
  set->SetTH2(h_R_ph_y  ,"Focal Plane #phi v.s Y"  ,"Y (m)"           ,"#phi (rad)");
  set->SetTH2(h_R_tgph_tgth,"Target #phi v.s #theta"  ,"#theta_{t} (rad)","#phi_{t} (rad)");

  h_R_beta      = new TH1D("h_R_beta"     ,"h_R_beta"     ,400,   0,  2); 
  h_R_m2        = new TH1D("h_R_m2"       ,"h_R_m2"       ,400,-0.5,2.5); 
  h_R_beta_p    = new TH2D("h_R_beta_p"   ,"h_R_beta_p"   ,200,   1,   3,200,   0,   2); 
  h_R_beta_m2   = new TH2D("h_R_beta_m2"  ,"h_R_beta_m2"  ,200,   1,   3,200,-0.5, 2.5); 
  h_R_dedx_p    = new TH2D("h_R_dedx_p"   ,"h_R_dedx_p"   ,200,   0,  10,200,   0,   2); 
  h_R_dedx_m2   = new TH2D("h_R_dedx_m2"  ,"h_R_dedx_m2"  ,200,   0,  10,200,-0.5, 2.5); 
  h_R_s0_dedx   = new TH1D("h_R_s0_dedx"  ,"h_R_s0_dedx"  ,400,   0,  10); 
  h_R_s0_beta_x = new TH2D("h_R_s0_beta_x","h_R_s0_beta_x",200,  -1,   1,200,   0,   2); 
  h_R_s0_dedx_x = new TH2D("h_R_s0_dedx_x","h_R_s0_dedx_x",200,  -1,   1,200,   0,  10); 
  h_R_s2_pad      = new TH1D("h_R_s2_pad"     ,"h_R_s2_pad"     , 16,   0, 16); 
  h_R_s2_dedx     = new TH1D("h_R_s2_dedx"    ,"h_R_s2_dedx"    ,400,   0, 10); 
  h_R_s2_beta_x   = new TH2D("h_R_s2_beta_x"  ,"h_R_s2_beta_x"  ,200,  -1,  1,200,   0,  2); 
  h_R_s2_dedx_x   = new TH2D("h_R_s2_dedx_x"  ,"h_R_s2_dedx_x"  ,200,  -1,  1,200,   0,10 ); 
  h_R_s2_beta_pad = new TH2D("h_R_s2_beta_pad","h_R_s2_beta_pad", 16,  0, 16,200,   0,10 ); 
  h_R_s2_dedx_pad = new TH2D("h_R_s2_dedx_pad","h_R_s2_dedx_pad", 16,  0, 16,200,   0, 2 ); 
  h_R_a1_sum    = new TH1D("h_R_a1_sum"   ,"h_R_a1_sum"   ,400,   0,1000);
  h_R_a1_sum_x  = new TH2D("h_R_a1_sum_x" ,"h_R_a1_sum_x" ,200,  -1,   1,200,   0,1000); 
  h_R_a1_sum_p  = new TH2D("h_R_a1_sum_p" ,"h_R_a1_sum_p" ,200,   1,   3,200,   0,1000); 
  h_R_a1_sum_m2 = new TH2D("h_R_a1_sum_m2","h_R_a1_sum_m2",200,-0.5, 2.5,200,   0,1000); 
  h_R_a2_sum    = new TH1D("h_R_a2_sum"   ,"h_R_a2_sum"   ,400,   0,1000);
  h_R_a2_sum_x  = new TH2D("h_R_a2_sum_x" ,"h_R_a2_sum_x" ,200,  -1,   1,200,   0,1000); 
  h_R_a2_sum_p  = new TH2D("h_R_a2_sum_p" ,"h_R_a2_sum_p" ,200,   1,   3,200,   0,1000); 
  h_R_a2_sum_m2 = new TH2D("h_R_a2_sum_m2","h_R_a2_sum_m2",200,-0.5, 2.5,200,   0,1000); 
  set->SetTH1(h_R_beta       ,"Track beta"                        ,"#beta"               ,"Counts");
  set->SetTH1(h_R_m2         ,"Mass Square"                       ,"M^{2} (GeV^{2}/c^{4}","Counts");
  set->SetTH2(h_R_beta_p     ,"#beta v.s Momentum"                ,"p (GeV/c)"           ,"#beta");
  set->SetTH2(h_R_beta_m2    ,"#beta v.s Mass Square"             ,"M^{2} (GeV^{2}/c^{4}","#beta");
  set->SetTH2(h_R_dedx_p     ,"Energy Deposit v.s Momentum"       ,"p (GeV/c)"           ,"dE/dx ()");
  set->SetTH2(h_R_dedx_m2    ,"Energy Deposit v.s Mass Square"    ,"M^{2} (GeV^{2}/c^{4}","dE/dx ()");
  set->SetTH1(h_R_s0_dedx    ,"Energy Deposit (S0)"               ,"dE/dx"               ,"Counts");
  set->SetTH2(h_R_s0_beta_x  ,"#beta v.s X-pos (S0)"              ,"X (m)"               ,"#beta");
  set->SetTH2(h_R_s0_dedx_x  ,"Energy Deposit (S0) v.s X-pos"     ,"X (m)"               ,"dE/dx ()");
  set->SetTH1(h_R_s2_pad     ,"Hit Paddle (S2)"                   ,"Paddle No."          ,"Counts");
  set->SetTH1(h_R_s2_dedx    ,"Energy Deposit (S2)"               ,"dE/dx"               ,"Counts");
  set->SetTH2(h_R_s2_beta_x  ,"#beta v.s X-pos (S2)"              ,"X (m)"               ,"#beta");
  set->SetTH2(h_R_s2_dedx_x  ,"Energy Deposit (S2) v.s X-pos"     ,"X (m)"               ,"dE/dx ()");
  set->SetTH2(h_R_s2_beta_pad,"#beta v.s Paddle (S2)"             ,"Paddle No."          ,"#beta");
  set->SetTH2(h_R_s2_dedx_pad,"Energy Deposit (S2) v.s Paddle"    ,"Paddle No."          ,"dE/dx ()");
  set->SetTH1(h_R_a1_sum     ,"Cherenkov SUM (A1)"                ,""                    ,"Counts");
  set->SetTH2(h_R_a1_sum_x   ,"Cherenkov SUM v.s X-pos (A1)"      ,"X (m)"               ,"");
  set->SetTH2(h_R_a1_sum_p   ,"Cherenkov SUM v.s Momentum (A1)"   ,"p (GeV/c)"           ,"");
  set->SetTH2(h_R_a1_sum_m2  ,"Cherenkov SUM v.s Mass Square (A1)","M^{2} (GeV^{2}/c^{4}","");
  set->SetTH1(h_R_a2_sum     ,"Cherenkov SUM (A2)"                ,""                    ,"Counts");
  set->SetTH2(h_R_a2_sum_x   ,"Cherenkov SUM v.s X-pos (A2)"      ,"X (m)"               ,"");
  set->SetTH2(h_R_a2_sum_p   ,"Cherenkov SUM v.s Momentum (A2)"   ,"p (GeV/c)"           ,"");
  set->SetTH2(h_R_a2_sum_m2  ,"Cherenkov SUM v.s Mass Square (A2)","M^{2} (GeV^{2}/c^{4}","");

  h_R_tgt       = new TH1D("h_R_tgt"      ,"h_R_tgt"       ,400,0,1000);
  h_R_s2pad_tgt = new TH2D("h_R_s2pad_tgt","h_R_s2pad_tgt" ,200,0,1000, 16,   0, 16);
  h_R_p_tgt     = new TH2D("h_R_p_tgt"    ,"h_R_p_tgt"     ,200,0,1000,200,   1,  3);
  h_R_pathl_tgt = new TH2D("h_R_pathl_tgt","h_R_pathl_tgt" ,200,0,1000,200,  20, 30);
  h_R_tgy_tgt   = new TH2D("h_R_tgy_tgt"  ,"h_R_tgy_tgt"   ,200,0,1000,200,-0.5,0.5);
  h_R_tgth_tgt  = new TH2D("h_R_tgth_tgt" ,"h_R_tgth_tgt"  ,200,0,1000,200,-0.2,0.2);
  h_R_tgph_tgt  = new TH2D("h_R_tgph_tgt" ,"h_R_tgph_tgt"  ,200,0,1000,200,-0.2,0.2);
  h_R_x_tgt     = new TH2D("h_R_x_tgt"    ,"h_R_x_tgt"     ,200,0,1000,200,  -1,  1);
  h_R_y_tgt     = new TH2D("h_R_y_tgt"    ,"h_R_y_tgt"     ,200,0,1000,200,-0.2,0.2);
  set->SetTH1(h_R_tgt      ,"Time at Target (S2-RF)","Time (ns)"         ,"Counts");
  set->SetTH2(h_R_s2pad_tgt,"S2 Paddle v.s Time at Target (S2-RF)"       ,"Time (ns)","Paddle No.");
  set->SetTH2(h_R_p_tgt    ,"Momentum v.s Time at Target (S2-RF)"        ,"Time (ns)","Momentum (GeV/c)");
  set->SetTH2(h_R_pathl_tgt,"Path Length v.s Time at Target (S2-RF)"     ,"Time (ns)","L (m)");
  set->SetTH2(h_R_tgy_tgt  ,"Y at Target v.s Time at Target (S2-RF)"     ,"Time (ns)","Y_{t}");
  set->SetTH2(h_R_tgth_tgt ,"#theta at Target v.s Time at Target (S2-RF)","Time (ns)","#theta_{t}");
  set->SetTH2(h_R_tgph_tgt ,"#phi at Target v.s Time at Target (S2-RF)"  ,"Time (ns)","#phi_{t}");
  set->SetTH2(h_R_x_tgt    ,"X at FP v.s Time at Target (S2-RF)"         ,"Time (ns)","X");
  set->SetTH2(h_R_y_tgt    ,"Y at FP v.s Time at Target (S2-RF)"         ,"Time (ns)","Y");


/////////////////////
//// Coincidence ////
/////////////////////
  h_ct       = new TH1D("h_ct"      ,"h_ct"      ,400, -100, 100); 
  h_Rs2x_ct  = new TH2D("h_Rs2x_ct" ,"h_Rs2x_ct" ,200, -100, 100,200,   -1,  1); 
  h_Ls2x_ct  = new TH2D("h_Ls2x_ct" ,"h_Ls2x_ct" ,200, -100, 100,200,   -1,  1); 
  h_a1sum_ct = new TH2D("h_a1sum_ct","h_a1sum_ct",200, -100, 100,200,    0,1000); 
  h_a2sum_ct = new TH2D("h_a2sum_ct","h_a2sum_ct",200, -100, 100,200,    0,1000); 
  h_mm       = new TH1D("h_mm"      ,"h_mm"      ,400,-0.05,0.25); 
  h_mmall    = new TH1D("h_mmall"   ,"h_mmall"   ,400,-0.05,0.25); 
  h_mmfoil   = new TH1D("h_mmfoil"  ,"h_mmfoil"  ,400,-0.05,0.25); 
  h_Rp_mm    = new TH2D("h_Rp_mm"   ,"h_Rp_mm"   ,200,-0.05,0.25,200,    1,   3); 
  h_Rl_mm    = new TH2D("h_Rl_mm"   ,"h_Rl_mm"   ,200,-0.05,0.25,200,   20,  30); 
  h_Rtgy_mm  = new TH2D("h_Rtgy_mm" ,"h_Rtgy_mm" ,200,-0.05,0.25,200,   20,  30); 
  h_Rtgth_mm = new TH2D("h_Rtgth_mm","h_Rtgth_mm",200,-0.05,0.25,200, -0.2, 0.2); 
  h_Rtgph_mm = new TH2D("h_Rtgph_mm","h_Rtgph_mm",200,-0.05,0.25,200, -0.2, 0.2); 
  h_Rvx_mm   = new TH2D("h_Rvx_mm"  ,"h_Rvx_mm"  ,200,-0.05,0.25,200,-0.01,0.01); 
  h_Rvy_mm   = new TH2D("h_Rvy_mm"  ,"h_Rvy_mm"  ,200,-0.05,0.25,200,-0.01,0.01); 
  h_Rvz_mm   = new TH2D("h_Rvz_mm"  ,"h_Rvz_mm"  ,200,-0.05,0.25,200, -0.5, 0.5); 
  h_Rx_mm    = new TH2D("h_Rx_mm"   ,"h_Rx_mm"   ,200,-0.05,0.25,200,   -1,   1); 
  h_Ry_mm    = new TH2D("h_Ry_mm"   ,"h_Ry_mm"   ,200,-0.05,0.25,200, -0.2, 0.2); 
  h_Rth_mm   = new TH2D("h_Rth_mm"  ,"h_Rth_mm"  ,200,-0.05,0.25,200, -0.3, 0.3); 
  h_Rph_mm   = new TH2D("h_Rph_mm"  ,"h_Rph_mm"  ,200,-0.05,0.25,200, -0.2, 0.2); 
  h_Lp_mm    = new TH2D("h_Lp_mm"   ,"h_Lp_mm"   ,200,-0.05,0.25,200,    1,   3); 
  h_Ll_mm    = new TH2D("h_Ll_mm"   ,"h_Ll_mm"   ,200,-0.05,0.25,200,   20,  30); 
  h_Ltgy_mm  = new TH2D("h_Ltgy_mm" ,"h_Ltgy_mm" ,200,-0.05,0.25,200,   20,  30); 
  h_Ltgth_mm = new TH2D("h_Ltgth_mm","h_Ltgth_mm",200,-0.05,0.25,200, -0.2, 0.2); 
  h_Ltgph_mm = new TH2D("h_Ltgph_mm","h_Ltgph_mm",200,-0.05,0.25,200, -0.2, 0.2); 
  h_Lvx_mm   = new TH2D("h_Lvx_mm"  ,"h_Lvx_mm"  ,200,-0.05,0.25,200,-0.01,0.01); 
  h_Lvy_mm   = new TH2D("h_Lvy_mm"  ,"h_Lvy_mm"  ,200,-0.05,0.25,200,-0.01,0.01); 
  h_Lvz_mm   = new TH2D("h_Lvz_mm"  ,"h_Lvz_mm"  ,200,-0.05,0.25,200, -0.5, 0.5); 
  h_Lx_mm    = new TH2D("h_Lx_mm"   ,"h_Lx_mm"   ,200,-0.05,0.25,200,   -1,   1); 
  h_Ly_mm    = new TH2D("h_Ly_mm"   ,"h_Ly_mm"   ,200,-0.05,0.25,200, -0.2, 0.2); 
  h_Lth_mm   = new TH2D("h_Lth_mm"  ,"h_Lth_mm"  ,200,-0.05,0.25,200, -0.3, 0.3); 
  h_Lph_mm   = new TH2D("h_Lph_mm"  ,"h_Lph_mm"  ,200,-0.05,0.25,200, -0.2, 0.2); 
  set->SetTH1(h_ct      ,"Coincidence Time"                      ,"Cointime (ns)"           ,"Counts");
  set->SetTH2(h_Rs2x_ct ,"RHRS S2 X-pos v.s Cointime"            ,"Cointime (ns)"           ,"X (m)");
  set->SetTH2(h_Ls2x_ct ,"LHRS S2 X-pos v.s Cointime"            ,"Cointime (ns)"           ,"X (m)");
  set->SetTH2(h_a1sum_ct,"RHRS A1 SUM v.s Cointime"              ,"Cointime (ns)"           ,"");
  set->SetTH2(h_a2sum_ct,"RHRS A1 SUM v.s Cointime"              ,"Cointime (ns)"           ,"");
  set->SetTH1(h_mm      ,"#Lambda Binding Energy (Gas)"          ,"-B_{#Lambda} (GeV/c^{2})","Counts");
  set->SetTH1(h_mmall   ,"#Lambda Binding Energy (w/o Z_{v} cut)","-B_{#Lambda} (GeV/c^{2})","Counts");
  set->SetTH1(h_mmfoil  ,"#Lambda Binding Energy (Al foil)"      ,"-B_{#Lambda} (GeV/c^{2})","Counts");
  set->SetTH2(h_Rp_mm   ,"RHRS Momentum v.s B_{#Lambda}"         ,"-B_{#Lambda} (GeV/c^{2})","Momentum (GeV/c)");
  set->SetTH2(h_Rl_mm   ,"RHRS Length v.s B_{#Lambda}"           ,"-B_{#Lambda} (GeV/c^{2})","Length (m)");
  set->SetTH2(h_Rtgy_mm ,"RHRS Target Y v.s B_{#Lambda}"         ,"-B_{#Lambda} (GeV/c^{2})","Y_{t} (m)");
  set->SetTH2(h_Rtgth_mm,"RHRS Target #theta v.s B_{#Lambda}"    ,"-B_{#Lambda} (GeV/c^{2})","#theta_{t} (rad)");
  set->SetTH2(h_Rtgph_mm,"RHRS Target #phi v.s B_{#Lambda}"      ,"-B_{#Lambda} (GeV/c^{2})","#phi_{t} (rad)");
  set->SetTH2(h_Rvx_mm  ,"RHRS Vertex X v.s B_{#Lambda}"         ,"-B_{#Lambda} (GeV/c^{2})","X_{v} (m)");
  set->SetTH2(h_Rvy_mm  ,"RHRS Vertex Y v.s B_{#Lambda}"         ,"-B_{#Lambda} (GeV/c^{2})","Y_{v} (m)");
  set->SetTH2(h_Rvz_mm  ,"RHRS Vertex Z v.s B_{#Lambda}"         ,"-B_{#Lambda} (GeV/c^{2})","Z_{v} (m)");
  set->SetTH2(h_Rx_mm   ,"RHRS FP X v.s B_{#Lambda}"             ,"-B_{#Lambda} (GeV/c^{2})","X_{FP} (m)");
  set->SetTH2(h_Ry_mm   ,"RHRS FP Y v.s B_{#Lambda}"             ,"-B_{#Lambda} (GeV/c^{2})","Y_{FP} (m)");
  set->SetTH2(h_Rth_mm  ,"RHRS FP #theta v.s B_{#Lambda}"        ,"-B_{#Lambda} (GeV/c^{2})","#theta_{FP} (rad)");
  set->SetTH2(h_Rph_mm  ,"RHRS FP #phi v.s B_{#Lambda}"          ,"-B_{#Lambda} (GeV/c^{2})","#phi_{FP} (rad)");
  set->SetTH2(h_Lp_mm   ,"LHRS Momentum v.s B_{#Lambda}"         ,"-B_{#Lambda} (GeV/c^{2})","Momentum (GeV/c)");
  set->SetTH2(h_Ll_mm   ,"LHRS Length v.s B_{#Lambda}"           ,"-B_{#Lambda} (GeV/c^{2})","Length (m)");
  set->SetTH2(h_Ltgy_mm ,"LHRS Target Y v.s B_{#Lambda}"         ,"-B_{#Lambda} (GeV/c^{2})","Y_{t} (m)");
  set->SetTH2(h_Ltgth_mm,"LHRS Target #theta v.s B_{#Lambda}"    ,"-B_{#Lambda} (GeV/c^{2})","#theta_{t} (rad)");
  set->SetTH2(h_Ltgph_mm,"LHRS Target #phi v.s B_{#Lambda}"      ,"-B_{#Lambda} (GeV/c^{2})","#phi_{t} (rad)");
  set->SetTH2(h_Lvx_mm  ,"LHRS Vertex X v.s B_{#Lambda}"         ,"-B_{#Lambda} (GeV/c^{2})","X_{v} (m)");
  set->SetTH2(h_Lvy_mm  ,"LHRS Vertex Y v.s B_{#Lambda}"         ,"-B_{#Lambda} (GeV/c^{2})","Y_{v} (m)");
  set->SetTH2(h_Lvz_mm  ,"LHRS Vertex Z v.s B_{#Lambda}"         ,"-B_{#Lambda} (GeV/c^{2})","Z_{v} (m)");
  set->SetTH2(h_Lx_mm   ,"LHRS FP X v.s B_{#Lambda}"             ,"-B_{#Lambda} (GeV/c^{2})","X_{FP} (m)");
  set->SetTH2(h_Ly_mm   ,"LHRS FP Y v.s B_{#Lambda}"             ,"-B_{#Lambda} (GeV/c^{2})","Y_{FP} (m)");
  set->SetTH2(h_Lth_mm  ,"LHRS FP #theta v.s B_{#Lambda}"        ,"-B_{#Lambda} (GeV/c^{2})","#theta_{FP} (rad)");
  set->SetTH2(h_Lph_mm  ,"LHRS FP #phi v.s B_{#Lambda}"          ,"-B_{#Lambda} (GeV/c^{2})","#phi_{FP} (rad)");

} // makehist()

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::DrawLine(double xmin, double ymin, double xmax, double ymax, int color, int style){
  line[LineID] = new TLine(xmin,ymin,xmax,ymax);
  line[LineID]->SetLineColor(color);
  line[LineID]->SetLineStyle(style);
  line[LineID]->Draw();
  LineID++;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::DrawText(double x, double y, string str, double size, int color){
  text[TextID] = new TLatex(x,y,str.c_str());
  text[TextID]->SetTextFont(42);
  text[TextID]->SetTextAlign(11);
  text[TextID]->SetTextColor(color);
  text[TextID]->SetTextSize(size);
  text[TextID]->Draw();
  TextID++;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
bool ana::Close(){

  ofp->Write();
  ofp->Close();

  return true;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::ReadParam(string name){

  L_s0l_toff = 10.;  L_s0r_toff = 0.;
  R_s0l_toff = 0.;  R_s0r_toff = 0.;
  for(int i=0;i<16;i++){
    L_s2l_toff[i] = 0.;
    L_s2r_toff[i] = 0.;
    R_s2l_toff[i] = 0.;
    R_s2r_toff[i] = 0.;
  }

}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
int main(int argc, char** argv){

  gErrorIgnoreLevel = kError;

  int ch;
  extern char *optarg;
  int MaxNum = 0;
  string ifname("test.list");
  string pname = "param/.default.param";

  while((ch=getopt(argc,argv,"hRLf:n:p:w:r:b"))!=-1){
    switch(ch){
    case 'R':
      RHRS= true;
      break;
    case 'L':
      LHRS= true;
      break;
    case 'f':
      ifname = optarg;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'p':
      pname = optarg;
      break;
    case 'w':
      ofname = optarg;
      pdf_out = true;
      batch = true;
      break;
    case 'r':
      ofroot = optarg;
      root_out = true;
      break;
    case 'b':
      gROOT->SetBatch(1);
      batch = true;
      break;
    case 'h':
      cout<<"-R: RHRS Analysis"<<endl;
      cout<<"-L: LHRS Analysis"<<endl;
      cout<<"-RL: Coin Analysis"<<endl;
      cout<<"-f (inputfile): input ROOT file name"<<endl;
      cout<<"-n (evnum): max no. of events"<<endl;
      cout<<"-p (input param): input parameter file name"<<endl;
      cout<<"-w (output pdffile): output PDF file name"<<endl;
      cout<<"-r (output rootfile): output ROOT file name"<<endl;
      cout<<"-b: batch mode"<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }

  ana *Ana = new ana();

  TApplication theApp("App", &argc, argv);

  Ana->SetMaxEvent(MaxNum);  // Set Maximum event roop
  Ana->MakeHist();           // Initialize histograms

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }

  Ana->ReadParam(Form("%s",pname.c_str()));

  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    Ana->ReadTree(Form("./rootfiles/%s",runname.c_str()));
  }

  Ana->Roop();
  Ana->Draw();               // save histograms to pdf file

  cout<<"Done!"<<endl;

  if( root_out ){ Ana->Close(); }
  if( batch ){ gSystem->Exit(1); }
//  gSystem->Exit(1);
  theApp.Run();
  return 0;

}

