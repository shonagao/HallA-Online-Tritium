using namespace std;

#include "ana.h"

bool batch = false;
bool RHRS = false;
bool LHRS = false;

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
ana::ana()
{
  set = new Setting();
  set->Initialize();

  tr = new Tree();

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

  for(int n=0;n<ENum;n++){
    GetEntry(n);
    int NLtr = (int)tr->L_tr_n;  if(NLtr>MAX) NLtr = MAX;
    int NRtr = (int)tr->R_tr_n;  if(NRtr>MAX) NRtr = MAX;

    for(int t=0;t<NLtr;t++){
      L_Tr = L_FP = false;
      // Cuts
      if( tr->L_tr_chi2[t]<0.01 ) L_Tr = true;
      if( tr->L_tr_th[t]<0.17*tr->L_tr_x[t]+0.025 && tr->L_tr_th[t]>0.17*tr->L_tr_x[t]-0.035 && tr->L_tr_th[t]<0.4*tr->L_tr_x[t]+0.13 ) L_FP = true;
    } // for NLtr

// Fill Hists

    if(n % 10000 == 0){
      cout<<n<<" / "<<ENum<<endl;
    }
  } // for ENum


  delete tr;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::Draw(){

//  TCanvas *c1 = new TCanvas("c1","c1",1500,900);

}


/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::MakeHist(){
// Hist name is defined by "h + LorR + variable" for TH1.
//                         "h + LorR + variableY + variableX" for TH2.
/////////////
//// BPM ////
/////////////
  h_rbay_rbax = new TH2D("h_rbay_rbax","h_rbay_rbax",200,-20,20,200,-20,20);
  h_rbby_rbbx = new TH2D("h_rbby_rbbx","h_rbby_rbbx",200,-20,20,200,-20,20);
  h_rby_rbx   = new TH2D("h_rby_rbx"  ,"h_rby_rbx"  ,200,-20,20,200,-20,20);
  set->SetTH2(h_rbay_rbax,"BPM A"         ,"X","Y");
  set->SetTH2(h_rbby_rbbx,"BPM B"         ,"X","Y");
  set->SetTH2(h_rby_rbx  ,"Raster Pattern","X","Y");

//////////////
//// LHRS ////
//////////////
  h_L_trig = new TH1D("h_L_trig","h_L_trig",10,0,10);
  set->SetTH1(h_L_trig,"Trigger Flag","Trig No.","Counts");

  h_L_tr_n   = new TH1D("h_L_tr_n"  ,"h_L_tr_n"  ,30 , 0    , 30);
  h_L_tr_ch2 = new TH1D("h_L_tr_ch2","h_L_tr_ch2",400, 0    , 10);
  h_L_p      = new TH1D("h_L_p"     ,"h_L_p"     ,400, 1    , 3);
  h_L_pathl  = new TH1D("h_L_pathl" ,"h_L_pathl" ,400, 20   , 30);
  h_L_px     = new TH1D("h_L_px"    ,"h_L_px"    ,400, 0    , 2);
  h_L_py     = new TH1D("h_L_py"    ,"h_L_px"    ,400, -1   , 1);
  h_L_pz     = new TH1D("h_L_pz"    ,"h_L_px"    ,400, 1    , 3);
  h_L_tgy    = new TH1D("h_L_tgy"   ,"h_L_tgy"   ,400, -1   , 1);
  h_L_tgth   = new TH1D("h_L_tgth"  ,"h_L_tgth"  ,400, -0.2 , 0.2);
  h_L_tgph   = new TH1D("h_L_tgph"  ,"h_L_tgph"  ,400, -0.2 , 0.2);
  h_L_vx     = new TH1D("h_L_vx"    ,"h_L_vx"    ,400, -0.01, 0.01);
  h_L_vy     = new TH1D("h_L_vy"    ,"h_L_vy"    ,400, -0.01, 0.01);
  h_L_vz     = new TH1D("h_L_vz"    ,"h_L_vz"    ,400, -1   , 1);
  h_L_y_x    = new TH2D("h_L_y_x"   ,"h_L_y_x"   ,200, -1   , 1  , 200, -0.2, 0.2);
  h_L_th_x   = new TH2D("h_L_th_x"  ,"h_L_th_x"  ,200, -1   , 1  , 200, -0.3, 0.3);
  h_L_ph_y   = new TH2D("h_L_ph_y"  ,"h_L_ph_y"  ,200, -0.2 , 0.2, 200, -0.2, 0.2);
  set->SetTH1(h_L_tr_n,"No. of Tracks","No. of Tracks","Counts");
  set->SetTH1(h_L_tr_ch2,"Tracking #chi^{2}","#chi^{2}","Counts");
  set->SetTH1(h_L_p     ,"Track Momentum","p (GeV/#it{c})","Counts");
  set->SetTH1(h_L_pathl ,"Track Path Length","l (m)","Counts");
  set->SetTH1(h_L_px    ,"Momentum X","px (GeV/#it{c})","Counts");
  set->SetTH1(h_L_py    ,"Momentum Y","py (GeV/#it{c})","Counts");
  set->SetTH1(h_L_pz    ,"Momentum Z","pz (GeV/#it{c})","Counts");
  set->SetTH1(h_L_tgy   ,"Target Plane Y","y_{t} (m)","Counts");
  set->SetTH1(h_L_tgth  ,"Target Plane #theta","#theta_{t} (rad)","Counts");
  set->SetTH1(h_L_tgph  ,"Target Plane #phi","#phi_{t} (rad)","Counts");
  set->SetTH1(h_L_vx    ,"Vertex X","x_{v} (m)","Counts");
  set->SetTH1(h_L_vy    ,"Vertex Y","y_{v} (m)","Counts");
  set->SetTH1(h_L_vz    ,"Vertex Z","z_{v} (m)","Counts");
  set->SetTH2(h_L_y_x   ,"Focal Plane Y v.s X","X (m)","Y (m)");
  set->SetTH2(h_L_th_x  ,"Focal Plane #theta v.s X","X (m)","#theta (rad)");
  set->SetTH2(h_L_ph_y  ,"Focal Plane #phi v.s Y","Y (m)","#phi (rad)");

  h_L_beta       = new TH1D("h_L_beta"      ,"h_L_beta"      ,400,0   ,2);
  h_L_m2         = new TH1D("h_L_m2"        ,"h_L_m2"        ,400,-0.5,2.5);
  h_L_beta_p     = new TH2D("h_L_beta_p"    ,"h_L_beta_p"    ,200,1,3 ,200,0,2);
  h_L_beta_m2    = new TH2D("h_L_beta_m2"   ,"h_L_beta_m2"   ,200,1,3 ,200,-0.5,2.5);
  h_L_dedx_p     = new TH2D("h_L_dedx_p"    ,"h_L_dedx_p"    ,200,0,10,200,0,2);
  h_L_dedx_m2    = new TH2D("h_L_dedx_m2"   ,"h_L_dedx_m2"   ,200,0,10,200,-0.5,2.5);
  h_L_s0_dedx    = new TH1D("h_L_s0_dedx"   ,"h_L_s0_dedx"   ,400,0,10);
  h_L_s0_beta_x  = new TH2D("h_L_s0_beta_x" ,"h_L_s0_beta_x" ,200,-1,1,200,0,2);
  h_L_s0_dedx_x  = new TH2D("h_L_s0_dedx_x" ,"h_L_s0_dedx_x" ,200,-1,1,200,0,10);
  h_L_s2_dedx    = new TH1D("h_L_s2_dedx"   ,"h_L_s2_dedx"   ,400,0,10);
  h_L_s2_beta_x  = new TH2D("h_L_s2_beta_x" ,"h_L_s2_beta_x" ,200,-1,1,200,0, 2);
  h_L_s2_dedx_x  = new TH2D("h_L_s2_dedx_x" ,"h_L_s2_dedx_x" ,200,-1,1,200,0, 10);
  set->SetTH1(h_L_beta      ,"Track beta","#beta","Counts");
  set->SetTH1(h_L_m2        ,"Mass Square","M^{2} (GeV^{2}/c^{4}","Counts");
  set->SetTH2(h_L_beta_p    ,"#beta v.s Momentum","p (GeV/c)","#beta");
  set->SetTH2(h_L_beta_m2   ,"#beta v.s Mass Square","M^{2} (GeV^{2}/c^{4}","#beta");
  set->SetTH2(h_L_dedx_p    ,"Energy Deposit v.s Momentum","p (GeV/c)","dE/dx ()");
  set->SetTH2(h_L_dedx_m2   ,"Energy Deposit v.s Mass Square","M^{2} (GeV^{2}/c^{4}","dE/dx ()");
  set->SetTH1(h_L_s0_dedx   ,"Energy Deposit (S0)","dE/dx","Counts");
  set->SetTH2(h_L_s0_beta_x ,"#beta v.s X-pos (S0)","X (m)","#beta");
  set->SetTH2(h_L_s0_dedx_x ,"Energy Deposit (S0) v.s X-pos","X (m)","dE/dx ()");
  set->SetTH1(h_L_s2_dedx   ,"Energy Deposit (S2)","dE/dx","Counts");
  set->SetTH2(h_L_s2_beta_x ,"#beta v.s X-pos (S2)","X (m)","#beta");
  set->SetTH2(h_L_s2_dedx_x ,"Energy Deposit (S2) v.s X-pos","X (m)","dE/dx ()");

  h_L_tgtime    = new TH1D("h_L_tgtime"   ,"h_L_tgtime" ,400,0,1000);
  set->SetTH1(h_L_tgtime,"Time at Target (S2-RF)","Time (ns)","Counts");

//////////////
//// RHRS ////
//////////////
  h_R_trig = new TH1D("h_R_trig","h_R_trig",10,0,10);
  set->SetTH1(h_R_trig,"Trigger Flag","Trig No.","Counts");

  h_R_tr_n   = new TH1D("h_R_tr_n"  ,"h_R_tr_n"  ,30 , 0    , 30);
  h_R_tr_ch2 = new TH1D("h_R_tr_ch2","h_R_tr_ch2",400, 0    , 10);
  h_R_p      = new TH1D("h_R_p"     ,"h_R_p"     ,400, 1    , 3);
  h_R_pathl  = new TH1D("h_R_pathl" ,"h_R_pathl" ,400, 20   , 30);
  h_R_px     = new TH1D("h_R_px"    ,"h_R_px"    ,400, 0    , 2);
  h_R_py     = new TH1D("h_R_py"    ,"h_R_px"    ,400, -1   , 1);
  h_R_pz     = new TH1D("h_R_pz"    ,"h_R_px"    ,400, 1    , 3);
  h_R_tgy    = new TH1D("h_R_tgy"   ,"h_R_tgy"   ,400, -1   , 1);
  h_R_tgth   = new TH1D("h_R_tgth"  ,"h_R_tgth"  ,400, -0.2 , 0.2);
  h_R_tgph   = new TH1D("h_R_tgph"  ,"h_R_tgph"  ,400, -0.2 , 0.2);
  h_R_vx     = new TH1D("h_R_vx"    ,"h_R_vx"    ,400, -0.01, 0.01);
  h_R_vy     = new TH1D("h_R_vy"    ,"h_R_vy"    ,400, -0.01, 0.01);
  h_R_vz     = new TH1D("h_R_vz"    ,"h_R_vz"    ,400, -1   , 1);
  h_R_y_x    = new TH2D("h_R_y_x"   ,"h_R_y_x"   ,200, -1   , 1  , 200, -0.2, 0.2);
  h_R_th_x   = new TH2D("h_R_th_x"  ,"h_R_th_x"  ,200, -1   , 1  , 200, -0.3, 0.3);
  h_R_ph_y   = new TH2D("h_R_ph_y"  ,"h_R_ph_y"  ,200, -0.2 , 0.2, 200, -0.2, 0.2);
  set->SetTH1(h_R_tr_n,"No. of Tracks","No. of Tracks","Counts");
  set->SetTH1(h_R_tr_ch2,"Tracking #chi^{2}","#chi^{2}","Counts");
  set->SetTH1(h_R_p     ,"Track Momentum","p (GeV/#it{c})","Counts");
  set->SetTH1(h_R_pathl ,"Track Path Length","l (m)","Counts");
  set->SetTH1(h_R_px    ,"Momentum X","px (GeV/#it{c})","Counts");
  set->SetTH1(h_R_py    ,"Momentum Y","py (GeV/#it{c})","Counts");
  set->SetTH1(h_R_pz    ,"Momentum Z","pz (GeV/#it{c})","Counts");
  set->SetTH1(h_R_tgy   ,"Target Plane Y","y_{t} (m)","Counts");
  set->SetTH1(h_R_tgth  ,"Target Plane #theta","#theta_{t} (rad)","Counts");
  set->SetTH1(h_R_tgph  ,"Target Plane #phi","#phi_{t} (rad)","Counts");
  set->SetTH1(h_R_vx    ,"Vertex X","x_{v} (m)","Counts");
  set->SetTH1(h_R_vy    ,"Vertex Y","y_{v} (m)","Counts");
  set->SetTH1(h_R_vz    ,"Vertex Z","z_{v} (m)","Counts");
  set->SetTH2(h_R_y_x   ,"Focal Plane Y v.s X","X (m)","Y (m)");
  set->SetTH2(h_R_th_x  ,"Focal Plane #theta v.s X","X (m)","#theta (rad)");
  set->SetTH2(h_R_ph_y  ,"Focal Plane #phi v.s Y","Y (m)","#phi (rad)");

  h_R_beta       = new TH1D("h_R_beta"      ,"h_R_beta"      ,400,0   ,2);
  h_R_m2         = new TH1D("h_R_m2"        ,"h_R_m2"        ,400,-0.5,2.5);
  h_R_beta_p     = new TH2D("h_R_beta_p"    ,"h_R_beta_p"    ,200,1,3 ,200,0,2);
  h_R_beta_m2    = new TH2D("h_R_beta_m2"   ,"h_R_beta_m2"   ,200,1,3 ,200,-0.5,2.5);
  h_R_dedx_p     = new TH2D("h_R_dedx_p"    ,"h_R_dedx_p"    ,200,0,10,200,0,2);
  h_R_dedx_m2    = new TH2D("h_R_dedx_m2"   ,"h_R_dedx_m2"   ,200,0,10,200,-0.5,2.5);
  h_R_s0_dedx    = new TH1D("h_R_s0_dedx"   ,"h_R_s0_dedx"   ,400,0,10);
  h_R_s0_beta_x  = new TH2D("h_R_s0_beta_x" ,"h_R_s0_beta_x" ,200,-1,1,200,0,2);
  h_R_s0_dedx_x  = new TH2D("h_R_s0_dedx_x" ,"h_R_s0_dedx_x" ,200,-1,1,200,0,10);
  h_R_s2_dedx    = new TH1D("h_R_s2_dedx"   ,"h_R_s2_dedx"   ,400,0,10);
  h_R_s2_beta_x  = new TH2D("h_R_s2_beta_x" ,"h_R_s2_beta_x" ,200,-1,1,200,0, 2);
  h_R_s2_dedx_x  = new TH2D("h_R_s2_dedx_x" ,"h_R_s2_dedx_x" ,200,-1,1,200,0, 10);
  h_R_a1_sum    = new TH1D("h_R_a1_sum"   ,"h_R_a1_sum"   ,400,0,1000);
  h_R_a1_sum_x  = new TH2D("h_R_a1_sum_x" ,"h_R_a1_sum_x" ,200,-1,1,200,0,1000);
  h_R_a1_sum_p  = new TH2D("h_R_a1_sum_p" ,"h_R_a1_sum_p" ,200,1,3,200,0,1000);
  h_R_a1_sum_m2 = new TH2D("h_R_a1_sum_m2","h_R_a1_sum_m2",200,-0.5,2.5,200,0,1000);
  h_R_a2_sum    = new TH1D("h_R_a2_sum"   ,"h_R_a2_sum"   ,400,0,1000);
  h_R_a2_sum_x  = new TH2D("h_R_a2_sum_x" ,"h_R_a2_sum_x" ,200,-1,1,200,0,1000);
  h_R_a2_sum_p  = new TH2D("h_R_a2_sum_p" ,"h_R_a2_sum_p" ,200,1,3,200,0,1000);
  h_R_a2_sum_m2 = new TH2D("h_R_a2_sum_m2","h_R_a2_sum_m2",200,-0.5,2.5,200,0,1000);
  set->SetTH1(h_R_beta      ,"Track beta","#beta","Counts");
  set->SetTH1(h_R_m2        ,"Mass Square","M^{2} (GeV^{2}/c^{4}","Counts");
  set->SetTH2(h_R_beta_p    ,"#beta v.s Momentum","p (GeV/c)","#beta");
  set->SetTH2(h_R_beta_m2   ,"#beta v.s Mass Square","M^{2} (GeV^{2}/c^{4}","#beta");
  set->SetTH2(h_R_dedx_p    ,"Energy Deposit v.s Momentum","p (GeV/c)","dE/dx ()");
  set->SetTH2(h_R_dedx_m2   ,"Energy Deposit v.s Mass Square","M^{2} (GeV^{2}/c^{4}","dE/dx ()");
  set->SetTH1(h_R_s0_dedx   ,"Energy Deposit (S0)","dE/dx","Counts");
  set->SetTH2(h_R_s0_beta_x ,"#beta v.s X-pos (S0)","X (m)","#beta");
  set->SetTH2(h_R_s0_dedx_x ,"Energy Deposit (S0) v.s X-pos","X (m)","dE/dx ()");
  set->SetTH1(h_R_s2_dedx   ,"Energy Deposit (S2)","dE/dx","Counts");
  set->SetTH2(h_R_s2_beta_x ,"#beta v.s X-pos (S2)","X (m)","#beta");
  set->SetTH2(h_R_s2_dedx_x ,"Energy Deposit (S2) v.s X-pos","X (m)","dE/dx ()");
  set->SetTH1(h_R_a1_sum   ,"Cherenkov SUM (A1)","","Counts");
  set->SetTH2(h_R_a1_sum_x ,"Cherenkov SUM v.s X-pos (A1)","X (m)","");
  set->SetTH2(h_R_a1_sum_p ,"Cherenkov SUM v.s Momentum (A1)","p (GeV/c)","");
  set->SetTH2(h_R_a1_sum_m2,"Cherenkov SUM v.s Mass Square (A1)","M^{2} (GeV^{2}/c^{4}","");
  set->SetTH1(h_R_a2_sum   ,"Cherenkov SUM (A2)","","Counts");
  set->SetTH2(h_R_a2_sum_x ,"Cherenkov SUM v.s X-pos (A2)","X (m)","");
  set->SetTH2(h_R_a2_sum_p ,"Cherenkov SUM v.s Momentum (A2)","p (GeV/c)","");
  set->SetTH2(h_R_a2_sum_m2,"Cherenkov SUM v.s Mass Square (A2)","M^{2} (GeV^{2}/c^{4}","");

  h_R_tgtime    = new TH1D("h_R_tgtime"   ,"h_R_tgtime" ,400,0,1000);
  set->SetTH1(h_R_tgtime,"Time at Target (S2-RF)","Time (ns)","Counts");

/////////////////////
//// Coincidence ////
/////////////////////
  h_ct       = new TH1D("h_ct","h_ct",400,-100,100);
  h_Rs2x_ct  = new TH2D("h_Rs2x_ct","h_Rs2x_ct",200,-100,100,200,-1,1);
  h_Ls2x_ct  = new TH2D("h_Ls2x_ct","h_Ls2x_ct",200,-100,100,200,-1,1);
  h_m2_ct    = new TH2D("h_m2_ct","h_m2_ct",200,-0.5,2.5,200,-1,1);
  h_beta_ct  = new TH2D("h_beta_ct","h_beta_ct",200,0,2,200,-1,1);
  h_dedx_ct  = new TH2D("h_dedx_ct","h_dedx_ct",200,0,10,200,-1,1);
  h_a1sum_ct = new TH2D("h_a1sum_ct","h_a1sum_ct",200,0,100,200,-1,1);
  h_a2sum_ct = new TH2D("h_a2sum_ct","h_a2sum_ct",200,0,100,200,-1,1);
  h_mm       = new TH1D("h_mm","h_mm",400,-0.05,0.25);
  h_mmall    = new TH1D("h_mmall","h_mmall",400,-0.05,0.25);
  h_mmfoil   = new TH1D("h_mmfoil","h_mmfoil",400,-0.05,0.25);
  h_Rp_mm    = new TH2D("h_Rp_mm","h_Rp_mm",200,-0.05,0.25,200,1,3);
  h_Rl_mm    = new TH2D("h_Rl_mm","h_Rl_mm",200,-0.05,0.25,200,20,30);
  h_Rtgy_mm  = new TH2D("h_Rtgy_mm","h_Rtgy_mm",200,-0.05,0.25,200,20,30);
  h_Rtgth_mm = new TH2D("h_Rtgth_mm","h_Rtgth_mm",200,-0.05,0.25,200,-0.2,0.2);
  h_Rtgph_mm = new TH2D("h_Rtgph_mm","h_Rtgph_mm",200,-0.05,0.25,200,-0.2,0.2);
  h_Rvx_mm   = new TH2D("h_Rvx_mm","h_Rvx_mm",200,-0.05,0.25,200,-0.01,0.01);
  h_Rvy_mm   = new TH2D("h_Rvy_mm","h_Rvy_mm",200,-0.05,0.25,200,-0.01,0.01);
  h_Rvz_mm   = new TH2D("h_Rvz_mm","h_Rvz_mm",200,-0.05,0.25,200,-1,1);
  h_Rx_mm    = new TH2D("h_Rx_mm","h_Rx_mm",200,-0.05,0.25,200,-1,1);
  h_Ry_mm    = new TH2D("h_Ry_mm","h_Ry_mm",200,-0.05,0.25,200,-0.2,0.2);
  h_Rth_mm   = new TH2D("h_Rth_mm","h_Rth_mm",200,-0.05,0.25,200,-0.3,0.3);
  h_Rph_mm   = new TH2D("h_Rph_mm","h_Rph_mm",200,-0.05,0.25,200,-0.2,0.2);
  h_Lp_mm    = new TH2D("h_Lp_mm","h_Lp_mm",200,-0.05,0.25,200,1,3);
  h_Ll_mm    = new TH2D("h_Ll_mm","h_Ll_mm",200,-0.05,0.25,200,20,30);
  h_Ltgy_mm  = new TH2D("h_Ltgy_mm","h_Ltgy_mm",200,-0.05,0.25,200,20,30);
  h_Ltgth_mm = new TH2D("h_Ltgth_mm","h_Ltgth_mm",200,-0.05,0.25,200,-0.2,0.2);
  h_Ltgph_mm = new TH2D("h_Ltgph_mm","h_Ltgph_mm",200,-0.05,0.25,200,-0.2,0.2);
  h_Lvx_mm   = new TH2D("h_Lvx_mm","h_Lvx_mm",200,-0.05,0.25,200,-0.01,0.01);
  h_Lvy_mm   = new TH2D("h_Lvy_mm","h_Lvy_mm",200,-0.05,0.25,200,-0.01,0.01);
  h_Lvz_mm   = new TH2D("h_Lvz_mm","h_Lvz_mm",200,-0.05,0.25,200,-1,1);
  h_Lx_mm    = new TH2D("h_Lx_mm","h_Lx_mm",200,-0.05,0.25,200,-1,1);
  h_Ly_mm    = new TH2D("h_Ly_mm","h_Ly_mm",200,-0.05,0.25,200,-0.2,0.2);
  h_Lth_mm   = new TH2D("h_Lth_mm","h_Lth_mm",200,-0.05,0.25,200,-0.3,0.3);
  h_Lph_mm   = new TH2D("h_Lph_mm","h_Lph_mm",200,-0.05,0.25,200,-0.2,0.2);
  set->SetTH1(h_ct      ,"Coincidence Time","Cointime (ns)","Counts");
  set->SetTH2(h_Rs2x_ct ,"RHRS S2 X-pos v.s Cointime","Cointime (ns)","X (m)");
  set->SetTH2(h_Ls2x_ct ,"LHRS S2 X-pos v.s Cointime","Cointime (ns)","X (m)");
  set->SetTH2(h_m2_ct   ,"RHRS Mass Square v.s Cointime","Cointime (ns)","M^{2} (GeV^{2}/c^{4}");
  set->SetTH2(h_beta_ct ,"RHRS beta v.s Cointime","Cointime (ns)","#beta");
  set->SetTH2(h_dedx_ct ,"RHRS Energy Deposit v.s Cointime","Cointime (ns)","dE/dx ()");
  set->SetTH2(h_a1sum_ct,"RHRS A1 SUM v.s Cointime","Cointime (ns)","");
  set->SetTH2(h_a2sum_ct,"RHRS A1 SUM v.s Cointime","Cointime (ns)","");
  set->SetTH1(h_mm      ,"#Lambda Binding Energy (Gas)","-B_{#Lambda} (GeV/c^{2})","Counts");
  set->SetTH1(h_mmall   ,"#Lambda Binding Energy (w/o Z_{v} cut)","-B_{#Lambda} (GeV/c^{2})","Counts");
  set->SetTH1(h_mmfoil  ,"#Lambda Binding Energy (Al foil)","-B_{#Lambda} (GeV/c^{2})","Counts");
  set->SetTH2(h_Rp_mm   ,"RHRS Momentum v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Momentum (GeV/c)");
  set->SetTH2(h_Rl_mm   ,"RHRS Length v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Length (m)");
  set->SetTH2(h_Rtgy_mm ,"RHRS Target Y v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Y_{t} (m)");
  set->SetTH2(h_Rtgth_mm,"RHRS Target #theta v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","#theta_{t} (rad)");
  set->SetTH2(h_Rtgph_mm,"RHRS Target #phi v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","#phi_{t} (rad)");
  set->SetTH2(h_Rvx_mm  ,"RHRS Vertex X v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","X_{v} (m)");
  set->SetTH2(h_Rvy_mm  ,"RHRS Vertex Y v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Y_{v} (m)");
  set->SetTH2(h_Rvz_mm  ,"RHRS Vertex Z v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Z_{v} (m)");
  set->SetTH2(h_Rx_mm   ,"RHRS FP X v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","X_{FP} (m)");
  set->SetTH2(h_Ry_mm   ,"RHRS FP Y v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Y_{FP} (m)");
  set->SetTH2(h_Rth_mm  ,"RHRS FP #theta v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","#theta_{FP} (rad)");
  set->SetTH2(h_Rph_mm  ,"RHRS FP #phi v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","#phi_{FP} (rad)");
  set->SetTH2(h_Lp_mm   ,"LHRS Momentum v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Momentum (GeV/c)");
  set->SetTH2(h_Ll_mm   ,"LHRS Length v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Length (m)");
  set->SetTH2(h_Ltgy_mm ,"LHRS Target Y v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Y_{t} (m)");
  set->SetTH2(h_Ltgth_mm,"LHRS Target #theta v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","#theta_{t} (rad)");
  set->SetTH2(h_Ltgph_mm,"LHRS Target #phi v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","#phi_{t} (rad)");
  set->SetTH2(h_Lvx_mm  ,"LHRS Vertex X v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","X_{v} (m)");
  set->SetTH2(h_Lvy_mm  ,"LHRS Vertex Y v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Y_{v} (m)");
  set->SetTH2(h_Lvz_mm  ,"LHRS Vertex Z v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Z_{v} (m)");
  set->SetTH2(h_Lx_mm   ,"LHRS FP X v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","X_{FP} (m)");
  set->SetTH2(h_Ly_mm   ,"LHRS FP Y v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","Y_{FP} (m)");
  set->SetTH2(h_Lth_mm  ,"LHRS FP #theta v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","#theta_{FP} (rad)");
  set->SetTH2(h_Lph_mm  ,"LHRS FP #phi v.s B_{#Lambda}","-B_{#Lambda} (GeV/c^{2})","#phi_{FP} (rad)");

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

  return true;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
int main(int argc, char** argv){

  gErrorIgnoreLevel = kError;

  int ch;
  extern char *optarg;
  int MaxNum = 0;
  string ifname("test.list");
  while((ch=getopt(argc,argv,"hRLf:n:b"))!=-1){
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
    case 'b':
      gROOT->SetBatch(1);
      batch = true;
      break;
    case 'h':
      cout<<"-f : input root file"<<endl;
      cout<<"-n : max no. of events"<<endl;
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

  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    Ana->ReadTree(Form("../rootfiles/%s",runname.c_str()));
  }

  Ana->Roop();
  Ana->Draw();               // save histograms to pdf file

  cout<<"Done!"<<endl;

  Ana->Close();
  if(batch){ gSystem->Exit(1); }
  gSystem->Exit(1);
  theApp.Run();
  return 0;

}

