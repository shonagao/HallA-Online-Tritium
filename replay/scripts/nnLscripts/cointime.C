void cointime(Int_t flag){
//  TFile *ifp = new TFile("hogehoge.root","read");

  gStyle->SetOptStat(0);

  TTree *tree = (TTree*)gDirectory->Get("T");

  string L_rftime = "LTDC.F1FirstHit[47] - LTDC.F1FirstHit[40] * 0.05623";
//  string L_s2 = Form("(L.s2.lt_c[%d] + L.s2.rt_c[%d])/2.",(int)L.s2.trpad[0],(int)L.s2.trpad[0]);
  string L_s2 = Form("(L.s2.lt_c[L.s2.trpad[0]] + L.s2.rt_c[L.s2.trpad[0]])/2.");
  string L_time = Form("%s - %s",L_s2.c_str(),L_rftime.c_str());
  string L_beta = "L.tr.p[0] / sqrt(0.000511*0.000511 + L.tr.p[0]*L.tr.p[0])";
  string L_tgtime = Form("%s - (L.tr.pathl[0] + L.s2.trpath[0])/%s",L_time.c_str(),L_beta.c_str());

  string R_rftime = "RTDC.F1FirstHit[15] - RTDC.F1FirstHit[9]  * 0.05623";
  string R_s2 = Form("(R.s2.lt_c[R.s2.trpad[0]] + R.s2.rt_c[R.s2.trpad[0]])/2.");
  string R_time = Form("%s - %s",R_s2.c_str(),R_rftime.c_str());
  string R_beta = "R.tr.p[0] / sqrt(0.5*0.5 + R.tr.p[0]*R.tr.p[0])";
  string R_tgtime = Form("%s - (R.tr.pathl[0] + R.s2.trpath[0])/%s",R_time.c_str(),R_beta.c_str());

  string ctime = Form("%s - %s",L_tgtime.c_str(),R_tgtime.c_str());

  if( flag == 1 ){
    TH1F *h_ctime = new TH1F("h_ctime","h_ctime",800,-40,40);
    h_ctime->SetTitle("Coincidece Time");
    h_ctime->GetXaxis()->SetTitle("L-R (ns)");
    h_ctime->GetYaxis()->SetTitle("Counts / 1ns");
    tree->Draw(Form("%s>>h_time",ctime.c_str()));
  }

  if( flag == 2 ){
    TH2F *h_ctime_Ls2 = new TH2F("h_ctime_Ls2","h_ctime_Ls2",16,0,16,800,-40,40);
    h_ctime_Ls2->SetTitle("Coincidece Time v.s LHRS S2 Paddle");
    h_ctime_Ls2->GetYaxis()->SetTitle("S2 Paddle No.");
    h_ctime_Ls2->GetYaxis()->SetTitle("L-R (ns)");
    tree->Draw(Form("%s:L.s2.trpad[0]>>h_time_Ls2",ctime.c_str()));
  }

  if( flag == 3 ){
    TH2F *h_ctime_Rs2 = new TH2F("h_ctime_Rs2","h_ctime_Rs2",16,0,16,800,-40,40);
    h_ctime_Rs2->SetTitle("Coincidece Time v.s RHRS S2 Paddle");
    h_ctime_Rs2->GetYaxis()->SetTitle("S2 Paddle No.");
    h_ctime_Rs2->GetYaxis()->SetTitle("L-R (ns)");
    tree->Draw(Form("%s:R.s2.trpad[0]>>h_time_Rs2",ctime.c_str()));
  }

}
