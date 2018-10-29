// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
void R_S2S0TOF(int runnum=93926, int seg=0){
  TChain *tree = new TChain("T");
  tree->Add(Form("rootfiles/tritium_online_%d.root",runnum));

  TH1D *h1 = new TH1D("h1","h1",1000,-90,10);
  TH2D *h2_s0t = new TH2D("h2_s0t","h2_s0t",200,-100,3900,200,-90,10);
  TH2D *h2_s0b = new TH2D("h2_s0b","h2_s0b",200,-100,3900,200,-90,10);
  TH2D *h2_s2l = new TH2D("h2_s2l","h2_s2l",200,-100, 900,200,-90,10);
  TH2D *h2_s2r = new TH2D("h2_s2r","h2_s2r",200,-100, 900,200,-90,10);
  SetTH1(h1,"TOF","S2 - S0 (ns)","Counts / 0.1 ns");
  SetTH2(h2_s0t,"TOF v.s S0t ADC(FastBus)","S0t ADC","S2 - S0 (ns)");
  SetTH2(h2_s0b,"TOF v.s S0b ADC(FastBus)","S0b ADC","S2 - S0 (ns)");
  SetTH2(h2_s2l,"TOF v.s S2l ADC(FastBus)","S2l ADC","S2 - S0 (ns)");
  SetTH2(h2_s2r,"TOF v.s S2r ADC(FastBus)","S2r ADC","S2 - S0 (ns)");

  string F1 = "RTDC.F1FirstHit";
  double fact = 0.05623;
  string S0t = Form("%s[%d]",F1.c_str(),43);
  string S0b = Form("%s[%d]",F1.c_str(),44);
  string S2l = Form("%s[%d]",F1.c_str(),16+seg);
  string S2r = Form("%s[%d]",F1.c_str(),48+seg);
  string L1  = Form("%s[%d]",F1.c_str(),9);
  string L2  = Form("%s[%d]",F1.c_str(),46);

  string S0at = "R.s0.la_c";
  string S0ab = "R.s0.ra_c";
  string S2al = Form("R.s2.la_c[%d]",seg);
  string S2ar = Form("R.s2.ra_c[%d]",seg);

  string cut = Form("%s>0 && %s>0 && %s>0 && %s>0 && %s>0 && %s>0"
            ,S0t.c_str(),S0b.c_str(),S2l.c_str(),S2r.c_str(),L1.c_str(),L2.c_str());

  string cut_adc = Form("abs(%s-1100)<100 && abs(%s-1100)<100 && abs(%s-150)<30 && abs(%s-150)<30"
            ,S0at.c_str(),S0ab.c_str(),S2al.c_str(),S2ar.c_str());

  string tof = Form("(((%s-%s)+(%s-%s))/2. - ((%s-%s)+(%s-%s))/2.)*%lf"
                    ,S2l.c_str(),L1.c_str(),S2r.c_str(),L2.c_str()
                    ,S0t.c_str(),L2.c_str(),S0b.c_str(),L2.c_str(),fact);

  tree->Project("h1",Form("%s",tof.c_str()),Form("%s",cut.c_str()));
//  tree->Project("h1",Form("%s",tof.c_str()),Form("%s && %s",cut.c_str(),cut_adc.c_str()));
  tree->Project("h2_s0t",Form("%s:%s",tof.c_str(),S0at.c_str()),Form("%s",cut.c_str()));
  tree->Project("h2_s0b",Form("%s:%s",tof.c_str(),S0ab.c_str()),Form("%s",cut.c_str()));
  tree->Project("h2_s2l",Form("%s:%s",tof.c_str(),S2al.c_str()),Form("%s",cut.c_str()));
  tree->Project("h2_s2r",Form("%s:%s",tof.c_str(),S2ar.c_str()),Form("%s",cut.c_str()));

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd()->SetMargin(0.15,0.15,0.15,0.10);
  gPad->SetLogy();
  h1->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c2->Divide(2,2,1E-5,1E-5);
  c2->cd(1)->SetMargin(0.15,0.15,0.15,0.10);  gPad->SetLogz();  h2_s0t->Draw("colz");
  c2->cd(2)->SetMargin(0.15,0.15,0.15,0.10);  gPad->SetLogz();  h2_s0b->Draw("colz");
  c2->cd(3)->SetMargin(0.15,0.15,0.15,0.10);  gPad->SetLogz();  h2_s2l->Draw("colz");
  c2->cd(4)->SetMargin(0.15,0.15,0.15,0.10);  gPad->SetLogz();  h2_s2r->Draw("colz");

}

