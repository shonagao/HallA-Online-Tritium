// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
void R_S2S0TOF(int seg=0){
  TChain *tree = new TChain("T");
  tree->Add("rootfiles/tritium_online_94000.root");
//  tree->Add("rootfiles/tritium_online_93926.root");

  TH1D *h1 = new TH1D("h1","h1",600,-20,10);
  string F1 = "RTDC.F1FirstHit";
  double fact = 0.005623;
  string S0t = Form("%s[%d]",F1.c_str(),43);
  string S0b = Form("%s[%d]",F1.c_str(),44);
  string S2l = Form("%s[%d]",F1.c_str(),16+seg);
  string S2r = Form("%s[%d]",F1.c_str(),48+seg);
  string L1  = Form("%s[%d]",F1.c_str(),9);
  string L2  = Form("%s[%d]",F1.c_str(),46);
  string cut = Form("%s>0 && %s>0 && %s>0 && %s>0 && %s>0 && %s>0"
            ,S0t.c_str(),S0b.c_str(),S2l.c_str(),S2r.c_str(),L1.c_str(),L2.c_str());
  tree->Project("h1",Form("(((%s-%s)+(%s-%s))/2. - ((%s-%s)+(%s-%s))/2.)*%lf"
                    ,S2l.c_str(),L1.c_str(),S2r.c_str(),L2.c_str()
                    ,S0t.c_str(),L2.c_str(),S0b.c_str(),L2.c_str(),fact)
                    ,Form("%s",cut.c_str()));

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  gPad->SetLogy();
  h1->Draw();
}

