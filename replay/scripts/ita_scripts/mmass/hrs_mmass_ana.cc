
void hrs_mmass_ana(){

  TFile *fin=new TFile("../rootfiles/Lambda_small1_mm_new.root","read");
  TH1D* hL=(TH1D*)fin->Get("h_peak_L");
  hL->SetName("hL");

  TF1* fL=new TF1("fL","gausn",1.1,1.3);
  fL->SetNpx(2000);
  TF1* fS=new TF1("fS","gausn",1.1,1.3);
  fS->SetNpx(2000);
  hL->Fit("fL","","",1.1,1.13);
  hL->Fit("fS","","",1.19,1.22);  

  TCanvas* c0=new TCanvas("c0","Lambda Missing Mass Hist");
  c0->cd();
  hL->Draw();
  fL->Draw("same");
  fS->Draw("same");
  //  fin->Close();
}
