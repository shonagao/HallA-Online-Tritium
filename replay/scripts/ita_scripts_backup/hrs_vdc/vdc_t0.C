#include "rootalias.h"
void vdc_t0(Int_t runnum){

  gStyle->SetOptStat(0);
  double nbin=500;
  double max =2000;
  double min =1000;
  TChain *t   = LoadRun(runnum);
  TCanvas *c1 = new TCanvas("c1","",800,600);
  TH1F *h1=new TH1F("h1","T1",nbin,min,max);
  TH1F *h2=new TH1F("h2","T2",nbin,min,max);
  TH1F *h3=new TH1F("h3","T3",nbin,min,max);
  
  t->Draw("L.vdc.u1.rawtime>>h1",track_L+"DL.bit1>0","");
  t->Draw("L.vdc.u1.rawtime>>h2",track_L+"DL.bit2>0","same");
  t->Draw("L.vdc.u1.rawtime>>h3",track_L+"DL.bit3>0","same");
  h1->SetLineColor(kRed);
  h2->SetLineColor(kGreen);
  h3->SetLineColor(kBlue);

  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetTextSize(0.035);
  leg->SetTextFont(60);
  leg->AddEntry(h1,"","l");
  leg->AddEntry(h2,"","l");
  leg->AddEntry(h3,"","l");
  leg->Draw();
}
