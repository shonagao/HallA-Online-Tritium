#include "Setting.cc"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus2_pol3bg(double *x, double *par) {
  //par[0]=area (gaus1)
  //par[1]=mean (gaus1)
  //par[2]=sigma (gaus1)
  //par[3]=area (gaus2)
  //par[4]=mean (gaus2)
  //par[5]=sigma (gaus2)
  //par[6]=const bg
  //par[7]=1st bg
  //par[8]=2st bg
  //par[9]=3st bg
  double val;
  double bg;
  double ga1,ga2,ga3;

  ga1=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  ga2=par[3]*TMath::Gaus(x[0],par[4],par[5],1);
  bg =  par[6] + par[7]*x[0] + par[8]*x[0]*x[0] + par[9]*x[0]*x[0]*x[0];
  val = ga1 + ga2 + bg;
  return val;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void fit_mm(){
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.07);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);
  

  TF1 *ga_L, *ga_S, *f_bg, *f_full;

  Setting *set = new Setting();
  ga_L  = new TF1("ga_L","gausn",-0.1,0.2);
  ga_S  = new TF1("ga_S","gausn",-0.1,0.2);
  f_bg  = new TF1("f_bg" ,"pol3",-0.1,0.2);
  f_full = new TF1("f_full",gaus2_pol3bg,-0.1,0.2,10);

  set->SetTF1(ga_L ,  1,2,2);
  set->SetTF1(ga_S ,  4,2,2);
  set->SetTF1(f_bg ,  6,2,2);
  set->SetTF1(f_full, 2,1,2);

  double param[10];



  TFile *ifp = new TFile("../nnLscripts/output.root");

  TH1F *h_mm   = (TH1F*)ifp->Get("h_mm");
  TH1F *h_mmbg = (TH1F*)ifp->Get("h_mmbg");
  h_mm  ->Rebin(2);
  h_mmbg->Rebin(2);
  h_mm  ->GetXaxis()->SetRangeUser(-0.1,0.2);
  h_mmbg->GetXaxis()->SetRangeUser(-0.1,0.2);
  double binw = h_mm->GetBinWidth(1);
  cout<<"bin width "<<binw<<endl;

  h_mm->GetYaxis()->SetTitle(Form("counts/%.3lfMeV",1000*binw));

  for(int i=0;i<5;i++){
    h_mmbg ->Fit(f_bg,"0QR","",-0.1,-0.2);
  }
  double bg[4];
  f_bg->GetParameters(&bg[0]);

  param[0] = 500.*binw; 
  param[1] = 0.;
  param[2] = 0.0010;
  param[3] = 200.*binw;
  param[4] = 0.084;
  param[5] = 0.0010;
  param[6] = bg[0]; 
  param[7] = bg[1]; 
  param[8] = bg[2]; 
  param[9] = bg[3];

  f_full -> SetParameters(&param[0]);
  //f_full -> FixParameter(4,param[4]);
  //f_full -> FixParameter(4,param[4]);
  //f_full -> FixParameter(4,param[4]);
  //f_full -> FixParameter(4,param[4]);
  //f_full -> SetParLimits(3,0.5*param[3],2.*param[3]);
  //f_full -> SetParLimits(4,0.5*param[4],2.*param[4]);
  //f_full -> SetParLimits(5,0.5*param[5],2.*param[5]);
  //f_full -> FixParameter(5,param[5]);
  
  
  h_mm ->Fit(f_full ,"0QR","",-0.1, 0.2);

  cout<<"Lambda peak param"   <<param[0] <<" "<<param[1]<<" "<<param[2]<<endl;
  cout<<"Sigma0 peak param"   <<param[3] <<" "<<param[4]<<" "<<param[5]<<endl;
  f_full -> GetParameters(&param[0]);
  cout<<"----after full range fit-------"<<endl;
  cout<<"Lambda peak param"   <<param[0] <<" "<<param[1]<<" "<<param[2]<<endl;
  cout<<"Sigma0 peak param"   <<param[3] <<" "<<param[4]<<" "<<param[5]<<endl;

  ga_L ->SetParameters(&param[0]);
  ga_S  ->SetParameters(&param[3]);  
  f_bg  ->SetParameters(&param[6]);
  

  double NL,NS;
  double eNL,eNS;
  NL = param[0]/binw;
  NS = param[3]/binw;
  
  eNL = f_full->GetParError(0)/binw;
  eNS = f_full->GetParError(3)/binw;

  cout<<"Num. of Lambda :"<<NL<<" +/- "<<eNL<<endl;
  cout<<"Num. of Sigma0 :"<<NS <<" +/- "<<eNS <<endl;

  //counting num of Kaon by integrating bin content - BG (for double check)
  /*
  double INT_pi = h->Integral(h->FindBin(-2.),h->FindBin(2.)) - h->Integral(h->FindBin(-10.),h->FindBin(-6.));
  double INT_K = h->Integral(h->FindBin(2.),h->FindBin(3.5))  - h->Integral(h->FindBin(-10.),h->FindBin(-8.5));

  cout<<"--------double check of num----------"<<endl;
  cout<<"Num. of pi :"<<INT_pi<<endl;
  cout<<"Num. of K :"<<INT_K  <<endl;
 */ 

  TCanvas *c[3];
  for(int i=0;i<3;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200,900);
  }

  c[0]->Clear();
  c[0]->cd(1);
  h_mm->Draw();
  ga_L  ->Draw("same");
  ga_S   ->Draw("same");  
  f_bg   ->Draw("same");
  f_full ->Draw("same");

  c[1]->Clear();
  c[1]->cd(1);
  h_mm->Draw();
  h_mmbg->Draw("same");
  f_full ->Draw("same");

  c[2]->Clear();
  c[2]->cd(1);
  h_mm->Draw();
  ga_L  ->Draw("same");
  ga_S   ->Draw("same");  
  f_full ->Draw("same");


//  c[0]->Print("pdf/mm_fit.pdf[");
//  c[0]->Print("pdf/mm_fit.pdf");
//  c[1]->Print("pdf/mm_fit.pdf");
//  c[2]->Print("pdf/mm_fit.pdf");
//  c[2]->Print("pdf/mm_fit.pdf]");
}
