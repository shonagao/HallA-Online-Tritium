#include "Setting.cc"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus3_pol0bg(double *x, double *par) {
  //par[0]=area (gaus1)
  //par[1]=mean (gaus1)
  //par[2]=sigma (gaus1)
  //par[3]=area (gaus2)
  //par[4]=mean (gaus2)
  //par[5]=sigma (gaus2)
  //par[6]=area (gaus3)
  //par[7]=mean (gaus3)
  //par[8]=sigma (gaus3)
  //par[9]=const bg
  double val;
  double bg;
  double ga1,ga2,ga3;

  ga1=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  ga2=par[3]*TMath::Gaus(x[0],par[4],par[5],1);
  ga3=par[6]*TMath::Gaus(x[0],par[7],par[8],1);
  bg =  par[9];
  val = ga1 + ga2 + ga3 + bg;
  return val;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void fit_coin(){
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);
  

  TF1 *ga_pi, *ga_k, *ga_p, *f_bg, *f_full;

  Setting *set = new Setting();
  ga_pi = new TF1("ga_pi","gausn(0)+pol0(3)",-10,10);
  ga_k  = new TF1("ga_k" ,"gausn(0)+pol0(3)",  0, 5);
  ga_p  = new TF1("ga_p" ,"gausn(0)+pol0(3)",  0,15);
  f_bg  = new TF1("f_bg" ,"pol0",-15,15);
  f_full = new TF1("f_full",gaus3_pol0bg,-15,15,10);

  set->SetTF1(ga_pi,  1,2,2);
  set->SetTF1(ga_k ,  4,2,2);
  set->SetTF1(ga_p ,  6,2,2);
  set->SetTF1(f_bg ,  1,2,2);
  set->SetTF1(f_full, 2,1,2);

  double param[10];



  TFile *ifp = new TFile("pdf/coin_all.root");

  TH1F *h = (TH1F*)ifp->Get("h_s2ccoin_f1acz");
  h->Rebin(4);
  double binw = h->GetBinWidth(1);
  cout<<"bin width "<<binw<<endl;

  h->GetYaxis()->SetTitle(Form("counts/%.0lfps",1000*binw));

  h ->Fit(f_bg,"0QR","",-10,-5);
  double bg=f_bg->GetParameter(0);

  //ga_pi -> SetParameter(0,100.);
  ga_pi -> SetParameter(1,0.);
  ga_pi -> SetParameter(2,0.4);
  ga_pi -> SetParameter(3,bg);
  //ga_pi -> FixParameter(3,bg);

  //ga_k  -> SetParameter(0,10.);
  ga_k  -> SetParameter(1,3.);
  ga_k  -> SetParameter(2,0.4);
  ga_k  -> SetParameter(3,bg);
  //ga_k  -> FixParameter(3,bg);

  //ga_p  -> SetParameter(0,10.);
  ga_p  -> SetParameter(1,10.);
  ga_p  -> SetParameter(2,1.0);
  ga_p  -> SetParameter(3,bg);
  //ga_p  -> FixParameter(3,bg);


  h ->Fit(ga_pi,"0QR","",-1.5, 1.0);
  h ->Fit(ga_k ,"0QR","", 2.3, 8.2);
  h ->Fit(ga_p ,"0QR","", 10., 15.);

  param[0] = ga_pi->GetParameter(0); 
  param[1] = ga_pi->GetParameter(1);
  param[2] = ga_pi->GetParameter(2);
  param[3] = ga_k ->GetParameter(0);
  param[4] = ga_k ->GetParameter(1);
  param[5] = ga_k ->GetParameter(2);
  param[6] = ga_p ->GetParameter(0); 
  param[7] = ga_p ->GetParameter(1); 
  param[8] = ga_p ->GetParameter(2); 
  param[9] = bg;

  f_full -> SetParameters(&param[0]);
  //f_full -> FixParameter(4,param[4]);
  //f_full -> SetParLimits(3,0.5*param[3],2.*param[3]);
  //f_full -> SetParLimits(4,0.5*param[4],2.*param[4]);
  //f_full -> SetParLimits(5,0.5*param[5],2.*param[5]);
  f_full -> FixParameter(5,param[5]);
  h ->Fit(f_full ,"0QR","",-10., 15.);

  cout<<"Pion peak param"   <<param[0] <<" "<<param[1]<<" "<<param[2]<<endl;
  cout<<"Kaon peak param"   <<param[3] <<" "<<param[4]<<" "<<param[5]<<endl;
  cout<<"Proton peak param" <<param[6] <<" "<<param[7]<<" "<<param[8]<<endl;
  f_full -> GetParameters(&param[0]);
  cout<<"-----------"<<endl;
  cout<<"Pion peak param"   <<param[0] <<" "<<param[1]<<" "<<param[2]<<endl;
  cout<<"Kaon peak param"   <<param[3] <<" "<<param[4]<<" "<<param[5]<<endl;
  cout<<"Proton peak param" <<param[6] <<" "<<param[7]<<" "<<param[8]<<endl;

  ga_pi ->SetParameters(&param[0]);
  ga_k  ->SetParameters(&param[3]);  
  ga_p  ->SetParameters(&param[6]);
  f_bg  ->SetParameters(&param[9]);
  ga_pi ->SetParameter(3,0.);
  ga_k  ->SetParameter(3,0.);  
  ga_p  ->SetParameter(3,0.);
  

  double Npi,NK,Np;
  double eNpi,eNK,eNp;
  Npi = param[0]/binw;
  NK  = param[3]/binw;
  Np  = param[6]/binw;
  
  eNpi = f_full->GetParError(0)/binw;
  eNK  = f_full->GetParError(3)/binw;
  eNp  = f_full->GetParError(6)/binw;

  cout<<"Num. of pi :"<<Npi<<" +/- "<<eNpi<<endl;
  cout<<"Num. of K+ :"<<NK <<" +/- "<<eNK <<endl;
  cout<<"Num. of p  :"<<Np <<" +/- "<<eNp <<endl;

  //counting num of Kaon by integrating bin content - BG (for double check)
  double INT_pi = h->Integral(h->FindBin(-2.),h->FindBin(2.)) - h->Integral(h->FindBin(-10.),h->FindBin(-6.));
  double INT_K = h->Integral(h->FindBin(2.),h->FindBin(3.5))  - h->Integral(h->FindBin(-10.),h->FindBin(-8.5));

  cout<<"--------double check of num----------"<<endl;
  cout<<"Num. of pi :"<<INT_pi<<endl;
  cout<<"Num. of K :"<<INT_K  <<endl;
  

  TCanvas *c[3];
  for(int i=0;i<3;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200,900);
  }

  c[0]->Clear();
  c[0]->cd(1);
  h->Draw();
  ga_pi  ->Draw("same");
  ga_k   ->Draw("same");  
  ga_p   ->Draw("same");
  f_bg   ->Draw("same");
  f_full ->Draw("same");

  c[1]->Clear();
  c[1]->cd(1);
  h->Draw();
  f_full ->Draw("same");

  c[2]->Clear();
  c[2]->cd(1);
  h->Draw();
  ga_pi  ->Draw("same");
  ga_k   ->Draw("same");  
  ga_p   ->Draw("same");
  f_full ->Draw("same");


  c[0]->Print("pdf/coin_fit.pdf[");
  c[0]->Print("pdf/coin_fit.pdf");
  c[1]->Print("pdf/coin_fit.pdf");
  c[2]->Print("pdf/coin_fit.pdf");
  c[2]->Print("pdf/coin_fit.pdf]");
}
