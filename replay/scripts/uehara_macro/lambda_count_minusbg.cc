
void lambda_count_minusbg(){

  TFile *fp = new TFile("/adaqfs/home/a-onl/tritium_work/nagao/HallA-Online-Tritium/replay/scripts/nnLscripts/output.root");
  TFile *fpbg = new TFile("/adaqfs/home/a-onl/tritium_work/uehara/tohoku_analysis/HallA-Online-Tritium/replay/scripts/nnLscripts_new/background.root");

  
  TH1F *h1 = (TH1F*)fp->Get("h_mm");
  h1->Rebin(4);
  TH1F *h2 = (TH1F*)fpbg->Get("h_mm");
  h2->Rebin(4);
  h2->Scale(0.016);
  h1->Add(h2,-1);
// -----------parameters for fitting-------------
  
  double L_start = 0.015;
  double L_end = 0.05;
  double sig_start = 0.09 ;
  double sig_end = 0.137;
  double bg_start = -0.1;
  double bg_end = 0.2;

  const int npar = 8;
  double par[npar]={};
  double g_par[6]={3.5,0.03,0.03,2,0.11,0.04};
  double fit_par[npar]={};
  double fit_parerr[npar]={};

// -----------fitting-------------

  TF1 *f_L = new TF1("f_L","gausn",L_start,L_end);
  SetTF1(f_L,3,2,1,500);
  TF1 *f_sig = new TF1("f_sig","gausn",sig_start,sig_end);
  SetTF1(f_sig,4,2,1,500);
  TF1 *f_bg = new TF1("f_bg","pol1",bg_start,bg_end);
  SetTF1(f_bg,5,2,1,500);

  f_L->SetParameter(0,g_par[0]);
  f_L->SetParLimits(0,g_par[0]-2,g_par[0]+2);
  f_L->SetParameter(1,g_par[1]);
  f_L->SetParLimits(1,g_par[1]-0.1,g_par[1]+0.1);
  f_L->SetParameter(2,g_par[2]);
  f_L->SetParLimits(2,g_par[2]-0.1,g_par[2]+0.1);

  f_sig->SetParameter(0,g_par[3]);
  f_sig->SetParLimits(0,g_par[3]-2,g_par[3]+2);
  f_sig->SetParameter(1,g_par[4]);
  f_sig->SetParLimits(1,g_par[4]-0.1,g_par[4]+0.1);
  f_sig->SetParameter(2,g_par[5]);
  f_sig->SetParLimits(2,g_par[5]-0.1,g_par[5]+0.1);

  h1->Fit(f_L,"0","",L_start,L_end);
  h1->Fit(f_sig,"0","",sig_start,sig_end);
  h1->Fit(f_bg,"0","",bg_start,bg_end);

  for(int i=0;i<npar;i++){
    if(i<3){
      par[i] = f_L->GetParameter(i);
    }else if(i<6){
      par[i] = f_sig->GetParameter(i-3);
    }else if(i<npar){
      par[i] = f_bg->GetParameter(i-6);
    }
  }

  TF1 *f1 = new TF1("f1","gausn(0)+gausn(3)+pol1(6)",bg_start,bg_end);
  SetTF1(f1,2,2,1,500);
  for(int i=0;i<npar;i++){
    f1->SetParameter(i,par[i]);
  }

  h1->Fit(f1,"0","",bg_start,bg_end);

  for(int i=0;i<npar;i++){
    fit_par[i] = f1->GetParameter(i);
    fit_parerr[i] = f1->GetParError(i);
  }

  TF1 *f_fitbg1 = new TF1("f_bg1","pol1",bg_start,bg_end);
  SetTF1(f_fitbg1,2,2,1,500);
  for(int i=6;i<8;i++){
    f_fitbg1->SetParameter(i-6,fit_par[i]);
  }


// -----------calculation-------------

  double Ndf;
  double Chisquare;
  double Reduced_Chisquare;

  Ndf = f1->GetNDF();
  Chisquare = f1->GetChisquare();
  Reduced_Chisquare = Chisquare / Ndf;
  
  cout<<Form("NDf : %lf\n",Ndf)<<Form("Chisquare : %lf\n",Chisquare)<<Form("Reduced_Chisquare : %lf",Reduced_Chisquare)<<endl;
  
  double BinWidth;
  BinWidth = h1->GetBinWidth(1);

  double nL,nLe;
  double nsig,nsige;
  
  nL = fit_par[0]/BinWidth;
  nLe = fit_parerr[0]/BinWidth;
  nsig = fit_par[3]/BinWidth;
  nsige = fit_parerr[3]/BinWidth;

  
  cout<<"Lambda : "<<nL<<" +/- "<<nLe<<endl;
  cout<<"Sigma  : "<<nsig<<" +/- "<<nsige<<endl;

// -----------drawing-------------

  h1->GetYaxis()->SetRangeUser(-10,80);
  h1->Draw();
  f1->Draw("same");
  // f_L->Draw("same");
  // f_sig->Draw("same");
  // f_bg->Draw("same");
  f_fitbg1->Draw("same");
 
}
