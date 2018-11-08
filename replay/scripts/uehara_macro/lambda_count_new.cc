
void lambda_count_new(){

  TFile *fp = new TFile("/adaqfs/home/a-onl/tritium_work/nagao/HallA-Online-Tritium/replay/scripts/nnLscripts/output.root");
  TFile *fpbg = new TFile("/adaqfs/home/a-onl/tritium_work/uehara/tohoku_analysis/HallA-Online-Tritium/replay/scripts/nnLscripts_new/background.root");

  
  TH1F *h1 = (TH1F*)fp->Get("h_mm");
  h1->Rebin(4);
  TH1F *h2 = (TH1F*)fpbg->Get("h_mm");
  h2->Rebin(4);

// -----------parameters for fitting-------------
  
  double L_start = 0.015;
  double L_end = 0.047;
  double sig_start = 0.09 ;
  double sig_end = 0.137;
  double bg1_start = -0.06;
  double bg1_end = 0.07;
  double bg2_start = 0.08;
  double bg2_end = 0.2;

  const int npar = 10;
  double par[npar]={};
  double g_par[6]={2.5,0.03,0.03,2,0.11,0.04};
  double fit_par[npar]={};
  double fit_parerr[npar]={};

// -----------fitting-------------

  TF1 *f_L = new TF1("f_L","gausn",L_start,L_end);
  SetTF1(f_L,3,2,1,500);
  TF1 *f_bg1 = new TF1("f_bg1","pol1",bg1_start,bg1_end);
  SetTF1(f_bg1,5,2,1,500);
  TF1 *f_sig = new TF1("f_sig","gausn",sig_start,sig_end);
  SetTF1(f_sig,4,2,1,500);
  TF1 *f_bg2 = new TF1("f_bg2","pol1",bg2_start,bg2_end);
  SetTF1(f_bg2,5,2,1,500);


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
  h1->Fit(f_bg1,"0","",bg1_start,bg1_end);
  h1->Fit(f_sig,"0","",sig_start,sig_end);
  h1->Fit(f_bg2,"0","",bg2_start,bg2_end);

  //koko
  f_bg2->SetParameter(0,58);
  f_bg2->SetParameter(1,-280);


  for(int i=0;i<npar;i++){
    if(i<3){
      par[i] = f_L->GetParameter(i);
    }else if(i<5){
      par[i] = f_bg1->GetParameter(i-3);
    }else if(i<8){
      par[i] = f_sig->GetParameter(i-5);
    }else if(i<npar){
      par[i] = f_bg2->GetParameter(i-8);
    }
  }


  TF1 *f1 = new TF1("f1","gausn(0)+pol1(3)",bg1_start,bg1_end);
  SetTF1(f1,2,2,1,500);
  for(int i=0;i<5;i++){
    f1->SetParameter(i,par[i]);
  }
  TF1 *f2 = new TF1("f2","gausn(0)+pol1(3)",bg2_start,bg2_end);
  SetTF1(f2,2,2,1,500);
  for(int i=5;i<npar;i++){
    f2->SetParameter(i-5,par[i]);
  }


  f2->SetParLimits(3,par[8]-0.01,par[8]+0.01);
  f2->SetParLimits(4,par[9]-0.01,par[9]+0.01);


  h1->Fit(f1,"0","",bg1_start,bg1_end);
  h1->Fit(f2,"0","",bg2_start,bg2_end);

  for(int i=0;i<npar;i++){
    if(i<5){
      fit_par[i] = f1->GetParameter(i);
      fit_parerr[i] = f1->GetParError(i);
    }else if(i<npar){
      fit_par[i] = f2->GetParameter(i-5);
      fit_parerr[i] = f2->GetParError(i-5);
    }
  }

  TF1 *f_fitbg1 = new TF1("f_bg1","pol1",bg1_start,bg1_end);
  SetTF1(f_fitbg1,2,2,1,500);
  TF1 *f_fitbg2 = new TF1("f_bg2","pol1",bg2_start,bg2_end);
  SetTF1(f_fitbg2,2,2,1,500);
  for(int i=3;i<5;i++){
    f_fitbg1->SetParameter(i-3,fit_par[i]);
  }
  for(int i=8;i<10;i++){
    f_fitbg2->SetParameter(i-8,fit_par[i]);
  }


// -----------calculation-------------
  
  double BinWidth;
  BinWidth = h1->GetBinWidth(1);

  double nL,nLe;
  double nsig,nsige;
  
  nL = fit_par[0]/BinWidth;
  nLe = fit_parerr[0]/BinWidth;
  nsig = fit_par[5]/BinWidth;
  nsige = fit_parerr[5]/BinWidth;

  
  cout<<"Lambda : "<<nL<<" +/- "<<nLe<<endl;
  cout<<"Sigma  : "<<nsig<<" +/- "<<nsige<<endl;

// -----------drawing-------------
  h2->Scale(0.015);

  h1->Draw();
  h2->Draw("same");
  f1->Draw("same");
  f2->Draw("same");
  // f_L->Draw("same");
  // f_sig->Draw("same");
  // f_bg1->Draw("same");
  // f_bg2->Draw("same");
  f_fitbg1->Draw("same");
  f_fitbg2->Draw("same");
 
}
