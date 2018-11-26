//make 3 plots for checking current z-position reconstruction
//fit rpl.z(z-position of reaction point reconstruction)
//get sigma and residual of each peak

//change below
//When you cannot fit z-reconstruction well, change the values of peak_range_start[i][j] and peak_range_end[i][j] that are fitting range of each peak. 
//When you compare the runs which are different meam momentum or theta, change peak_range_start[i][j] and peak_range_end[i][j] for each i. 


void multifoil(){

// -----------calling rootfiles & settings-------------

  int read_nf; 

  cout<<"Enter the number of rootfiles you want to analyze:";
  cin>>read_nf;

  const int nf = read_nf;
  int runnum[nf];

  for(int i=0;i<read_nf;i++){
    runnum[i]=0;

    if(runnum[i]==0){
      cout<<"Enter run number:";
      cin>>runnum[i];

      if(runnum[i]==0){
	i--;
      }
    }
  }

  TFile *fp[nf];
  for(int i=0;i<nf;i++){
    // fp[i] = new TFile(Form(" /adaqfs/home/a-onl/tritium_work/uehara/tohoku_analysis/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",runnum[i]));
    fp[i] = new TFile(Form("/chafs1/work1/tritium/Rootfiles/tritium_%d.root",runnum[i]));
  }

  const int np=10; //the number of peaks


// -----------making histgram of rpl.z-------------

  TTree *tree[nf];
  TH1D *h[nf];
  for(int i=0;i<nf;i++){
    tree[i]= (TTree*)fp[i]->Get("T");
    h[i]= new TH1D(Form("h[%d]",i),Form("f[%d]",i),400,-20,20);
    SetTH1(h[i],"rpl.z","Z on beam line [cm]","Counts",1.,3002,i+2);
    tree[i]->Project(Form("h[%d]",i),"rpl.z*100","");
    
  }


// -----------parameters for fitting-------------

  double peak_range_start[nf][np+1];
  double peak_range_end[nf][np+1];
  
  for(int i=0;i<nf;i++){
    for(int j=0;j<np+1;j++){
      peak_range_start[i][j]=-2000;
      peak_range_end[i][j]=-2000;
    }
    peak_range_start[i][0]=-20;   //fit range
    peak_range_start[i][1]=-10.4;
    peak_range_start[i][2]=-8.1;
    peak_range_start[i][3]=-6.2;
    peak_range_start[i][4]=-4.4;
    peak_range_start[i][5]=-2.4;
    peak_range_start[i][6]=-0.7;
    peak_range_start[i][7]=1.4;
    peak_range_start[i][8]=3.2;
    peak_range_start[i][9]=7;
    peak_range_start[i][10]=9;

    peak_range_end[i][0]=20;
    peak_range_end[i][1]=-9.1;
    peak_range_end[i][2]=-7;
    peak_range_end[i][3]=-5;
    peak_range_end[i][4]=-3.3;
    peak_range_end[i][5]=-1.2;
    peak_range_end[i][6]=0.8;
    peak_range_end[i][7]=2.4;
    peak_range_end[i][8]=4.8;
    peak_range_end[i][9]=8.4;
    peak_range_end[i][10]=10.6;

  }

  double par[nf][np*3+2];
  double fit_par[nf][np*3+2];
  double init_pos = -11.5;
  double init_dist = 2.5;

  double bg[nf][2];
  
  for(int i=0;i<nf;i++){
    for(int j=0;j<np*3+2;j++){
      par[i][j] = 0;
      fit_par[i][j]=0;
    }
  }
  for(int i=0;i<nf;i++){
    for(int j=0;j<np*3+2;j++){
      if(j>1){
	if(j%3==2){ //const
	  par[i][j] = 50000;
	}
	if(j%3==0){ //mean
	  if(j<26){
	    par[i][j] = init_dist*j/3+init_pos;
	  }else{
	    par[i][j] = init_dist*(j+3)/3+init_pos;
	  }
	}
	if(j%3==1){ //sigma
	  par[i][j] = 0.3;
	}
      }
    }
  }


// -----------fitting-------------


  TF1 *f[nf][np+4];
  for(int i=0;i<nf;i++){
    for(int j=0;j<np+2;j++){
      if(j==0){
	f[i][j] = new TF1(Form("f%d%d",i,j),"[0]*x+[1]",peak_range_start[i][j],peak_range_end[i][j]); //background
	SetTF1(f[i][j],i+2,2,1,500);
	h[i]->Fit(f[i][j],"0","",peak_range_start[i][j],peak_range_end[i][j]);
	fit_par[i][0] = f[i][j]->GetParameter(0);
	fit_par[i][1] = f[i][j]->GetParameter(1);
      }else{
	f[i][j] = new TF1(Form("f%d%d",i,j),"gaus",peak_range_start[i][j],peak_range_end[i][j]); //each peak
	SetTF1(f[i][j],i+2,2,1,500);
	f[i][j]->SetParLimits(0,par[i][j*3-1]-50000,par[i][j*3-1]+50000);  //const
	f[i][j]->SetParameter(1,par[i][j*3]);                              //mean
	f[i][j]->SetParLimits(2,0.1,0.6);                                  //sigma
	f[i][j]->SetParameter(2,par[i][j*3+1]);                            //sigma
	h[i]->Fit(f[i][j],"0","",peak_range_start[i][j],peak_range_end[i][j]);
	fit_par[i][j*3-1] = f[i][j]->GetParameter(0);    //const
	fit_par[i][j*3] = f[i][j]->GetParameter(1);      //mean
	fit_par[i][j*3+1] = f[i][j]->GetParameter(2);    //sigma
      }
    }

    f[i][11] = new TF1(Form("f[%d]11",i),"[0]*x+[1] + gaus(2) + gaus(5) + gaus(8) + gaus(11) + gaus(14) + gaus(17) + gaus(20) + gaus(23) + gaus(26) + gaus(29)",-20,20);
    SetTF1(f[i][11],i+2,2,1,500);

    for(int j=0;j<np*3+2;j++){
      
      f[i][11]->SetParameter(j,fit_par[i][j]);

      if(j>=2){
	if(j%3==2){ //const
	  f[i][11]->SetParLimits(j,0,100000);
	}
	if(j%3==0){ //mean
	  f[i][11]->SetParLimits(j,fit_par[i][j]-1,fit_par[i][j]+1);
	}
	if(j%3==1){ //sigma
	  f[i][11]->SetParLimits(j,0.1,0.6);
	}
      }
    }

    h[i]->Fit(f[i][11],"0","",-20,20);

    bg[i][0] = f[i][11]->GetParameter(0);
    bg[i][1] = f[i][11]->GetParameter(1);

    f[i][12] = new TF1(Form("f%d12",i),"[0]*x+[1]",peak_range_start[i][0],peak_range_end[i][0]); //background
    SetTF1(f[i][12],i+2,2,2,500);
    f[i][12]->SetParameter(0,bg[i][0]);
    f[i][12]->SetParameter(1,bg[i][1]);

  }

// -----------getting parameter(bg & sigma & mean & reduced-chisquare)-------------
  
  double sigma[nf][np];
  double FWHM[nf][np];
  double mean[nf][np];
  double esigma[nf][np];
  double eFWHM[nf][np];
  double emean[nf][np];
  double real_pos[nf][np];
  double def_pos[nf][np];
  double Ndf[nf];
  double Chisquare[nf];
  double Reduced_Chisquare[nf];
  
  int s[nf];
  int m[nf];

  for(int i=0;i<nf;i++){
    
    for(int j=0;j<np;j++){
      s[i] = j*3+4;
      m[i] = j*3+3;
      sigma[i][j] = f[i][11]->GetParameter(s[i]);
      FWHM[i][j] = 2*sqrt(2*log(2))*sigma[i][j];
      mean[i][j] = f[i][11]->GetParameter(m[i]);
      esigma[i][j] = f[i][11]->GetParError(s[i]);
      eFWHM[i][j] = 2*sqrt(2*log(2))*esigma[i][j];
      emean[i][j] = f[i][11]->GetParError(m[i]);
      if(j<8){
	real_pos[i][j] = j*2.5-12.5;
      }else{
	real_pos[i][j] = (j+1)*2.5-12.5;
      }
      def_pos[i][j] = real_pos[i][j]-mean[i][j];
    }

    Ndf[i] = f[i][11]->GetNDF();
    Chisquare[i] = f[i][11]->GetChisquare();
    Reduced_Chisquare[i] = Chisquare[i]/Ndf[i];
    cout<<Form("\nNDf%d:%lf\n",i,Ndf[i])<<Form("Chisquare%d:%lf\n",i,Chisquare[i])<<Form("Reduced_Chisquare%d:%lf\n",i,Reduced_Chisquare[i])<<endl;
  }


// -----------z-position depedence of FWHM-------------

  TGraphErrors *g_FWHM[nf];
  for(int i=0;i<nf;i++){
    g_FWHM[i] = new TGraphErrors(10,mean[i],FWHM[i],emean[i],eFWHM[i]);
    SetGr(g_FWHM[i],"z-position dependence of FWHM","Z on beam line [cm]","FWHM [cm]",i+2,1,1,i+2,i+20,1);
  }


// -----------z-position depedence of residual-------------

  TGraphErrors *g_res[nf];
  for(int i=0;i<nf;i++){
    g_res[i] = new TGraphErrors(10,real_pos[i],def_pos[i],0,emean[i]);
    SetGr(g_res[i],"z-position dependence of residual","Z on beam line [cm]","real - exp [cm]",i+2,1,1,i+2,i+20,1);
  }


// -----------legend-------------

  TLegend *leg1 = new TLegend(0.75,0.9-nf*0.05,0.87,0.9,"");
  for(int i=0;i<nf;i++){
    leg1->AddEntry(f[i][11],Form("%d",runnum[i]),"l");
  }
  leg1->SetBorderSize(1);
  leg1->SetFillColor(0);

  TLegend *leg2 = new TLegend(0.75,0.9-nf*0.05,0.87,0.9,"");
  for(int i=0;i<nf;i++){
    leg2->AddEntry(g_FWHM[i],Form("%d",runnum[i]),"p");
  }
  leg2->SetBorderSize(1);
  leg2->SetFillColor(0);

  TLegend *leg3 = new TLegend(0.75,0.9-nf*0.05,0.87,0.9,"");
  for(int i=0;i<nf;i++){
    leg3->AddEntry(g_res[i],Form("%d",runnum[i]),"p");
  }
  leg3->SetBorderSize(1);
  leg3->SetFillColor(0);


// -----------drawing-------------

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd(1)->SetMargin(0.15,0.10,0.15,0.07);
  for(int i=0;i<nf;i++){
    if(i==0){
      h[i]->Draw();
      f[i][11]->Draw("same");
      f[i][12]->Draw("same");
    }else{
      h[i]->Draw("same");
      f[i][11]->Draw("same");
      f[i][12]->Draw("same");
    }
  }
  leg1->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",1600,600);
  c2->Divide(2,1,1E-5,1E-5);
  c2->cd(1)->SetMargin(0.15,0.10,0.15,0.07);
  // g_FWHM[0]->SetMinimum(0.4);
  for(int i=0;i<nf;i++){
     if(i==0){
      g_FWHM[i]->Draw("AP");
    }else{
      g_FWHM[i]->Draw("P");
    }
  }
  leg2->Draw();

  c2->cd(2)->SetMargin(0.15,0.10,0.15,0.07);
  // g_res[0]->SetMinimum(1);
  // g_res[0]->SetMaximum(-1);
  for(int i=0;i<nf;i++){

    if(i==0){
      g_res[i]->Draw("AP");
    }else{
      g_res[i]->Draw("P");
    }
  }
  leg3->Draw();
}
