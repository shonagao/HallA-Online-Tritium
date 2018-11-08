

void cointime_accut(){

  gStyle->SetOptStat(0);
  TTree *T = (TTree*)gDirectory->Get("T");

  TString Rarm;
  TString Larm;
  
  Double_t Rs2time[16],Ls2time[16],Rs2lt[16],Rs2rt[16],Ls2lt[16],Ls2rt[16];
  Double_t Ra1asumc,Ra2asumc;

  T->SetBranchAddress("R.s2.time",Rs2time);
  T->SetBranchAddress("L.s2.time",Ls2time);
  T->SetBranchAddress("R.s2.lt",Rs2lt);
  T->SetBranchAddress("R.s2.rt",Rs2rt);
  T->SetBranchAddress("L.s2.lt",Ls2lt);
  T->SetBranchAddress("L.s2.rt",Ls2rt);
  T->SetBranchAddress("R.a1.asum_c",&Ra1asumc);
  T->SetBranchAddress("R.a2.asum_c",&Ra2asumc);

  Double_t bin,b1,b2;
  bin=200;
  b1=-25;
  b2=15;
  // bin=400;
  // b1=-0.25/1000000;
  // b2=0;

  TH1F *cointime = new TH1F("cointime","Coincidence time(Rs2time-Ls2time)",bin,b1,b2);
  TH1F *cointime_a1cut = new TH1F("cointime_a1cut","Coincidence time(Rs2time-Ls2time)",bin,b1,b2);
  TH1F *cointime_a2cut = new TH1F("cointime_a2cut","Coincidence time(Rs2time-Ls2time)",bin,b1,b2);
  TH1F *cointime_accut = new TH1F("cointime_accut","Coincidence time accut(Rs2time-Ls2time)",bin,b1,b2);

  //----------------Get Entry---------------------------

  Int_t nentries = T->GetEntries();
  Bool_t a1cut=false;
  Bool_t a2cut=false;
  Double_t data;
  cout<<"Total Number of Events = "<<nentries<<endl;

  for(Int_t i=0;i<nentries;i++){

    if(i%100000==0) cout << " events processed = " << i << endl;
    T->GetEntry(i);

    if(Ra1asumc>130){
      a1cut = true;
    }else{
      a1cut = false;
    }
    if(Ra2asumc>3000){
      a2cut = true;
    }else{
      a2cut = false;
    }

    if(!a1cut){
      for(Int_t k=0;k<16;k++){
	for(Int_t j=0;j<16;j++){
	  if(Rs2lt[k]>0 && Rs2rt[k]>0 && Ls2lt[j]>0 && Ls2rt[j]>0){
	    data = (Rs2time[k]-Ls2time[j])*1000000000+150.7;
	    cointime->Fill(data);
	    cointime_a1cut->Fill(data);
	  }
	}
      }
    }
    if(a2cut){
      for(Int_t k=0;k<16;k++){
	for(Int_t j=0;j<16;j++){
	  if(Rs2lt[k]>0 && Rs2rt[k]>0 && Ls2lt[j]>0 && Ls2rt[j]>0){
	    data = (Rs2time[k]-Ls2time[j])*1000000000+150.7;
	    cointime->Fill(data);
	    cointime_a2cut->Fill(data);
	  }
	}
      }
    }
    if(!a1cut && a2cut){
      for(Int_t k=0;k<16;k++){
	for(Int_t j=0;j<16;j++){
	  if(Rs2lt[k]>0 && Rs2rt[k]>0 && Ls2lt[j]>0 && Ls2rt[j]>0){
	    data = (Rs2time[k]-Ls2time[j])*1000000000+150.7;
	    cointime->Fill(data);
	    cointime_accut->Fill(data);
	  }
	}
      }
    }
    
  }

  // cointime_a1cut->SetLineColor(2);
  // cointime_a2cut->SetLineColor(3);
  // cointime_accut->SetLineColor(4);

  // cointime->Draw();
  // cointime_a1cut->Draw("same");
  // cointime_a2cut->Draw("same");
  cointime_accut->Draw();
}
