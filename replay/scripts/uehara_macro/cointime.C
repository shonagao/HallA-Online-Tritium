

void cointime(){

  gStyle->SetOptStat(0);
  TTree *T = (TTree*)gDirectory->Get("T");

  TString Rarm;
  TString Larm;
  
  Double_t Rs2time[16],Ls2time[16],Rs2lt[16],Rs2rt[16],Ls2lt[16],Ls2rt[16];

  T->SetBranchAddress("R.s2.time",Rs2time);
  T->SetBranchAddress("L.s2.time",Ls2time);
  T->SetBranchAddress("R.s2.lt",Rs2lt);
  T->SetBranchAddress("R.s2.rt",Rs2rt);
  T->SetBranchAddress("L.s2.lt",Ls2lt);
  T->SetBranchAddress("L.s2.rt",Ls2rt);

  Double_t bin,b1,b2;
  // bin=400;
  // b1=-0.25/1000000;
  // b2=0;
  bin=200;
  b1=-25;
  b2=15;
  TH1F *cointime = new TH1F("cointime","Coincidence time [ns] (Rs2time-Ls2time)",bin,b1,b2);

  //----------------Get Entry---------------------------

  Int_t nentries = T->GetEntries();
  Double_t data;
  cout<<"Total Number of Events = "<<nentries<<endl;
  for(Int_t i=0;i<nentries;i++){

    if(i%100000==0) cout << " events processed = " << i << endl;
    T->GetEntry(i);

    for(Int_t k=0;k<16;k++){
      for(Int_t j=0;j<16;j++){
	if(Rs2lt[k]>0 && Rs2rt[k]>0 && Ls2lt[j]>0 && Ls2rt[j]>0){
	  data = (Rs2time[k]-Ls2time[j])*1000000000+150.7;
	  cointime->Fill(data);
	}
      }
    }
  }

  cointime->Draw();
}
