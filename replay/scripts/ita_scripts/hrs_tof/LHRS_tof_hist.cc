//Author K. Itabashi
//HRS-R TOF Histgram macro

double rrange_para(int i,int j);
double lrange_para(int i,int j);
  int chamx=17;
double tdcTime=5.0e-11;

void LHRS_tof_hist(){

//-------- TTree data input ---------------//
  int nrun; //run number
  char* arm;
  cout<<"=======< HRS TOF Hist macro >========"<<ednl;
  cout<<"HRS arm : ">>arm<<endl;
  cout<<"Run number : ">>nrun<<endl; 
  // TFile* f1 = new TFile(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d.root",run));

 TFile* f1 = new TFile(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",nrun));
  //      TFile* f1 = new TFile("tritium_100028.root");
  TTree* t1 = (TTree*)f1->Get("T");

  int ent=T->GetEntries();
  cout<<"Get Entries : "<<ent;
 int max=10000; 
 //================== Right arm ==========================================//
 int Ndata_Rs0la_c,Ndata_Rs0ra_c,Ndata_Rs0rt_c,Ndata_Rs0lt_c;//S0 Ndata
 int Ndata_Rs2la_c,Ndata_Rs2ra_c,Ndata_Rs2rt_c,Ndata_Rs2lt_c;//S2 Ndata
 double DR_evtypebits; // trigger condition  
 double Rs0la_c[max],Rs0ra_c[max],Rs0lt_c[max],Rs0rt_c[max]; //S0 TDC & ADC  
 double Rs2la_c[max],Rs2ra_c[max],Rs2lt_c[max],Rs2rt_c[max]; //S2 TDC & ADC 
 //================== Left arm ===========================================//
 int Ndata_Ls0la_c,Ndata_Ls0ra_c,Ndata_Ls0rt_c,Ndata_Ls0lt_c;//S0 Ndata
 int Ndata_Ls2la_c,Ndata_Ls2ra_c,Ndata_Ls2rt_c,Ndata_Ls2lt_c;//S2 Ndata
 //double DR_evtypebits; // trigger condition  
 double Ls0la_c[max],Ls0ra_c[max],Ls0lt_c[max],Ls0rt_c[max]; //S0 TDC & ADC  
 double Ls2la_c[max],ls2ra_c[max],Ls2lt_c[max],Ls2rt_c[max]; //S2 TDC & ADC 



  if(arm="right"||arm="Right" || arm="R" || arm="r"){

 //============= Set Branch Status ==================//
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("Ndata.R.s0.la_c",1);
  T->SetBranchStatus("Ndata.R.s0.ra_c",1);
  T->SetBranchStatus("Ndata.R.s0.lt_c",1);
  T->SetBranchStatus("Ndata.R.s0.rt_c",1);
  T->SetBranchStatus("Ndata.R.s2.la_c",1);
  T->SetBranchStatus("Ndata.R.s2.ra_c",1);
  T->SetBranchStatus("Ndata.R.s2.lt_c",1);
  T->SetBranchStatus("R.s0.la_c",1);
  T->SetBranchStatus("R.s0.ra_c",1);
  T->SetBranchStatus("R.s0.lt_c",1);
  T->SetBranchStatus("R.s0.rt_c",1);
  T->SetBranchStatus("R.s2.la_c",1);
  T->SetBranchStatus("R.s2.ra_c",1);
  T->SetBranchStatus("R.s2.lt_c",1);
  T->SetBranchStatus("R.s2.rt_c",1);
  T->SetBranchStatus("DR.evtypebits",1);   
  T->SetBranchStatus(nHRS_trig,1);
  T->SetBranchStatus(nRs2_pads,1);
  T->SetBranchStatus(nRs2_nthit,1);
  T->SetBranchStatus(nRs2_tdchit,1);
  T->SetBranchStatus(nR_trx,1);
  T->SetBranchStatus(nR_trth,1);
  T->SetBranchStatus(nRs2_lac,1);
  T->SetBranchStatus(nR_trbeta,1);

  // Right Arm S0 Ndata// 
  T->SetBranchAddress("Ndata.R.s0.la_c",&Ndata_Rs0la_c); // Right arm S0-Top(B) ADC
  T->SetBranchAddress("Ndata.R.s0.ra_c",&Ndata_Rs0ra_c); // Right arm S0-Bottom(A) ADC
  T->SetBranchAddress("Ndata.R.s2.la_c",&Ndata_Rs2la_c); // Right arm S2-Top(B) ADC
  T->SetBranchAddress("Ndata.R.s2.ra_c",&Ndata_Rs2ra_c); // Right arm S2-Bottom(A) ADC
  // Right Arm S2 Ndata//
  T->SetBranchAddress("Ndata.R.s0.lt_c",&Ndata_Rs0lt_c); // Right arm S0-Top(B) TDC
  T->SetBranchAddress("Ndata.R.s0.rt_c",&Ndata_Rs0rt_c); // Right arm S0-Bottom(A) TDC
  T->SetBranchAddress("Ndata.R.s2.lt_c",&Ndata_Rs2lt_c); // Right arm S2-Top(B) TDC
  T->SetBranchAddress("Ndata.R.s2.rt_c",&Ndata_Rs2rt_c); // Right arm S2-Bottom(A) TDC
  // Right Arm S0 // 
  T->SetBranchAddress("R.s0.la_c",Rs0la_c); // Right arm S0-Top(B) ADC
  T->SetBranchAddress("R.s0.ra_c",Rs0ra_c); // Right arm S0-Bottom(A) ADC
  T->SetBranchAddress("R.s0.lt_c",Rs0lt_c); // Right arm S0-Top(B) TDC
  T->SetBranchAddress("R.s0.rt_c",Rs0rt_c); // Right arm S0-Bottom(A) TDC
  // Right Arm S2 // 
  T->SetBranchAddress("R.s2.la_c",Rs2la_c); // Right arm S2-Top(B) ADC
  T->SetBranchAddress("R.s2.ra_c",Rs2ra_c); // Right arm S2-Bottom(A) ADC
  T->SetBranchAddress("R.s2.lt_c",Rs2lt_c); // Right arm S2-Top(B) TDC
  T->SetBranchAddress("R.s2.rt_c",Rs2rt_c); // Right arm S2-Bottom(A) TDC
  // Trigger condtion //
  T->SetBranchAddress("DR.evtypebits",&DR_evtypebits); 

  //=== Coincidence F1TDC Branch ================//
	  T->SetBranchAddress(nHRS_trig, &HRS_trig);
          T->SetBranchAddress(nRs2_pads, &Rs2_pads);
          T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
          T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
          T->SetBranchAddress(nR_trx, &R_trx);
          T->SetBranchAddress(nR_trth, &R_trth);
          T->SetBranchAddress(nRs2_lac, &Rs2_lac);
          T->SetBranchAddress(nR_trbeta, &R_trbeta);
  //============================================//

  }


  if(arm="left"||arm="Left" || arm="L" || arm="l"){

 //============= Set Branch Status ==================//
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("Ndata.L.s0.la_c",1);
  T->SetBranchStatus("Ndata.L.s0.ra_c",1);
  T->SetBranchStatus("Ndata.L.s0.lt_c",1);
  T->SetBranchStatus("Ndata.L.s0.rt_c",1);
  T->SetBranchStatus("Ndata.L.s2.la_c",1);
  T->SetBranchStatus("Ndata.L.s2.ra_c",1);
  T->SetBranchStatus("Ndata.L.s2.lt_c",1); 
  T->SetBranchStatus("L.s0.la_c",1);
  T->SetBranchStatus("L.s0.ra_c",1);
  T->SetBranchStatus("L.s0.lt_c",1);
  T->SetBranchStatus("L.s0.rt_c",1);
  T->SetBranchStatus("L.s2.la_c",1);
  T->SetBranchStatus("L.s2.ra_c",1);
  T->SetBranchStatus("L.s2.lt_c",1);
  T->SetBranchStatus("L.s2.rt_c",1);
  T->SetBranchStatus("DR.evtypebits",1);   
  T->SetBranchStatus(nHRS_trig,1);
  T->SetBranchStatus(nLs2_pads,1);
  T->SetBranchStatus(nLs2_nthit,1);
  T->SetBranchStatus(nLs2_tdchit,1);
  T->SetBranchStatus(nL_trx,1);
  T->SetBranchStatus(nL_trth,1);
  T->SetBranchStatus(nLs2_lac,1);
  T->SetBranchStatus(nL_trbeta,1);

  // Right Arm S0 Ndata// 
  T->SetBranchAddress("Ndata.L.s0.la_c",&Ndata_Ls0la_c); // Right arm S0-Top(B) ADC
  T->SetBranchAddress("Ndata.L.s0.ra_c",&Ndata_Ls0ra_c); // Right arm S0-Bottom(A) ADC
  T->SetBranchAddress("Ndata.L.s2.la_c",&Ndata_Ls2la_c); // Right arm S2-Top(B) ADC
  T->SetBranchAddress("Ndata.L.s2.ra_c",&Ndata_Ls2ra_c); // Right arm S2-Bottom(A) ADC
  // Right Arm S2 Ndata//
  T->SetBranchAddress("Ndata.L.s0.lt_c",&Ndata_Ls0lt_c); // Right arm S0-Top(B) TDC
  T->SetBranchAddress("Ndata.L.s0.rt_c",&Ndata_Ls0rt_c); // Right arm S0-Bottom(A) TDC
  T->SetBranchAddress("Ndata.L.s2.lt_c",&Ndata_Ls2lt_c); // Right arm S2-Top(B) TDC
  T->SetBranchAddress("Ndata.L.s2.rt_c",&Ndata_Ls2rt_c); // Right arm S2-Bottom(A) TDC
  // Right Arm S0 // 
  T->SetBranchAddress("L.s0.la_c",Ls0la_c); // Right arm S0-Top(B) ADC
  T->SetBranchAddress("L.s0.ra_c",Ls0ra_c); // Right arm S0-Bottom(A) ADC
  T->SetBranchAddress("L.s0.lt_c",Ls0lt_c); // Right arm S0-Top(B) TDC
  T->SetBranchAddress("L.s0.rt_c",Ls0rt_c); // Right arm S0-Bottom(A) TDC
  // Right Arm S2 // 
  T->SetBranchAddress("L.s2.la_c",Ls2la_c); // Right arm S2-Top(B) ADC
  T->SetBranchAddress("L.s2.ra_c",Ls2ra_c); // Right arm S2-Bottom(A) ADC
  T->SetBranchAddress("L.s2.lt_c",Ls2lt_c); // Right arm S2-Top(B) TDC
  T->SetBranchAddress("L.s2.rt_c",Ls2rt_c); // Right arm S2-Bottom(A) TDC
  // Trigger condtion //
  T->SetBranchAddress("DR.evtypebits",&DR_evtypebits); 

  //=== Coincidence F1TDC Branch ================//
	  T->SetBranchAddress(nHRS_trig, &HRS_trig);
          T->SetBranchAddress(nRs2_pads, &Rs2_pads);
          T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
          T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
          T->SetBranchAddress(nR_trx, &R_trx);
          T->SetBranchAddress(nR_trth, &R_trth);
          T->SetBranchAddress(nRs2_lac, &Rs2_lac);
          T->SetBranchAddress(nR_trbeta, &R_trbeta);
  //============================================//

  }

  //============ Definition Hist =====================//

  double tof_tdc_min,tof_tdc_max;
  double bin_tof_tdc[chmax];
  double min_adc,max_adc;
  double tof[chmax];
  min_adc=100.; max_adc=2000.;
  int bin_adc=max-adc-min_adc;
  tof_tdc_min=-1.0e-6;
  tof_tdc_max=1.0e-6;
  bin_tof_tdc=(tof_tdc_max-tof_tdc_min)/tdcTime;
  bin_tof_tdc=(int)bin_tof_tdc;
  TH1F* htof_tdc[chmax];
  TH2F* htof_adc_t[chmax][2];

  for(int i=0;i<chmax;i++){
    htof_tdc[i]=new TH1F(Form("htof_tdc[%d]",i),Form("LHRS TOF with S2[%d]-S0 tdc hist",i),bin_tof_tdc,tof_tdc_min,tof_tdc_max);
   htof_adc_t[i][0]=new TH2F(Form("htof_adc_t[%d][0]",i),Form("LHRS TOF vs S2[%d]-R-PMT ADC Hist ",i),bin_tof_tdc,tof_tdc_min,tof_tdc_max,bin_adc,min_adc,max_adc);
   htof_adc_t[i][1]=new TH2F(Form("htof_adc_t[%d][1]",i),Form("LHRS TOF vs S2[%d]-L-PMT ADC Hist ",i),bin_tof_tdc,tof_tdc_min,tof_tdc_max,bin_adc,min_adc,max_adc);
  
  for(int k=0;k<ent;k++){
    T->GetEntry(k);


    tof[i]=(Ls2rt_c[i]+Ls2lt_c[i])/2.0-(Ls0rt_c[0]+Ls0lt_c[0])/2.0;
    htof_tdc[i]->Fill(tof[i]);
    htof_adc_t[i][0]->Fill(tof[i],Ls2ra_c[i]);
    htof_adc_t[i][1]->Fill(tof[i],Ls2la_c[i]);
}


  }



}




double rrange_para(int i,int j){
  int npara=6;// number of parameters
  double par[chmax][npara];
  double param;
 
  //=== Inital parameters========//
  for(int k=0;k<chmax;k++){
   

   par[k][0]=-4e-9, par[k][1]=4e-9, par[k][2]=-0.1e-7, par[k][3]=0.1e-7;//TOF


   if(k==16){//S0
  par[k][4]=500.0, par[k][5]=5000.;//ADC
    }else{ //S2
      par[k][4]=100.0, par[k][5]=500.;//ADC
    }
   }

  //===== Set Parameters ========//
  //  par[8][0]=0.18e-6, par[8][1]=0.22e-6, par[8][2]=-0.18e-6,par[8][3]=0.18e-6;
  //  par[16][0]=0.180e-6, par[16][1]=0.23e-6, par[16][2]=-0.18e-6, par[16][3]=0.25e-6;

  //============================//
  for(int k=0;k<chmax;k++){
   for(int l=0;l<npara;l++){
    if(k==i && l==j){
      param=par[k][l];
     
    }

    }
  }
  return param;
}


double lrange_para(int i,int j){
  int npara=6;// number of parameters
  double par[chmax][npara];
  double param;
 
  //=== Inital parameters========//
  for(int k=0;k<chmax;k++){
   

   par[k][0]=-4e-9, par[k][1]=4e-9, par[k][2]=-0.1e-7, par[k][3]=0.1e-7;//TOF


   if(k==16){//S0
  par[k][4]=500.0, par[k][5]=5000.;//ADC
    }else{ //S2
      par[k][4]=100.0, par[k][5]=500.;//ADC
    }
   }

  //===== Set Parameters ========//
  //  par[8][0]=0.18e-6, par[8][1]=0.22e-6, par[8][2]=-0.18e-6,par[8][3]=0.18e-6;
  //  par[16][0]=0.180e-6, par[16][1]=0.23e-6, par[16][2]=-0.18e-6, par[16][3]=0.25e-6;

  //============================//
  for(int k=0;k<chmax;k++){
   for(int l=0;l<npara;l++){
    if(k==i && l==j){
      param=par[k][l];
     
    }

    }
  }
  return param;
}

