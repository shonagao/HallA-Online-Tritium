//Author K. Itabashi
//HRS-R TOF Histgram macro
//void hrs_tof_hist(){
void hrs_tof_hist(){
  // gROOT->SetBatch();


//-------- TTree data input ---------------//
  TChain*  T=new TChain("T");
  int num;
  int run=94003;
  cout<<"==========< TOF Hist >========"<<endl;
  cout<<" run number : "<<run<<endl;
  //  cin>>run<<endl;
  //  T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d*.root",i));//a-online

  T->Add(Form("/u/home/itabashi/nnL/Tohoku_github/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d*.root",run));// ifarm



  /*
 for(Int_t i = 93500; i<93504;i++) /// these run goes 93495 to 93504
 {
   // T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));//single arm
  T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d*.root",i));//For short range run

}
  */



int max=10000; 
 int Ndata_Rs0la_c,Ndata_Rs0ra_c,Ndata_Rs0rt_c,Ndata_Rs0lt_c;//S0 Ndata
 int Ndata_Rs2la_c,Ndata_Rs2ra_c,Ndata_Rs2rt_c,Ndata_Rs2lt_c;//S2 Ndata
 double DR_evtypebits; // trigger condition  
 double Rs0la_c[max],Rs0ra_c[max],Rs0lt_c[max],Rs0rt_c[max]; //S0 TDC & ADC  
 double Rs2la_c[max],Rs2ra_c[max],Rs2lt_c[max],Rs2rt_c[max]; //S2 TDC & ADC 
 double Rpathl[max]; 
  // Variable to be used event by event 
 	 
   // Define Tree leaf variable to hold the values
          TString nHRS_trig = "DR.rrawt4";
	  //TString nHRS_trig = "DR.evtypebits"; // for beta =1 particle
          TString nRs2_pads = "R.s2.t_pads";
          TString nRs2_nthit = "R.s2.nthit";
	  TString nRs2_tdchit = "RTDC.F1FirstHit"; // for coincidence run
	  // TString nRs2_tdchit = "F1FirstHit"; // for marathan run
	  TString nR_trx = "R.tr.x";
          TString nR_trth = "R.tr.th";
          TString nRs2_lac = "R.s2.la_c";
          TString nR_trbeta="R.tr.beta";

    // Define Tree leaf variable
	  Double_t HRS_trig;
	  Double_t Rs2_pads[100];
	  Double_t Rs2_nthit;
	  Double_t Rs2_tdchit[100];
	  Double_t R_trx[100];
	  Double_t R_trth[100];
	  Double_t Rs2_lac[100];
          Double_t R_trbeta;

 
 //============= Set Branch Status ==================//
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("Ndata.R.s0.la_c",1);
  T->SetBranchStatus("Ndata.R.s0.ra_c",1);
  T->SetBranchStatus("Ndata.R.s0.lt_c",1);
  T->SetBranchStatus("Ndata.R.s0.rt_c",1);
  T->SetBranchStatus("Ndata.R.s2.la_c",1);
  T->SetBranchStatus("Ndata.R.s2.ra_c",1);
  T->SetBranchStatus("Ndata.R.s2.lt_c",1);
  T->SetBranchStatus("R.s2.rt_c",1); 
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
  T->SetBranchStatus("R.tr.pathl",1);


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
  // Path Length // 
 T->SetBranchAddress("R.tr.pathl",Rpathl); 
    // Set The Branch addres
	  T->SetBranchAddress(nHRS_trig, &HRS_trig);
          T->SetBranchAddress(nRs2_pads, &Rs2_pads);
          T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
          T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
          T->SetBranchAddress(nR_trx, &R_trx);
          T->SetBranchAddress(nR_trth, &R_trth);
          T->SetBranchAddress(nRs2_lac, &Rs2_lac);
          T->SetBranchAddress(nR_trbeta, &R_trbeta);



	  // original values
  Double_t corr_x[14] = { 4.91939e-09, 5.41034e-09, 6.87688e-09, 9.22121e-09, 0.795016e-8, 0.802636e-08, 0.787479e-08, 7.54862e-09, 7.5127e-9, 8.48865e-9, 5.42156e-9, 6.27864e-9, 5.78027e-9, 8.605971e-10};
	 
  Double_t corr_th[14] = {-4.46911e-08, -4.79507e-08, -3.59540e-08, -3.04303e-08, -3.25733e-08, -3.10881e-08, -3.18107e-08, -3.5057e-08,  -3.59703e-08, -3.76206e-08, -3.68166e-08,  -3.51979e-08, -3.54868e-08, -1.56998e-08};
 	 
  Double_t corr_adc[14] = {-8.66369e-13, -3.84824e-13, -1.45016e-12, -12.08217e-13,-6.95305e-13, -5.37148e-13, -5.95287e-13, -1.01789e-12, -1.02612e-12, -1.87664e-12, -3.19282e-12, -3.37812e-12, -7.80362e-13,-7.65658e-13};

  Double_t alignment[14] = {-1.915e-9, -1.917e-9, 0.85e-9, 1.90e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};

  //====== Definition ===========//
  int chmax=17;//s2 channel 
  double tdcTime = 56.23e-12; // this is the conversion per chanel ie 56.23 picosec
   TString name;
  TH1F* htof[chmax];
  TH2F* htof_adc_r[chmax];
  TH2F* htof_adc_l[chmax];
  TH2F* htof_pathl[chmax];
  TH2F* htof_x[chmax];
  TH2F* htof_th[chmax];
  // TCanvas* ctof[chmax];
  //TCanvas* ctof_adc[chmax];
  // TCanvas* ctof_adc_r[chmax];
  // TCanvas* ctof_pathl[chmax];
  double min_tof,max_tof,min_adc,max_adc,min_pathl,max_pathl;
  double bin_tof,bin_adc,bin_pathl;
  min_tof=-0.65e-6; max_tof=0.65e-6;
   // min_tof=-0.65; max_tof=0.65;
  bin_tof=(max_tof-min_tof)/tdcTime;
  bin_tof=(int)bin_tof;
  min_adc=100.;
  max_adc=1000.;
  bin_adc=max_adc-min_adc;
  bin_adc=(int)bin_adc;
  min_pathl=25.;
  max_pathl=29.;
  bin_pathl=(max_pathl-min_pathl)*100;
  bin_pathl=(int)bin_pathl;
  cout<<"bin_tof :"<<bin_tof<<endl;
  cout<<"bin_adc :"<<bin_adc<<endl;
  cout<<"bin_pathl :"<<bin_pathl<<endl;



  for(int i=0;i<chmax;i++){
    //=========== Definition Hist =========================//
    htof[i]=new TH1F(Form("htof[%d]",i),Form("F1TDC TOF Histgram RF-S2ch%d",i),bin_tof,min_tof,max_tof);
    htof_adc_r[i]=new TH2F(Form("htof_adc_r[%d]",i),Form("F1TDC TOF vs S2ch%d R-PMT ADC Hist",i),bin_tof,min_tof,max_tof,bin_adc,min_adc,max_adc);
 htof_adc_l[i]=new TH2F(Form("htof_adc_l[%d]",i),Form("F1TDC TOF vs S2ch%d L-PMT ADC Hist",i),bin_tof,min_tof,max_tof,bin_adc,min_adc,max_adc);
    htof_pathl[i]=new TH2F(Form("htof_pathl[%d]",i),Form("F1TDC TOF (RF-S2ch%d) vs Path Length  Hist",i),bin_tof,min_tof,max_tof,bin_pathl,min_pathl,max_pathl);

  }

 
   //======= Definition Canvas ==========================//
    TCanvas *ctof=new TCanvas("ctof","ctof");
    //    TCanvas* ctof_adc_r=new TCanvas("ctof_adc_r","ctof_adc_r");
    //  TCanvas* ctof_adc_l=new TCanvas("ctof_adc_l","ctof_adc_l");
    //  TCanvas* ctof_pathl=new TCanvas("ctof_pathl","ctof_pathl");
    // ctof->Divide(2,8);
     ctof->Divide(2,2);
    //   ctof_adc_r->Divide(2,8);
    //  ctof_adc_l->Divide(2,8);
    // ctof_pathl->Divide(2,8);

  double RF_s2_mean[chmax];
  int nentries = T->GetEntries();
  cout<<"Get Event : "<<nentries<<endl;

TCanvas* c1= new TCanvas("c1","c1");
   cout<<"make canvas done "<<endl;
 
  //====== Fill Hist ============//


for( int i=0;i<4;i++){
for(int j=0;j<nentries;j++){
	      T->GetEntry(j); 	 
	      
   if( HRS_trig>1 && Rs2_nthit==1 && Rs2_pads[0]==i &&  R_trbeta>0.76 && R_trbeta<1.01  ) // R.tr.beta>0.76 && R.tr.beta<80
	{
	  RF_s2_mean[i] = (Rs2_tdchit[8] - Rs2_tdchit[1]) * tdcTime  - ((Rs2_tdchit[8]-Rs2_tdchit[i+16])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[i+48])*tdcTime)/2.0;	  
           htof[i]->Fill(RF_s2_mean[i]);
	}

	     

 } // end entry loop ie irun=90854 loop

   
 }

 cout<<"Fill done "<<endl;

//ctof->cd(i);
 c1->cd(1);
 cout<<"move canvas done "<<endl;
  htof[7]->Draw();
  
  /* 
for( int i=0;i<16;i++){
for(int j=0;j<nentries;j++){
	      T->GetEntry(j); 	 
	      
   if( HRS_trig>1 && Rs2_nthit==1 && Rs2_pads[0]==i &&  R_trbeta>0.76 && R_trbeta<1.01  ) // R.tr.beta>0.76 && R.tr.beta<80
	{
	  //Calc RF_S2mean depending on i
	  //  RF_s2_mean[i] = (Rs2_tdchit[8] - Rs2_tdchit[1]) * tdcTime  - (((Rs2_tdchit[8]-Rs2_tdchit[i+16])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[i+48])*tdcTime)/2.0 + corr_x[i-2]*R_trx[0] + corr_th[i-2]*R_trth[0] + corr_adc[i-2]*Rs2_lac[i] + alignment[i-2]);
	  RF_s2_mean[i] = (Rs2_tdchit[8] - Rs2_tdchit[1]) * tdcTime  - ((Rs2_tdchit[8]-Rs2_tdchit[i+16])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[i+48])*tdcTime)/2.0;
	  
           htof[i]->Fill(RF_s2_mean[i]);
	   htof_adc_r[i]->Fill(RF_s2_mean[i],Rs2rt_c[i]);
	   htof_adc_l[i]->Fill(RF_s2_mean[i],Rs2lt_c[i]);
	   htof_pathl[i]->Fill(RF_s2_mean[i],Rpathl[0]);
	}
	     

 } // end entry loop ie irun=90854 loop

 ctof->cd(1+i);
 htof[i]->Draw();
 ctof_adc_r->cd(1+i);
 htof_adc_r[i]->Draw();
 ctof_adc_l->cd(1+i);
 htof_adc_l[i]->Draw();
 ctof_pathl->cd(1+i);
 htof_pathl[i]->Draw();
 
 }


 // name.Form("./pdf/tof_%d.pdf",run);
 name.Form("./pdf/tof_94003.pdf");
 ctof->Print(name+"[","pdf");
 ctof->Print(name,"pdf");
 ctof_adc_r->Print(name,"pdf");
 ctof_adc_l->Print(name,"pdf");
 ctof_pathl->Print(name,"pdf");
 ctof_pathl->Print(name +"]","pdf");
cout<<"========================"<<endl;
cout<<"Have Finished crating PDF File !!"<<endl;
  */ 


}
