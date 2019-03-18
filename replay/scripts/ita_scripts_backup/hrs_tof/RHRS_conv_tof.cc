// 10/4/2018 Author Itabshi
// Recreate ROOTfile 





void RHRS_conv_tof(){
//-------- TTree data input ---------------//
      //TChain *T = new TChain("T");
   TChain chain("T");
 for(Int_t i = 93495; i<93504;i++) /// these run goes 93495 to 93504
 {
    chain.Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));
 }

 chain.Merge("chain.root");
 TFile *fin=new TFile("chain.root");
 TTree* T=(TTree*)fin->Get("T");  

  //======== Recreate TOF_twcorr.root =============//
 TFile* fnew=new TFile("newROOTfiles/TOF_twcorr.root","recreate");
 TTree* tree=new TTree("tree","TOF_twcorr.root");



 int max=10000; 
 int Ndata_Rs0la_c,Ndata_Rs0ra_c,Ndata_Rs0rt_c,Ndata_Rs0lt_c;//S0 Ndata
 int Ndata_Rs2la_c,Ndata_Rs2ra_c,Ndata_Rs2rt_c,Ndata_Rs2lt_c;//S2 Ndata
 double DR_evtypebits; // trigger condition  
 double Rs0la_c[max],Rs0ra_c[max],Rs0lt_c[max],Rs0rt_c[max]; //S0 TDC & ADC  
 double Rs2la_c[max],Rs2ra_c[max],Rs2lt_c[max],Rs2rt_c[max]; //S2 TDC & ADC 



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


  T->Fill();



  //======= Recreate TOF_twcorr.root Ndata Branch Definition ===========//

  //S0 Ndata //
  tree->Branch("Ndata_Rs0ra_c",&Ndata_Rs0ra_c,"Ndata_Rs0ra_c/I");
  tree->Branch("Ndata_Rs0la_c",&Ndata_Rs0la_c,"Ndata_Rs0la_c/I");
  tree->Branch("Ndata_Rs0rt_c",&Ndata_Rs0rt_c,"Ndata_Rs0rt_c/I");
  tree->Branch("Ndata_Rs0lt_c",&Ndata_Rs0lt_c,"Ndata_Rs0lt_c/I");
  //S0 Ndata //
  tree->Branch("Ndata_Rs2ra_c",&Ndata_Rs2ra_c,"Ndata_Rs2ra_c/I");
  tree->Branch("Ndata_Rs2la_c",&Ndata_Rs2la_c,"Ndata_Rs2la_c/I");
  tree->Branch("Ndata_Rs2rt_c",&Ndata_Rs2rt_c,"Ndata_Rs2rt_c/I");
  tree->Branch("Ndata_Rs2lt_c",&Ndata_Rs2lt_c,"Ndata_Rs2lt_c/I");

  //===== Fill Ndata TOF_twcorr.root =======//
 
 int T_ent;
  T_ent=T->GetEntries();
  cout<<"T Entries : "<<T_ent<<endl;
  for(int i=0;i<T_ent;i++){
    T->GetEntry(i);   
    tree->Fill();
  }



  //===== Recreate TOF_twcorr.root Branch Definition =========//
 //S0 Branch //
  tree->Branch("Rs0ra_c",Rs0ra_c,"Rs0ra_c[Ndata_Rs0ra_c]/D");
  tree->Branch("Rs0la_c",Rs0la_c,"Rs0la_c[Ndata_Rs0la_c]/D");
  tree->Branch("Rs0rt_c",Rs0rt_c,"Rs0rt_c[Ndata_Rs0rt_c]/D");
  tree->Branch("Rs0lt_c",Rs0lt_c,"Rs0lt_c[Ndata_Rs0lt_c]/D");
  //S0 Branch //
  tree->Branch("Rs2ra_c",Rs2ra_c,"Rs2ra_c[Ndata_Rs2ra_c]/D");
  tree->Branch("Rs2la_c",Rs2la_c,"Rs2la_c[Ndata_Rs2la_c]/D");
  tree->Branch("Rs2rt_c",Rs2rt_c,"Rs2rt_c[Ndata_Rs2rt_c]/D");
  tree->Branch("Rs2lt_c",Rs2lt_c,"Rs2lt_c[Ndata_Rs2lt_c]/D");
  // Trigger condition //
  tree->Branch("DR_evtypebits",&DR_evtypebits,"DR_evtypebits/D");

  //=========================================================//
  

  double Rtof[16];
  tree->Branch("Rtof",Rtof,"Rtof/D");
  // std::string cut_trig("DR_evtypebits==112");
 char* cut_trig="DR_evtypebits==112";
  for(int i=0;i<T_ent;i++){
    T->GetEntry(i);
     
   //TOF Fill Branch //
 for(int j=0;j<16;j++){   
        Rtof[j]=(Rs2lt_c[j]+Rs2rt_c[j])/2.0-(Rs0lt_c[0]+Rs0rt_c[0])/2.0;	

 }

  

 	tree->Fill(); 
    

   tree->Write();
   fnew->Close();

  }

}
