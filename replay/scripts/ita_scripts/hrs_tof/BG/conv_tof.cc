/*
  conv.cc
  
Auther K. Itabashi
2018 Sep. 19th
Convert tritium_rootfile to TOF hist

*/

const double ch2ns = 0.05; // ch --> ns 

void conv_tof(){

  int run,evn;
  TChain *t1 = new TChain("T");
 for(Int_t i = 100670; i<100672;i++) /// these run goes 100200 to 100685
 {
 t1->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

 }

 TFile* fnew = new TFile("tritium_new.root","recreate"); 
 TTree* tnew = new TTree("T","Tritium test");






  // Right Arm Pion Rehector 1,2 //
  tnew->Branch("NRa1a",&NRa1a,"NRa1a");
  tnew->Branch("NRa2a",&NRa2a,"NRa2a");
  tnew->Branch("NRa1t",&NRa1t,"NRa1t");
  tnew->Branch("NRa2t",&NRa2t,"NRa2t");
  tnew->Branch("Rsha",Rsha,"Rsha/D"); // Right arm PionRehector1 ADC
  tnew->Branch("Rpsa",Rpsa,"Rpsa/D"); // Right arm Preshower ADC 
  // Right Arm S0 // 
  tnew->Branch("Rs0la",Rs0la,"Rs0la/D"); // Right arm S0-Top(B) ADC
  tnew->Branch("Rs0ra",Rs0ra,"Rs0ra/D"); // Right arm S0-Bottom(A) ADC
  tnew->Branch("Rs0lt",Rs0lt,"Rs0lt/D"); // Right arm S0-Top(B) TDC
  tnew->Branch("Rs0rt",Rs0rt,"Rs0rt/D"); // Right arm S0-Bottom(A) TDC
  tnew->Branch("Rs0la_p",Rs0la_p,"Rs01a_p/D"); // Right arm S0-Topo(B) ADC
  tnew->Branch("Rs0ra_p",Rs0ra_p,"Rs0ra_p/D"); // Right arm S2-Bottom(A) ADC
  // Right Arm S2 // 
  tnew->Branch("Rs2la",Rs2la,"Rs2la/D"); // Right arm S2-Top(B) ADC
  tnew->Branch("Rs2ra",Rs2ra,"Rs2ra/D"); // Right arm S2-Bottom(A) ADC
  tnew->Branch("Rs2lt",Rs2lt,"Rs2lt/D"); // Right arm S2-Top(B) TDC
  tnew->Branch("Rs2rt",Rs2rt,"Rs2rt/D"); // Right arm S2-Bottom(A) TDC
  tnew->Branch("Rs2la_p",Rs2la_p,"Rs21a_p/D"); // Right arm S2-Topo(B) ADC
  tnew->Branch("Rs2ra_p",Rs2ra_p,"Rs2ra_p/D"); // Right arm S2-Bottom(A) ADC


  // Left Arm Pion Rehector 1,2 //
  tnew->Branch("Lsha",Lsha,"Lsha/D"); // Left arm PionRehector1 ADC
  tnew->Branch("Lpsa",Lpsa,"Lpsa/D"); // Left arm Preshower ADC 
  // Left Arm S0 // 
  tnew->Branch("Ls0la",Ls0la,"Ls0la/D"); // Left arm S0-Top(B) ADC
  tnew->Branch("Ls0ra",Ls0ra,"Ls0ra/D"); // Left arm S0-Bottom(A) ADC
  tnew->Branch("Ls0lt",Ls0lt,"Ls0lt/D"); // Left arm S0-Top(B) TDC
  tnew->Branch("Ls0rt",Ls0rt,"Ls0rt/D"); // Left arm S0-Bottom(A) TDC
  tnew->Branch("Ls0la_p",Ls0la_p,"Ls01a_p/D"); // Left arm S0-Topo(B) ADC
  tnew->Branch("Ls0ra_p",Ls0ra_p,"Ls0ra_p/D"); // Leftt arm S2-Bottom(A) ADC
  // Left Arm S2 // 
  tnew->Branch("Ls2la",Ls2la,"Ls2la/D"); // Left arm S2-Top(B) ADC
  tnew->Branch("Ls2ra",Ls2ra,"Ls2ra/D"); // Left arm S2-Bottom(A) ADC
  tnew->Branch("Ls2lt",Ls2lt,"Ls2lt/D"); // Left arm S2-Top(B) TDChttps://github.com/itabasi/AC_HRS_at_ESB.git
  tnew->Branch("Ls2rt",Ls2rt,"Ls2rt/D"); // Left arm S2-Bottom(A) TDC
  tnew->Branch("Ls2la_p",Ls2la_p,"Ls21a_p/D"); // Left arm S2-Topo(B) ADC
  tnew->Branch("Ls2ra_p",Ls2ra_p,"Ls2ra_p/D"); // Left arm S2-Bottom(A) ADC
  // Right Arm A1 //                                                                   
  tnew->Branch("Ra1a",Ra1a,"Ra1a[NRala]/D"); // Right arm AC1 ADC                             
  tnew->Branch("Ra1t",Ra1t,"Ra1t[NRa1t]/D"); // Right arm AC1 TDC                             
  tnew->Branch("Ra1sum",&Ra1sum,"Ra1sum/D"); // Right arm AC1 FADC sum                 
  tnew->Branch("Ra1a_p",Ra1a_p,"Ra1a_p/D"); // Right arm AC1 FADC (Pedestal)PMT        
  tnew->Branch("Ra1a_c",Ra1a_c,"Ra1a_c/D"); // Right arm AC1 FADC (correction) PMT     
                             
  // Right Arm AC2 //                                                                  
  tnew->Branch("Ra2a",Ra2a,"Ra2a/D"); // Right arm AC2 ADC                             
  tnew->Branch("Ra2t",Ra2t,"Ra2t[NRa2t]/D"); // Right arm AC2 TDC                             
  tnew->Branch("Ra2sum",&Ra2sum,"Ra2sum/D"); // Right arm AC2 FADC sum                 
  tnew->Branch("Ra2a_p",Ra2a_p,"Ra2a_p/D"); // Right arm AC2 FADC (Pedestal)PMT        
  tnew->Branch("Ra2a_c",Ra2a_c,"Ra2a_c/D"); // Right arm AC2 FADC (correction) PMT     

 double ent = t1->GetEntries();
 cout<<"GetEntries"<<ent<<ent<<endl; 
   tnew->Write();
   fnew->Close();
 

}

