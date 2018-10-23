/*
  conv.cc
  
Auther K. Itabashi
2018 Sep. 19th
Convert tritium_rootfile to AC hist

*/

const double ch2ns = 0.05; // ch --> ns 

void conv(){

  int run,evn;

  cout<<"run number (defolt(-1):93267) "<<endl;
  cin>>run;
  cout<<"Event number (Max: -1)"<<endl;
  cin>>evn;
  if(run<0){
    run=93276;
}
 

  double E,M,me,mm;
  me=0.511; //electron mass [MeV/c^2]
  M=938.;// target mass [MeV/c^2]
  E=4000.; //Beam energy [MeV]

  TFile* f1 = new TFile(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d.root",run));
  //      TFile* f1 = new TFile("tritium_100028.root");
  TTree* t1 = (TTree*)f1->Get("T");
   TFile* fnew = new TFile(Form("/adaqfs/home/a-onl/tritium_work/itabashi/ita_macro/rootfiles/tritium_new%d.root",run),"recreate");


  //  TFile* fnew = new TFile("rootfiles/tritium_new.root","recreate"); 
 TTree* tnew = new TTree("tree","Tritium test");


 double Rsha[10000],Rpsa[10000],Rs0ra[10000],Rs0la[10000],Rs0lt[10000],Rs0rt[10000],Rs0la_p[10000],Rs0ra_p[10000],Rs2ra[10000],Rs2la[10000],Rs2lt[10000],Rs2rt[10000],Rs2la_p[10000],Rs2ra_p[10000];
 double Lsha[10000],Lpsa[10000],Ls0ra[10000],Ls0la[10000],Ls0lt[10000],Ls0rt[10000],Ls0la_p[10000],Ls0ra_p[10000],Ls2ra[10000],Ls2la[10000],Ls2lt[10000],Ls2rt[10000],Ls2la_p[10000],Ls2ra_p[10000];
 double Ra2t[10000],Ra2a[10000],Ra2sum,Ra2a_p[10000],Ra2a_c[10000];
 double Ra1t[10000],Ra1a[10000],Ra1sum,Ra1a_p[10000],Ra1a_c[10000],R_tr_beta,L_tr_beta,EKRNU,EKLNU;
 int NRa1a,NRa2a,NRa1t,NRa2t,EKLQ2,EKRQ2;
 //NRa1a=24;
 //NRa2a=26;
 // (AC1)Aerogel Chrenkov Right ARM ADC //                                            
  t1->SetBranchAddress("Ndata.R.a1.a",&NRa1a);
 t1->SetBranchAddress("Ndata.R.a2.a",&NRa2a);
 t1->SetBranchAddress("Ndata.R.a1.t",&NRa1t);
 t1->SetBranchAddress("Ndata.R.a2.t",&NRa2t);


 t1->SetBranchAddress("R.a1.t",Ra1t);
 t1->SetBranchAddress("R.a1.a",Ra1a);
 t1->SetBranchAddress("R.a1.asum_c",&Ra1sum);
 t1->SetBranchAddress("R.a1.a_p",Ra1a_p);
 t1->SetBranchAddress("R.a1.a_c",Ra1a_c);

 // (AC2)Aerogel Chrenkov Right ARM ADC //                                            
 t1->SetBranchAddress("R.a2.t",Ra2t);
 t1->SetBranchAddress("R.a2.a",Ra2a);
 t1->SetBranchAddress("R.a2.asum_c",&Ra2sum);
 t1->SetBranchAddress("R.a2.a_p",Ra2a_p);
 t1->SetBranchAddress("R.a2.a_c",Ra2a_c);


  // Cherenkov Right arm ADC //
  t1->SetBranchAddress("R.sh.a",Rsha);
  t1->SetBranchAddress("R.ps.a",Rpsa); 
 // S0 Right arm ADC//
  t1->SetBranchAddress("R.s0.la",Rs0la);
  t1->SetBranchAddress("R.s0.ra",Rs0ra);
  t1->SetBranchAddress("R.s0.lt",Rs0lt);
  t1->SetBranchAddress("R.s0.rt",Rs0rt);
  t1->SetBranchAddress("R.s0.la_p",Rs0ra_p);
  t1->SetBranchAddress("R.s0.ra_p",Rs0la_p);
  // S2 Right arm ADC//
  t1->SetBranchAddress("R.s2.la",Rs2la);
  t1->SetBranchAddress("R.s2.ra",Rs2ra);
  t1->SetBranchAddress("R.s2.lt",Rs2lt);
  t1->SetBranchAddress("R.s2.rt",Rs2rt);
  t1->SetBranchAddress("R.s2.la_p",Rs2ra_p);
  t1->SetBranchAddress("R.s2.ra_p",Rs2la_p);
 
  // Cherenkov Left arm ADC //
  t1->SetBranchAddress("L.sh.a",Lsha);
  t1->SetBranchAddress("L.ps.a",Lpsa);
   
 // S0 Left arm ADC//
  t1->SetBranchAddress("L.s0.la",Ls0la);
  t1->SetBranchAddress("L.s0.ra",Ls0ra);
  t1->SetBranchAddress("L.s0.lt",Ls0lt);
  t1->SetBranchAddress("L.s0.rt",Ls0rt);
  t1->SetBranchAddress("L.s0.la_p",Ls0ra_p);
  t1->SetBranchAddress("L.s0.ra_p",Ls0la_p);
  // S2 Left arm ADC//
  t1->SetBranchAddress("L.s2.la",Ls2la);
  t1->SetBranchAddress("L.s2.ra",Ls2ra);
  t1->SetBranchAddress("L.s2.lt",Ls2lt);
  t1->SetBranchAddress("L.s2.rt",Ls2rt);
  t1->SetBranchAddress("L.s2.la_p",Ls2ra_p);
  t1->SetBranchAddress("L.s2.ra_p",Ls2la_p);

  // beta information //
  t1->SetBranchAddress("R.tr.beta",&R_tr_beta);
  t1->SetBranchAddress("L.tr.beta",&L_tr_beta);
 
 //Energy inforamtion //
  t1->SetBranchAddress("EKL.omega",&EKLNU);                                            
  t1->SetBranchAddress("EKR.omega",&EKRNU);                             

  
  ///======= particle information =================//
  // beta //
  tnew->Branch("L_tr_beta",&L_tr_beta,"L_tr_beta/D"); // Left arm beta
  tnew->Branch("R_tr_beta",&R_tr_beta,"R_tr_beta/D"); // Right arm beta 
  // momentum //
  tnew->Branch("EKLQ2",&EKLQ2,"EKLQ2/D"); // Left arm beta
  tnew->Branch("EKRQ2",&EKRQ2,"EKRQ2/D"); // Right arm beta 
  //  Energy //
  tnew->Branch("EKLNU",&EKLNU,"EKLNU/D"); // (Ee-Ee') energy transfer
  tnew->Branch("EKRNU",&EKRNU,"EKRNU/D"); // (Q-Ep) energy transfer 

  

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


  // missing mass //
  tnew->Branch("mm",&mm,"mm/D");
 
 double ent = t1->GetEntries();
 if(ent>evn && evn>0){
   ent=evn;
 }
  cout<<"event counts "<<ent<<endl;  

  for(int i=0 ; i<ent ; i++){
       t1->GetEntry(i);
       }
       mm=sqrt(pow(E+M-(E-EKLNU)-EKRNU,2)-pow(sqrt(pow(E,2)-pow(me,2))-L_tr_beta*(E-EKLNU)-R_tr_beta*EKRNU,2));
   tnew->Fill();


 
   tnew->Write();
   fnew->Close();
 
}

