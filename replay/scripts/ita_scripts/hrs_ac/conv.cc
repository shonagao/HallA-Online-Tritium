/*
  conv.cc
  
Auther K. Itabashi
2018 Sep. 19th
Convert tritium_rootfile to AC hist

*/

const double ch2ns = 0.05; // ch --> ns 

void conv(){

  int run,evn,Nbeta;
  Nbeta=0;
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
 double Ra2t[10000],Ra2a[10000],Ra2sum,Ra2a_p[10000],Ra2a_c[10000],Ra2npe[10000];
 double Ra1t[10000],Ra1a[10000],Ra1sum,Ra1a_p[10000],Ra1a_c[10000],Ra1npe[10000],R_tr_beta[100],L_tr_beta[100],EKRNU,EKLNU,EKLQ2,EKRQ2;

 //==== AC 1PE position ===========//
//====== Pedestal and 1PE position ===============//
double ped1[24],gain1[24],ped2[26],gain2[26];

gain1[0]=379.;
gain1[1]=345.;
gain1[2]=341.;
gain1[3]=341.;
gain1[4]=364.;
gain1[5]=459.;
gain1[6]=375.;
gain1[7]=365.;
gain1[8]=370.;
gain1[9]=370.;
gain1[10]=362.;
gain1[11]=357.;
gain1[12]=320.;
gain1[13]=455.;
gain1[14]=341.;
gain1[15]=312.;
gain1[16]=340.;
gain1[17]=350.;
gain1[18]=356.;
gain1[19]=330.;
gain1[20]=347.;
gain1[21]=400.;
gain1[22]=372.;
gain1[23]=300.;

gain2[0]=397.;
gain2[1]=413.;
gain2[2]=378.;
gain2[3]=400.;
gain2[4]=389.;
gain2[5]=410.;
gain2[6]=380.;
gain2[7]=378.;
gain2[8]=410.;
gain2[9]=437.;
gain2[10]=433.;
gain2[11]=480.;
gain2[12]=402.;
gain2[13]=390.;
gain2[14]=378.;
gain2[15]=376.;
gain2[16]=375.;
gain2[17]=400.;
gain2[18]=409.;
gain2[19]=405.;
gain2[20]=388.;
gain2[21]=404.;
gain2[22]=391.;
gain2[23]=437.;
gain2[24]=412.;
gain2[25]=368.;


 int NRa1a,NRa2a,NRa1t,NRa2t,NRa1ac,NRa2ac; // A1 and A2 Ndata
 int NRs0lt,NRs0rt,NRs0la,NRs0ra;//s0 Right arm Ndata 
 int  NRs2lt,NRs2rt,NRs2la,NRs2ra; //s2 Right arm Ndata
 int NLs0lt,NLs0rt,NLs0la,NLs0ra;//s0 Left arm Ndata 
 int  NLs2lt,NLs2rt,NLs2la,NLs2ra; //s2 Left arm Ndata
 int NR_tr_beta,NL_tr_beta; // beta Ndata
 //================================//
 //====== Fill Ndata ==============//
 //================================//



 //================== HRS Right Arm ===============================//
 //--- Right arm AC1 and AC2 --------/
 t1->SetBranchAddress("Ndata.R.a1.a",&NRa1a);
 t1->SetBranchAddress("Ndata.R.a2.a",&NRa2a);
 t1->SetBranchAddress("Ndata.R.a1.t",&NRa1t);
 t1->SetBranchAddress("Ndata.R.a2.t",&NRa2t);
 t1->SetBranchAddress("Ndata.R.a1.a_c",&NRa1ac);
 t1->SetBranchAddress("Ndata.R.a2.a_c",&NRa2ac);

 tnew->Branch("NRa1a",&NRa1a,"NRa1a/I"); // Right arm AC1 ADC              
 tnew->Branch("NRa1t",&NRa1t,"NRa1t/I"); // Right arm AC1 TDC              
 tnew->Branch("NRa2a",&NRa1a,"NRa2a/I"); // Right arm AC2 ADC              
 tnew->Branch("NRa2t",&NRa1t,"NRa2t/I"); // Right arm AC2 TDC
              
 tnew->Branch("NRa1ac",&NRa1ac,"NRa1ac/I"); // Right arm AC1 pedestal sub ADC  
 tnew->Branch("NRa2ac",&NRa2ac,"NRa2ac/I"); // Right arm AC2 pedestal sub ADC
 //---- Right arm S0 ----------// 
 t1->SetBranchAddress("Ndata.R.s0.la",&NRs0la);
 t1->SetBranchAddress("Ndata.R.s0.ra",&NRs0ra);
 t1->SetBranchAddress("Ndata.R.s0.lt",&NRs0lt);
 t1->SetBranchAddress("Ndata.R.s0.rt",&NRs0rt);
 tnew->Branch("NRs0ra",&NRs0ra,"NRs0ra/I"); // Right arm s0 ADC               
 tnew->Branch("NRs0la",&NRs0la,"NRs0la/I"); // Right arm s0 ADC              
 tnew->Branch("NRs0rt",&NRs0rt,"NRs0rt/I"); // Right arm s0 TDC              
 tnew->Branch("NRs0lt",&NRs0lt,"NRs0lt/I"); // Right arm s0 TDC              
 //---- Right arm S2 ----------// 
 t1->SetBranchAddress("Ndata.R.s2.la",&NRs2la);
 t1->SetBranchAddress("Ndata.R.s2.ra",&NRs2ra);
 t1->SetBranchAddress("Ndata.R.s2.lt",&NRs2lt);
 t1->SetBranchAddress("Ndata.R.s2.rt",&NRs2rt);
 tnew->Branch("NRs2ra",&NRs2ra,"NRs2ra/I"); // Right arm s2 ADC               
 tnew->Branch("NRs2la",&NRs2la,"NRs2la/I"); // Right arm s2 ADC              
 tnew->Branch("NRs2rt",&NRs2rt,"NRs2rt/I"); // Right arm s2 TDC              
 tnew->Branch("NRs2lt",&NRs2lt,"NRs2lt/I"); // Right arm s2 TDC      

 //============== HRS Left Arm ==========================================//
 //---- Left arm S0 ----------// 
 t1->SetBranchAddress("Ndata.L.s0.la",&NLs0la);
 t1->SetBranchAddress("Ndata.L.s0.ra",&NLs0ra);
 t1->SetBranchAddress("Ndata.L.s0.lt",&NLs0lt);
 t1->SetBranchAddress("Ndata.L.s0.rt",&NLs0rt);
 tnew->Branch("NLs0ra",&NLs0ra,"NLs0ra/I"); // Left arm s0 ADC               
 tnew->Branch("NLs0la",&NLs0la,"NLs0la/I"); // Left arm s0 ADC              
 tnew->Branch("NLs0rt",&NLs0rt,"NLs0rt/I"); // Left arm s0 TDC              
 tnew->Branch("NLs0lt",&NLs0lt,"NLs0lt/I"); // Left arm s0 TDC              
 //---- Left arm S2 ----------// 
 t1->SetBranchAddress("Ndata.L.s2.la",&NLs2la);
 t1->SetBranchAddress("Ndata.L.s2.ra",&NLs2ra);
 t1->SetBranchAddress("Ndata.L.s2.lt",&NLs2lt);
 t1->SetBranchAddress("Ndata.L.s2.rt",&NLs2rt);
 tnew->Branch("NLs2ra",&NLs2ra,"NLs2ra/I"); // Left arm s2 ADC               
 tnew->Branch("NLs2la",&NLs2la,"NLs2la/I"); // Left arm s2 ADC              
 tnew->Branch("NLs2rt",&NLs2rt,"NLs2rt/I"); // Left arm s2 TDC              
 tnew->Branch("NLs2lt",&NLs2lt,"NLs2lt/I"); // Left arm s2 TDC              


 //============= Particle Ndata ===================================//        
 //----- beta --------------//
 t1->SetBranchAddress("Ndata.R.tr.beta",&NR_tr_beta);
 t1->SetBranchAddress("Ndata.L.tr.beta",&NL_tr_beta);
 tnew->Branch("NR_tr_beta",&NR_tr_beta,"NR_tr_beta/I"); // Right arm beta              
 tnew->Branch("NR_tr_beta",&NL_tr_beta,"NL_tr_beta/I"); // Left arm beta              


 //=========== Fill Ndata ===================================//
 double ent = t1->GetEntries();
 if(ent>evn && evn>0){
   ent=evn;
 }
  cout<<"event counts "<<ent<<endl;  
  for(int i=0 ; i<ent ; i++){
       t1->GetEntry(i);
       tnew->Fill();
  }
 

  //==========================================================//



 // (AC1)Aerogel Chrenkov Right ARM ADC //                                      
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
  t1->SetBranchAddress("R.tr.beta",R_tr_beta);
  t1->SetBranchAddress("L.tr.beta",L_tr_beta);
 
 //Energy inforamtion //
  t1->SetBranchAddress("EKL.omega",&EKLNU);                                            
  t1->SetBranchAddress("EKR.omega",&EKRNU);                             

  
  ///======= particle information =================//
  // beta //
  tnew->Branch("L_tr_beta",L_tr_beta,"L_tr_beta[NL_tr_beta]/D"); // Left arm beta
  tnew->Branch("R_tr_beta",R_tr_beta,"R_tr_beta[NR_tr_beta]/D"); // Right arm beta 
  // momentum //
  tnew->Branch("EKLQ2",&EKLQ2,"EKLQ2/D"); // Left arm beta
  tnew->Branch("EKRQ2",&EKRQ2,"EKRQ2/D"); // Right arm beta 
  //  Energy //
  tnew->Branch("EKLNU",&EKLNU,"EKLNU/D"); // (Ee-Ee') energy transfer
  tnew->Branch("EKRNU",&EKRNU,"EKRNU/D"); // (Q-Ep) energy transfer 

 

 // Right Arm A1 //                                                                   
  tnew->Branch("Ra1a",Ra1a,"Ra1a[NRa1a]/D"); // Right arm AC1 ADC                        
  tnew->Branch("Ra1t",Ra1t,"Ra1t[NRa1t]/D"); // Right arm AC1 TDC                        
  tnew->Branch("Ra1sum",&Ra1sum,"Ra1sum/D"); // Right arm AC1 FADC sum                 
  tnew->Branch("Ra1a_p",Ra1a_p,"Ra1a_p/D"); // Right arm AC1 FADC (Pedestal)PMT        
  tnew->Branch("Ra1a_c",Ra1a_c,"Ra1a_c[NRa1ac]/D"); // Right arm AC1 FADC (correction) PMT
  tnew->Branch("Ra1npe",Ra1npe,"Ra1npe[NRa1ac]/D"); // Right arm AC1 NPE             
                             
  // Right Arm AC2 //                                                                  
  tnew->Branch("Ra2a",Ra2a,"Ra2a[NRa2a]/D"); // Right arm AC2 ADC                         
  tnew->Branch("Ra2t",Ra2t,"Ra2t[NRa2t]/D"); // Right arm AC2 TDC                         
  tnew->Branch("Ra2sum",&Ra2sum,"Ra2sum/D"); // Right arm AC2 FADC sum                 
  tnew->Branch("Ra2a_p",Ra2a_p,"Ra2a_p/D"); // Right arm AC2 FADC (Pedestal)PMT        
  tnew->Branch("Ra2a_c",Ra2a_c,"Ra2a_c[NRa2ac]/D"); // Right arm AC2 FADC (correction) PMT
  tnew->Branch("Ra2npe",Ra1npe,"Ra1npe[NRa1ac]/D"); // Right arm AC2 NPE

 // Right Arm S0 // 
  tnew->Branch("Rs0la",Rs0la,"Rs0la[NRs0la]/D"); // Right arm S0-Top(B) ADC
  tnew->Branch("Rs0ra",Rs0ra,"Rs0ra[NRs0ra]/D"); // Right arm S0-Bottom(A) ADC
  tnew->Branch("Rs0lt",Rs0lt,"Rs0lt[NRs0lt]/D"); // Right arm S0-Top(B) TDC
  tnew->Branch("Rs0rt",Rs0rt,"Rs0rt[NRs0rt]/D"); // Right arm S0-Bottom(A) TDC
  tnew->Branch("Rs0la_p",Rs0la_p,"Rs01a_p/D"); // Right arm S0-Topo(B) ADC
  tnew->Branch("Rs0ra_p",Rs0ra_p,"Rs0ra_p/D"); // Right arm S2-Bottom(A) ADC
  // Right Arm S2 // 
  tnew->Branch("Rs2la",Rs2la,"Rs2la[NRs2la]/D"); // Right arm S2-Top(B) ADC
  tnew->Branch("Rs2ra",Rs2ra,"Rs2ra[NRs2ra]/D"); // Right arm S2-Bottom(A) ADC
  tnew->Branch("Rs2lt",Rs2lt,"Rs2lt[NRs2lt]/D"); // Right arm S2-Top(B) TDC
  tnew->Branch("Rs2rt",Rs2rt,"Rs2rt[NRs2rt]/D"); // Right arm S2-Bottom(A) TDC
  tnew->Branch("Rs2la_p",Rs2la_p,"Rs21a_p/D"); // Right arm S2-Topo(B) ADC
  tnew->Branch("Rs2ra_p",Rs2ra_p,"Rs2ra_p/D"); // Right arm S2-Bottom(A) ADC


  // Left Arm S0 // 
  tnew->Branch("Ls0la",Ls0la,"Ls0la[NLs0la]/D"); // Left arm S0-Top(B) ADC
  tnew->Branch("Ls0ra",Ls0ra,"Ls0ra[NLs0ra]/D"); // Left arm S0-Bottom(A) ADC
  tnew->Branch("Ls0lt",Ls0lt,"Ls0lt[NLs0lt]/D"); // Left arm S0-Top(B) TDC
  tnew->Branch("Ls0rt",Ls0rt,"Ls0rt[NLs0rt]/D"); // Left arm S0-Bottom(A) TDC
  tnew->Branch("Ls0la_p",Ls0la_p,"Ls01a_p/D"); // Left arm S0-Topo(B) ADC
  tnew->Branch("Ls0ra_p",Ls0ra_p,"Ls0ra_p/D"); // Leftt arm S2-Bottom(A) ADC
  // Left Arm S2 // 
  tnew->Branch("Ls2la",Ls2la,"Ls2la[NLs2la]/D"); // Left arm S2-Top(B) ADC
  tnew->Branch("Ls2ra",Ls2ra,"Ls2ra[NLs2ra]/D"); // Left arm S2-Bottom(A) ADC
  tnew->Branch("Ls2lt",Ls2lt,"Ls2lt[NLs2lt]/D"); // Left arm S2-Top(B) TDC
  tnew->Branch("Ls2rt",Ls2rt,"Ls2rt[NLs2rt]/D"); // Left arm S2-Bottom(A) TDC
  tnew->Branch("Ls2la_p",Ls2la_p,"Ls21a_p/D"); // Left arm S2-Topo(B) ADC
  tnew->Branch("Ls2ra_p",Ls2ra_p,"Ls2ra_p/D"); // Left arm S2-Bottom(A) ADC
  // missing mass //
  tnew->Branch("mm",&mm,"mm/D");
 



 if(ent>evn && evn>0){
   ent=evn;
 }
  cout<<"event counts "<<ent<<endl;  

  for(int i=0 ; i<ent ; i++){
       t1->GetEntry(i);
       
       if(NR_tr_beta>=NL_tr_beta){
	 Nbeta=NR_tr_beta;
       }else{
	 Nbeta=NL_tr_beta;
}
       for(int j=0;j<Nbeta+1;j++){
       mm=sqrt(pow(E+M-(E-EKLNU)-EKRNU,2)-pow(sqrt(pow(E,2)-pow(me,2))-L_tr_beta[j]*(E-EKLNU)-R_tr_beta[j]*EKRNU,2));
       }  

       for(int j=0;j<NRa1ac;j++){
	 Ra1npe[j]=;      
 }
 tnew->Fill();

  }
 
   tnew->Write();
   fnew->Close();
 
}

