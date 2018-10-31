// Author K. Itabashi Aug. 30th
// HRS nnL experiment missing mass analysis

const double c=299792458;// [m/s]
const double mk=493.7;// Kaon mass [MeV/c^2]
const double me=0.511;// electron mass [MeV/c^2] 
const double ml=1115.7;//Lambda mass [MeV/c^2]
void hrs_mmass(){

  TChain*  T=new TChain("T");
  T->Add(Form("/w/halla-scifs17exp/triton/itabashi/Tohoku_github/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",nrun));


 //============= Set Branch Status ==================//
  int max=10000; 
  double RF1[max],LF1[chmax];
  double Rs0r_ac[max],Rs0l_ac[max],Ls0r_ac[max],Ls0l_ac[max];
  double Rs2r_ac[max],Rs2l_ac[max],Ls2r_ac[max],Ls2l_ac[max];
  double Ra1t[max],Ra1a[max],Ra1a_p[max],Ra1sum;
  double Ra2t[max],Ra2a[max],Ra2a_p[max],Ra2sum;
  double Rp[max],Rpx[max],Rpy[max],Lp[max],Lpx[max],Lpy[max];
  double Rth[max],Rph[max],Lth[max],Lph[max];
  int trig;
  double hallap;


  T->SetBranchStatus("*",0);  
  T->SetBranchStatus("HALLA_p",1);
  T->SetBranchAddress("HALLA_p",hallap); 
 //------ Right Arm -------------//
  T->SetBranchStatus("RTDC.F1FirstHit",1);
  T->SetBranchAddress("RTDC.F1FirstHit",RF1); 
  T->SetBranchStatus("R.s0.ra_c",1);        // Right arm S0 R-PMT  ADC
  T->SetBranchAddress("R.s0.ra_c",&s0radc); // Right arm S0 R-PMT  ADC
  T->SetBranchStatus("R.s0.la_c",1);        // Right arm S0 L-PMT  ADC
  T->SetBranchAddress("R.s0.la_c",&s0ladc); // Right arm S0 L-PMT  ADC
  T->SetBranchStatus("R.s2.ra_c",1);        // Right arm S2 R-PMT  ADC
  T->SetBranchAddress("R.s2.ra_c",s2radc);  // Right arm S2 R-PMT  ADC
  T->SetBranchStatus("R.s2.la_c",1);        // Right arm S2 L-PMT  ADC
  T->SetBranchAddress("R.s2.la_c",s2ladc);  // Right arm S2 L-PMT  ADC

 // (AC1)Aerogel Chrenkov Right ARM ADC //         
 T->SetBranchStatus("R.a1.t",1);
 T->SetBranchStatus("R.a1.a",1);
 T->SetBranchStatus("R.a1.asum_c",1);
 T->SetBranchStatus("R.a1.a_p",1);
 T->SetBranchStatus("R.a1.a_c",1);
 T->SetBranchAddress("R.a1.t",Ra1t);
 T->SetBranchAddress("R.a1.a",Ra1a);
 T->SetBranchAddress("R.a1.asum_c",&Ra1sum);
 T->SetBranchAddress("R.a1.a_p",Ra1a_p);
 T->SetBranchAddress("R.a1.a_c",Ra1a_c);
 // (AC2)Aerogel Chrenkov Right ARM ADC //                                   
 T->SetBranchStatus("R.a2.t",1);
 T->SetBranchStatus("R.a2.a",1);
 T->SetBranchStatus("R.a2.asum_c",1);
 T->SetBranchStatus("R.a2.a_p",1);
 T->SetBranchStatus("R.a2.a_c",1);
 T->SetBranchAddress("R.a2.t",Ra2t);
 T->SetBranchAddress("R.a2.a",Ra2a);
 T->SetBranchAddress("R.a2.asum_c",&Ra2sum);
 T->SetBranchAddress("R.a2.a_p",Ra2a_p);
 T->SetBranchAddress("R.a2.a_c",Ra2a_c);

 // Target positon information //
 T->SetBranchAddress("R.tr.p",1);
 T->SetBranchAddress("R.tr.p",Rp);
 T->SetBranchAddress("R.tr.px",1);
 T->SetBranchAddress("R.tr.px",Rpx);
 T->SetBranchAddress("R.tr.py",1);
 T->SetBranchAddress("R.tr.py",Rpy);
 T->SetBranchAddress("R.tr.ph",1);
 T->SetBranchAddress("R.tr.ph",Rph);
 T->SetBranchAddress("R.tr.th",1);
 T->SetBranchAddress("R.tr.th",Rth);

 //------ Left Arm ---------------//
  T->SetBranchStatus("LTDC.F1FirstHit",1);
  T->SetBranchAddress("LTDC.F1FirstHit",LF1); 
  T->SetBranchStatus("L.s0.ra_c",1);        // Left arm S0 R-PMT  ADC
  T->SetBranchAddress("L.s0.ra_c",&s0radc); // Left arm S0 R-PMT  ADC
  T->SetBranchStatus("L.s0.la_c",1);        // Left arm S0 L-PMT  ADC
  T->SetBranchAddress("L.s0.la_c",&s0ladc); // Left arm S0 L-PMT  ADC
  T->SetBranchStatus("L.s2.ra_c",1);        // Left arm S2 R-PMT  ADC
  T->SetBranchAddress("L.s2.ra_c",s2radc);  // Left arm S2 R-PMT  ADC
  T->SetBranchStatus("L.s2.la_c",1);        // Left arm S2 L-PMT  ADC
  T->SetBranchAddress("L.s2.la_c",s2ladc);  // Left arm S2 L-PMT  ADC

 T->SetBranchAddress("L.tr.p",1);
 T->SetBranchAddress("L.tr.p",Lp);  
 T->SetBranchAddress("L.tr.px",1);
 T->SetBranchAddress("L.tr.px",Lpx);
 T->SetBranchAddress("L.tr.py",1);
 T->SetBranchAddress("L.tr.py",Lpy);
 T->SetBranchAddress("L.tr.ph",1);
 T->SetBranchAddress("L.tr.ph",Lph);
 T->SetBranchAddress("L.tr.th",1);
 T->SetBranchAddress("L.tr.th",Lth);

 //==================================================//

 TH1F* hmm=new TH1F("hmm","Missing mass -mass of Lambda Hist [MeV]",5000,-1000,1000);
 TCanvas* c0=new TCanvas("c0","c0");
 int evnt=T->GetEntries();
 cout<<"Get Entries: "<<evnt<<endl;
 double mtr;
 mtr=938.27;// proton mass [MeV/c^2]
 double mh;
 double Ee,Ee_,Ek;
 double pe,pe_,pk;

 for(int i=0;k<evnt;k++){
 T->GetEntries(k)
 pe_=Lp[0]*sqrt(1+pow(Lpx[0],2)+pow(Lpy[0],2));
 pk=Rp[0]*sqrt(1+pow(Rpx[0],2)+pow(Rpy[0],2));
 pe=hallap;
 Ee=sqrt(pow(pe,2)+pow(me,2));
 Ee_=sqrt(pow(pe_,2)+pow(me,2));
 Ek=sqrt(pow(pk,2)+pow(mk,2));
 mh=sqrt(pow(Ee+mtr-Ee_-Ek,2)-pow(pe-pe_-pk,2))-ml;
 hmm->Fill(mh);

}

 hmm->Draw();

}
