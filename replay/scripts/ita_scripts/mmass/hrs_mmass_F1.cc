// Author K. Itabashi Aug. 30th
// HRS nnL experiment missing mass analysis

const double c=299792458;// [m/s]
const double mk=493.7;// Kaon mass [MeV/c^2]
const double me=0.511;// electron mass [MeV/c^2] 
const double ml=1115.7;//Lambda mass [MeV/c^2]
void hrs_mmass_F1(){

  TChain*  T=new TChain("T");
  int nrun=111145;
  //  T->Add(Form("/w/halla-scifs17exp/triton/itabashi/Tohoku_github/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",nrun));
  T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",nrun));

 //============= Set Branch Status ==================//
  int max=10000; 
  double RF1[max],LF1[max];
  double Rs0r_ac[max],Rs0l_ac[max],Ls0r_ac[max],Ls0l_ac[max];
  double Rs2r_ac[max],Rs2l_ac[max],Ls2r_ac[max],Ls2l_ac[max];
  double Rs0r_tc,Rs0l_tc,Ls0r_tc,Ls0l_tc;
  double Rs2r_tc[max],Rs2l_tc[max],Ls2r_tc[max],Ls2l_tc[max];
  double Ra1t[max],Ra1a[max],Ra1a_p[max],Ra1a_c[max],Ra1sum;
  double Ra2t[max],Ra2a[max],Ra2a_p[max],Ra2a_c[max],Ra2sum;
  double La1t[max],La1a[max],La1a_p[max],La1a_c[max],La1sum;
  double La2t[max],La2a[max],La2a_p[max],La2a_c[max],La2sum;
  double Rp[max],Rpx[max],Rpy[max],Lp[max],Lpx[max],Lpy[max];
  double Rth[max],Rph[max],Lth[max],Lph[max];
  int trig;
  double hallap;


  T->SetBranchStatus("*",0);  
  T->SetBranchStatus("HALLA_p",1);
  T->SetBranchAddress("HALLA_p",&hallap); 
 //------ Right Arm -------------//
  T->SetBranchStatus("RTDC.F1FirstHit",1);
  T->SetBranchAddress("RTDC.F1FirstHit",RF1); 
  T->SetBranchStatus("R.s0.ra_c",1);        // Right arm S0 R-PMT  ADC
  T->SetBranchAddress("R.s0.ra_c",&Rs0r_ac); // Right arm S0 R-PMT  ADC
  T->SetBranchStatus("R.s0.la_c",1);        // Right arm S0 L-PMT  ADC
  T->SetBranchAddress("R.s0.la_c",&Rs0l_ac); // Right arm S0 L-PMT  ADC
  T->SetBranchStatus("R.s2.ra_c",1);        // Right arm S2 R-PMT  ADC
  T->SetBranchAddress("R.s2.ra_c",Rs2r_ac);  // Right arm S2 R-PMT  ADC
  T->SetBranchStatus("R.s2.la_c",1);        // Right arm S2 L-PMT  ADC
  T->SetBranchAddress("R.s2.la_c",Rs2l_ac);  // Right arm S2 L-PMT  ADC

  T->SetBranchStatus("R.s0.rt_c",1);        // Right arm S0 R-PMT  ADC
  T->SetBranchAddress("R.s0.rt_c",&Rs0r_tc); // Right arm S0 R-PMT  ADC
  T->SetBranchStatus("R.s0.lt_c",1);        // Right arm S0 L-PMT  ADC
  T->SetBranchAddress("R.s0.lt_c",&Rs0l_tc); // Right arm S0 L-PMT  ADC
  T->SetBranchStatus("R.s2.rt_c",1);        // Right arm S2 R-PMT  ADC
  T->SetBranchAddress("R.s2.rt_c",Rs2r_tc);  // Right arm S2 R-PMT  ADC
  T->SetBranchStatus("R.s2.lt_c",1);        // Right arm S2 L-PMT  ADC
  T->SetBranchAddress("R.s2.lt_c",Rs2l_tc);  // Right arm S2 L-PMT  ADC

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
 T->SetBranchStatus("R.tr.p",1);
 T->SetBranchAddress("R.tr.p",Rp);
 T->SetBranchStatus("R.tr.px",1);
 T->SetBranchAddress("R.tr.px",Rpx);
 T->SetBranchStatus("R.tr.py",1);
 T->SetBranchAddress("R.tr.py",Rpy);
 T->SetBranchStatus("R.tr.ph",1);
 T->SetBranchAddress("R.tr.ph",Rph);
 T->SetBranchStatus("R.tr.th",1);
 T->SetBranchAddress("R.tr.th",Rth);

 //------ Left Arm ---------------//
  T->SetBranchStatus("LTDC.F1FirstHit",1);
  T->SetBranchAddress("LTDC.F1FirstHit",LF1); 
  T->SetBranchStatus("L.s0.ra_c",1);        // Left arm S0 R-PMT  ADC
  T->SetBranchAddress("L.s0.ra_c",&Ls0r_ac); // Left arm S0 R-PMT  ADC
  T->SetBranchStatus("L.s0.la_c",1);        // Left arm S0 L-PMT  ADC
  T->SetBranchAddress("L.s0.la_c",&Ls0l_ac); // Left arm S0 L-PMT  ADC
  T->SetBranchStatus("L.s2.ra_c",1);        // Left arm S2 R-PMT  ADC
  T->SetBranchAddress("L.s2.ra_c",Ls2r_ac);  // Left arm S2 R-PMT  ADC
  T->SetBranchStatus("L.s2.la_c",1);        // Left arm S2 L-PMT  ADC
  T->SetBranchAddress("L.s2.la_c",Ls2l_ac);  // Left arm S2 L-PMT  ADC

  T->SetBranchStatus("L.s0.rt_c",1);        // Left arm S0 R-PMT  ADC
  T->SetBranchAddress("L.s0.rt_c",&Ls0r_tc); // Left arm S0 R-PMT  ADC
  T->SetBranchStatus("L.s0.lt_c",1);        // Left arm S0 L-PMT  ADC
  T->SetBranchAddress("L.s0.lt_c",&Ls0l_tc); // Left arm S0 L-PMT  ADC
  T->SetBranchStatus("L.s2.rt_c",1);        // Left arm S2 R-PMT  ADC
  T->SetBranchAddress("L.s2.rt_c",Ls2r_tc);  // Left arm S2 R-PMT  ADC
  T->SetBranchStatus("L.s2.lt_c",1);        // Left arm S2 L-PMT  ADC
  T->SetBranchAddress("L.s2.lt_c",Ls2l_tc);  // Left arm S2 L-PMT  ADC

 T->SetBranchStatus("L.tr.p",1);
 T->SetBranchAddress("L.tr.p",Lp);  
 T->SetBranchStatus("L.tr.px",1);
 T->SetBranchAddress("L.tr.px",Lpx);
 T->SetBranchStatus("L.tr.py",1);
 T->SetBranchAddress("L.tr.py",Lpy);
 T->SetBranchStatus("L.tr.ph",1);
 T->SetBranchAddress("L.tr.ph",Lph);
 T->SetBranchStatus("L.tr.th",1);
 T->SetBranchAddress("L.tr.th",Lth);

 //==================================================//

 TH1F* hmm=new TH1F("hmm","Missing mass -mass of Lambda Hist [MeV]",5000,-1000,1000);
 TCanvas* c0=new TCanvas("c0","c0");
 TH1F* hcoin_t=new TH1F("hcoin_t","Coincidence time S2R-S2L[sec] ",5000,-1.0e-6,1.0e-6);
 TCanvas* c1=new TCanvas("c1","c1");
 int evnt=T->GetEntries();
 cout<<"Get Entries: "<<evnt<<endl;
 double mtr;
 mtr=938.27;// proton mass [MeV/c^2]
 double mh;
 double Ee,Ee_,Ek;
 double pe,pe_,pk;
 double coin_t;
 int i=8;
 for(int k=0;k<evnt;k++){
   T->GetEntry(k);
   pe_=Lp[0]*sqrt(1+pow(Lpx[0],2)+pow(Lpy[0],2))*1.0e6;
   pk=Rp[0]*sqrt(1+pow(Rpx[0],2)+pow(Rpy[0],2))*1.0e4;
 pe=hallap;
 Ee=sqrt(pow(pe,2)+pow(me,2));
 Ee_=sqrt(pow(pe_,2)+pow(me,2));
 Ek=sqrt(pow(pk,2)+pow(mk,2));
 mh=sqrt(pow(Ee+mtr-Ee_-Ek,2)-pow(pe-pe_-pk,2))-ml;
 coin_t=(Rs0r_tc+Rs0l_tc)/2.0-(Ls0r_tc+Ls0l_tc)/2.0;

 if((Rs0r_tc>0 && Rs0l_tc>0) && ( Ls0r_tc>0 && Ls0l_tc>0)) hcoin_t->Fill(coin_t);
 if((Rs0r_tc>0 && Rs0l_tc>0) && ( Ls0r_tc>0 && Ls0l_tc>0) && coin_t>-0.15e-6 && coin_t<-0.14e-6)hmm->Fill(mh);


	 /*
if(500<k && k<600){
   cout<<" k : "<<k<<endl;
   cout<<"Lp :"<<Lp[0]<<endl;
   cout<<"Rp :"<<Rp[0]<<endl;
   cout<<"Lpx: "<<Lpx[0]<<endl;
   cout<<"Rpx: "<<Rpx[0]<<endl;
   cout<<"Lpy: "<<Lpy[0]<<endl;
   cout<<"Rpy: "<<Rpy[0]<<endl;
   cout<<"pe :"<<pe<<endl;
   cout<<"pe_ :"<<pe_<<endl;
   cout<<"pk :"<<pk<<endl;
   cout<<"Ee :"<<Ee<<endl;
   cout<<"Ee_ :"<<Ee_<<endl;
   cout<<"Ek :"<<Ek<<endl;
   cout<<"Rs2r_tc :"<<Rs2r_tc[i]<<endl;
   cout<<"Rs2l_tc :"<<Rs2l_tc[i]<<endl;
   cout<<"Ls2r_tc :"<<Ls2r_tc[i]<<endl;
   cout<<"Ls2l_tc :"<<Ls2l_tc[i]<<endl;
 }
	 */

 

 }

   c0->cd();
 hmm->Draw();
 c1->cd();
 hcoin_t->Draw();
}
