/*
  conv.cc
Auther K. Itabashi
2018 Nov. 1st
Convert tritium_rootfile to missing mass hist
*/


void conv(){

  int target;
  target=0;
  int run,evn,Nbeta;
  double E,M,me,mm,mk,mtr;
  me=0.511e-3; //electron mass [GeV/c^2]
  mk=493.7e-3; //Kaon mass [GeV/c^2]
  mtr=938.27e-3;// target mass [GeV/c^2] 
  if(target==1){
    mtr=3.016;// target mass [GeV/c^2]
  }else{
  mtr=938.27e-3;// target mass [GeV/c^2] 
  }



 TChain* t1=new TChain("T");
 int nrun_st=111160;
 int nrun_end=111161;
 cout<<"=== Recreate new ROOTFiles ====="<<endl;
 cout<<"Start run : ";
 cin>>nrun_st;
 cout<<"End run : ";
 cin>>nrun_end; 

 for(int i=nrun_st;i<nrun_end+1;i++){
   t1->Add(Form("/adaqfs/home/a-onl/tritium/replay/t2root/itabashi/Rootfiles/tritium_%d.root",i)); 
   for(int j=0;j<5;j++){
     t1->Add(Form("/adaqfs/home/a-onl/tritium/replay/t2root/itabashi/Rootfiles/tritium_%d_%d.root",i,j)); 
}

}
 
  TTree* tnew = new TTree("T","recreate");
  TFile* fnew;
  fnew =new TFile(Form("/adaqfs/home/a-onl/tritium/replay/t2root/itabashi/ita_rootfiles/tritium_ita%d_%d.root",nrun_st,nrun_end),"recreate");


 double Rsha[100],Rpsa[100],Rs0ra[100],Rs0la[100],Rs0lt[100],Rs0rt[100],Rs0lt_c[100],Rs0rt_c[100],Rs0la_p[100],Rs0ra_p[100],Rs2ra[100],Rs2la[100],Rs2lt[100],Rs2rt[100],Rs2lt_c[100],Rs2rt_c[100],Rs2la_p[100],Rs2ra_p[100];
 double Lsha[100],Lpsa[10000],Ls0ra[10000],Ls0la[10000],Ls0lt[10000],Ls0rt[10000],Ls0lt_c[10000],Ls0rt_c[10000],Ls0la_p[10000],Ls0ra_p[10000],Ls2ra[10000],Ls2la[10000],Ls2lt[10000],Ls2rt[10000],Ls2lt_c[10000],Ls2rt_c[10000],Ls2la_p[10000],Ls2ra_p[100];
 double Rs0ra_c[1000],Rs0la_c[10000],Rs2ra_c[1000],Rs2la_c[1000],Ls0ra_c[1000],Ls0la_c[10000],Ls2ra_c[1000],Ls2la_c[1000];
 double Ls2pads[100],Rs2pads[100];
 double Ls2rnhits[100],Rs2rnhits[100];
 double Ls2lnhits[100],Rs2lnhits[100];
 double Ls0rnhits[100],Rs0rnhits[100];
 double Ls0lnhits[100],Rs0lnhits[100];
 double Ra2t[100],Ra2a[100],Ra2sum,Ra2a_p[100],Ra2a_c[100],Ra2npe[100];
 double Ra1t[100],Ra1a[100],Ra1sum,Ra1a_p[100],Ra1a_c[100],Ra1npe[100],R_tr_beta[100],L_tr_beta[100];
 double RF1[100],LF1[100];
 double Beam_p;
 double trig;
 double Rp[100],Rpx[10000],Rpy[10000];
 double Lp[100],Lpx[10000],Lpy[10000];
 double Rph[10000],Rth[10000];
 double Lph[10000],Lth[10000];
 double Rtr_pathl[10000],Ltr_pathl[10000];
 double Rs2_pathl[10000],Rs0_pathl[10000],Ls2_pathl[10000],Ls0_pathl[10000]; 
 double coin_time[100][100];
 double DRT1[100],DRT2[100],DRT3[100],DRT4[100],DRT5[100],DRT6[100]; 
 double DLT1[100],DLT2[100],DLT3[100],DLT4[100],DLT5[100],DLT6[100]; 
 

 //================================//
 //====== Fill Ndata ==============//
 //================================//
 
 
  // Beam branch //
 t1->SetBranchAddress("HALLA_p",&Beam_p);
 // Trigger branch //
 t1->SetBranchAddress("DR.evtypebits",&trig);
 t1->SetBranchAddress("DR.T1",DRT1);
 t1->SetBranchAddress("DR.T2",DRT2);
 t1->SetBranchAddress("DR.T3",DRT3);
 t1->SetBranchAddress("DR.T4",DRT4);
 t1->SetBranchAddress("DR.T5",DRT5);
 t1->SetBranchAddress("DR.T6",DRT6);
 t1->SetBranchAddress("DL.T1",DLT1);
 t1->SetBranchAddress("DL.T2",DLT2);
 t1->SetBranchAddress("DL.T3",DLT3);
 t1->SetBranchAddress("DL.T4",DLT4);
 t1->SetBranchAddress("DL.T5",DLT5);
 t1->SetBranchAddress("DL.T6",DLT6);
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
 // S0 Right arm //
  t1->SetBranchAddress("R.s0.la",Rs0la);
  t1->SetBranchAddress("R.s0.ra",Rs0ra);
  t1->SetBranchAddress("R.s0.la_c",Rs0la_c);
  t1->SetBranchAddress("R.s0.ra_c",Rs0ra_c);
  t1->SetBranchAddress("R.s0.lt",Rs0lt);
  t1->SetBranchAddress("R.s0.rt",Rs0rt);
  t1->SetBranchAddress("R.s0.lt_c",Rs0lt_c);
  t1->SetBranchAddress("R.s0.rt_c",Rs0rt_c);
  t1->SetBranchAddress("R.s0.la_p",Rs0ra_p);
  t1->SetBranchAddress("R.s0.ra_p",Rs0la_p);
  t1->SetBranchAddress("R.s0.rnhits",Rs0rnhits);
  t1->SetBranchAddress("R.s0.lnhits",Rs0lnhits);
  // S2 Right arm ADC//
  t1->SetBranchAddress("R.s2.la",Rs2la);
  t1->SetBranchAddress("R.s2.ra",Rs2ra);
  t1->SetBranchAddress("R.s2.la_c",Rs2la_c);
  t1->SetBranchAddress("R.s2.ra_c",Rs2ra_c);
  t1->SetBranchAddress("R.s2.lt",Rs2lt);
  t1->SetBranchAddress("R.s2.rt",Rs2rt);
  t1->SetBranchAddress("R.s2.lt_c",Rs2lt_c);
  t1->SetBranchAddress("R.s2.rt_c",Rs2rt_c);
  t1->SetBranchAddress("R.s2.la_p",Rs2ra_p);
  t1->SetBranchAddress("R.s2.ra_p",Rs2la_p);
  t1->SetBranchAddress("R.s2.rnhits",Rs2rnhits);
  t1->SetBranchAddress("R.s2.lnhits",Rs2lnhits);   
  t1->SetBranchAddress("R.s2.t_pads",Rs2pads); 
  // S0 Left arm ADC//
  t1->SetBranchAddress("L.s0.la",Ls0la);
  t1->SetBranchAddress("L.s0.ra",Ls0ra);
  t1->SetBranchAddress("L.s0.la_c",Ls0la_c);
  t1->SetBranchAddress("L.s0.ra_c",Ls0ra_c);
  t1->SetBranchAddress("L.s0.lt",Ls0lt);
  t1->SetBranchAddress("L.s0.rt",Ls0rt);
  t1->SetBranchAddress("L.s0.lt_c",Ls0lt_c);
  t1->SetBranchAddress("L.s0.rt_c",Ls0rt_c);
  t1->SetBranchAddress("L.s0.la_p",Ls0ra_p);
  t1->SetBranchAddress("L.s0.ra_p",Ls0la_p);
  t1->SetBranchAddress("L.s0.rnhits",Ls0rnhits);
  t1->SetBranchAddress("L.s0.lnhits",Ls0lnhits);
  // S2 Left arm ADC//
  t1->SetBranchAddress("L.s2.la",Ls2la);
  t1->SetBranchAddress("L.s2.ra",Ls2ra);
  t1->SetBranchAddress("L.s2.la_c",Ls2la_c);
  t1->SetBranchAddress("L.s2.ra_c",Ls2ra_c);
  t1->SetBranchAddress("L.s2.lt",Ls2lt);
  t1->SetBranchAddress("L.s2.rt",Ls2rt);
  t1->SetBranchAddress("L.s2.lt_c",Ls2lt_c);
  t1->SetBranchAddress("L.s2.rt_c",Ls2rt_c);
  t1->SetBranchAddress("L.s2.la_p",Ls2ra_p);
  t1->SetBranchAddress("L.s2.ra_p",Ls2la_p);
  t1->SetBranchAddress("L.s2.rnhits",Ls2rnhits);
  t1->SetBranchAddress("L.s2.lnhits",Ls2lnhits);   
  t1->SetBranchAddress("L.s2.t_pads",Ls2pads);
  // beta information //
  t1->SetBranchAddress("R.tr.beta",R_tr_beta);
  t1->SetBranchAddress("L.tr.beta",L_tr_beta);
  // F1TDC //
  t1->SetBranchAddress("RTDC.F1FirstHit",RF1);                                  
  t1->SetBranchAddress("LTDC.F1FirstHit",LF1);     
  // Path Lenght //
  t1->SetBranchAddress("R.tr.pathl",Rtr_pathl);                                
  t1->SetBranchAddress("R.s0.trpath",Rs0_pathl);                                
  t1->SetBranchAddress("R.s2.trpath",Rs2_pathl);                                
  t1->SetBranchAddress("L.tr.pathl",Ltr_pathl);   
  t1->SetBranchAddress("L.s0.trpath",Ls0_pathl);     
  t1->SetBranchAddress("L.s2.trpath",Ls2_pathl);
  // at Target angle //                                
  t1->SetBranchAddress("R.tr.ph",Rph);               
  t1->SetBranchAddress("R.tr.th",Rth);                      
  t1->SetBranchAddress("L.tr.ph",Lph);               
  t1->SetBranchAddress("L.tr.th",Lth);                       
  t1->SetBranchAddress("R.tr.p",Rp);                       
  t1->SetBranchAddress("R.tr.px",Rpx);                       
  t1->SetBranchAddress("R.tr.py",Rpy);                       
  t1->SetBranchAddress("L.tr.p",Lp);                       
  t1->SetBranchAddress("L.tr.px",Lpx);                       
  t1->SetBranchAddress("L.tr.py",Lpy);                       
 
 ///======= particle information =================//
  // Beam //
  tnew->Branch("HALLA_p",&Beam_p,"Beam_p/D"); // Electron beam momentom
  tnew->Branch("DR.evtypebits",&trig,"trig/D"); // Trigger conditon
  tnew->Branch("DR.T1",DRT1,"DRT1[10]/D"); // Trigger conditon
  tnew->Branch("DR.T2",DRT2,"DRT2[10]/D"); // Trigger conditon
  tnew->Branch("DR.T3",DRT3,"DRT3[10]/D"); // Trigger conditon
  tnew->Branch("DR.T4",DRT4,"DRT4[10]/D"); // Trigger conditon
  tnew->Branch("DR.T5",DRT5,"DRT5[10]/D"); // Trigger conditon
  tnew->Branch("DR.T6",DRT6,"DRT6[10]/D"); // Trigger conditon
  // beta //
  tnew->Branch("L.tr.beta",L_tr_beta,"L_tr_beta[5]/D"); // Left arm beta
  tnew->Branch("R.tr.beta",R_tr_beta,"R_tr_beta[5]/D"); // Right arm beta 
  // momentum //
  tnew->Branch("R.tr.p",Rp,"Rp[5]/D"); // Right arm momentom
  tnew->Branch("R.tr.px",Rpx,"Rpx[5]/D"); // Right arm x-momentom 
  tnew->Branch("R.tr.py",Rpy,"Rpx[5]/D"); // Right arm y-momentom 
  tnew->Branch("L.tr.p",Lp,"Lp[5]/D"); // Left arm momentom
  tnew->Branch("L.tr.px",Lpx,"Lpx[5]/D"); // Left arm x-momentom 
  tnew->Branch("L.tr.py",Lpy,"Lpx[5]/D"); // Left arm y-momentom 
  // Path Length //
  tnew->Branch("R.tr.pathl",Rtr_pathl,"Rtr_pathl[16]/D"); 
  tnew->Branch("R.s2.trpath",Rs2_pathl,"Rs2_pathl[16]/D"); 
  tnew->Branch("R.s0.trpath",Rs0_pathl,"Rs0_pathl[1]/D"); 
  tnew->Branch("L.tr.pathl",Ltr_pathl,"Ltr_pathl[16]/D"); 
  tnew->Branch("L.s2.trpath",Ls2_pathl,"Ls2_pathl[16]/D"); 
  tnew->Branch("L.s0.trpath",Ls0_pathl,"Ls0_pathl[1]/D"); 
 // F1TDC //
  tnew->Branch("RTDC.F1FirstHit",RF1,"RF1[100]/D");                                       
  tnew->Branch("LTDC.F1FirstHit",LF1,"LF1[100]/D");  
  // at Target angle //
  tnew->Branch("R.tr.th",Rth,"Rth[5]/D");  
  tnew->Branch("R.tr.ph",Rph,"Rph[5]/D");
  tnew->Branch("L.tr.th",Lth,"Lth[5]/D");  
  tnew->Branch("L.tr.ph",Lph,"Lph[5]/D");
 // Right Arm A1 //    
  tnew->Branch("R.a1.a",Ra1a,"Ra1a[23]/D"); // Right arm AC1 ADC 
  tnew->Branch("R.a1.t",Ra1t,"Ra1t[23]/D"); // Right arm AC1 TDC  
  tnew->Branch("R.a1.asum_c",&Ra1sum,"Ra1sum/D"); // Right arm AC1 FADC sum    
  tnew->Branch("R.a1.a_p",Ra1a_p,"Ra1a_p[23]/D"); // Right arm AC1 FADC (Pedestal)
  tnew->Branch("R.a1.a_c",Ra1a_c,"Ra1a_c[23]/D"); // Right arm AC1 FADC (correction) PMT
  // Right Arm AC2 //                                                                 
  tnew->Branch("R.a2.a",Ra2a,"Ra2a[26]/D"); // Right arm AC2 ADC             
  tnew->Branch("R.a2.t",Ra2t,"Ra2t[26]/D"); // Right arm AC2 TDC             
  tnew->Branch("R.a2.asum_c",&Ra2sum,"Ra2sum/D"); // Right arm AC2 FADC sum
  tnew->Branch("R.a2.a_p",Ra2a_p,"Ra2a_p[26]/D"); // Right arm AC2 FADC (Pedestal)PM
  tnew->Branch("R.a2.a_c",Ra2a_c,"Ra2a_c[26]/D"); // Right arm AC2 FADC (correction) PMT
 // Right Arm S0 // 
  tnew->Branch("R.s0.la",Rs0la,"Rs0la[5]/D"); // Right arm S0-Top(B) ADC
  tnew->Branch("R.s0.ra",Rs0ra,"Rs0ra[5]/D"); // Right arm S0-Bottom(A) ADC
  tnew->Branch("R.s0.la_c",Rs0la_c,"Rs0la_c[5]/D"); // Right arm S0-Top(B) ADC
  tnew->Branch("R.s0.ra_c",Rs0ra_c,"Rs0ra_c[5]/D"); // Right arm S0-Bottom(A) ADC
  tnew->Branch("R.s0.lt",Rs0lt,"Rs0lt[5]/D"); // Right arm S0-Top(B) TDC
  tnew->Branch("R.s0.rt",Rs0rt,"Rs0rt[5]/D"); // Right arm S0-Bottom(A) TDC
  tnew->Branch("R.s0.lt_c",Rs0lt_c,"Rs0lt_c[5]/D"); // Right arm S0-Top(B) TDC
  tnew->Branch("R.s0.rt_c",Rs0rt_c,"Rs0rt_c[5]/D"); // Right arm S0-Bottom(A) TDC
  tnew->Branch("R.s0.la_p",Rs0la_p,"Rs01a_p/D"); // Right arm S0-Topo(B) ADC
  tnew->Branch("R.s0.ra_p",Rs0ra_p,"Rs0ra_p/D"); // Right arm S2-Bottom(A) ADC
  tnew->Branch("R.s0.rnhits",Rs0rnhits,"Rs0rnhits[10]/D"); //
  tnew->Branch("R.s0.lnhits",Rs0lnhits,"Rs0lnhits[10]/D"); //
  // Right Arm S2 // 
  tnew->Branch("R.s2.la",Rs2la,"Rs2la[54]/D"); // Right arm S2-Top(B) ADC
  tnew->Branch("R.s2.ra",Rs2ra,"Rs2ra[54]/D"); // Right arm S2-Bottom(A) ADC
tnew->Branch("R.s2.la_c",Rs2la_c,"Rs2la[54]/D"); // Right arm S2-Top(B) ADC
  tnew->Branch("R.s2.ra_c",Rs2ra_c,"Rs2ra[54]/D"); // Right arm S2-Bottom(A) ADC
  tnew->Branch("R.s2.lt",Rs2lt,"Rs2lt[54]/D"); // Right arm S2-Top(B) TDC
  tnew->Branch("R.s2.rt",Rs2rt,"Rs2rt[54]/D"); // Right arm S2-Bottom(A) TDC
  tnew->Branch("R.s2.lt_c",Rs2lt_c,"Rs2lt_c[54]/D"); // Right arm S2-Top(B) TDC
  tnew->Branch("R.s2.rt_c",Rs2rt_c,"Rs2rt_c[54]/D"); // Right arm S2-Bottom(A) TDC
  tnew->Branch("R.s2.la_p",Rs2la_p,"Rs21a_p/D"); // Right arm S2-Topo(B) ADC
  tnew->Branch("R.s2.ra_p",Rs2ra_p,"Rs2ra_p/D"); // Right arm S2-Bottom(A) ADC
  tnew->Branch("R.s2.rnhits",Rs2rnhits,"Rs2rnhits[10]/D"); //
  tnew->Branch("R.s2.lnhits",Rs2lnhits,"Rs2lnhits[10]/D"); //
  tnew->Branch("R.s2.t_pads",Rs2pads,"Rs2pads[16]/D"); //
  

  // Left Arm S0 // 
  tnew->Branch("L.s0.la",Ls0la,"Ls0la[3]/D"); // Left arm S0-Top(B) ADC
  tnew->Branch("L.s0.ra_c",Ls0ra,"Ls0ra[3]/D"); // Left arm S0-Bottom(A) ADC
  tnew->Branch("L.s0.la_c",Ls0la_c,"Ls0la_c[3]/D"); // Left arm S0-Top(B) ADC
  tnew->Branch("L.s0.ra",Ls0ra_c,"Ls0ra_c[3]/D"); // Left arm S0-Bottom(A) AD
  tnew->Branch("L.s0.lt",Ls0lt,"Ls0lt[3]/D"); // Left arm S0-Top(B) TDC
  tnew->Branch("L.s0.rt",Ls0rt,"Ls0rt[3]/D"); // Left arm S0-Bottom(A) TDC
  tnew->Branch("L.s0.lt_c",Ls0lt_c,"Ls0lt_c[3]/D"); // Left arm S0-Top(B) TDC
  tnew->Branch("L.s0.rt_c",Ls0rt_c,"Ls0rt_c[3]/D"); // Left arm S0-Bottom(A) TDC
  tnew->Branch("L.s0.la_p",Ls0la_p,"Ls01a_p/D"); // Left arm S0-Topo(B) ADC
  tnew->Branch("L.s0.ra_p",Ls0ra_p,"Ls0ra_p/D"); // Leftt arm S2-Bottom(A) ADC
  tnew->Branch("L.s0.rnhits",Ls0rnhits,"Ls0rnhits[10]/D"); //
  tnew->Branch("L.s0.lnhits",Ls0lnhits,"Ls0lnhits[10]/D"); //
  
  // Left Arm S2 // 
  tnew->Branch("L.s2.la",Ls2la,"Ls2la[54]/D"); // Left arm S2-Top(B) ADC
  tnew->Branch("L.s2.ra",Ls2ra,"Ls2ra[54]/D"); // Left arm S2-Bottom(A) ADC
  tnew->Branch("L.s2.la_c",Ls2la_c,"Ls2la[54]/D"); // Left arm S2-Top(B) ADC
  tnew->Branch("L.s2.ra_c",Ls2ra_c,"Ls2ra[54]/D"); // Left arm S2-Bottom(A) ADC
  tnew->Branch("L.s2.lt",Ls2lt,"Ls2lt[54]/D"); // Left arm S2-Top(B) TDC
  tnew->Branch("L.s2.rt",Ls2rt,"Ls2rt[54]/D"); // Left arm S2-Bottom(A) TDC
  tnew->Branch("L.s2.lt_c",Ls2lt_c,"Ls2lt_c[54]/D"); // Left arm S2-Top(B) TDC
  tnew->Branch("L.s2.rt_c",Ls2rt_c,"Ls2rt_c[54]/D"); // Left arm S2-Bottom(A) TDC
  tnew->Branch("L.s2.la_p",Ls2la_p,"Ls21a_p/D"); // Left arm S2-Topo(B) ADC
  tnew->Branch("L.s2.ra_p",Ls2ra_p,"Ls2ra_p/D"); // Left arm S2-Bottom(A) ADC
  tnew->Branch("L.s2.rnhits",Ls2rnhits,"Ls2rnhits[10]/D"); //
  tnew->Branch("L.s2.lnhits",Ls2lnhits,"Ls2lnhits[10]/D"); //
  tnew->Branch("L.s2.t_pads",Ls2pads,"Ls2pads[16]/D"); // 
  // missing mass //
  tnew->Branch("mm",&mm,"mm/D");
  // Coincidence Time R-L //
  tnew->Branch("coin_time",coin_time,"coin_time[16][16]/D");
  double pe,pe_,pk,Ee,Ee_,Ek;


 int ent=t1->GetEntries();
 cout<<"Get Event: "<<ent<<endl;
  for(int i=0 ; i<ent ; i++){
       t1->GetEntry(i);
       pe=Beam_p;
       pe_=Lp[0]*sqrt(1-pow(Lth[0],2)+pow(Lph[0],2));
       pk=Rp[0]*sqrt(1-pow(Rth[0],2)+pow(Rph[0],2));
       Ee=sqrt(pow(pe,2)+pow(me,2)); 
       Ee_=sqrt(pow(pe_,2)+pow(me,2));
       Ek=sqrt(pow(pk,2)+pow(mk,2));
       mm=sqrt(pow(Ee+mtr-Ee_-Ek,2)-pow(pe-pe_-pk,2));
       for(int r=0;r<16;r++){
	 for(int l=0;l<16;l++){
	   coin_time[r][l]=(Rs2rt_c[r]+Rs2lt_c[r])/2.0-(Ls2rt_c[l]+Ls2lt_c[l])/2.0;
    }
       }
   tnew->Fill();
  }
 
   tnew->Write();
   fnew->Close();

}
