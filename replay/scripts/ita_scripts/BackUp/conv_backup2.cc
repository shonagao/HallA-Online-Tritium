/*
  conv.cc
Auther K. Itabashi
2018 Nov. 1st
Convert tritium_rootfile to missing mass hist
*/


void conv_backup2(){

  int target;
  cout<<"Target(H:0,3H:1,Other:-1):  "<<endl;
  cin>>target;
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
 int nrun_st=111166;
 int nrun_end=111166;
 
 for(int i=nrun_st;i<nrun_end+1;i++){
   t1->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",i)); //a-onl
 }
 
  TTree* tnew = new TTree("T","recreate");
  TFile* fnew;
  /*
if(target==1){
  fnew= new TFile(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/ita_Rootfiles/ita%d-%d_nnLcoin.root",nrun_st,nrun_end),"recreate");
 }else if(target==0){
  fnew = new TFile(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/ita_Rootfiles/ita%d-%d_Hcoin.root",nrun_st,nrun_end),"recreate");
 }else{
  fnew = new TFile(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/ita_Rootfiles/ita%d-%d_Optcoin.root",nrun_st,nrun_end),"recreate");
}
  */
  //fnew = new TFile("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/ita_Rootfiles/ita_trititum_calib.root","recreate");
fnew = new TFile("./ita_trititum_calib.root","recreate");

 double Rsha[100],Rpsa[100],Rs0ra[100],Rs0la[100],Rs0lt[100],Rs0rt[100],Rs0lt_c[100],Rs0rt_c[100],Rs0la_p[100],Rs0ra_p[100],Rs2ra[100],Rs2la[100],Rs2lt[100],Rs2rt[100],Rs2lt_c[100],Rs2rt_c[100],Rs2la_p[100],Rs2ra_p[100];
 double Lsha[10000],Lpsa[10000],Ls0ra[10000],Ls0la[10000],Ls0lt[10000],Ls0rt[10000],Ls0lt_c[10000],Ls0rt_c[10000],Ls0la_p[10000],Ls0ra_p[10000],Ls2ra[10000],Ls2la[10000],Ls2lt[10000],Ls2rt[10000],Ls2lt_c[10000],Ls2rt_c[10000],Ls2la_p[10000],Ls2ra_p[10000];
 double Rs0ra_c[1000],Rs0la_c[10000],Rs2ra_c[1000],Rs2la_c[1000],Ls0ra_c[1000],Ls0la_c[10000],Ls2ra_c[1000],Ls2la_c[1000];
 double Ra2t[10000],Ra2a[10000],Ra2sum,Ra2a_p[10000],Ra2a_c[10000],Ra2npe[10000];
 double Ra1t[10000],Ra1a[10000],Ra1sum,Ra1a_p[10000],Ra1a_c[10000],Ra1npe[10000],R_tr_beta[100],L_tr_beta[100];
 double RF1[100000],LF1[100000];
 double Beam_p;
 double trig;
 double Rp[10000],Rpx[10000],Rpy[10000];
 double Lp[10000],Lpx[10000],Lpy[10000];
 double Rph[10000],Rth[10000];
 double Lph[10000],Lth[10000];
 double Rtr_pathl[10000],Ltr_pathl[10000];
 double Rs2_pathl[10000],Rs0_pathl[10000],Ls2_pathl[10000],Ls0_pathl[10000]; 
 double coin_time[100][100];
 int NRa1a,NRa2a,NRa1t,NRa2t,NRa1ac,NRa2ac; // A1 and A2 Ndata
 int NRs0lt,NRs0rt,NRs0la,NRs0ra,NRs0lt_c,NRs0rt_c;//s0 Right arm Ndata 
 int NRs2lt,NRs2rt,NRs2la,NRs2ra,NRs2lt_c,NRs2rt_c; //s2 Right arm Ndata
 int NLs0lt,NLs0rt,NLs0lt_c,NLs0rt_c,NLs0la,NLs0ra;//s0 Left arm Ndata 
 int NLs2lt,NLs2rt,NLs2lt_c,NLs2rt_c,NLs2la,NLs2ra; //s2 Left arm Ndata
 int NR_tr_beta,NL_tr_beta; // beta Ndata
 int NRF1,NLF1;
 int NRp,NRpx,NRpy,NLp,NLpx,NLpy;
 int NRph,NRth,NLph,NLth;
 int NRtr_pathl,NLtr_pathl,NRs2_pathl,NRs0_pathl,NLs2_pathl,NLs0_pathl;  
 int NRs0ra_c,NRs0la_c,NRs2ra_c,NRs2la_c,NLs0ra_c,NLs0la_c,NLs2ra_c,NLs2la_c;
 //================================//
 //====== Fill Ndata ==============//
 //================================//
 // Momentum //
 t1->SetBranchAddress("Ndata.R.tr.p",&NRp);
 t1->SetBranchAddress("Ndata.R.tr.px",&NRpx); 
 t1->SetBranchAddress("Ndata.R.tr.py",&NRpy); 
 t1->SetBranchAddress("Ndata.L.tr.p",&NLp);
 t1->SetBranchAddress("Ndata.L.tr.px",&NLpx); 
 t1->SetBranchAddress("Ndata.L.tr.py",&NLpy);
 tnew->Branch("Ndata.R.tr.p",&NRp,"NRp/I"); 
 tnew->Branch("Ndata.R.tr.px",&NRpx,"NRpx/I");
 tnew->Branch("Ndata.R.tr.py",&NRpy,"NRpy/I");
 tnew->Branch("Ndata.L.tr.p",&NLp,"NLp/I"); 
 tnew->Branch("Ndata.L.tr.px",&NLpx,"NLpx/I");
 tnew->Branch("Ndata.tr.L.py",&NLpy,"NLpy/I");
 // Path Lenght //
 t1->SetBranchAddress("Ndata.R.tr.pathl",&NRtr_pathl);
 t1->SetBranchAddress("Ndata.R.s2.trpath",&NRs2_pathl); 
 t1->SetBranchAddress("Ndata.R.s0.trpath",&NRs0_pathl); 
 t1->SetBranchAddress("Ndata.L.tr.pathl",&NLtr_pathl);
 t1->SetBranchAddress("Ndata.L.s2.trpath",&NLs2_pathl); 
 t1->SetBranchAddress("Ndata.L.s0.trpath",&NLs0_pathl); 
 tnew->Branch("Ndata.R.tr.path",&NRtr_pathl,"NRtr_pathl/I"); 
 tnew->Branch("Ndata.R.s2.trpath",&NRs2_pathl,"NRs2_pathl/I");
 tnew->Branch("Ndata.R.s0.trpath",&NRs0_pathl,"NRs0_pathl/I");
 tnew->Branch("Ndata.L.tr.pathl",&NLtr_pathl,"NLtr_pathl/I"); 
 tnew->Branch("Ndata.L.s2.trpath",&NLs2_pathl,"NLs2_pathl/I");
 tnew->Branch("Ndata.L.s0.trpath",&NLs0_pathl,"NLs0_pathl/I");
 // at Target angle //
 t1->SetBranchAddress("Ndata.L.tr.ph",&NLph); 
 t1->SetBranchAddress("Ndata.L.tr.th",&NLth);
 t1->SetBranchAddress("Ndata.R.tr.ph",&NRph); 
 t1->SetBranchAddress("Ndata.R.tr.th",&NRth);
 tnew->Branch("Ndata.R.ph",&NRph,"NRph/I"); 
 tnew->Branch("Ndata.R.th",&NRth,"NRth/I");
 tnew->Branch("Ndata.L.ph",&NLph,"NLph/I"); 
 tnew->Branch("Ndata.L.th",&NLth,"NLth/I");
 //================== HRS Right Arm ===============================//
 //--- Right arm AC1 and AC2 --------/
 t1->SetBranchAddress("Ndata.R.a1.a",&NRa1a);
 t1->SetBranchAddress("Ndata.R.a2.a",&NRa2a);
 t1->SetBranchAddress("Ndata.R.a1.t",&NRa1t);
 t1->SetBranchAddress("Ndata.R.a2.t",&NRa2t);
 t1->SetBranchAddress("Ndata.R.a1.a_c",&NRa1ac);
 t1->SetBranchAddress("Ndata.R.a2.a_c",&NRa2ac);

 tnew->Branch("Ndata.R.a1.a",&NRa1a,"NRa1a/I"); // Right arm AC1 ADC            
 tnew->Branch("Ndata.R.a1.t",&NRa1t,"NRa1t/I"); // Right arm AC1 TDC           
 tnew->Branch("Ndata.R.a2.a",&NRa1a,"NRa2a/I"); // Right arm AC2 ADC            
 tnew->Branch("Ndata.R.a2.t",&NRa1t,"NRa2t/I"); // Right arm AC2 TDC
 tnew->Branch("Ndata.R.a1.a_c",&NRa1ac,"NRa1ac/I"); // Right arm AC1 pedestal sub ADC  
 tnew->Branch("Ndata.R.a2.a_c",&NRa2ac,"NRa2ac/I"); // Right arm AC2 pedestal sub ADC
 //---- Right arm S0 ----------// 
 t1->SetBranchAddress("Ndata.R.s0.la",&NRs0la);
 t1->SetBranchAddress("Ndata.R.s0.ra",&NRs0ra);
 t1->SetBranchAddress("Ndata.R.s0.la_c",&NRs0la_c);
 t1->SetBranchAddress("Ndata.R.s0.ra_c",&NRs0ra_c);
 t1->SetBranchAddress("Ndata.R.s0.lt",&NRs0lt);
 t1->SetBranchAddress("Ndata.R.s0.rt",&NRs0rt);
 t1->SetBranchAddress("Ndata.R.s0.lt_c",&NRs0lt_c);
 t1->SetBranchAddress("Ndata.R.s0.rt_c",&NRs0rt_c);
 tnew->Branch("Ndata.R.s0.ra_c",&NRs0ra,"NRs0ra_c/I"); // Right arm s0 ADC
 tnew->Branch("Ndata.R.s0.la_c",&NRs0la,"NRs0la_c/I"); // Right arm s0 ADC
 tnew->Branch("Ndata.R.s0.ra",&NRs0ra,"NRs0ra/I"); // Right arm s0 ADC          
 tnew->Branch("Ndata.R.s0.la",&NRs0la,"NRs0la/I"); // Right arm s0 ADC         
 tnew->Branch("Ndata.R.s0.rt",&NRs0rt,"NRs0rt/I"); // Right arm s0 TDC          
 tnew->Branch("Ndata.R.s0.lt",&NRs0lt,"NRs0lt/I"); // Right arm s0 TDC          
 tnew->Branch("Ndata.R.s0.rt_c",&NRs0rt_c,"NRs0rt_c/I"); // Right arm s0 TDC    
 tnew->Branch("Ndata.R.s0.lt_c",&NRs0lt_c,"NRs0lt_c/I"); // Right arm s0 TDC              
 //---- Right arm S2 ----------// 
 t1->SetBranchAddress("Ndata.R.s2.la",&NRs2la);
 t1->SetBranchAddress("Ndata.R.s2.ra",&NRs2ra);
 t1->SetBranchAddress("Ndata.R.s2.la_c",&NRs2la_c);
 t1->SetBranchAddress("Ndata.R.s2.ra_c",&NRs2ra_c);
 t1->SetBranchAddress("Ndata.R.s2.lt",&NRs2lt);
 t1->SetBranchAddress("Ndata.R.s2.rt",&NRs2rt);
 t1->SetBranchAddress("Ndata.R.s2.lt_c",&NRs2lt_c);
 t1->SetBranchAddress("Ndata.R.s2.rt_c",&NRs2rt_c);
 tnew->Branch("Ndata.R.s2.ra",&NRs2ra,"NRs2ra/I"); // Right arm s2 ADC          
 tnew->Branch("Ndata.R.s2.la",&NRs2la,"NRs2la/I"); // Right arm s2 ADC          
 tnew->Branch("Ndata.R.s2.ra_c",&NRs2ra_c,"NRs2ra_c/I"); // Right arm s2 ADC
 tnew->Branch("Ndata.R.s2.la_c",&NRs2la_c,"NRs2la_c/I"); // Right arm s2 ADC   
 tnew->Branch("Ndata.R.s2.rt",&NRs2rt,"NRs2rt/I"); // Right arm s2 TDC          
 tnew->Branch("Ndata.R.s2.lt",&NRs2lt,"NRs2lt/I"); // Right arm s2 TDC
 tnew->Branch("Ndata.R.s2.rt_c",&NRs2rt_c,"NRs2rt_c/I"); // Right arm s2 TDC    
 tnew->Branch("Ndata.R.s2.lt_c",&NRs2lt_c,"NRs2lt_c/I"); // Right arm s2 TDC
 t1->SetBranchAddress("Ndata.RTDC.F1FirstHit",&NRF1);
 tnew->Branch("Ndata.RTDC.F1FirstHit",&NLF1,"NRF1/I"); // Right arm F1tdc

 //============== HRS Left Arm ==========================================//
 //---- Left arm S0 ----------// 
 t1->SetBranchAddress("Ndata.L.s0.la",&NLs0la);
 t1->SetBranchAddress("Ndata.L.s0.ra",&NLs0ra);
 t1->SetBranchAddress("Ndata.L.s0.lt",&NLs0lt);
 t1->SetBranchAddress("Ndata.L.s0.rt",&NLs0rt);
 t1->SetBranchAddress("Ndata.L.s0.lt_c",&NLs0lt_c);
 t1->SetBranchAddress("Ndata.L.s0.rt_c",&NLs0rt_c);
 tnew->Branch("Ndata.L.s0.ra",&NLs0ra,"NLs0ra/I"); // Left arm s0 ADC           
 tnew->Branch("Ndata.L.s0.la",&NLs0la,"NLs0la/I"); // Left arm s0 ADC           
 tnew->Branch("Ndata.L.s0.ra_c",&NLs0ra_c,"NLs0ra_c/I"); // Left arm s0 ADC
 tnew->Branch("Ndata.L.s0.la_c",&NLs0la_c,"NLs0la_c/I"); // Left arm s0 ADC
 tnew->Branch("Ndata.L.s0.rt",&NLs0rt,"NLs0rt/I"); // Left arm s0 TDC          
 tnew->Branch("Ndata.L.s0.lt",&NLs0lt,"NLs0lt/I"); // Left arm s0 TDC  
 tnew->Branch("Ndata.L.s0.rt_c",&NLs0rt_c,"NLs0rt_c/I"); // Left arm s0 TDC    
 tnew->Branch("Ndata.L.s0.lt_c",&NLs0lt_c,"NLs0lt_c/I"); // Left arm s0 TDC
 //---- Left arm S2 ----------// 
 t1->SetBranchAddress("Ndata.L.s2.la",&NLs2la);
 t1->SetBranchAddress("Ndata.L.s2.ra",&NLs2ra);
 t1->SetBranchAddress("Ndata.L.s2.lt",&NLs2lt);
 t1->SetBranchAddress("Ndata.L.s2.rt",&NLs2rt);
 t1->SetBranchAddress("Ndata.L.s2.lt_c",&NLs2lt_c);
 t1->SetBranchAddress("Ndata.L.s2.rt_c",&NLs2rt_c);
 tnew->Branch("Ndata.L.s2.ra",&NLs2ra,"NLs2ra/I"); // Left arm s2 ADC    
 tnew->Branch("Ndata.L.s2.la",&NLs2la,"NLs2la/I"); // Left arm s2 ADC
 tnew->Branch("Ndata.L.s2.ra_c",&NLs2ra_c,"NLs2ra_c/I"); // Left arm s2 ADC    
 tnew->Branch("Ndata.L.s2.la_c",&NLs2la_c,"NLs2la_c/I"); // Left arm s2 ADC
 tnew->Branch("Ndata.L.s2.rt",&NLs2rt,"NLs2rt/I"); // Left arm s2 TDC           
 tnew->Branch("Ndata.L.s2.lt",&NLs2lt,"NLs2lt/I"); // Left arm s2 TDC
 tnew->Branch("Ndata.L.s2.rt_c",&NLs2rt_c,"NLs2rt_c/I"); // Left arm s2 TDC  
 tnew->Branch("Ndata.L.s2.lt_c",&NLs2lt_c,"NLs2lt_c/I"); // Left arm s2 TDC 
 t1->SetBranchAddress("Ndata.LTDC.F1FirstHit",&NLF1);
 tnew->Branch("Ndata.LTDCF1FirstHit",&NLF1,"NLF1/I"); // Left arm F1tdc              
 //============= Particle Ndata ===================================//        
 //----- beta --------------//
 t1->SetBranchAddress("Ndata.R.tr.beta",&NR_tr_beta);
 t1->SetBranchAddress("Ndata.L.tr.beta",&NL_tr_beta);
 tnew->Branch("Ndata.R.tr.beta",&NR_tr_beta,"NR_tr_beta/I"); // Right arm beta
 tnew->Branch("Ndata.L.tr.beta",&NL_tr_beta,"NL_tr_beta/I"); // Left arm beta              

 //=========== Fill Ndata ===================================//
 double ent = t1->GetEntries();
  cout<<"event counts "<<ent<<endl;  
  for(int i=0 ; i<ent ; i++){
       t1->GetEntry(i);
       tnew->Fill();
  }
 

  //==========================================================//

  // Beam branch //
 t1->SetBranchAddress("HALLA_p",&Beam_p);
 // Trigger branch //
 t1->SetBranchAddress("DR.evtypebits",&trig);
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
  t1->SetBranchAddress("R.s0.la_p",Rs0ra_p);
  t1->SetBranchAddress("R.s0.ra_p",Rs0la_p);

  // S2 Right arm ADC//
  t1->SetBranchAddress("R.s2.la",Rs2la);
  t1->SetBranchAddress("R.s2.ra",Rs2ra);
  t1->SetBranchAddress("R.s2.la_c",Rs2la_c);
  t1->SetBranchAddress("R.s2.ra_c",Rs2ra_c);
  t1->SetBranchAddress("R.s2.lt",Rs2lt);
  t1->SetBranchAddress("R.s2.rt",Rs2rt);
  t1->SetBranchAddress("R.s2.la_p",Rs2ra_p);
  t1->SetBranchAddress("R.s2.ra_p",Rs2la_p);
   
  // S0 Left arm ADC//
  t1->SetBranchAddress("L.s0.la",Ls0la);
  t1->SetBranchAddress("L.s0.ra",Ls0ra);
  t1->SetBranchAddress("L.s0.la_c",Ls0la_c);
  t1->SetBranchAddress("L.s0.ra_c",Ls0ra_c);
  t1->SetBranchAddress("L.s0.lt",Ls0lt);
  t1->SetBranchAddress("L.s0.rt",Ls0rt);
  t1->SetBranchAddress("L.s0.la_p",Ls0ra_p);
  t1->SetBranchAddress("L.s0.ra_p",Ls0la_p);
  // S2 Left arm ADC//
  t1->SetBranchAddress("L.s2.la",Ls2la);
  t1->SetBranchAddress("L.s2.ra",Ls2ra);
  t1->SetBranchAddress("L.s2.la_c",Ls2la_c);
  t1->SetBranchAddress("L.s2.ra_c",Ls2ra_c);
  t1->SetBranchAddress("L.s2.lt",Ls2lt);
  t1->SetBranchAddress("L.s2.rt",Ls2rt);
  t1->SetBranchAddress("L.s2.la_p",Ls2ra_p);
  t1->SetBranchAddress("L.s2.ra_p",Ls2la_p);
  // beta information //
  t1->SetBranchAddress("R.tr.beta",R_tr_beta);
  t1->SetBranchAddress("L.tr.beta",L_tr_beta);
  // F1TDC //
  t1->SetBranchAddress("RTDC.F1FirstHit",RF1);                                    t1->SetBranchAddress("LTDC.F1FirstHit",LF1);     
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
  tnew->Branch("DR.ev.typebits",&trig,"trig/D"); // Trigger conditon
  // beta //
  tnew->Branch("L.tr.beta",L_tr_beta,"L_tr_beta[NL_tr_beta]/D"); // Left arm beta
  tnew->Branch("R.tr.beta",R_tr_beta,"R_tr_beta[NR_tr_beta]/D"); // Right arm beta 
  // momentum //
  tnew->Branch("R.tr.p",Rp,"Rp[NRp]/D"); // Right arm momentom
  tnew->Branch("R.tr.px",Rpx,"Rpx[NRpx]/D"); // Right arm x-momentom 
  tnew->Branch("R.tr.py",Rpy,"Rpx[NRpy]/D"); // Right arm y-momentom 
  tnew->Branch("L.tr.p",Lp,"Lp[NLp]/D"); // Left arm momentom
  tnew->Branch("L.tr.px",Lpx,"Lpx[NLpx]/D"); // Left arm x-momentom 
  tnew->Branch("L.tr.py",Lpy,"Lpx[NLpy]/D"); // Left arm y-momentom 

  // Path Length //
  tnew->Branch("R.tr.pathl",Rtr_pathl,"Rtr_pathl[NRtr_pathl]/D"); 
  tnew->Branch("R.s2.trpath",Rs2_pathl,"Rs2_pathl[NRs2_pathl]/D"); 
  tnew->Branch("R.s0.trpath",Rs0_pathl,"Rs0_pathl[NRs0_pathl]/D"); 
  tnew->Branch("L.tr.pathl",Ltr_pathl,"Ltr_pathl[NLtr_pathl]/D"); 
  tnew->Branch("L.s2.trpath",Ls2_pathl,"Ls2_pathl[NLs2_pathl]/D"); 
  tnew->Branch("L.s0.trpath",Ls0_pathl,"Ls0_pathl[NLs0_pathl]/D"); 
 // F1TDC //
  tnew->Branch("RTDC.F1FirstHit",RF1,"RF1[NRF1]/D");                                       
  tnew->Branch("LTDC.F1FirstHit",LF1,"LF1[NLF1]/D");  
  // at Target angle //
  tnew->Branch("R.tr.th",Rth,"Rth[NRth]/D");  
  tnew->Branch("R.tr.ph",Rph,"Rph[NRph]/D");
  tnew->Branch("L.tr.th",Lth,"Lth[NLth]/D");  
  tnew->Branch("L.tr.ph",Lph,"Lph[NLph]/D");
 // Right Arm A1 //    
  tnew->Branch("R.a1.a",Ra1a,"Ra1a[NRa1a]/D"); // Right arm AC1 ADC 
  tnew->Branch("R.a1.t",Ra1t,"Ra1t[NRa1t]/D"); // Right arm AC1 TDC  
  tnew->Branch("R.a1.asum_c",&Ra1sum,"Ra1sum/D"); // Right arm AC1 FADC sum    
  tnew->Branch("R.a1.a_p",Ra1a_p,"Ra1a_p/D"); // Right arm AC1 FADC (Pedestal)
  tnew->Branch("R.a1.a_c",Ra1a_c,"Ra1a_c[NRa1ac]/D"); // Right arm AC1 FADC (correction) PMT

  // Right Arm AC2 //                                                                 
  tnew->Branch("R.a2.a",Ra2a,"Ra2a[NRa2a]/D"); // Right arm AC2 ADC             
  tnew->Branch("R.a2.t",Ra2t,"Ra2t[NRa2t]/D"); // Right arm AC2 TDC             
  tnew->Branch("R.a2.asum_c",&Ra2sum,"Ra2sum/D"); // Right arm AC2 FADC sum
  tnew->Branch("R.a2.a_p",Ra2a_p,"Ra2a_p/D"); // Right arm AC2 FADC (Pedestal)PM
  tnew->Branch("R.a2.a_c",Ra2a_c,"Ra2a_c[NRa2ac]/D"); // Right arm AC2 FADC (correction) PMT
 // Right Arm S0 // 
  tnew->Branch("R.s0.la",Rs0la,"Rs0la[NRs0la]/D"); // Right arm S0-Top(B) ADC
  tnew->Branch("R.s0.ra",Rs0ra,"Rs0ra[NRs0ra]/D"); // Right arm S0-Bottom(A) ADC
tnew->Branch("R.s0.la_c",Rs0la_c,"Rs0la_c[NRs0la_c]/D"); // Right arm S0-Top(B) ADC
  tnew->Branch("R.s0.ra_c",Rs0ra_c,"Rs0ra_c[NRs0ra_c]/D"); // Right arm S0-Bottom(A) ADC
  tnew->Branch("R.s0.lt",Rs0lt,"Rs0lt[NRs0lt]/D"); // Right arm S0-Top(B) TDC
  tnew->Branch("R.s0.rt",Rs0rt,"Rs0rt[NRs0rt]/D"); // Right arm S0-Bottom(A) TDC
  tnew->Branch("R.s0.lt_c",Rs0lt_c,"Rs0lt_c[NRs0lt_c]/D"); // Right arm S0-Top(B) TDC
  tnew->Branch("R.s0.rt_c",Rs0rt_c,"Rs0rt_c[NRs0rt_c]/D"); // Right arm S0-Bottom(A) TDC
  tnew->Branch("R.s0.la_p",Rs0la_p,"Rs01a_p/D"); // Right arm S0-Topo(B) ADC
  tnew->Branch("R.s0.ra_p",Rs0ra_p,"Rs0ra_p/D"); // Right arm S2-Bottom(A) ADC
  // Right Arm S2 // 
  tnew->Branch("R.s2.la",Rs2la,"Rs2la[NRs2la]/D"); // Right arm S2-Top(B) ADC
  tnew->Branch("R.s2.ra",Rs2ra,"Rs2ra[NRs2ra]/D"); // Right arm S2-Bottom(A) ADC
tnew->Branch("R.s2.la_c",Rs2la_c,"Rs2la[NRs2la_c]/D"); // Right arm S2-Top(B) ADC
  tnew->Branch("R.s2.ra_c",Rs2ra_c,"Rs2ra[NRs2ra_c]/D"); // Right arm S2-Bottom(A) ADC
  tnew->Branch("R.s2.lt",Rs2lt,"Rs2lt[NRs2lt]/D"); // Right arm S2-Top(B) TDC
  tnew->Branch("R.s2.rt",Rs2rt,"Rs2rt[NRs2rt]/D"); // Right arm S2-Bottom(A) TDC
  tnew->Branch("R.s2.lt_c",Rs2lt_c,"Rs2lt_c[NRs2lt_c]/D"); // Right arm S2-Top(B) TDC
  tnew->Branch("R.s2.rt_c",Rs2rt_c,"Rs2rt_c[NRs2rt_c]/D"); // Right arm S2-Bottom(A) TDC
  tnew->Branch("R.s2.la_p",Rs2la_p,"Rs21a_p/D"); // Right arm S2-Topo(B) ADC
  tnew->Branch("R.s2.ra_p",Rs2ra_p,"Rs2ra_p/D"); // Right arm S2-Bottom(A) ADC
  // Left Arm S0 // 
  tnew->Branch("L.s0.la",Ls0la,"Ls0la[NLs0la]/D"); // Left arm S0-Top(B) ADC
  tnew->Branch("L.s0.ra_c",Ls0ra,"Ls0ra[NLs0ra]/D"); // Left arm S0-Bottom(A) ADC
  tnew->Branch("L.s0.la_c",Ls0la_c,"Ls0la_c[NLs0la_c]/D"); // Left arm S0-Top(B) ADC
  tnew->Branch("L.s0.ra",Ls0ra_c,"Ls0ra_c[NLs0ra_c]/D"); // Left arm S0-Bottom(A) AD
  tnew->Branch("L.s0.lt",Ls0lt,"Ls0lt[NLs0lt]/D"); // Left arm S0-Top(B) TDC
  tnew->Branch("L.s0.rt",Ls0rt,"Ls0rt[NLs0rt]/D"); // Left arm S0-Bottom(A) TDC
  tnew->Branch("L.s0.lt_c",Ls0lt_c,"Ls0lt_c[NLs0lt_c]/D"); // Left arm S0-Top(B) TDC
  tnew->Branch("L.s0.rt_c",Ls0rt_c,"Ls0rt_c[NLs0rt_c]/D"); // Left arm S0-Bottom(A) TDC
  tnew->Branch("L.s0.la_p",Ls0la_p,"Ls01a_p/D"); // Left arm S0-Topo(B) ADC
  tnew->Branch("L.s0.ra_p",Ls0ra_p,"Ls0ra_p/D"); // Leftt arm S2-Bottom(A) ADC
  // Left Arm S2 // 
  tnew->Branch("L.s2.la",Ls2la,"Ls2la[NLs2la]/D"); // Left arm S2-Top(B) ADC
  tnew->Branch("L.s2.ra",Ls2ra,"Ls2ra[NLs2ra]/D"); // Left arm S2-Bottom(A) ADC
  tnew->Branch("L.s2.la_c",Ls2la_c,"Ls2la[NLs2la_c]/D"); // Left arm S2-Top(B) ADC
  tnew->Branch("L.s2.ra_c",Ls2ra_c,"Ls2ra[NLs2ra_c]/D"); // Left arm S2-Bottom(A) ADC
  tnew->Branch("L.s2.lt",Ls2lt,"Ls2lt[NLs2lt]/D"); // Left arm S2-Top(B) TDC
  tnew->Branch("L.s2.rt",Ls2rt,"Ls2rt[NLs2rt]/D"); // Left arm S2-Bottom(A) TDC
  tnew->Branch("L.s2.lt_c",Ls2lt_c,"Ls2lt_c[NLs2lt_c]/D"); // Left arm S2-Top(B) TDC
  tnew->Branch("L.s2.rt_c",Ls2rt_c,"Ls2rt_c[NLs2rt_c]/D"); // Left arm S2-Bottom(A) TDC
  tnew->Branch("L.s2.la_p",Ls2la_p,"Ls21a_p/D"); // Left arm S2-Topo(B) ADC
  tnew->Branch("L.s2.ra_p",Ls2ra_p,"Ls2ra_p/D"); // Left arm S2-Bottom(A) ADC
  // missing mass //
  tnew->Branch("mm",&mm,"mm/D");
  // Coincidence Time R-L //
  tnew->Branch("coin_time",coin_time,"coin_time[16][16]/D");
  double pe,pe_,pk,Ee,Ee_,Ek;

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
