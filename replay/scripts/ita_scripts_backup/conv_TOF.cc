/*
  conv.cc
Auther K. Itabashi
2018 Nov. 1st
Convert tritium_rootfile to missing mass hist
*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
using namespace std;
#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TRandom.h"


int main(int argc, char** argv){
  TApplication *theApp =new TApplication("App",&argc,argv);
  int target;
  target=0;


  int run,evn,Nbeta;
  double E,M,me,mm,mk,mtr;
 int nrun_st=111160;
 int nrun_end=111161;
 cout<<"========================================="<<endl;
 cout<<"======== Recreate new ROOTFiles ========="<<endl;
 cout<<"========================================="<<endl;
 cout<<"Start run : ";
 cin>>nrun_st;
 cout<<"End run : ";
 cin>>nrun_end; 
 cout<<"Target (H;0,3H;1,aother;-1) : ";
 cin>>target;

  me=0.511e-3; //electron mass [GeV/c^2]
  mk=493.7e-3; //Kaon mass [GeV/c^2]
  mtr=938.27e-3;// target mass [GeV/c^2] 
  if(target==1){
    mtr=3.016;// target mass [GeV/c^2]
  }else{
  mtr=938.27e-3;// target mass [GeV/c^2] 
  }



 TChain* t1=new TChain("T");
 

 for(int i=nrun_st;i<nrun_end+1;i++){
   t1->Add(Form("/adaqfs/home/a-onl/tritium/replay/t2root/itabashi/Rootfiles/tritium_%d.root",i)); 
   for(int j=0;j<5;j++){
     t1->Add(Form("/adaqfs/home/a-onl/tritium/replay/t2root/itabashi/Rootfiles/tritium_%d_%d.root",i,j)); 
}

}
 char* tar;
 if(target==0){tar="Hydrogen target run";}else if(target==1){tar="Tritium target run";}else{tar="Optics run ";}
  TTree* tnew = new TTree("T",Form("TChain %s: %d- %d ",tar,nrun_st,nrun_end));
  TFile* fnew;
  fnew =new TFile(Form("/adaqfs/home/a-onl/tritium/replay/t2root/itabashi/ita_rootfiles/tritium_coin%d_%d.root",nrun_st,nrun_end),"recreate");


  double Rsha[100],Rpsa[100],Rs0ra[100],Rs0la[100],Rs0lt[100],Rs0rt[100],Rs0lt_c[100],Rs0rt_c[100],Rs0la_p[100],Rs0ra_p[100],Rs2ra[100],Rs2la[100],Rs2lt[100],Rs2rt[100],Rs2lt_c[100],Rs2rt_c[100],Rs2la_p[100],Rs2ra_p[100];
 double Lsha[100],Lpsa[10000],Ls0ra[10000],Ls0la[10000],Ls0lt[10000],Ls0rt[10000],Ls0lt_c[10000],Ls0rt_c[10000],Ls0la_p[10000],Ls0ra_p[10000],Ls2ra[10000],Ls2la[10000],Ls2lt[10000],Ls2rt[10000],Ls2lt_c[10000],Ls2rt_c[10000],Ls2la_p[10000],Ls2ra_p[100];
 double Rs0ra_c[1000],Rs0la_c[10000],Rs2ra_c[1000],Rs2la_c[1000],Ls0ra_c[1000],Ls0la_c[100],Ls2ra_c[100],Ls2la_c[100];
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
 double Rp[100],Rpx[100],Rpy[100],Rvz[100];
 double Lp[100],Lpx[100],Lpy[100],Lvz[100];
 double Rph[100],Rth[100],Rx[100];
 double Lph[100],Lth[100],Lx[100];
 double Rtr_pathl[100],Ltr_pathl[100];
 double Rs2_pathl[100],Rs0_pathl[100],Ls2_pathl[100],Ls0_pathl[100]; 
 double DRT1[100],DRT2[100],DRT3[100],DRT4[100],DRT5[100],DRT6[100]; 
 

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
 // (AC1)Aerogel Chrenkov Right ARM ADC //                                    
 t1->SetBranchAddress("R.a1.asum_c",&Ra1sum);
 t1->SetBranchAddress("R.a2.asum_c",&Ra2sum);
 
 
  // S2 Right arm ADC//
  t1->SetBranchAddress("R.s2.rnhits",Rs2rnhits);
  t1->SetBranchAddress("R.s2.lnhits",Rs2lnhits);   
  t1->SetBranchAddress("R.s2.t_pads",Rs2pads); 
  // S2 Left arm ADC//
  t1->SetBranchAddress("L.s2.rnhits",Ls2rnhits);
  t1->SetBranchAddress("L.s2.lnhits",Ls2lnhits);   
  t1->SetBranchAddress("L.s2.t_pads",Ls2pads);


  // F1TDC //
  t1->SetBranchAddress("RTDC.F1FirstHit",RF1);                                 
  t1->SetBranchAddress("LTDC.F1FirstHit",LF1);     
  // Path Lenght //
  t1->SetBranchAddress("R.tr.pathl",Rtr_pathl);                                
  t1->SetBranchAddress("R.s2.trpath",Rs2_pathl);                               
  t1->SetBranchAddress("L.tr.pathl",Ltr_pathl);   
  t1->SetBranchAddress("L.s2.trpath",Ls2_pathl);
  // at Target angle //                                
  t1->SetBranchAddress("R.tr.th",Rth);                      
  t1->SetBranchAddress("L.tr.th",Lth);                       
  t1->SetBranchAddress("R.tr.p",Rp);                       
  t1->SetBranchAddress("R.tr.px",Rpx);                       
  t1->SetBranchAddress("R.tr.py",Rpy);                       
  t1->SetBranchAddress("R.tr.vz",Rvz);                       
  t1->SetBranchAddress("L.tr.p",Lp);                       
  t1->SetBranchAddress("L.tr.px",Lpx);                       
  t1->SetBranchAddress("L.tr.py",Lpy);                       
  t1->SetBranchAddress("R.tr.x",Rx);                       
  t1->SetBranchAddress("L.tr.x",Lx);          
  t1->SetBranchAddress("L.tr.vz",Lvz);
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
  tnew->Branch("L.tr.pathl",Ltr_pathl,"Ltr_pathl[16]/D"); 
  tnew->Branch("L.s2.trpath",Ls2_pathl,"Ls2_pathl[16]/D"); 
  // F1TDC //
  tnew->Branch("RTDC.F1FirstHit",RF1,"RF1[100]/D");                            
  tnew->Branch("LTDC.F1FirstHit",LF1,"LF1[100]/D");  
  // at Target angle and Z potion //
  tnew->Branch("R.tr.th",Rth,"Rth[5]/D");  
  tnew->Branch("R.tr.ph",Rph,"Rph[5]/D");
  tnew->Branch("L.tr.th",Lth,"Lth[5]/D");  
  tnew->Branch("L.tr.ph",Lph,"Lph[5]/D");
  tnew->Branch("R.tr.x",Rx,"Rx[5]/D");  
  tnew->Branch("L.tr.x",Lx,"Lx[5]/D");  
  tnew->Branch("R.tr.vz",Rx,"Rvz[5]/D");  
  tnew->Branch("L.tr.vz",Rx,"Lvz[5]/D");   
 // Right Arm A1 //    
  tnew->Branch("R.a1.asum_c",&Ra1sum,"Ra1sum/D"); // Right arm AC1 FADC sum    
  // Right Arm AC2 //                                                       
  tnew->Branch("R.a2.asum_c",&Ra2sum,"Ra2sum/D"); // Right arm AC2 FADC sum
  /*
  // Right Arm S2 // 
  tnew->Branch("R.s2.lt_c",Rs2lt_c,"Rs2lt_c[54]/D"); // Right arm S2-Top(B) TDC
  tnew->Branch("R.s2.rt_c",Rs2rt_c,"Rs2rt_c[54]/D"); // Right arm S2-Bottom(A) TDC
  
  // Left Arm S2 // 
  tnew->Branch("L.s2.lt_c",Ls2lt_c,"Ls2lt_c[54]/D"); // Left arm S2-Top(B) TDC
  tnew->Branch("L.s2.rt_c",Ls2rt_c,"Ls2rt_c[54]/D"); // Left arm S2-Bottom(A) TDC
  */
  tnew->Branch("R.s2.rnhits",Rs2rnhits,"Rs2rnhits[10]/D"); //
  tnew->Branch("R.s2.lnhits",Rs2lnhits,"Rs2lnhits[10]/D"); //
  tnew->Branch("R.s2.t_pads",Rs2pads,"Rs2pads[16]/D"); //  
  tnew->Branch("L.s2.rnhits",Ls2rnhits,"Ls2rnhits[10]/D"); //
  tnew->Branch("L.s2.lnhits",Ls2lnhits,"Ls2lnhits[10]/D"); //
  tnew->Branch("L.s2.t_pads",Ls2pads,"Ls2pads[16]/D"); // 
 


 int ent=t1->GetEntries();
 cout<<"Get Event: "<<ent<<endl;
  for(int i=0 ; i<ent ; i++){
       t1->GetEntry(i);
       tnew->Fill();
  }
 
   tnew->Write();
   fnew->Close();

 theApp->Run();
 return 0;

}

