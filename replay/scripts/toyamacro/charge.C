#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>

void charge(Int_t run){
string LR="R";
double current_cut = 2; //current cut of 2 uA (minimum value of current)
if(run==111193)exit(1);

//============  Reading the Rootfile =======================//

  const TString rootfilePath = "/adaqfs/home/a-onl/tritium/replay/t2root/Rootfiles/";//a-onl machine
  std::ostringstream str;
  str << rootfilePath<<"tritium_"<<run;
  TString basename = str.str().c_str();
  TString rootfile = basename + ".root";
  TChain* T;
  TChain*ev;
  bool Left_flag=false, Right_flag=true;

  if(LR=="L")     Left_flag = true;
  else if(LR=="R")Right_flag = true;
  else {cout<<"unknown arm info."<<endl;exit(1);}

  if(Left_flag){
    T = new TChain("TSLeft");
    ev = new TChain("evLeft");
  }
  else if(Right_flag){
    T = new TChain("TSRight");
    ev = new TChain("evRight");
  }

  Long_t split=0;

  char* file = 0;

  //====adding splits rootfiles =======================//

  Long_t i=0;
  while ( !gSystem->AccessPathName(rootfile.Data()) ) {
     T->Add(rootfile.Data());
     ev->Add(rootfile.Data());
     cout << "ROOT file " << rootfile << " added to TChain." << endl;
     i++;
     rootfile = basename + "_" + i + ".root";
   }

  if(!T->GetEntries()){
     cerr<< "Not root file was found" << endl;
     return;
  }

    //==finish adding splits rootfiles=====================//  

  Double_t d3c, d10c, unewc, dnewc;
  Double_t d3r, d10r, unewr, dnewr;
  Double_t clk;


  if(Left_flag){
    ev->SetBranchAddress("evLeftd3"    ,&d3c   );
    ev->SetBranchAddress("evLeftd10"   ,&d10c  );
    ev->SetBranchAddress("evLeftunew"  ,&unewc );
    ev->SetBranchAddress("evLeftdnew"  ,&dnewc );
    ev->SetBranchAddress("evLeftd3_r"  ,&d3r   );
    ev->SetBranchAddress("evLeftd10_r" ,&d10r  );
    ev->SetBranchAddress("evLeftunew_r",&unewr );
    ev->SetBranchAddress("evLeftdnew_r",&dnewr );
    ev->SetBranchAddress("evLeftLclock",&clk   );
  }
  else if(Right_flag){
    ev->SetBranchAddress("evRightd3"    ,&d3c   );
    ev->SetBranchAddress("evRightd10"   ,&d10c  );
    ev->SetBranchAddress("evRightunew"  ,&unewc );
    ev->SetBranchAddress("evRightdnew"  ,&dnewc );
    ev->SetBranchAddress("evRightd3_r"  ,&d3r   );
    ev->SetBranchAddress("evRightd10_r" ,&d10r  );
    ev->SetBranchAddress("evRightunew_r",&unewr );
    ev->SetBranchAddress("evRightdnew_r",&dnewr );
    ev->SetBranchAddress("evRightLclock",&clk   );
  }
  
  //========Calibration Values of the bcm accorfing to the bcm[] array
  //d3, d10, unew and dnew
   Double_t g[4] = {0.0001076,3.72e-05,0.0004183, 0.0003297};
   Double_t of[4] = { 0.1694,0.05958,0.08409,0.09209};

  //============== Variables ====================================================//
  
   Int_t evnentries = ev->GetEntries();

   //==========preliminary charge with the final and initial values ======//
  
  Double_t cd3_i, cd10_i, cunew_i, cdnew_i, clk_i, cd3_f, cd10_f, cunew_f, cdnew_f, clk_f;
  Double_t Qd3_tot, Qd10_tot, Qunew_tot, Qdnew_tot;

  ev->GetEntry(0);
  cd3_i = d3c; cd10_i = d10c; cunew_i = unewc; cdnew_i = dnewc;
  clk_i = clk;

  ev->GetEntry(evnentries-1);
  cd3_f = d3c; cd10_f = d10c; cunew_f = unewc; cdnew_f = dnewc;
  clk_f = clk;
  cout<<cd3_i<<" "<<cd3_f<<endl;
  
  Qd3_tot   = (g[0]*(cd3_f-cd3_i))    + ((of[0]*(clk_f-clk_i))/103700.0 );
  Qd10_tot  = (g[1]*(cd10_f-cd10_i))  + ((of[1]*(clk_f-clk_i))/103700.0 );
  Qunew_tot = (g[2]*(cunew_f-cunew_i))+ ((of[2]*(clk_f-clk_i))/103700.0 );
  Qdnew_tot = (g[3]*(cdnew_f-cdnew_i))+ ((of[3]*(clk_f-clk_i))/103700.0 );

  cout << "The preliminary charge for d10: "  << Qd3_tot << " uC" << endl;
  cout << "The preliminary charge for d10: "  << Qd10_tot << " uC" << endl;
  cout << "The preliminary charge for unew: " << Qunew_tot << " uC" << endl;
  cout << "The preliminary charge for dnew: " << Qdnew_tot << " uC" << endl;
   
  cout << "time " << (clk_f-clk_i)/103700.0 << endl;
   //==========Current Calculation==============================================

  Double_t I_d3[evnentries], I_d10[evnentries], I_unew[evnentries], I_dnew[evnentries];
  Double_t I_d3er[evnentries], I_d10er[evnentries], I_unewer[evnentries], I_dnewer[evnentries];
  
  for (Int_t i = 0; i<evnentries ; i++) {
    ev->GetEntry(i);
    //========current
    I_d3[i]   = g[0]*d3r   + of[0];
    I_d10[i]  = g[1]*d10r  + of[1];
    I_unew[i] = g[2]*unewr + of[2];
    I_dnew[i] = g[3]*dnewr + of[3];
  }

   //============   Current Average  ================================================

   Double_t curd3, curd10, curdnew, curunew, co3, co10, counew, codnew;
   Double_t curd3er, curd10er, curdnewer, curunewer;
   Double_t curd3er1, curd10er1, curdnewer1, curunewer1;
   curd3 = curd10 = curdnew = curunew = co3 = co10 = counew = codnew = 0;

   for(Int_t i=0; i<evnentries-1; i++){
    if (I_d3[i]>current_cut && I_d3[i]<30 ){curd3+=I_d3[i]; co3+=1; }
    if (I_d10[i]>current_cut && I_d10[i]<30 ){curd10+=I_d10[i]; co10+=1; }
    if (I_unew[i]>current_cut && I_unew[i]<30 ){curunew+=I_unew[i]; counew+=1;}
    if (I_dnew[i]>current_cut && I_dnew[i]<30 ){curdnew+=I_dnew[i]; codnew+=1; }
   }
  
  curd3/=co3; curd10/=co10; curunew/=counew; curdnew/=codnew;
  
  cout << "Current d3: "   << curd3   << endl;
  cout << "Current d10: "  << curd10  << endl;
  cout << "Current unew: " << curunew << endl;
  cout << "Current dnew: " << curdnew << endl;

   //============   Charge Calculation  ================================================
   
   Double_t Time[evnentries], d3c_t[evnentries], d10c_t[evnentries], unc_t[evnentries], dnc_t[evnentries];
   Double_t d3c_T, d10c_T, unc_T, dnc_T, test;
   d3c_T = d10c_T = unc_T = dnc_T = test=0;

   for (Int_t i = 0; i < evnentries; i++) {
    ev->GetEntry(i);
    Time[i] = clk - test;
    d3c_t[i] = d3c - d3c_T; d10c_t[i] = d10c -  d10c_T;
    unc_t[i] = unewc - unc_T; dnc_t[i] = dnewc - dnc_T;    
    d3c_T = d3c; d10c_T = d10c; unc_T = unewc; dnc_T = dnewc;
    test = clk;
  }


   Double_t  Q_d3[evnentries], Q_d10[evnentries], Q_un[evnentries], Q_dn[evnentries];
   Double_t  Q_d3er[evnentries], Q_d10er[evnentries], Q_uner[evnentries], Q_dner[evnentries];
   
  for (Int_t i = 0; i<evnentries ; i++) {
  //=========charge
    Q_d3[i]  = (g[0]*d3c_t[i]  + (of[0]*Time[i]/103700.0 ));
    Q_d10[i] = (g[1]*d10c_t[i] + (of[1]*Time[i]/103700.0 ));
    Q_un[i]  = (g[2]*unc_t[i]  + (of[2]*Time[i]/103700.0 ));
    Q_dn[i]  = (g[3]*dnc_t[i]  + (of[3]*Time[i]/103700.0 ));
  }


  Double_t Charge_d3, Charge_d10, Charge_un, Charge_dn;
  Charge_d3 = Charge_d10 = Charge_un = Charge_dn = 0;
  Double_t x[evnentries];

  Double_t df[4];

  for(Int_t j=0; j<4; j++){
    df[j] = (current_cut - of[i])/g[i];
  }

  Int_t tester=0;

  for(Int_t i = 0; i<evnentries; i++){
    ev->GetEntry(i);
    x[i] = i; 
    if(d3r >df[0])  Charge_d3  += Q_d3[i];  
    if(d10r>df[1])  Charge_d10 += Q_d10[i]; 
    if(unewr>df[2]) Charge_un  += Q_un[i]; 
    if(dnewr>df[3]) Charge_dn  += + Q_dn[i];  
  }

  cout << "Total Charge d3:   " << Charge_d3  << endl;
  //cout << "Total Charge d10:  " << Charge_d10 << endl;
  cout << "Total Charge unew: " << Charge_un  << endl;
  cout << "Total Charge dnew: " << Charge_dn  << endl;
 
  ofstream ofp("scaler_charge.txt",ios::app);
  ofp<<Form("%d %.2lf %.2lf %.2lf",run,1.e-3*Charge_d3,1.e-3*Charge_un, 1.e-3*Charge_dn)<<endl;

}
