const double c=299792458e-9;// [m/ns]
const double mk=493.7e-3;// Kaon mass [GeV/c^2]
const double me=0.511e-3;// electron mass [GeV/c^2] 
const double ml=1115.7e-3;//Lambda mass [GeV/c^2]
const double mn=939.6e-3; // neutron mass [GeV/c^2]
const double mpi=139.6e-3;// pion mass [GeV/c^2]
const  int nth=3; //th num
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




//===================================================================//
//============================= Main ================================//
//===================================================================//

int main(int argc, char** argv){

  int ch; char* mode;
  string ifname = "/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/scripts/ita_scripts/run_list/Lambda_test.list";
  string ofname = "/pdf/hydro1_AC_eff_test.pdf";
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"h:f:w:n:bcop:GHT"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'w':
      output_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
  
    case 'G':
    mode="G";
      break;
  
    case 'H':
    mode="H";
      break;

    case 'T':
      mode="T";    
	break;

    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }


  TApplication *theApp =new TApplication("App",&argc,argv);
 if(draw_flag==0)gROOT->SetBatch(1);



  //=============== ROOT File Mode ================//
  /*if(ifname.c_str()=="../run_list/coin_H2_1.root")mode="G";
  else if(ifname.c_str()=="../run_list/Lambda_small.list" ||ifname.c_str()=="../run_list/Lambda_test.list")mode="H_1";
  else{cout<<"false to read mode types ";};
  */


 //double tdc_time=56.23e-3;//[ns]
 double tdc_time=58.e-3;//[ns] 

/*
 if(mode=="G"){tdc_time=56.23e-3;}
 else if(mode=="H"){tdc_time=56.23e-3;}
 else if(mode=="T"){tdc_time=58e-3;}
 */


 //TChain //
  TChain* T;
  if(mode=="G"){T=new TChain("tree"); }
  else {T=new TChain("T");}


  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());
    //  cout<<buf<<endl;
  }


  cout<<"mode :"<<mode<<endl;
  cout<<"tdc_time[ns]: "<<tdc_time<<endl;
  int evnt=T->GetEntries();
  cout<<"Get Entries: "<<evnt<<endl;


  double Ra1a_p[100],Ra2a_p[100],Ra1sum,Ra2sum;
 
 T->SetBranchStatus("*",0);
 T->SetBranchStatus("R.a1.a_p",1);
 T->SetBranchAddress("R.a1.a_p",Ra1a_p);
 T->SetBranchStatus("R.a2.a_p",1);
 T->SetBranchAddress("R.a2.a_p",Ra2a_p);
 T->SetBranchStatus("R.a1.asum_c",1);
 T->SetBranchAddress("R.a1.asum_c",&Ra1sum);
 T->SetBranchStatus("R.a2.asum_c",1);
 T->SetBranchAddress("R.a2.asum_c",&Ra2sum);
 TH1F *ha1_adc[24];
 TH1F *ha2_adc[26];
 TH1F *ha1_npe[24];
 TH1F *ha2_npe[26];
 TF1 *fconv_a1[24];
 TF1 *fconv_a2[26];
 
 double bin_ac1,min_ac1,max_ac1,bin_ac2,min_ac2,max_ac2;
 min_ac1=-100.0;
 max_ac1=10000.;
 min_ac2=-100;
 max_ac2=10000.;
 bin_ac1=(max_ac1-min_ac1);
 bin_ac1=(int)bin_ac1;
 bin_ac2=(max_ac2-min_ac2);
 bin_ac2=(int)bin_ac2;

 double bin_npe1,min_npe1,max_npe1,bin_npe2,min_npe2,max_npe2;
 min_npe1=-1.0;
 max_npe1=30.0;
 bin_npe1=max_npe1-min_npe1;
 bin_npe1=(int)bin_npe1;
 min_npe2=-1.0;
 max_npe2=30.0;
 bin_npe2=max_npe2-min_npe2;
 bin_npe2=(int)bin_npe2;

 TH1F *ha1_adc_sum=new TH1F("ha1_adc_sum","AC1 ADC SUM HIST",bin_ac1,min_ac1,max_ac1);
 TH1F *ha2_adc_sum=new TH1F("ha2_adc_sum","AC2 ADC SUM HIST",bin_ac2,min_ac2,max_ac2);
 TH1F *ha1_npe_sum=new TH1F("ha1_npe_sum","AC1 NPE SUM HIST",bin_npe1,min_npe1,max_npe1);
 TH1F *ha2_npe_sum=new TH1F("ha2_npe_sum","AC2 NPE SUM HIST",bin_npe2,min_npe2,max_npe2);
 
 for(int i=0;i<24;i++){
   ha1_adc[i]=new TH1F(Form("ha1_adc[%d]",i),"AC1 ADC HIST",bin_ac1,min_ac1,max_ac1);
   ha1_npe[i]=new TH1F(Form("ha1_npe[%d]",i),"AC1 NPE HIST",bin_npe1,min_npe1,max_npe1);
   fconv_a1[i]=new TF1(Form("fconv_a1[%d]",i),"gaus",min_ac1,max_ac1);
 }

for(int i=0;i<26;i++){
   ha2_adc[i]=new TH1F(Form("ha2_adc[%d]",i),"AC2 ADC HIST",bin_ac2,min_ac2,max_ac2);
   ha2_npe[i]=new TH1F(Form("ha2_npe[%d]",i),"AC2 NPE HIST",bin_npe2,min_npe2,max_npe2);   
   fconv_a2[i]=new TF1(Form("fconv_a2[%d]",i),"gaus",min_ac2,max_ac2);
 }
 
 
 for(int k=0;k<evnt;k++){
   T->GetEntry(k);
   for(int i=0;i<24;i++){
     if(Ra1a_p[i]>150.)ha1_adc[i]->Fill(Ra1a_p[i]);}
   for(int i=0;i<26;i++){
     if(Ra2a_p[i]>150.)ha2_adc[i]->Fill(Ra2a_p[i]);}
   ha1_adc_sum->Fill(Ra1sum);
   ha2_adc_sum->Fill(Ra2sum);
 }

 double pe1_a1[24],pe1_a2[26],n_a1[24],n_a2[26];


 for(int i=0;i<24;i++){
 pe1_a1[i]=ha1_adc[i]->GetBinContent(ha1_adc[i]->GetMaximumBin());
 //n_a1[i]=ha1_adc[i]->GetYaxis()->GetBinContent(ha1_adc[i]->GetMaximumBin());
 // fconv_a1->SetParameter(0,n_a1[i]);
 fconv_a1[i]->SetParameter(1,pe1_a1[i]);
 ha1_adc[i]->Fit(Form("fconv_a1[%d]",i),"Rq","",150,1000.);
 pe1_a1[i]=fconv_a1[i]->GetParameter(1);
 }


 for(int i=0;i<26;i++){
 pe1_a2[i]=ha2_adc[i]->GetBinContent(ha2_adc[i]->GetMaximumBin());
 //n_a2[i]=ha2_adc[i]->GetYaxis()->GetBinContent(ha2_adc[i]->GetMaximumBin());
 //fconv_a2->SetParameter(0,n_a2[i]);
 fconv_a2[i]->SetParameter(1,pe1_a2[i]);
 ha2_adc[i]->Fit(Form("fconv_a2[%d]",i),"Rq","",150.,1000.);
 pe1_a2[i]=fconv_a2[i]->GetParameter(1);
 }


for(int k=0;k<evnt;k++){
   T->GetEntry(k);
   for(int i=0;i<24;i++){
     ha1_npe[i]->Fill(Ra1a_p[i]/pe1_a1[i]);
     ha1_npe_sum->Fill(Ra1a_p[i]/pe1_a1[i]);}
   for(int i=0;i<26;i++){
     ha2_npe[i]->Fill(Ra2a_p[i]/pe1_a2[i]);
     ha2_npe_sum->Fill(Ra2a_p[i]/pe1_a2[i]);}
 }



 TCanvas* c0=new TCanvas("c0","AC1 ADC HIST");
 c0->Divide(4,6);
 TCanvas* c1=new TCanvas("c1","AC1 NPE HIST");
 c1->Divide(4,6);
 TCanvas* c2=new TCanvas("c2","AC2 ADC HIST");
 c2->Divide(4,7);
 TCanvas* c3=new TCanvas("c3","AC2 NPE HIST");
 c3->Divide(4,7);
 TCanvas* c4=new TCanvas("c4","AC SUM HIST");
 c4->Divide(2,2);
 
 for(int i=0;i<24;i++){
   c0->cd(i+1);
   ha1_adc[i]->Draw();
   fconv_a1[i]->Draw("same");
   c1->cd(i+1);
   ha1_npe[i]->Draw();
   
 }

 for(int i=0;i<26;i++){
   c2->cd(i+1);
   ha2_adc[i]->Draw();
   fconv_a2[i]->Draw("same");
   c3->cd(i+1);
   ha1_npe[i]->Draw();
 }

     c4->cd(1);
   ha1_adc_sum->Draw();
   c4->cd(2);
   ha2_adc_sum->Draw();
   c4->cd(3);
   ha1_npe_sum->Draw();
   c4->cd(4);
   ha2_npe_sum->Draw();

   
TString name;
 if(output_flag){
 
 name.Form(ofname.c_str());
 c0->Print(name+"[","pdf");
 c0->Print(name,"pdf");
 c1->Print(name,"pdf");
 c2->Print(name,"pdf");
 c3->Print(name,"pdf");
 c4->Print(name,"pdf");
 c4->Print(name+"]","pdf");


 }     
    
 cout<<"Print is done "<<endl;
   
 if(draw_flag==0)gSystem->Exit(1);
 theApp->Run();
 return 0;

		    
}
