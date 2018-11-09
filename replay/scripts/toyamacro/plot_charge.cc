//////////////////////////////////
// plot beam charge for each run//
//////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
#include <vector>
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
#include "TMath.h"
#include "TGaxis.h"
#include "TRandom.h"

#include "Setting.h"

int main(int argc, char** argv){
  string ifname = argv[1];


  TApplication *theApp = new TApplication("App", &argc, argv);
  Setting *set = new Setting();


  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  int runnum;
  double Charge_d3, Charge_un, Charge_dn, Charge_av;
  double d3_total, un_total, dn_total, av_total;

  d3_total = un_total = dn_total = av_total = 0.;
  TGraph *tg_charge = new TGraph();

  int n=0;

  while(1){

    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runnum;
    sbuf >> Charge_d3;
    sbuf >> Charge_un;
    sbuf >> Charge_dn;
    Charge_av = (Charge_d3 + Charge_un + Charge_dn)/3.;
    //cout<<buf<<endl;

    tg_charge -> SetPoint(n,runnum,Charge_av);
    d3_total += Charge_d3;
    un_total += Charge_un;
    dn_total += Charge_dn;
    av_total += Charge_av;
    n++;
  }


  cout<<"----total charge----"<<endl;
  cout<<   av_total   <<" mC(average)"<<endl;
  TCanvas *c1 = new TCanvas("c1","c1",1800,900);
  set->SetGr(tg_charge, "","Run number", "charge [mC]",1,4,22);
  tg_charge->SetMarkerSize(0.8);
  tg_charge -> Draw("APL");
  //gSystem->Exit(1);
  theApp->Run();
return 0;
}
