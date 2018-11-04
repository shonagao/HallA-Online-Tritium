///////////////////////
// calc. total charge//
///////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
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

int main(int argc, char** argv){
  string ifname = argv[1];


  TApplication *theApp = new TApplication("App", &argc, argv);

  TChain *E = new TChain("E");

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    E->Add(runname.c_str());
    cout<<buf<<endl;
  }


  int ENum = E->GetEntries();

  double time_interval=4.;//sec
  double bcm_av;
  double bcm_total=0.;
  E->SetBranchAddress("hac_bcm_average",&bcm_av);
  
  for(int n=0;n<ENum;n++){
    E->GetEntry(n);

    //cout<<bcm_av<<endl;
    bcm_total += bcm_av;
  }

  cout<<bcm_total<<endl;
  cout<<"==total charge=="<<endl;
  cout<<bcm_total * 1.e-3 *time_interval<<" [mC]"<<endl;

  gSystem->Exit(1);
  theApp->Run();
return 0;
}
