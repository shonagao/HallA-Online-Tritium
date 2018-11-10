///////////////////////////////////////////////////
// Remake tree for S2 time walk correction(F1TDC)//
// by Y. Toyama Nov. 2018                        //
///////////////////////////////////////////////////

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

#include "Tree.h"
#include "ParamMan.h"
#include "define.h"

#define Calibration


const int NCanvas = 9;//num of canvas

int main(int argc, char** argv){
  string ifname = "list/run.txt";
  string ofname = "rootfiles/Rs2_seg0.root";
  string paramname = "param/f1_tuned_Lambda.param";
  int ch;
  int lr=0;
  int MaxNum = 1000000;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:p:"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input run list: "<<ifname<<endl;
      break;
    case 'w':
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'p':
      paramname = optarg;
      break;
    case 'h':
      cout<<"-f : input runlist filename"<<endl;
      cout<<"-w : output root filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-p : input param filename"<<endl;
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


  TApplication *theApp = new TApplication("App", &argc, argv);
  Tree *tr = new Tree();
  ParamMan *param = new ParamMan(paramname.c_str());
  if(param->SetVal())cout<<"initial parameter setted"<<endl;

  TFile *ofp = new TFile(ofname.c_str(),"recreate");
  TTree *tree_out = new TTree("tree","tree");

  
  double ctime;

  int s2lseg;
  double s2ltime,s2lat,s2lab;
  double s2ltt,s2ltb;
  double ltof;

  int s2rseg;
  double s2rtime,s2rat,s2rab;
  double s2rtt,s2rtb;
  double rtof;
  tree_out->Branch("ctime"   ,&ctime  ,"ctime/D"  );//coin time w/ pathl cor
  tree_out->Branch("s2lseg"  ,&s2lseg ,"s2lseg/I" );
  tree_out->Branch("s2ltime" ,&s2ltime,"s2ltime/D");
  tree_out->Branch("s2ltt"   ,&s2ltt  ,"s2ltt/D"  );
  tree_out->Branch("s2ltb"   ,&s2ltb  ,"s2ltb/D"  );
  tree_out->Branch("s2lat"   ,&s2lat  ,"s2lat/D"  );
  tree_out->Branch("s2lab"   ,&s2lab  ,"s2lab/D"  );
  tree_out->Branch("ltof"    ,&ltof   ,"ltof/D"   );
  tree_out->Branch("s2rseg"  ,&s2rseg ,"s2rseg/I" );
  tree_out->Branch("s2rtime" ,&s2rtime,"s2rtime/D");
  tree_out->Branch("s2rtt"   ,&s2rtt  ,"s2rtt/D"  );
  tree_out->Branch("s2rtb"   ,&s2rtb  ,"s2rtb/D"  );
  tree_out->Branch("s2rat"   ,&s2rat  ,"s2rat/D"  );
  tree_out->Branch("s2rab"   ,&s2rab  ,"s2rab/D"  );
  tree_out->Branch("rtof"    ,&rtof   ,"rtof/D"   );

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    tr->add_tree(runname);

  }
  tr->pack_tree();

  tr->readtreeTrackL();
  tr->readtreeTrackR();
  tr->readtreeF1TDCL();
  tr->readtreeF1TDCR();
  tr->readtreeS2L();
  tr->readtreeS2R();


  
  int EN=tr->ENum;
  int fill_ev=0;

  bool Rtr_flag, Ltr_flag;
  bool RZ_flag, LZ_flag, Z_flag;
  double beta_k[MAX], Vk[MAX];

  for(int n=0;n<EN;n++){
    if(n%10000==0)cout<<n <<" / "<<EN<<"  "<<fill_ev<<" filled"<<endl;
    if(fill_ev>MaxNum)break;

    Rtr_flag=Ltr_flag=RZ_flag=LZ_flag=false;

    tr->tree->GetEntry(n);
    if(abs(tr->R_evtype-5.)>0.5)continue;
    tr->convertF1TDCR(param);
    tr->convertF1TDCL(param);

    for(int rtr=0;rtr<tr->R_tr_n;rtr++){
      for(int ltr=0;ltr<tr->L_tr_n;ltr++){
        if(tr->R_tr_pathl[rtr] + tr->R_s2_trpath[rtr]>28.7 && tr->R_tr_pathl[rtr] + tr->R_s2_trpath[rtr]<29.4)Rtr_flag=true;
        if(tr->L_tr_pathl[ltr] + tr->L_s2_trpath[ltr]>28.6 && tr->L_tr_pathl[ltr] + tr->L_s2_trpath[ltr]<29.2)Ltr_flag=true;
        if(tr->R_tr_vz[rtr]>-0.1 && tr->R_tr_vz[rtr]< 0.1)RZ_flag=true;
        if(tr->L_tr_vz[ltr]>-0.1 && tr->L_tr_vz[ltr]< 0.1)LZ_flag=true;

        beta_k[rtr] = tr->R_tr_p[rtr]/sqrt(tr->R_tr_p[rtr]*tr->R_tr_p[rtr] + MK*MK);
        Vk[rtr] = beta_k[rtr]*LightVelocity;
        s2lseg = (int)tr->L_s2_trpad[ltr];
        s2rseg = (int)tr->R_s2_trpad[rtr];

        ltof = (tr->L_tr_pathl[ltr] + tr->L_s2_trpath[ltr])/LightVelocity;
        rtof = (tr->R_tr_pathl[rtr] + tr->R_s2_trpath[rtr])/Vk[rtr];

        double L_tgt = tr->LS2_F1time[s2lseg] - ltof;
        double R_tgt = tr->RS2_F1time[s2rseg] - rtof;
        ctime = L_tgt - R_tgt;

        s2ltime = tr->LS2_F1time[s2lseg];
        s2lat   = tr->L_s2_ra_p[s2lseg];
        s2lab   = tr->L_s2_la_p[s2lseg];
        s2ltt   = tr->LS2T_F1time[s2lseg];
        s2ltb   = tr->LS2B_F1time[s2lseg];

        s2rtime = tr->RS2_F1time[s2rseg];
        s2rat   = tr->R_s2_ra_p[s2rseg];
        s2rab   = tr->R_s2_la_p[s2rseg];
        s2rtt   = tr->RS2T_F1time[s2rseg];
        s2rtb   = tr->RS2B_F1time[s2rseg];

        // cout<<ctime<<endl;
        if(Rtr_flag && Ltr_flag && RZ_flag && LZ_flag && abs(ctime-3)<3.){ // 
          tree_out->Fill();
          fill_ev++;
        }
      }
    }

  
  }

  ofp->Write();
  ofp->Close();
  return 0;
}

