#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
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
//#include "TSignalHandler.h"

#include "TRandom.h"
#include "Tree.h"


////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "output0000.root";
  string ofname = "hoge.pdf";
  int ch;
  int MaxNum = 0;
  bool print_flag = false;
  bool draw_flag = true;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:pb"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'w':
      ofname = optarg;
      cout<<"onput filename : "<<ofname<<endl;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'p':
      print_flag = true;
//      pngname = optarg;
      cout<<"print flag = TRUE"<<endl;
      break;
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
    case 'h':
      cout<<"-f : input filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-p : print png file"<<endl;
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

  TApplication theApp("App", &argc, argv);

  time_t start, end;
  start = time(NULL);
  time(&start);

  if(!draw_flag) gROOT->SetBatch(1);
//  gROOT->SetBatch(1);
  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");

  gStyle->SetOptDate(0);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetOptDate(0);
  gStyle->SetOptStat("m");
  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.10);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  gStyle->SetGridWidth(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first

  TChain *tree = new TChain("tree");
  tree->Add(ifname.c_str());

////////////////////////////////////////////////////////////////////////
//////////////        hist
////////////////////////////////////////////////////////////////////
  TH2F *h2_adc_tof_gt, *h2_adc_tof_gt_dec;
  TH2F *h2_adc_tof_gb, *h2_adc_tof_gb_dec;
  TH2F *h2_adc_tof_ht, *h2_adc_tof_ht_dec;
  TH2F *h2_adc_tof_hb, *h2_adc_tof_hb_dec;
  TH1F *h1_tof, *h1_tof_dec;

  h2_adc_tof_gt     = new TH2F("h2_adc_tof_gt"    ,"h2_adc_tof_gt"    ,400,0,8000,350, 0,1400);
  h2_adc_tof_gb     = new TH2F("h2_adc_tof_gb"    ,"h2_adc_tof_gb"    ,400,0,8000,350, 0,1400);
  h2_adc_tof_ht     = new TH2F("h2_adc_tof_ht"    ,"h2_adc_tof_ht"    ,400,0,8000,350, 0,1400);
  h2_adc_tof_hb     = new TH2F("h2_adc_tof_hb"    ,"h2_adc_tof_hb"    ,400,0,8000,350, 0,1400);
  h2_adc_tof_gt_dec = new TH2F("h2_adc_tof_gt_dec","h2_adc_tof_gt_dec",400,0,8000,350, 0,1400);
  h2_adc_tof_gb_dec = new TH2F("h2_adc_tof_gb_dec","h2_adc_tof_gb_dec",400,0,8000,350, 0,1400);
  h2_adc_tof_ht_dec = new TH2F("h2_adc_tof_ht_dec","h2_adc_tof_ht_dec",400,0,8000,350, 0,1400);
  h2_adc_tof_hb_dec = new TH2F("h2_adc_tof_hb_dec","h2_adc_tof_hb_dec",400,0,8000,350, 0,1400);
  h2_adc_tof_gt    ->SetMinimum(1.); h2_adc_tof_gb    ->SetMinimum(1.);
  h2_adc_tof_ht    ->SetMinimum(1.); h2_adc_tof_hb    ->SetMinimum(1.);
  h2_adc_tof_gt_dec->SetMinimum(1.); h2_adc_tof_gb_dec->SetMinimum(1.);
  h2_adc_tof_ht_dec->SetMinimum(1.); h2_adc_tof_hb_dec->SetMinimum(1.);

  h1_tof     = new TH1F("h1_tof"    ,"h1_tof"    ,800,0,40);
  h1_tof_dec = new TH1F("h1_tof_dec","h1_tof_dec",800,0,40);

////////////////////////////////////////////////////////////////////////
//////////////        Fill
////////////////////////////////////////////////////////////////////
  int ENum;
  ENum = tree->GetEntries();
  if(MaxNum>0 && MaxNum<ENum) ENum = MaxNum;
//  int acadct, acadcb, ac2adct, ac2adcb;
  int evID, runID;
  double dec_gadct, dec_gtdct, dec_gadcb, dec_gtdcb;
  double dec_iadct, dec_itdct, dec_iadcb, dec_itdcb;
  double dec_hadct, dec_htdct, dec_hadcb, dec_htdcb;
  double tmp_gt[5], tmp_gb[5];
  double tmp_ht[5], tmp_hb[5];
  double chi2_xz;

  tree->SetBranchAddress("evID",     &evID);
  tree->SetBranchAddress("runID",    &runID);
  tree->SetBranchAddress("d_gadct",  &dec_gadct);
  tree->SetBranchAddress("d_gadcb",  &dec_gadcb);
  tree->SetBranchAddress("d_iadct",  &dec_iadct);
  tree->SetBranchAddress("d_iadcb",  &dec_iadcb);
  tree->SetBranchAddress("d_hadct",  &dec_hadct);
  tree->SetBranchAddress("d_hadcb",  &dec_hadcb);
  tree->SetBranchAddress("d_gtdct",  &dec_gtdct);
  tree->SetBranchAddress("d_gtdcb",  &dec_gtdcb);
  tree->SetBranchAddress("d_itdct",  &dec_itdct);
  tree->SetBranchAddress("d_itdcb",  &dec_itdcb);
  tree->SetBranchAddress("d_htdct",  &dec_htdct);
  tree->SetBranchAddress("d_htdcb",  &dec_htdcb);
  tree->SetBranchAddress("chi2_xz",  &chi2_xz);

// timewalk
  double p_ht[5][10]  = {{  },
                         {  },
                         {  },
                         {  },
                         {  }};
  double p_hb[5][10]  = {{  },
                         {  },
                         {  },
                         {  },
                         {  }};
  double p_gt[5][10]  = {{  },
                         {  },
                         {  },
                         {  },
                         {  }};
  double p_gb[5][10]  = {{  },
                         {  },
                         {  },
                         {  },
                         {  }};

  for(int i=0;i<8;i++){
    p_ht[0][i]  = -1.5*p_ht[0][i]; p_ht[1][i] = -1.5*p_ht[1][i]; p_ht[2][i] = -1.5*p_ht[2][i]; p_ht[3][i] = -1.5*p_ht[3][i]; p_ht[4][i] = -1.5*p_ht[4][i];
    p_hb[0][i]  = -1.8*p_hb[0][i]; p_hb[1][i] = -1.5*p_hb[1][i]; p_hb[2][i] = -1.5*p_hb[2][i]; p_hb[3][i] = -1.5*p_hb[3][i]; p_hb[4][i] = -1.5*p_hb[4][i];
    p_gt[0][i]  = +1.0*p_gt[0][i]; p_gt[1][i] = +1.0*p_gt[1][i]; p_gt[2][i] = +1.0*p_gt[2][i]; p_gt[3][i] = +1.0*p_gt[3][i]; p_gt[4][i] = +1.0*p_gt[4][i];
    p_gb[0][i]  = +1.0*p_gb[0][i]; p_gb[1][i] = +1.0*p_gb[1][i]; p_gb[2][i] = +1.0*p_gb[2][i]; p_gb[3][i] = +1.0*p_gb[3][i]; p_gb[4][i] = +1.0*p_gb[4][i];
  }

  for(int i=0;i<5;i++){
    tmp_gt[i] = - (p_gt[i][0] + p_gt[i][1]*pow(p_gt[i][9],1) + p_gt[i][2]*pow(p_gt[i][9],2) + p_gt[i][3]*pow(p_gt[i][9],3)
                           + p_gt[i][4]*pow(p_gt[i][9],4) + p_gt[i][5]*pow(p_gt[i][9],5) + p_gt[i][6]*pow(p_gt[i][9],6)
                           + p_gt[i][7]*pow(p_gt[i][9],7) );
    tmp_gb[i] = - (p_gb[i][0] + p_gb[i][1]*pow(p_gb[i][9],1) + p_gb[i][2]*pow(p_gb[i][9],2) + p_gb[i][3]*pow(p_gb[i][9],3)
                           + p_gb[i][4]*pow(p_gb[i][9],4) + p_gb[i][5]*pow(p_gb[i][9],5) + p_gb[i][6]*pow(p_gb[i][9],6)
                           + p_gb[i][7]*pow(p_gb[i][9],7) );

    tmp_ht[i] = - (p_ht[i][0] + p_ht[i][1]*pow(p_ht[i][9],1) + p_ht[i][2]*pow(p_ht[i][9],2) + p_ht[i][3]*pow(p_ht[i][9],3)
                           + p_ht[i][4]*pow(p_ht[i][9],4) + p_ht[i][5]*pow(p_ht[i][9],5) + p_ht[i][6]*pow(p_ht[i][9],6)
                           + p_ht[i][7]*pow(p_ht[i][9],7) );

    tmp_hb[i] = - (p_hb[i][0] + p_hb[i][1]*pow(p_hb[i][9],1) + p_hb[i][2]*pow(p_hb[i][9],2) + p_hb[i][3]*pow(p_hb[i][9],3)
                              + p_hb[i][4]*pow(p_hb[i][9],4) + p_hb[i][5]*pow(p_hb[i][9],5) + p_hb[i][6]*pow(p_hb[i][9],6)
                              + p_hb[i][7]*pow(p_hb[i][9],7) );

  }

  double tof =0.;
  double tdc_gt, tdc_gb, tdc_ht, tdc_hb;
  double adc_gt, adc_gb, adc_ht, adc_hb;

  for(int n=0;n<ENum;n++){
    tree->GetEntry(n);

    tdc_gt = dec_gtdct;  tdc_gb = dec_gtdcb; tdc_ht = dec_htdct;  tdc_hb = dec_htdcb;
    adc_gt = dec_gadct;  adc_gb = dec_gadcb; adc_ht = dec_hadct;  adc_hb = dec_hadcb;

    tof = ((tdc_gt+tdc_gb)/2. - (tdc_ht+tdc_hb)/2.);

    if(   fabs(tdc_gt-4000)<4000 && fabs(tdc_gb-4000)<4000
       && fabs(tdc_ht-4000)<4000 && fabs(tdc_hb-4000)<4000){
      h2_adc_tof_gt->Fill(adc_gt,tof);
      h2_adc_tof_gb->Fill(adc_gb,tof);
      h2_adc_tof_ht->Fill(adc_ht,tof);
      h2_adc_tof_hb->Fill(adc_hb,tof);
      h1_tof->Fill(tof*0.025); 
    }

    for(int i=0;i<5;i++){
      if(adc_gt>p_gt[i][8] && adc_gt<p_gt[i][9]){
        tdc_gt = tdc_gt - (p_gt[i][0] + p_gt[i][1]*pow(adc_gt,1) + p_gt[i][2]*pow(adc_gt,2) + p_gt[i][3]*pow(adc_gt,3)
                                      + p_gt[i][4]*pow(adc_gt,4) + p_gt[i][5]*pow(adc_gt,5) + p_gt[i][6]*pow(adc_gt,6)
                                      + p_gt[i][7]*pow(adc_gt,7) );
        tdc_gt = tdc_gt - tmp_gt[i];
      }

      if(adc_gb>p_gb[i][8] && adc_gb<p_gb[i][9]){
        tdc_gb = tdc_gb - (p_gb[i][0] + p_gb[i][1]*pow(adc_gb,1) + p_gb[i][2]*pow(adc_gb,2) + p_gb[i][3]*pow(adc_gb,3)
                                      + p_gb[i][4]*pow(adc_gb,4) + p_gb[i][5]*pow(adc_gb,5) + p_gb[i][6]*pow(adc_gb,6)
                                      + p_gb[i][7]*pow(adc_gb,7) );
        tdc_gb = tdc_gb - tmp_gb[i];
      }

      if(adc_ht>p_ht[i][8] && adc_ht<p_ht[i][9]){
        tdc_ht = tdc_ht - (p_ht[i][0] + p_ht[i][1]*pow(adc_ht,1) + p_ht[i][2]*pow(adc_ht,2) + p_ht[i][3]*pow(adc_ht,3)
                                      + p_ht[i][4]*pow(adc_ht,4) + p_ht[i][5]*pow(adc_ht,5) + p_ht[i][6]*pow(adc_ht,6)
                                      + p_ht[i][7]*pow(adc_ht,7) );
        tdc_ht = tdc_ht - tmp_ht[i];
      }

      if(adc_hb>p_hb[i][8] && adc_hb<p_hb[i][9]){
        tdc_hb = tdc_hb - (p_hb[i][0] + p_hb[i][1]*pow(adc_hb,1) + p_hb[i][2]*pow(adc_hb,2) + p_hb[i][3]*pow(adc_hb,3)
                                      + p_hb[i][4]*pow(adc_hb,4) + p_hb[i][5]*pow(adc_hb,5) + p_hb[i][6]*pow(adc_hb,6)
                                      + p_hb[i][7]*pow(adc_hb,7) );
        tdc_hb = tdc_hb - tmp_hb[i];
      }

    }

    //tofcut
    if(   fabs(tdc_gt-2000)<1500 && fabs(tdc_gb-2000)<1500
       && fabs(tdc_ht-2000)<1500 && fabs(tdc_hb-2000)<1500
       && tof>0 && tof<370
//       && runID<11061700
//       && adc_gt>1000 && adc_gt<1200
//       && adc_gb>800 && adc_gb<1000
//       && adc_ht>1000 && adc_ht<1800
//       && adc_hb>1000 && adc_hb<8000
       && chi2_xz<2.0
      ){
      tof = ((tdc_gt+tdc_gb)/2. - (tdc_ht+tdc_hb)/2.);
      h2_adc_tof_gt_dec->Fill(adc_gt,tof);
      h2_adc_tof_gb_dec->Fill(adc_gb,tof);
      h2_adc_tof_ht_dec->Fill(adc_ht,tof);
      h2_adc_tof_hb_dec->Fill(adc_hb,tof);
      h1_tof_dec->Fill(tof*0.025);
    }

    if(n%100000==0){
     char clout[100];
     end = time(NULL);
     time(&end);
     sprintf(clout,"%.0f sec",difftime(end,start));
     cout<<n<<" / "<<ENum<<" : "<<clout<<endl;
    } // if n%100000==0 std::cout

  } // for ENum

//  double mean_gt[500], mean_gb[500], mean_it[500], mean_ib[500];
//  TF1 *f;
//  TH1F *h1;
//  int cut_bin = 
//  for(int i=0;i<50;i++){
//    double max = h2_time_de->ProjectionX("_px",0,);
//  }
//  h2_time_de->ProjectionX("_px",0,cut_bin)->Fit(f,"0QR+","",p_fit[1]-p_fit[2]*1.0,p_fit[1]+p_fit[2]*1.0);

  TObjArray aSlices;
  double mean_gt = h2_adc_tof_gt_dec->ProjectionY()->GetBinCenter(h2_adc_tof_gt_dec->ProjectionY()->GetMaximumBin());
  double mean_gb = h2_adc_tof_gb_dec->ProjectionY()->GetBinCenter(h2_adc_tof_gb_dec->ProjectionY()->GetMaximumBin());
  double mean_ht = h2_adc_tof_ht_dec->ProjectionY()->GetBinCenter(h2_adc_tof_ht_dec->ProjectionY()->GetMaximumBin());
  double mean_hb = h2_adc_tof_hb_dec->ProjectionY()->GetBinCenter(h2_adc_tof_hb_dec->ProjectionY()->GetMaximumBin());
  TF1 *fit_gt = new TF1("fit_gt","gaus",mean_gt-30,mean_gt+20);
  TF1 *fit_gb = new TF1("fit_gb","gaus",mean_gb-30,mean_gb+30);
  TF1 *fit_ht = new TF1("fit_ht","gaus",mean_ht-30,mean_ht+30);
  TF1 *fit_hb = new TF1("fit_hb","gaus",mean_hb-30,mean_hb+30);
  h2_adc_tof_gt_dec->FitSlicesY(fit_gt,0,500,0,"QRG3");
  h2_adc_tof_gb_dec->FitSlicesY(fit_gb,0,500,0,"QRG3");
  h2_adc_tof_ht_dec->FitSlicesY(fit_ht,0,500,0,"QRG3");
  h2_adc_tof_hb_dec->FitSlicesY(fit_hb,0,500,0,"QRG3");
  TH1D *h_mean_gt = (TH1D*)gROOT->FindObject("h2_adc_tof_gt_dec_1");
  TH1D *h_mean_gb = (TH1D*)gROOT->FindObject("h2_adc_tof_gb_dec_1");
  TH1D *h_mean_ht = (TH1D*)gROOT->FindObject("h2_adc_tof_ht_dec_1");
  TH1D *h_mean_hb = (TH1D*)gROOT->FindObject("h2_adc_tof_hb_dec_1");
  h_mean_gt->SetLineColor(1); h_mean_gt->SetLineWidth(2);
  h_mean_gb->SetLineColor(1); h_mean_gb->SetLineWidth(2);
  h_mean_ht->SetLineColor(1); h_mean_ht->SetLineWidth(2);
  h_mean_hb->SetLineColor(1); h_mean_hb->SetLineWidth(2);

  double min_gt = 400.,  max_gt = 1800;
  double min_gb = 400.,  max_gb = 1800;
  double min_ht = 400.,  max_ht = 2000;
  double min_hb = 400.,  max_hb = 2500;
  TF1 *fit_mean_gt = new TF1("fit_mean_gt","pol7");
  TF1 *fit_mean_gb = new TF1("fit_mean_gb","pol7");
  TF1 *fit_mean_ht = new TF1("fit_mean_ht","pol7");
  TF1 *fit_mean_hb = new TF1("fit_mean_hb","pol7");
  fit_mean_gt->SetLineWidth(1);
  fit_mean_gb->SetLineWidth(1);
  fit_mean_ht->SetLineWidth(1);
  fit_mean_hb->SetLineWidth(1);
  h_mean_gt->Fit(fit_mean_gt,"QR+","",min_gt,max_gt);
  h_mean_gb->Fit(fit_mean_gb,"QR+","",min_gb,max_gb);
  h_mean_ht->Fit(fit_mean_ht,"QR+","",min_ht,max_ht);
  h_mean_hb->Fit(fit_mean_hb,"QR+","",min_hb,max_hb);
  double p_mean_gt[10], p_mean_gb[10];
  double p_mean_ht[10], p_mean_hb[10];
  fit_mean_gt->GetParameters(p_mean_gt);
  fit_mean_gb->GetParameters(p_mean_gb);
  fit_mean_ht->GetParameters(p_mean_ht);
  fit_mean_hb->GetParameters(p_mean_hb);
  cout<<"fitting result HT"<<endl;
  for(int i=0;i<8;i++){ cout<<Form("%.6e, ",p_mean_ht[i])<<flush; }
  cout<<Form("%.4e, %.4e",min_ht,max_ht-100)<<endl;
  cout<<"fitting result HB"<<endl;
  for(int i=0;i<8;i++){ cout<<Form("%.6e, ",p_mean_hb[i])<<flush; }
  cout<<Form("%.4e, %.4e",min_hb,max_hb-100)<<endl;
  cout<<"fitting result GT"<<endl;
  for(int i=0;i<8;i++){ cout<<Form("%.6e, ",p_mean_gt[i])<<flush; }
  cout<<Form("%.4e, %.4e",min_gt,max_gt-100)<<endl;
  cout<<"fitting result GB"<<endl;
  for(int i=0;i<8;i++){ cout<<Form("%.6e, ",p_mean_gb[i])<<flush; }
  cout<<Form("%.4e, %.4e",min_gb,max_gb-100)<<endl;
  cout<<endl;

  TF1 *fit, *fit_dec;
  double fit_min, fit_max;
  fit_min = h1_tof->GetBinCenter(h1_tof->GetMaximumBin()) - 1.;
  fit_max = h1_tof->GetBinCenter(h1_tof->GetMaximumBin()) + 1.;
  double sigma1, sigma2;
  double mean1, mean2;
  for(int i=0;i<10;i++){
    fit = new TF1("fit","gaus");
    h1_tof->Fit(fit,"0QR","",fit_min,fit_max);
    fit_min = fit->GetParameter(1) - fit->GetParameter(2)*1.0;
    fit_max = fit->GetParameter(1) + fit->GetParameter(2)*1.0;
    fit->Clear();
  }
  fit = new TF1("fit","gaus");
  fit->SetLineColor(1);
  h1_tof->Fit(fit,"QR","",fit_min,fit_max);
  sigma1 = fit->GetParameter(2);
  mean1 = fit->GetParameter(1);

  fit_min = 10.; fit_max = 30.;
  fit_min = h1_tof_dec->GetBinCenter(h1_tof_dec->GetMaximumBin()) - 1.;
  fit_max = h1_tof_dec->GetBinCenter(h1_tof_dec->GetMaximumBin()) + 1.;
  for(int i=0;i<10;i++){
    fit_dec = new TF1("fit_dec","gaus");
    h1_tof_dec->Fit(fit_dec,"0QR","",fit_min,fit_max);
    fit_min = fit_dec->GetParameter(1) - fit_dec->GetParameter(2)*1.0; 
    fit_max = fit_dec->GetParameter(1) + fit_dec->GetParameter(2)*1.0;
    fit_dec->Clear();
  }
  fit_dec = new TF1("fit_dec","gaus");
  fit_dec->SetLineColor(1);
  h1_tof_dec->Fit(fit_dec,"QR","",fit_min,fit_max);
  sigma2 = fit_dec->GetParameter(2);
  mean2 = fit_dec->GetParameter(1);

  TText *txt1, *txt2;
  char str_sigma[100];
  double max_hist;
  sprintf(str_sigma,"sigma = %.3f ns",sigma1);
  max_hist = h1_tof->GetMaximum();
  txt1 = new TLatex(28.,max_hist*0.98,str_sigma);
  txt1->SetTextSize(0.04);
  txt1->SetTextColor(4);
  sprintf(str_sigma,"sigma = %.3f ns",sigma2);
  max_hist = h1_tof_dec->GetMaximum();
  txt2 = new TLatex(28.,max_hist*0.98,str_sigma);
  txt2->SetTextSize(0.04);
  txt2->SetTextColor(4);

  char str_mean1[100], str_mean2[100], str_meandiff[100];
  sprintf(str_mean1,"%.4e",mean1/0.025);
  sprintf(str_mean2,"%.4e",mean2/0.025);
//  sprintf(str_meandiff,"%.4e",mean1/0.025-mean2/0.025); // Hwall
  sprintf(str_meandiff,"%.4e",mean2/0.025-mean1/0.025); // Gwall
  cout<<"TOF diff : "<<str_mean1<<" - "<<str_mean2<<" = "<<str_meandiff<<endl;

  char str[10][100];
  cout<<"time walk parameters"<<endl;

  for(int j=0;j<5;j++){
    cout<<Form("H TOP %d",j)<<endl;
    for(int i=0;i<10;i++){
      if(i==0) sprintf(str[i],"%.6e",p_ht[j][i]+tmp_ht[j]);
      else sprintf(str[i],"%.6e",p_ht[j][i]);
    }
    cout<<str[0]<<", "<<str[1]<<", "<<str[2]<<", "<<str[3]<<", "<<str[4]<<", "
        <<str[5]<<", "<<str[6]<<", "<<str[7]<<", "<<str[8]<<", "<<str[9]<<endl;
  }
  cout<<endl;

  for(int j=0;j<5;j++){
    cout<<Form("H BOT %d",j)<<endl;
    for(int i=0;i<10;i++){
      if(i==0) sprintf(str[i],"%.6e",p_hb[j][i]+tmp_hb[j]);
      else sprintf(str[i],"%.6e",p_hb[j][i]);
    }
    cout<<str[0]<<", "<<str[1]<<", "<<str[2]<<", "<<str[3]<<", "<<str[4]<<", "
        <<str[5]<<", "<<str[6]<<", "<<str[7]<<", "<<str[8]<<", "<<str[9]<<endl;
  }
  cout<<endl;

  for(int j=0;j<5;j++){
    cout<<Form("G TOP %d",j)<<endl;
    for(int i=0;i<10;i++){
      if(i==0) sprintf(str[i],"%.6e",p_gt[j][i]+tmp_gt[j]);
      else sprintf(str[i],"%.6e",p_gt[j][i]);
    }
    cout<<str[0]<<", "<<str[1]<<", "<<str[2]<<", "<<str[3]<<", "<<str[4]<<", "
        <<str[5]<<", "<<str[6]<<", "<<str[7]<<", "<<str[8]<<", "<<str[9]<<endl;
  }
  cout<<endl;

  for(int j=0;j<5;j++){
    cout<<Form("G BOT %d",j)<<endl;
    for(int i=0;i<10;i++){
      if(i==0) sprintf(str[i],"%.6e",p_gb[j][i]+tmp_gb[j]);
      else sprintf(str[i],"%.6e",p_gb[j][i]);
    }
    cout<<str[0]<<", "<<str[1]<<", "<<str[2]<<", "<<str[3]<<", "<<str[4]<<", "
        <<str[5]<<", "<<str[6]<<", "<<str[7]<<", "<<str[8]<<", "<<str[9]<<endl;
  }
  cout<<endl;

////////////////////////////////////////////////////////////////////////
///////////////////     draw
///////////////////////////////////////////////////////////////////////


  TCanvas *c1 = new TCanvas("c1","c1",1700,850);
  TCanvas *c2 = new TCanvas("c2","c2",1700,850);
  TCanvas *c3 = new TCanvas("c3","c3",1700,850);

  c1->Divide(2,2);
  c1->cd(1); gPad->SetLogz(); h2_adc_tof_gt->Draw("colz");
  c1->cd(3); gPad->SetLogz(); h2_adc_tof_gb->Draw("colz");
  c1->cd(2); gPad->SetLogz(); h2_adc_tof_gt_dec->Draw("colz"); h_mean_gt->Draw("same");
  c1->cd(4); gPad->SetLogz(); h2_adc_tof_gb_dec->Draw("colz"); h_mean_gb->Draw("same");

  c2->Divide(2,2);
  c2->cd(1); gPad->SetLogz(); h2_adc_tof_ht->Draw("colz");
  c2->cd(3); gPad->SetLogz(); h2_adc_tof_hb->Draw("colz");
  c2->cd(2); gPad->SetLogz(); h2_adc_tof_ht_dec->Draw("colz"); h_mean_ht->Draw("same");
  c2->cd(4); gPad->SetLogz(); h2_adc_tof_hb_dec->Draw("colz"); h_mean_hb->Draw("same");

  c3->Divide(2,1);
  c3->cd(1);
  h1_tof->Draw(); txt1->Draw("same");
  c3->cd(2);
  h1_tof_dec->Draw(); txt2->Draw("same");

  if(print_flag==1){
    c1->Print("c1.png");
//    c4->Print("c4.png");
//    c5->Print("c5.png");
    gSystem->Exit(1);
    exit(1);
  }

  cout<<"fihish draw!!"<<endl;

  if(draw_flag==0) gSystem->Exit(1);
//  exit(1);

  theApp.Run();
  return 0;

} 






