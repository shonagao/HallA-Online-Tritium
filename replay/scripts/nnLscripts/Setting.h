#ifndef Setting_h
#define Setting_h 1

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <algorithm>

#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TArrow.h"
#include "TDatime.h"
#include "TMarker.h"

class Setting
{
public:
  Setting();
//  ~Setting();
  virtual ~Setting(){};

private:

public:
  void SetTH1(TH1 *h1, TString hname, TString xname, TString yname, int LColor=1, int FStyle=0, int FColor=0);
  void SetTH2(TH2 *h, TString name, TString xname, TString yname, double min=0.8);
  void SetGr(TGraph *gr, TString name, TString xname, TString yname, int LColor=1, int LWidth=1, int LStyle=1, int MColor=1, int MStyle=20, double MSize=1.);
  void SetGr(TGraphErrors *gr, TString name, TString xname, TString yname, int LColor=1, int LWidth=1, int LStyle=1, int MColor=1, int MStyle=20, double MSize=1.);
  void SetGr(TGraphAsymmErrors *gr, TString name, TString xname, TString yname, int LColor=1, int LWidth=1, int LStyle=1, int MColor=1, int MStyle=20, double MSize=1.);
  void SetGr(TGraph2D *gr, TString name, TString xname, TString yname, TString zname, int LColor=1, int LWidth=1, int LStyle=1, int MColor=1, int MStyle=20, double MSize=1.);
  void SetTF1(TF1 *f, int LColor=2, int LWidth=0, int LStyle=1, int Npx=1000);
  void BinLogX(TH2 *h);
  void BinLogY(TH2 *h);
  double landaugaus(double *x, double *par);
  void Initialize();

};

#endif
