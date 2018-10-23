using namespace std;

#include "ana.h"

bool batch = false;
bool RHRS = false;
bool LHRS = false;

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
ana::ana()
{
  set = new Setting();
  set->Initialize();

  tr = new Tree();

}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
ana::~ana(){
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::ReadTree(string name){
  tr->merge(name);
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
bool ana::GetEntry(int n){
  tr->tree->GetEntry(n);

  return true;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::Roop(){
  tr->readtree_COMN();
  if(RHRS)
    tr->readtree_RHRS();
  if(LHRS)
    tr->readtree_LHRS();

  ENum = tr->tree->GetEntries();
  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();

// L.tr.chi2<0.01 && L.tr.th<0.17*L.tr.x+0.025 && L.tr.th>0.17*L.tr.x-0.035 && L.tr.th<0.4*L.tr.x+0.13 (LHRS plane cut)

  for(int n=0;n<ENum;n++){
    GetEntry(n);

    if(n % 10000 == 0){
      cout<<n<<" / "<<ENum<<endl;
    }
  }


  delete tr;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::Draw(){

//  TCanvas *c1 = new TCanvas("c1","c1",1500,900);

}


/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::DrawLine(double xmin, double ymin, double xmax, double ymax, int color, int style){
  line[LineID] = new TLine(xmin,ymin,xmax,ymax);
  line[LineID]->SetLineColor(color);
  line[LineID]->SetLineStyle(style);
  line[LineID]->Draw();
  LineID++;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::DrawText(double x, double y, string str, double size, int color){
  text[TextID] = new TLatex(x,y,str.c_str());
  text[TextID]->SetTextFont(42);
  text[TextID]->SetTextAlign(11);
  text[TextID]->SetTextColor(color);
  text[TextID]->SetTextSize(size);
  text[TextID]->Draw();
  TextID++;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
void ana::MakeHist(){
// Hist name is defined by "h + LorR + variable" for TH1.
//                         "h + LorR + variableY + variableX" for TH2.
/////////////
//// BPM ////
/////////////
  h_rbay_rbax = new TH2D("h_rbay_rbax","h_rbay_rbax",200,-20,20,200,-20,20);
  h_rbby_rbbx = new TH2D("h_rbby_rbbx","h_rbby_rbbx",200,-20,20,200,-20,20);
  h_rby_rbx   = new TH2D("h_rby_rbx"  ,"h_rby_rbx"  ,200,-20,20,200,-20,20);
  set->SetTH2(h_rbay_rbax,"BPM A"         ,"X","Y");
  set->SetTH2(h_rbby_rbbx,"BPM B"         ,"X","Y");
  set->SetTH2(h_rby_rbx  ,"Raster Pattern","X","Y");

//////////////
//// RHRS ////
//////////////
  h_R_trig = new TH1D("h_R_trig","h_R_trig",10,0,10);
  set->SetTH1(h_R_trig,"Trigger Flag","Trig No.","Counts");

  h_R_tr_n   = new TH1D("h_R_tr_n"  ,"h_R_tr_n"  ,30 , 0    , 30);
  h_R_tr_ch2 = new TH1D("h_R_tr_ch2","h_R_tr_ch2",400, 0    , 10);
  h_R_p      = new TH1D("h_R_p"     ,"h_R_p"     ,400, 1    , 3);
  h_R_pathl  = new TH1D("h_R_pathl" ,"h_R_pathl" ,400, 20   , 30);
  h_R_px     = new TH1D("h_R_px"    ,"h_R_px"    ,400, 0    , 2);
  h_R_py     = new TH1D("h_R_py"    ,"h_R_px"    ,400, -1   , 1);
  h_R_pz     = new TH1D("h_R_pz"    ,"h_R_px"    ,400, 1    , 3);
  h_R_tgy    = new TH1D("h_R_tgy"   ,"h_R_tgy"   ,400, -1   , 1);
  h_R_tgth   = new TH1D("h_R_tgth"  ,"h_R_tgth"  ,400, -0.2 , 0.2);
  h_R_tgph   = new TH1D("h_R_tgph"  ,"h_R_tgph"  ,400, -0.2 , 0.2);
  h_R_vx     = new TH1D("h_R_vx"    ,"h_R_vx"    ,400, -0.01, 0.01);
  h_R_vy     = new TH1D("h_R_vy"    ,"h_R_vy"    ,400, -0.01, 0.01);
  h_R_vz     = new TH1D("h_R_vz"    ,"h_R_vz"    ,400, -1   , 1);
  h_R_y_x    = new TH2D("h_R_y_x"   ,"h_R_y_x"   ,200, -1   , 1  , 200, -0.2, 0.2);
  h_R_th_x   = new TH2D("h_R_th_x"  ,"h_R_th_x"  ,200, -1   , 1  , 200, -0.3, 0.3);
  h_R_ph_y   = new TH2D("h_R_ph_y"  ,"h_R_ph_y"  ,200, -0.2 , 0.2, 200, -0.2, 0.2);
  set->SetTH1(h_R_tr_n,"No. of Tracks","No. of Tracks","Counts");
  set->SetTH1(h_R_tr_ch2,"Tracking #chi^{2}","#chi^{2}","Counts");
  set->SetTH1(h_R_p     ,"Track Momentum","p (GeV/#it{c})","Counts");
  set->SetTH1(h_R_pathl ,"Track Path Length","l (m)","Counts");
  set->SetTH1(h_R_px    ,"Momentum X","px (GeV/#it{c})","Counts");
  set->SetTH1(h_R_py    ,"Momentum Y","py (GeV/#it{c})","Counts");
  set->SetTH1(h_R_pz    ,"Momentum Z","pz (GeV/#it{c})","Counts");
  set->SetTH1(h_R_tgy   ,"Target Plane Y","y_{t} (m)","Counts");
  set->SetTH1(h_R_tgth  ,"Target Plane #theta","#theta_{t} (rad)","Counts");
  set->SetTH1(h_R_tgph  ,"Target Plane #phi","#phi_{t} (rad)","Counts");
  set->SetTH1(h_R_vx    ,"Vertex X","x_{v} (m)","Counts");
  set->SetTH1(h_R_vy    ,"Vertex Y","y_{v} (m)","Counts");
  set->SetTH1(h_R_vz    ,"Vertex Z","z_{v} (m)","Counts");
  set->SetTH2(h_R_y_x   ,"Focal Plane Y v.s X","X (m)","Y (m)");
  set->SetTH2(h_R_th_x  ,"Focal Plane #theta v.s X","X (m)","#theta (rad)");
  set->SetTH2(h_R_ph_y  ,"Focal Plane #phi v.s Y","Y (m)","#phi (rad)");

  h_R_beta       = new TH1D("h_R_beta"      ,"h_R_beta"      ,400,0   ,2);
  h_R_m2         = new TH1D("h_R_m2"        ,"h_R_m2"        ,400,-0.5,2.5);
  h_R_beta_p     = new TH2D("h_R_beta_p"    ,"h_R_beta_p"    ,200,1,3 ,200,0,2);
  h_R_beta_m2    = new TH2D("h_R_beta_m2"   ,"h_R_beta_m2"   ,200,1,3 ,200,-0.5,2.5);
  h_R_dedx_p     = new TH2D("h_R_dedx_p"    ,"h_R_dedx_p"    ,200,0,10,200,0,2);
  h_R_dedx_m2    = new TH2D("h_R_dedx_m2"   ,"h_R_dedx_m2"   ,200,0,10,200,-0.5,2.5);
  h_R_s0_dedx    = new TH1D("h_R_s0_dedx"   ,"h_R_s0_dedx"   ,400,0,10);
  h_R_s0_beta_x  = new TH2D("h_R_s0_beta_x" ,"h_R_s0_beta_x" ,200,-1,1,200,0,2);
  h_R_s0_dedx_x  = new TH2D("h_R_s0_dedx_x" ,"h_R_s0_dedx_x" ,200,-1,1,200,0,10);
  h_R_s2_dedx    = new TH1D("h_R_s2_dedx"   ,"h_R_s2_dedx"   ,400,0,10);
  h_R_s2_beta_x  = new TH2D("h_R_s2_beta_x" ,"h_R_s2_beta_x" ,200,-1,1,200,0, 2);
  h_R_s2_dedx_x  = new TH2D("h_R_s2_dedx_x" ,"h_R_s2_dedx_x" ,200,-1,1,200,0, 10);
  h_R_a1_asum_x  = new TH2D("h_R_a1_asum_x" ,"h_R_a1_asum_x" ,200,-1,1,200,0, 100);
  h_R_a1_asum_p  = new TH2D("h_R_a1_asum_p" ,"h_R_a1_asum_p" ,200,1,3,200,0,100);
  h_R_a1_asum_m2 = new TH2D("h_R_a1_asum_m2","h_R_a1_asum_m2",200,-0.5,2.5,200,0,100);
  h_R_a2_asum_x  = new TH2D("h_R_a2_asum_x" ,"h_R_a2_asum_x" ,200,-1,1,200,0,100);
  h_R_a2_asum_p  = new TH2D("h_R_a2_asum_p" ,"h_R_a2_asum_p" ,200,1,3,200,0,100);
  h_R_a2_asum_m2 = new TH2D("h_R_a2_asum_m2","h_R_a2_asum_m2",200,-0.5,2.5,200,0,100);

  h_R_tgtime    = new TH1D("h_R_tgtime"   ,"h_R_tgtime" ,400,0,1000);

//////////////
//// LHRS ////
//////////////

/////////////////////
//// Coincidence ////
/////////////////////


} // makehist()

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
bool ana::Close(){

  return true;
}

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
int main(int argc, char** argv){

  gErrorIgnoreLevel = kError;

  int ch;
  extern char *optarg;
  int MaxNum = 0;
  string ifname("test.list");
  while((ch=getopt(argc,argv,"hRLf:n:b"))!=-1){
    switch(ch){
    case 'R':
      RHRS= true;
      break;
    case 'L':
      LHRS= true;
      break;
    case 'f':
      ifname = optarg;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'b':
      gROOT->SetBatch(1);
      batch = true;
      break;
    case 'h':
      cout<<"-f : input root file"<<endl;
      cout<<"-n : max no. of events"<<endl;
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

  ana *Ana = new ana();

  TApplication theApp("App", &argc, argv);

  Ana->SetMaxEvent(MaxNum);  // Set Maximum event roop
  Ana->MakeHist();           // Initialize histograms

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }

  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    Ana->ReadTree(Form("../rootfiles/%s",runname.c_str()));
  }

  Ana->Roop();
  Ana->Draw();               // save histograms to pdf file

  cout<<"Done!"<<endl;

  Ana->Close();
  if(batch){ gSystem->Exit(1); }
  gSystem->Exit(1);
  theApp.Run();
  return 0;

}

