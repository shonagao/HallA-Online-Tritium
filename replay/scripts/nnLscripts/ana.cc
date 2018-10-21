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
//// BPM ////
  h_rbay_rbax = new TH2D("h_rbay_rbax","h_rbay_rbax",200,-20,20,200,-20,20);
  set->SetTH2(h_rbay_rbax,"BPM A","X","Y");

//////////////
//// RHRS ////
//////////////

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

