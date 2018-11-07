#include "Setting.h"

using namespace std;
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  int ch;
  extern char *optarg;
  string ifname("input.root");
  string ofname("output.root");

  while((ch=getopt(argc,argv,"hf:w:"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      break;
    case 'w':
      ofname = optarg;
      break;
    case 'h':
      std::cout<<"-f (inputfile): input ROOT file name"<<std::endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }

  Setting *set = new Setting();
  set->Initialize();

  TApplication theApp("App", &argc, argv);

  TChain *oldtree = new TChain("T");
  oldtree->Add(Form("%s",ifname.c_str()));

  double L_tr_n, L_tr_x[20], L_tr_th[20], L_tr_p[20];
  double R_tr_n, R_tr_x[20], R_tr_th[20], R_tr_p[20];

  oldtree->SetBranchStatus("*",0);
  oldtree->SetBranchStatus("rbax"   ,1);
  oldtree->SetBranchStatus("rbay"   ,1);
  oldtree->SetBranchStatus("rbbx"   ,1);
  oldtree->SetBranchStatus("rbby"   ,1);
  oldtree->SetBranchStatus("rbx"    ,1);
  oldtree->SetBranchStatus("rby"    ,1);
  oldtree->SetBranchStatus("bpmaws" ,1);
  oldtree->SetBranchStatus("bpmbws" ,1);
  oldtree->SetBranchStatus("DL.evtype"       ,1);

  oldtree->SetBranchStatus("L.s0.la"         ,1);
  oldtree->SetBranchStatus("L.s0.la_c"       ,1);
  oldtree->SetBranchStatus("L.s0.la_p"       ,1);
  oldtree->SetBranchStatus("L.s0.lt"         ,1);
  oldtree->SetBranchStatus("L.s0.lt_c"       ,1);
  oldtree->SetBranchStatus("L.s0.ra"         ,1);
  oldtree->SetBranchStatus("L.s0.ra_c"       ,1);
  oldtree->SetBranchStatus("L.s0.ra_p"       ,1);
  oldtree->SetBranchStatus("L.s0.rt"         ,1);
  oldtree->SetBranchStatus("L.s0.rt_c"       ,1);
  oldtree->SetBranchStatus("L.s0.t_pads"     ,1);
  oldtree->SetBranchStatus("L.s0.time"       ,1);
  oldtree->SetBranchStatus("L.s0.dedx"       ,1);
  oldtree->SetBranchStatus("L.s0.trdy"       ,1);
  oldtree->SetBranchStatus("L.s0.troff"      ,1);
  oldtree->SetBranchStatus("L.s0.trpad"      ,1);
  oldtree->SetBranchStatus("L.s0.trpath"     ,1);
  oldtree->SetBranchStatus("L.s0.trx"        ,1);
  oldtree->SetBranchStatus("L.s0.try"        ,1);
  oldtree->SetBranchStatus("L.s0.lnhits"     ,1);
  oldtree->SetBranchStatus("L.s0.lpeak"      ,1);
  oldtree->SetBranchStatus("L.s0.lt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s0.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("L.s0.rnhits"     ,1);
  oldtree->SetBranchStatus("L.s0.rpeak"      ,1);
  oldtree->SetBranchStatus("L.s0.rt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s0.rtc_fadc"   ,1);

  oldtree->SetBranchStatus("L.s2.la"         ,1);
  oldtree->SetBranchStatus("L.s2.la_c"       ,1);
  oldtree->SetBranchStatus("L.s2.la_p"       ,1);
  oldtree->SetBranchStatus("L.s2.lt"         ,1);
  oldtree->SetBranchStatus("L.s2.lt_c"       ,1);
  oldtree->SetBranchStatus("L.s2.ra"         ,1);
  oldtree->SetBranchStatus("L.s2.ra_c"       ,1);
  oldtree->SetBranchStatus("L.s2.ra_p"       ,1);
  oldtree->SetBranchStatus("L.s2.rt"         ,1);
  oldtree->SetBranchStatus("L.s2.rt_c"       ,1);
  oldtree->SetBranchStatus("L.s2.t_pads"     ,1);
  oldtree->SetBranchStatus("L.s2.time"       ,1);
  oldtree->SetBranchStatus("L.s2.dedx"       ,1);
  oldtree->SetBranchStatus("L.s2.trdx"       ,1);
  oldtree->SetBranchStatus("L.s2.troff"      ,1);
  oldtree->SetBranchStatus("L.s2.trpad"      ,1);
  oldtree->SetBranchStatus("L.s2.trpath"     ,1);
  oldtree->SetBranchStatus("L.s2.trx"        ,1);
  oldtree->SetBranchStatus("L.s2.try"        ,1);
  oldtree->SetBranchStatus("L.s2.lnhits"     ,1);
  oldtree->SetBranchStatus("L.s2.lpeak"      ,1);
  oldtree->SetBranchStatus("L.s2.lt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s2.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("L.s2.rnhits"     ,1);
  oldtree->SetBranchStatus("L.s2.rpeak"      ,1);
  oldtree->SetBranchStatus("L.s2.rt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s2.rtc_fadc"   ,1);

  oldtree->SetBranchStatus("L.tr.n"          ,1);  oldtree->SetBranchAddress("L.tr.n" ,&L_tr_n );
  oldtree->SetBranchStatus("L.tr.flag"       ,1);
  oldtree->SetBranchStatus("L.tr.ndof"       ,1);
  oldtree->SetBranchStatus("L.tr.chi2"       ,1);
  oldtree->SetBranchStatus("L.tr.beta"       ,1);
  oldtree->SetBranchStatus("L.tr.d_x"        ,1);
  oldtree->SetBranchStatus("L.tr.d_y"        ,1);
  oldtree->SetBranchStatus("L.tr.d_th"       ,1);
  oldtree->SetBranchStatus("L.tr.d_ph"       ,1);
  oldtree->SetBranchStatus("L.tr.r_x"        ,1);
  oldtree->SetBranchStatus("L.tr.r_y"        ,1);
  oldtree->SetBranchStatus("L.tr.r_th"       ,1);
  oldtree->SetBranchStatus("L.tr.r_ph"       ,1);
  oldtree->SetBranchStatus("L.tr.x"          ,1);  oldtree->SetBranchAddress("L.tr.x" ,L_tr_x );
  oldtree->SetBranchStatus("L.tr.y"          ,1);
  oldtree->SetBranchStatus("L.tr.th"         ,1);  oldtree->SetBranchAddress("L.tr.th",L_tr_th);
  oldtree->SetBranchStatus("L.tr.ph"         ,1);
  oldtree->SetBranchStatus("L.tr.time"       ,1);
  oldtree->SetBranchStatus("L.tr.p"          ,1);  oldtree->SetBranchAddress("L.tr.p" ,L_tr_p );
  oldtree->SetBranchStatus("L.tr.pathl"      ,1);
  oldtree->SetBranchStatus("L.tr.px"         ,1);
  oldtree->SetBranchStatus("L.tr.py"         ,1);
  oldtree->SetBranchStatus("L.tr.pz"         ,1);
  oldtree->SetBranchStatus("L.tr.tg_dp"      ,1);
  oldtree->SetBranchStatus("L.tr.tg_y"       ,1);
  oldtree->SetBranchStatus("L.tr.tg_th"      ,1);
  oldtree->SetBranchStatus("L.tr.tg_ph"      ,1);
  oldtree->SetBranchStatus("L.tr.vx"         ,1);
  oldtree->SetBranchStatus("L.tr.vy"         ,1);
  oldtree->SetBranchStatus("L.tr.vz"         ,1);

  oldtree->SetBranchStatus("DL.ltRFtime"     ,1);
  oldtree->SetBranchStatus("LTDC.F1FirstHit" ,1);

  oldtree->SetBranchStatus("Ndata.DR.T1"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T2"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T3"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T4"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T5"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T6"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T7"     ,1);
  oldtree->SetBranchStatus("Ndata.DR.T8"     ,1);
  oldtree->SetBranchStatus("DR.evtype"       ,1);

  oldtree->SetBranchStatus("R.s0.la"         ,1);
  oldtree->SetBranchStatus("R.s0.la_c"       ,1);
  oldtree->SetBranchStatus("R.s0.la_p"       ,1);
  oldtree->SetBranchStatus("R.s0.lt"         ,1);
  oldtree->SetBranchStatus("R.s0.lt_c"       ,1);
  oldtree->SetBranchStatus("R.s0.ra"         ,1);
  oldtree->SetBranchStatus("R.s0.ra_c"       ,1);
  oldtree->SetBranchStatus("R.s0.ra_p"       ,1);
  oldtree->SetBranchStatus("R.s0.rt"         ,1);
  oldtree->SetBranchStatus("R.s0.rt_c"       ,1);
  oldtree->SetBranchStatus("R.s0.t_pads"     ,1);
  oldtree->SetBranchStatus("R.s0.time"       ,1);
  oldtree->SetBranchStatus("R.s0.dedx"       ,1);
  oldtree->SetBranchStatus("R.s0.trdy"       ,1);
  oldtree->SetBranchStatus("R.s0.troff"      ,1);
  oldtree->SetBranchStatus("R.s0.trpad"      ,1);
  oldtree->SetBranchStatus("R.s0.trpath"     ,1);
  oldtree->SetBranchStatus("R.s0.trx"        ,1);
  oldtree->SetBranchStatus("R.s0.try"        ,1);
  oldtree->SetBranchStatus("R.s0.lnhits"     ,1);
  oldtree->SetBranchStatus("R.s0.lpeak"      ,1);
  oldtree->SetBranchStatus("R.s0.lt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s0.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("R.s0.rnhits"     ,1);
  oldtree->SetBranchStatus("R.s0.rpeak"      ,1);
  oldtree->SetBranchStatus("R.s0.rt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s0.rtc_fadc"   ,1);

  oldtree->SetBranchStatus("R.s2.la"         ,1);
  oldtree->SetBranchStatus("R.s2.la_c"       ,1);
  oldtree->SetBranchStatus("R.s2.la_p"       ,1);
  oldtree->SetBranchStatus("R.s2.lt"         ,1);
  oldtree->SetBranchStatus("R.s2.lt_c"       ,1);
  oldtree->SetBranchStatus("R.s2.ra"         ,1);
  oldtree->SetBranchStatus("R.s2.ra_c"       ,1);
  oldtree->SetBranchStatus("R.s2.ra_p"       ,1);
  oldtree->SetBranchStatus("R.s2.rt"         ,1);
  oldtree->SetBranchStatus("R.s2.rt_c"       ,1);
  oldtree->SetBranchStatus("R.s2.t_pads"     ,1);
  oldtree->SetBranchStatus("R.s2.time"       ,1);
  oldtree->SetBranchStatus("R.s2.dedx"       ,1);
  oldtree->SetBranchStatus("R.s2.trdx"       ,1);
  oldtree->SetBranchStatus("R.s2.troff"      ,1);
  oldtree->SetBranchStatus("R.s2.trpad"      ,1);
  oldtree->SetBranchStatus("R.s2.trpath"     ,1);
  oldtree->SetBranchStatus("R.s2.trx"        ,1);
  oldtree->SetBranchStatus("R.s2.try"        ,1);
  oldtree->SetBranchStatus("R.s2.lnhits"     ,1);
  oldtree->SetBranchStatus("R.s2.lpeak"      ,1);
  oldtree->SetBranchStatus("R.s2.lt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s2.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("R.s2.rnhits"     ,1);
  oldtree->SetBranchStatus("R.s2.rpeak"      ,1);
  oldtree->SetBranchStatus("R.s2.rt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s2.rtc_fadc"   ,1);

  oldtree->SetBranchStatus("R.a1.a"          ,1);
  oldtree->SetBranchStatus("R.a1.a_c"        ,1);
  oldtree->SetBranchStatus("R.a1.a_p"        ,1);
  oldtree->SetBranchStatus("R.a1.t"          ,1);
  oldtree->SetBranchStatus("R.a1.t_c"        ,1);
  oldtree->SetBranchStatus("R.a1.asum_p"     ,1);
  oldtree->SetBranchStatus("R.a1.asum_c"     ,1);
  oldtree->SetBranchStatus("R.a1.trpath"     ,1);
  oldtree->SetBranchStatus("R.a1.trx"        ,1);
  oldtree->SetBranchStatus("R.a1.try"        ,1);
  oldtree->SetBranchStatus("R.a1.nhits"      ,1);
  oldtree->SetBranchStatus("R.a1.peak"       ,1);
  oldtree->SetBranchStatus("R.a1.t_fadc"     ,1);
  oldtree->SetBranchStatus("R.a1.tc_fadc"    ,1);

  oldtree->SetBranchStatus("R.a2.a"          ,1);
  oldtree->SetBranchStatus("R.a2.a_c"        ,1);
  oldtree->SetBranchStatus("R.a2.a_p"        ,1);
  oldtree->SetBranchStatus("R.a2.t"          ,1);
  oldtree->SetBranchStatus("R.a2.t_c"        ,1);
  oldtree->SetBranchStatus("R.a2.asum_p"     ,1);
  oldtree->SetBranchStatus("R.a2.asum_c"     ,1);
  oldtree->SetBranchStatus("R.a2.trpath"     ,1);
  oldtree->SetBranchStatus("R.a2.trx"        ,1);
  oldtree->SetBranchStatus("R.a2.try"        ,1);
  oldtree->SetBranchStatus("R.a2.nhits"      ,1);
  oldtree->SetBranchStatus("R.a2.peak"       ,1);
  oldtree->SetBranchStatus("R.a2.t_fadc"     ,1);
  oldtree->SetBranchStatus("R.a2.tc_fadc"    ,1);

  oldtree->SetBranchStatus("R.cer.a"         ,1);
  oldtree->SetBranchStatus("R.cer.a_c"       ,1);
  oldtree->SetBranchStatus("R.cer.a_p"       ,1);
  oldtree->SetBranchStatus("R.cer.t"         ,1);
  oldtree->SetBranchStatus("R.cer.t_c"       ,1);
  oldtree->SetBranchStatus("R.cer.asum_p"    ,1);
  oldtree->SetBranchStatus("R.cer.asum_c"    ,1);
  oldtree->SetBranchStatus("R.cer.trpath"    ,1);
  oldtree->SetBranchStatus("R.cer.trx"       ,1);
  oldtree->SetBranchStatus("R.cer.try"       ,1);
  oldtree->SetBranchStatus("R.cer.nhits"     ,1);
  oldtree->SetBranchStatus("R.cer.peak"      ,1);
  oldtree->SetBranchStatus("R.cer.t_fadc"    ,1);
  oldtree->SetBranchStatus("R.cer.tc_fadc"   ,1);

  oldtree->SetBranchStatus("R.ps.asum_p"     ,1);
  oldtree->SetBranchStatus("R.ps.asum_c"     ,1);

  oldtree->SetBranchStatus("R.sh.asum_p"     ,1);
  oldtree->SetBranchStatus("R.sh.asum_c"     ,1);

  oldtree->SetBranchStatus("R.tr.n"          ,1);  oldtree->SetBranchAddress("R.tr.n" ,&R_tr_n );
  oldtree->SetBranchStatus("R.tr.flag"       ,1);
  oldtree->SetBranchStatus("R.tr.ndof"       ,1);
  oldtree->SetBranchStatus("R.tr.chi2"       ,1);
  oldtree->SetBranchStatus("R.tr.beta"       ,1);
  oldtree->SetBranchStatus("R.tr.d_x"        ,1);
  oldtree->SetBranchStatus("R.tr.d_y"        ,1);
  oldtree->SetBranchStatus("R.tr.d_th"       ,1);
  oldtree->SetBranchStatus("R.tr.d_ph"       ,1);
  oldtree->SetBranchStatus("R.tr.r_x"        ,1);  oldtree->SetBranchAddress("R.tr.x" ,R_tr_x );
  oldtree->SetBranchStatus("R.tr.r_y"        ,1);
  oldtree->SetBranchStatus("R.tr.r_th"       ,1);  oldtree->SetBranchAddress("R.tr.th",R_tr_th);
  oldtree->SetBranchStatus("R.tr.r_ph"       ,1);
  oldtree->SetBranchStatus("R.tr.x"          ,1);
  oldtree->SetBranchStatus("R.tr.y"          ,1);  oldtree->SetBranchAddress("R.tr.p" ,R_tr_p );
  oldtree->SetBranchStatus("R.tr.th"         ,1);
  oldtree->SetBranchStatus("R.tr.ph"         ,1);
  oldtree->SetBranchStatus("R.tr.time"       ,1);
  oldtree->SetBranchStatus("R.tr.p"          ,1);
  oldtree->SetBranchStatus("R.tr.pathl"      ,1);
  oldtree->SetBranchStatus("R.tr.px"         ,1);
  oldtree->SetBranchStatus("R.tr.py"         ,1);
  oldtree->SetBranchStatus("R.tr.pz"         ,1);
  oldtree->SetBranchStatus("R.tr.tg_dp"      ,1);
  oldtree->SetBranchStatus("R.tr.tg_y"       ,1);
  oldtree->SetBranchStatus("R.tr.tg_th"      ,1);
  oldtree->SetBranchStatus("R.tr.tg_ph"      ,1);
  oldtree->SetBranchStatus("R.tr.vx"         ,1);
  oldtree->SetBranchStatus("R.tr.vy"         ,1);
  oldtree->SetBranchStatus("R.tr.vz"         ,1);

  oldtree->SetBranchStatus("DR.rtRFtime"     ,1);
  oldtree->SetBranchStatus("RTDC.F1FirstHit" ,1);


  TFile *ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  TTree *newtree = oldtree->CloneTree(0);

  int ENum = oldtree->GetEntries();

  for(int n=0;n<ENum;n++){
    oldtree->GetEntry(n);

    int NLtr = (int)L_tr_n;
    int NRtr = (int)R_tr_n;
    bool LOK = false;
    bool ROK = false;

    for(int lt=0;lt<NLtr;lt++){
      if(   L_tr_th[lt]<0.17*L_tr_x[lt]+0.050
         && L_tr_th[lt]>0.17*L_tr_x[lt]-0.060
         && L_tr_p[lt]>1. && L_tr_p[lt]<10.) LOK = true;
    }

    for(int rt=0;rt<NRtr;rt++){
      if(   R_tr_th[rt]<0.17*R_tr_x[rt]+0.050
         && R_tr_th[rt]>0.17*R_tr_x[rt]-0.060
         && R_tr_p[rt]>1. && R_tr_p[rt]<10. ) ROK = true;
    }

    if( LOK && ROK ){
      newtree->Fill();
    }

    if(n % 100000 == 0){ cout<<n<<" / "<<ENum<<endl; }

  } // for ENum

  newtree->AutoSave();
  delete ofp;

  gSystem->Exit(1);
  theApp.Run();
  return 0;

} 



