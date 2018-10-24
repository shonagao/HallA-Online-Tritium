// HRS S0 Time-Walk Correction
// Author K. Itabahis )ct. 16th 2018 @ Tohoku Univ.

double param(int i);

// Grobal Function //
 int chmax=17; // channel of S2 PMT + S0 PMT (16) 
 int ch=8;// defolt ch 8;
 double tdc=0.5e-9;//TDC converse ch->sec [ss/ch]


void hrs_s0_tdc(){


   gStyle->SetOptStat(111111111);


//-------- TTree data input ---------------//
 TChain *T = new TChain("T");
 for(Int_t i = 93500; i<93501;i++) /// these run goes 93495 to 93504
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));
 }
  

 //----- Defined Functions ----------//
  TCanvas* c0;
  TH1F* htdc;
  TH1F* htdc_c;
  TF1* ftdc[2];
  int ital=1; 


  char* s0m="(R.s0.rt_c+R.s0.lt_c)/2.0";
  char* s0m_c=Form("(R.s0.rt_c-(%lf)*1.0e-7(1./sqrt(R.s0.ra_c-1./sqrt(%lf)))+R.s0.lt_c-(%lf)*1.0e-7(1./sqrt(R.s0.la_c-1./sqrt(%lf))))",param[0]*1.0e7,param[1],param[2]*1.0e7,param[3]); 







}



double param(int i){

  double param[4];

  param[0]=-2.25e-7;
  param[1]=1.61;
  param[2]=-7.4e-7;
  param[3]=1.0e2;



  return param[i];

}
