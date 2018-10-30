// 10/4/2018 Author Itabshi
// time walk correction R-HRS S2


//====== Define Function ===========//
double range_para(int i,int j);
double range_adc(int i,int j);
long double fit_ini(int i,int j,int k);

//=================================//

#include<string>

// Grobal Function //
 int chmax=17; // channel of S2 PMT 
 int ch=8;// defolt ch 8;
 double tdc=0.5e-9;//TDC converse ch->sec [ss/ch]


void hrs_tof_twcorr(){

//-------- TTree data input ---------------//


  TChain*  T=new TChain("T");;
 for(Int_t i = 93500; i<93501;i++) /// these run goes 93495 to 93504
 {
    T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));
 }



 int max=10000; 
 int Ndata_Rs0la_c,Ndata_Rs0ra_c,Ndata_Rs0rt_c,Ndata_Rs0lt_c;//S0 Ndata
 int Ndata_Rs2la_c,Ndata_Rs2ra_c,Ndata_Rs2rt_c,Ndata_Rs2lt_c;//S2 Ndata
 double DR_evtypebits; // trigger condition  
 double Rs0la_c[max],Rs0ra_c[max],Rs0lt_c[max],Rs0rt_c[max]; //S0 TDC & ADC  
 double Rs2la_c[max],Rs2ra_c[max],Rs2lt_c[max],Rs2rt_c[max]; //S2 TDC & ADC 



 //============= Set Branch Status ==================//
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("Ndata.R.s0.la_c",1);
  T->SetBranchStatus("Ndata.R.s0.ra_c",1);
  T->SetBranchStatus("Ndata.R.s0.lt_c",1);
  T->SetBranchStatus("Ndata.R.s0.rt_c",1);
  T->SetBranchStatus("Ndata.R.s2.la_c",1);
  T->SetBranchStatus("Ndata.R.s2.ra_c",1);
  T->SetBranchStatus("Ndata.R.s2.lt_c",1);
  T->SetBranchStatus("R.s2.rt_c",1); 
  T->SetBranchStatus("R.s0.la_c",1);
  T->SetBranchStatus("R.s0.ra_c",1);
  T->SetBranchStatus("R.s0.lt_c",1);
  T->SetBranchStatus("R.s0.rt_c",1);
  T->SetBranchStatus("R.s2.la_c",1);
  T->SetBranchStatus("R.s2.ra_c",1);
  T->SetBranchStatus("R.s2.lt_c",1);
  T->SetBranchStatus("R.s2.rt_c",1);
  T->SetBranchStatus("DR.evtypebits",1);   


  // Right Arm S0 Ndata// 
  T->SetBranchAddress("Ndata.R.s0.la_c",&Ndata_Rs0la_c); // Right arm S0-Top(B) ADC
  T->SetBranchAddress("Ndata.R.s0.ra_c",&Ndata_Rs0ra_c); // Right arm S0-Bottom(A) ADC
  T->SetBranchAddress("Ndata.R.s2.la_c",&Ndata_Rs2la_c); // Right arm S2-Top(B) ADC
  T->SetBranchAddress("Ndata.R.s2.ra_c",&Ndata_Rs2ra_c); // Right arm S2-Bottom(A) ADC
  // Right Arm S2 Ndata//
  T->SetBranchAddress("Ndata.R.s0.lt_c",&Ndata_Rs0lt_c); // Right arm S0-Top(B) TDC
  T->SetBranchAddress("Ndata.R.s0.rt_c",&Ndata_Rs0rt_c); // Right arm S0-Bottom(A) TDC
  T->SetBranchAddress("Ndata.R.s2.lt_c",&Ndata_Rs2lt_c); // Right arm S2-Top(B) TDC
  T->SetBranchAddress("Ndata.R.s2.rt_c",&Ndata_Rs2rt_c); // Right arm S2-Bottom(A) TDC
  // Right Arm S0 // 
  T->SetBranchAddress("R.s0.la_c",Rs0la_c); // Right arm S0-Top(B) ADC
  T->SetBranchAddress("R.s0.ra_c",Rs0ra_c); // Right arm S0-Bottom(A) ADC
  T->SetBranchAddress("R.s0.lt_c",Rs0lt_c); // Right arm S0-Top(B) TDC
  T->SetBranchAddress("R.s0.rt_c",Rs0rt_c); // Right arm S0-Bottom(A) TDC
  // Right Arm S2 // 
  T->SetBranchAddress("R.s2.la_c",Rs2la_c); // Right arm S2-Top(B) ADC
  T->SetBranchAddress("R.s2.ra_c",Rs2ra_c); // Right arm S2-Bottom(A) ADC
  T->SetBranchAddress("R.s2.lt_c",Rs2lt_c); // Right arm S2-Top(B) TDC
  T->SetBranchAddress("R.s2.rt_c",Rs2rt_c); // Right arm S2-Bottom(A) TDC
 
  // Trigger condtion //
  T->SetBranchAddress("DR.evtypebits",&DR_evtypebits); 

  T->Fill();




  
 //========== Defined Parameters =======================//
 //========= Time walk parameters =================//
 //TH2F* htof_adc[chmax][2];
 TCanvas* cw[chmax];
 TCanvas* cwo[chmax];
 TCanvas* ctdc[chmax];
 TCanvas* ct[chmax];
 TProfile* tp[chmax][2];
 TF1* ffit[chmax][2];
 TH2F* htof_adc[chmax][2];
 TH2F* htof_adc_c[chmax][2];
 TH2F* htdc_adc[chmax][2];
 TH2F* htdc_adc_c[chmax][2];
 TH1F* htdc[chmax][2];
 TH1F* htdc_c[chmax][2];
 TF1* ftdc[chmax][2];
 TF1* ftdc_c[chmax][2];
 // TH2F* FS_tof[chmax][2];
 double min_tof,max_tof,min_adc,max_adc;
 double bin_tof;
 int bin_adc;

 
 // double min_tdc,max_tdc;
 int bin_tdc;
 long double p0[chmax][2],p1[chmax][2];
 // long double sig[chmax][2],sig_c[chmax][2];


 char* cut_time[chmax];
 char *adc[chmax][2];
 char *tdc[chmax][2];
 TCut cut_trig;
 cut_trig="DR.evtypebits==112";
 std::string str("DR.evtypebits==112");


 //================== TOF parameters ========================// 
 TCanvas* ctof_c[chmax];
 TH1F* htof[chmax];
 TH1F* htof_c[chmax];
 TF1* ftof[chmax];
 TF1* ftof_c[chmax];
 TH1F* htdc_m[chmax];
 TH1D* hmean_tof[chmax][2];
double tmin,tmax,tcmin,tcmax;
 // double mean_tof[chmax][2];
 
 double  tofp[chmax],tofp_c[chmax];
 char* ctdc_c[chmax][2];
 char* ctdc_mc[chmax];
 char* tof[chmax];

 double min_tof_c[chmax],max_tof_c[chmax],min_adc_c[chmax],max_adc_c[chmax];
 double bin_tof_c[chmax];
 int bin_adc_c[chmax];


 //============================================================//
 //========== Time Walk correction ============================//
 //============================================================//
 

 for(int i=0;i<chmax;i++){  
      if(i==ch || i==16){


       // Definition TCanvas //
     cwo[i]=new TCanvas(Form("cwo[%d]",i),Form("cwo[%d]",i));
                         cwo[i]->Divide(1,2);

   for(int j=0;j<2;j++){


     //======== Parameter Setting ==================//

 min_tof=range_para(i,0);
 max_tof=range_para(i,1);
 bin_tof=(range_para(i,1)-range_para(i,0))/0.5e-9;
 bin_tof=(int)bin_tof;

 min_adc=range_para(i,4);
 max_adc=range_para(i,5); 
 bin_adc=(-range_para(i,4)+range_para(i,5));

       //Cut conditon //
 cut_time[i]=Form("Rs2lt_c[%d]>0 && Rs2rt_c[%d]>0",ch,ch);
 cut_time[16]="Rs0lt_c>0 && Rs0rt_c>0";
 cut_time[i]=Form("Rs2lt_c[%d]>0 && Rs2rt_c[%d]>0",i,i);
 cut_time[16]="Rs0lt_c>0 && Rs0rt_c>0";



    // ADC & TDC inf//
 adc[i][0]=Form("Rs2ra_c[%d]",i);
 adc[i][1]=Form("Rs2la_c[%d]",i);
 adc[16][0]="Rs0ra_c";
 adc[16][1]="Rs0la_c";
 tdc[i][0]=Form("Rs2rt_c[%d]",i);
 tdc[i][1]=Form("Rs2lt_c[%d]",i);
 tdc[16][0]="Rs0rt_c";
 tdc[16][1]="Rs0lt_c";

 tof[i]=Form("(Rs2rt_c[%d]+Rs2lt_c[%d])/2.0-(Rs0rt_c+Rs0lt_c)/2.0",i,i);
 tof[16]=Form("(R.s2.rt_c[%d]+R.s2.lt_c[%d])/2.0-(R.s0.rt_c+R.s0.lt_c)/2.0",ch,ch);

 //======================================================//


     //--- TOF vs ADC 2D Hist --------//
   cwo[i]->cd(j+1);
 htof_adc[i][j]=new TH2F(Form("htof_adc[%d][%d]",i,j),Form("TOF vs Scinti CH %d,[%d]-PMT adc w/o correlation ",i,j),bin_adc,min_adc,max_adc,bin_tof,min_tof,max_tof);
   


       T->Project(Form("htof_adc[%d][%d]",i,j),Form("%s:%s",tof[i],adc[i][j]),cut_trig && cut_time[i] && cut_time[16]);
  
       /*
 if(cut_trig && cut_time[i] && cut_time[16]){

   htof_adc[i][j]->Fill(adc[i][j],tof[i]);

 }
       */

   htof_adc[i][j]->Draw();
   
      int xmin,xmax;
      //xmin=htof_adc[i][j]->GetXaxis()->FindBin(min_adc);
      //xmax=htof_adc[i][j]->GetXaxis()->FindBin(max_adc);
      xmin=htof_adc[i][j]->GetYaxis()->FindBin(min_tof);
      xmax=htof_adc[i][j]->GetYaxis()->FindBin(max_tof);
      //  cout<<"xmin :"<<xmin<<"xmax :"<<xmax<<endl;

      TObjArray aSlices;
      TF1 *test_gaus=new TF1("test_gaus","gaus",min_tof,max_tof);
 
      //FitSlicesY 
      //   double mean_tof[i][j]=htof_adc[i][j]->ProjectionY()->GetBinCenter(htof_adc[i][j]->ProjectionY()->GetMaximumBin());
      htof_adc[i][j]->FitSlicesY(test_gaus,0,-1,0,"QRG3");
 hmean_tof[i][j]=(TH1D*)gROOT->FindObject(Form("htof_adc[%d][%d]_1",i,j));





ffit[i][j]=new TF1(Form("ffit[%d][%d]",i,j),"[0]*(1./sqrt(x)-1./sqrt([1]))",min_adc,max_adc);

 ffit[i][j]->SetParameter(0,fit_ini(i,j,0));   
 ffit[i][j]->SetParameter(1,fit_ini(i,j,1));   

 // ffit[i][j]->FixParameter(0,fit_ini(i,j,0));   
 // ffit[i][j]->FixParameter(1,1.0e8);   

 // ffit[i][j]->SetParLimits(0,0.5*fit_ini(i,j,0),1.5*fit_ini(i,j,0));   
 //ffit[i][j]->SetParLimits(1,0.5*fit_ini(i,j,1),1.5*fit_ini(i,j,1));   


 hmean_tof[i][j]->Fit(Form("ffit[%d][%d]",i,j),"","");


         htof_adc[i][j]->Draw();
      hmean_tof[i][j]->Draw("same");


      p0[i][j]=ffit[i][j]->GetParameter(0);
      p1[i][j]=ffit[i][j]->GetParameter(1);
  
   p0[i][j]=p0[i][j]*1.0e9;
   p1[i][j]=p1[i][j]*1.0e5;

   }



 //======================================//
 //======= TOF ANALYSIS =================//
 //======================================//
 
 
  if(i==16){
 

       //Cut conditon //
 cut_time[ch]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",ch,ch);
 cut_time[16]="R.s0.lt_c>0 && R.s0.rt_c>0";
 TCut cut_trig="DR.evtypebits==112";
 cut_time[ch]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",ch,ch);
 cut_time[16]="R.s0.lt_c>0 && R.s0.rt_c>0";

    // ADC & TDC inf//
 adc[ch][0]=Form("R.s2.ra_c[%d]",ch);
 adc[ch][1]=Form("R.s2.la_c[%d]",ch);
 adc[16][0]="R.s0.ra_c";
 adc[16][1]="R.s0.la_c";

 tdc[ch][0]=Form("R.s2.rt_c[%d]",ch);
 tdc[ch][1]=Form("R.s2.lt_c[%d]",ch);
 tdc[16][0]="R.s0.rt_c";
 tdc[16][1]="R.s0.lt_c";

 //--- range parameters ------//



 tmin=-1.0e-8;
 tmax=1.0e-8;
 tcmin=-1.0e-8;
 tcmax=1.0e-8;
 int tbin=400;
 int tcbin=400;





 ct[ch]=new TCanvas(Form("ct[%d]",ch),Form("ct[%d]",ch),600,600);
 ct[ch]->Divide(1,2);
 ct[ch]->cd(1);
 


 //------- TOF w/o TW correction -------------//
htof[ch]=new TH1F(Form("htof[%d]",ch),Form("S2 Ch %d -S0  w/o TW correction TOF hist ",ch),tbin,tmin,tmax);
 T->Project(Form("htof[%d]",ch),Form("((%s)+(%s))/2.0-((%s)+(%s))/2.0",tdc[ch][0],tdc[ch][1],tdc[16][0],tdc[16][1]),cut_trig && cut_time[ch]&&cut_time[16]);
 ftof[ch]=new TF1(Form("ftof[%d]",ch),"gaus",tmin,tmax);
 htof[ch]->Fit(Form("ftof[%d]",ch),"","");
   htof[ch]->Draw();
   ftof[ch]->Draw("same");
 //------ TOF w/ TW correction-----------------//
   ct[ch]->cd(2); 
 htof_c[ch]=new TH1F(Form("htof_c[%d]",ch),Form("S2 Ch %d -S0  w/ TW correction TOF hist ",ch),tcbin,tcmin,tcmax);
 T->Project(Form("htof_c[%d]",ch),Form("(((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)))+((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))))/2.0-(((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)))+((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))))/2.0",tdc[ch][0],p0[ch][0],adc[ch][0],p1[ch][0],tdc[ch][1],p0[ch][1],adc[ch][1],p1[ch][1],tdc[16][0],p0[16][0],adc[16][0],p1[16][0],tdc[16][1],p0[16][1],adc[16][1],p1[16][1]),cut_trig && cut_time[ch] && cut_time[16]);



 ftof_c[ch]=new TF1(Form("ftof_c[%d]",ch),"gaus",tcmin,tcmax);
 htof_c[ch]->Fit(Form("ftof_c[%d]",ch),"","");
   htof_c[ch]->Draw();
   ftof_c[ch]->Draw("same");

   tofp[ch]=ftof[ch]->GetParameter(2);
   tofp_c[ch]=ftof_c[ch]->GetParameter(2);
   double tof_c_sig=ftof_c[ch]->GetParameter(3);



 //------- TOF vs ADC w/ TW Correction----------//


   min_tof_c[ch]=tofp[ch]-5*tof_c_sig;
   //range_para(ch,2);
 max_tof_c[ch]=tofp[ch]+5*tof_c_sig;
   //range_para(ch,3);
 bin_tof_c[ch]=(max_tof_c[ch]-min_tof_c[ch])/0.5e-9;
 bin_tof_c[ch]=(int)bin_tof_c[ch];
 


 min_adc_c[ch]=range_para(ch,4);
 max_adc_c[ch]=range_para(ch,5);
 bin_adc_c[ch]=(-range_para(ch,4)+range_para(ch,5));
 
 htof_adc_c[ch][0]=new TH2F(Form("htof_adc_c[%d][%d]",ch,0),Form("TOF vs ADC ch%d R-PMT ",ch),bin_adc_c[ch],min_adc_c[ch],max_adc_c[ch],bin_tof_c[ch],min_tof_c[ch],max_tof_c[ch]);
 htof_adc_c[ch][1]=new TH2F(Form("htof_adc_c[%d][%d]",ch,1),Form("TOF vs ADC ch%d L-PMT ",ch),bin_adc_c[ch],min_adc_c[ch],max_adc_c[ch],bin_tof_c[ch],min_tof_c[ch],max_tof_c[ch]);

 min_adc_c[16]=range_para(ch,4);
 max_adc_c[16]=range_para(ch,5);
 bin_adc_c[16]=(-range_para(ch,4)+range_para(ch,5));


 


 ctof_c[ch]=new TCanvas(Form("ctof_c[%d]",ch),Form("ctof_c[%d]",ch));
 ctof_c[16]=new TCanvas(Form("ctof_c[%d]",16),Form("ctof_c[%d]",16));
 //ctof_c[ch]->Divide(2,1);
 //ctof_c[16]->Divide(2,1);


 htof_adc_c[16][0]=new TH2F(Form("htof_adc_c[%d][%d]",16,0),Form("TOF vs ADC ch%d R-PMT ",16),bin_adc_c[16],min_adc_c[16],max_adc_c[16],bin_tof_c[ch],min_tof_c[ch],max_tof_c[ch]);
 htof_adc_c[16][1]=new TH2F(Form("htof_adc_c[%d][%d]",ch,1),Form("TOF vs ADC ch%d L-PMT ",ch),bin_adc_c[16],min_adc_c[16],max_adc_c[16],bin_tof_c[ch],min_tof_c[ch],max_tof_c[ch]);





       //Cut conditon //
 cut_time[i]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",ch,ch);
 cut_time[16]="R.s0.lt_c>0 && R.s0.rt_c>0";
 cut_trig="DR.evtypebits==112";
 cut_time[i]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",i,i);
 cut_time[16]="R.s0.lt_c>0 && R.s0.rt_c>0";



    // ADC & TDC inf//
 adc[ch][0]=Form("R.s2.ra_c[%d]",ch);
 adc[ch][1]=Form("R.s2.la_c[%d]",ch);
 adc[16][0]="R.s0.ra_c";
 adc[16][1]="R.s0.la_c";

 tdc[ch][0]=Form("R.s2.rt_c[%d]",ch);
 tdc[ch][1]=Form("R.s2.lt_c[%d]",ch);
 tdc[16][0]="R.s0.rt_c";
 tdc[16][1]="R.s0.lt_c";


 // ctof_c[ch]->cd(1);



 T->Project(Form("htof_adc_c[%d][%d]",ch,0),Form("(((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)))+((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))))/2.0-(((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)))+((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))))/2.0:%s",tdc[ch][0],p0[ch][0],adc[ch][0],p1[ch][0],tdc[ch][1],p0[ch][1],adc[ch][1],p1[ch][1],tdc[16][0],p0[16][0],adc[16][0],p1[16][0],tdc[16][1],p0[16][1],adc[16][1],p1[16][1],adc[ch][0]));
	    //,cut_trig && cut_time[ch]&& cut_time[16]);

 htof_adc_c[ch][0]->Draw();
 // ctof_c[ch]->cd(2);

 T->Project(Form("htof_adc_c[%d][%d]",ch,1),Form("(((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)))+((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))))/2.0-(((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)))+((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))))/2.0:%s",tdc[ch][0],p0[ch][0],adc[ch][0],p1[ch][0],tdc[ch][1],p0[ch][1],adc[ch][1],p1[ch][1],tdc[16][0],p0[16][0],adc[16][0],p1[16][0],tdc[16][1],p0[16][1],adc[16][1],p1[16][1],adc[ch][1]));
	    //,cut_trig && cut_time[ch]&& cut_time[16]);


 //htof_adc_c[ch][1]->Draw();
 //ctof_c[16]->cd(1);

 T->Project(Form("htof_adc_c[%d][%d]",16,0),Form("(((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)))+((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))))/2.0-(((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)))+((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))))/2.0:%s",tdc[ch][0],p0[ch][0],adc[ch][0],p1[ch][0],tdc[ch][1],p0[ch][1],adc[ch][1],p1[ch][1],tdc[16][0],p0[16][0],adc[16][0],p1[16][0],tdc[16][1],p0[16][1],adc[16][1],p1[16][1],adc[16][0]));
	    //,cut_trig && cut_time[ch]&& cut_time[16]);

 //htof_adc_c[16][0]->Draw();
 // ctof_c[16]->cd(2);

 T->Project(Form("htof_adc_c[%d][%d]",16,1),Form("(((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)))+((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))))/2.0-(((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)))+((%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))))/2.0:%s",tdc[ch][0],p0[ch][0],adc[ch][0],p1[ch][0],tdc[ch][1],p0[ch][1],adc[ch][1],p1[ch][1],tdc[16][0],p0[16][0],adc[16][0],p1[16][0],tdc[16][1],p0[16][1],adc[16][1],p1[16][1],adc[16][1]));
	    //,cut_trig && cut_time[ch]&& cut_time[16]);

 //htof_adc_c[16][1]->Draw();



 


 
  }
   }




 }






 
 for(int i=0;i<chmax;i++){
     if(i==ch || i==16){
//----------------- Comment out -----------------//
   cout<<"======== S2 R-PMT ch"<<i<<"parameters ========"<<endl;
   cout<<"p0 :"<<p0[i][0]*1.0e-9<<"  p1 :"<<p1[i][0]*1.0e-5<<endl;    
   cout<<"======== S2 L-PMT ch"<<i<<"parameters ========"<<endl;
   cout<<"p0 :"<<p0[i][1]*1.0e-9<<"  p1 :"<<p1[i][1]*1.0e-5<<endl;    


     }
 }

//---- Comment out TOF correction get fit parameters --------//

 cout<<"========== TOF get fit parameters ==========="<<endl;
  cout<<"TOF w/o correction sigma "<<tofp[ch]<<endl;
 cout<<"TOF w/ correction sigma "<<tofp_c[ch]<<endl;

 
}


double range_para(int i,int j){

  int npara=6;// number of parameters
  double par[chmax][npara];
  double param;
 
  //=== Inital parameters========//
  for(int k=0;k<chmax;k++){
   

   par[k][0]=-4e-9, par[k][1]=4e-9, par[k][2]=-0.1e-7, par[k][3]=0.1e-7;//TOF


   if(k==16){//S0
  par[k][4]=500.0, par[k][5]=5000.;//ADC
    }else{ //S2
      par[k][4]=100.0, par[k][5]=500.;//ADC
    }
   }

  //===== Set Parameters ========//
  //  par[8][0]=0.18e-6, par[8][1]=0.22e-6, par[8][2]=-0.18e-6,par[8][3]=0.18e-6;
  //  par[16][0]=0.180e-6, par[16][1]=0.23e-6, par[16][2]=-0.18e-6, par[16][3]=0.25e-6;

  //============================//
  for(int k=0;k<chmax;k++){
   for(int l=0;l<npara;l++){
    if(k==i && l==j){
      param=par[k][l];
     
    }

    }
  }
  return param;
}






long double fit_ini(int i,int j,int k){
  // COMMENT  i: ch num, j: R-PMT(0),L-PMT(1), k: para num 
  int npara=2;
  double fpar[chmax][2][npara];
  double f_ini_param;

  //======== Initial parameters =========//

  for(int k=0;k<chmax;k++){
  fpar[k][0][0]= -6.0e-9, fpar[k][0][1]=1.0e-3,  fpar[k][1][0]= -1.25e-8, fpar[k][1][1]=3.67e-3;
   }



  fpar[8][0][0]= -5.7e-7, fpar[8][0][1]=0.81e5,  fpar[8][1][0]= -1.25e-7, fpar[8][1][1]=3.67e5;
  fpar[8][0][0]=-2.0e-8;

  //2018 Oct. 17th
 fpar[8][0][0]= -1.67e-8, fpar[8][0][1]=3.59e2,  fpar[8][1][0]= -1.85e-8, fpar[8][1][1]=3.34e2;
 

  //fpar[16][0][0]= -2.55e-7, fpar[16][0][1]=1.61,  fpar[16][1][0]= -7.40e-7, fpar[16][1][1]=1.0e2;
fpar[16][0][0]= -2.55e-9, fpar[16][0][1]=5.0e3,  fpar[16][1][0]= -3.40e-8, fpar[16][1][1]=5.0e3;
//2018 Oct. 17th
fpar[16][0][0]= -1.08e-8, fpar[16][0][1]=8.17e13,  fpar[16][1][0]= -1.11e-8, fpar[16][1][1]=1.92e15;

  for(int x=0;x<chmax;x++){
    for(int y=0;y<2;y++){
      for(int z=0;z<npara;z++){
       if(x==i && y==j && z==k){
        f_ini_param=fpar[i][j][k];
      }
     }
    }
  }
  return f_ini_param;
}

