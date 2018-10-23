// 10/4/2018 Author Itabshi
// time walk correction R-HRS S0

//====== Define Function ===========//
double range_para(int i,int j);
long double fit_ini(int i,int j,int k);

//=================================//



// Grobal Function //
 int chmax=17; // channel of S2 PMT + S0 PMT (16) 
 int ch=8;// defolt ch 8;
 double tdc=0.5e-9;//TDC converse ch->sec [ss/ch]


void hrs_s0_wt_correction(){
   gStyle->SetOptStat(111111111);


//-------- TTree data input ---------------//
 TChain *T = new TChain("T");
 for(Int_t i = 93500; i<93504;i++) /// these run goes 93495 to 93504
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));
 }
  
  
 //========== Defined Parameters =======================//

 TH2F* htof_adc[chmax][2];
 TCanvas* cw[chmax];
 TCanvas* cwo[chmax];
 TCanvas* c0[chmax];
 TCanvas* ctdc[chmax];
 TCanvas* ctof[chmax];
 TCanvas* ctdcm_adc[chmax];
 TProfile* tp[chmax][2];
 TF1* ffit[chmax][2];
 TH2F* htdc_adc[chmax][2];
 TH2F* htdc_adc_c[chmax][2];
 TH1F* htdc[chmax][2];
 TH1F* htdc_c[chmax][2];
 TH1F* htof[chmax];
 TH1F* htof_c[chmax];
 TF1* ftdc[chmax][2];
 TF1* ftdc_c[chmax][2];
 TH2F* htdcm_adc[chmax][2];
 TH2F* htdcm_adc_c[chmax][2];
 double min_tof,max_tof,min_adc,max_adc;
 int bin_tof,bin_adc;
 min_tof=range_para(ch,0);
 max_tof=range_para(ch,1);


 //S0 condition 
 min_adc=500.;
 max_adc=3000.;
 bin_adc=1000.;

 double min_tof_c,max_tof_c,min_adc_c,max_adc_c;
 int bin_tof_c,bin_adc_c;
 min_tof_c=range_para(ch,2);
 max_tof_c=range_para(ch,3);
 bin_tof_c=200;
 min_adc_c=500.;
 max_adc_c=5000.;
 bin_adc_c=1000.;
 
 double min_tdc,max_tdc;
 min_tdc=-0.10e-6;
 max_tdc=0.22e-6;
 int bin_tdc=500;
 long double p0[chmax][2],p1[chmax][2];
 long double sig[chmax][2],sig_c[chmax][2];


 char* cut_time[chmax];
 cut_time[ch]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",ch,ch);
 cut_time[16]="R.s0.lt_c>0 && R.s0.rt_c>0";
 //cut_s0time="R.s0.lt_c>0 && R.s0.rt_c>0";
 TCut cut_trig="DR.evtypebits==112";


 //========================================================//
 //======== S2 time walk correction========================//
 //========================================================//

 for(int i=0;i<chmax;i++){
     if(i==16){

   for(int j=0;j<2;j++){
 
  htdc_adc[16][j]=new TH2F(Form("htdc_adc[16][%d]",j),Form("TDC vs S0 ,[%d] PMT adc w/o correlation ",j),bin_adc,min_adc,max_adc,bin_tof,range_para(16,0),range_para(16,1));
 htdc_adc_c[16][j]=new TH2F(Form("htdc_adc_c[16][%d]",j),Form("TDC vs S2 ,[%d] PMT adc w/ correlation ",j),bin_adc_c,min_adc,max_adc_c,bin_tof_c,range_para(16,2),range_para(16,3));

   }



 cut_time[16]="R.s0.lt_c>0 && R.s0.rt_c>0";
   


  //------ Project hist TDC vs ADC w/o Correction ----------------//

     T->Project("htdc_adc[16][0]","R.s0.rt_c:R.s0.ra_c",cut_trig && cut_time[16]);//Right PMT

  T->Project("htdc_adc[16][1]","R.s0.lt_c:R.s0.la_c",cut_trig && cut_time[16]);//Left PMT





     


			 



 
     //-------------- TProfile Function-----------------//

   tp[16][0]=new TProfile("tp[16][0]","TProfile S0 R-PMT",bin_adc,min_adc,max_adc,range_para(16,0),range_para(16,1));
   tp[16][1]=new TProfile("tp[16][1]","TProfile S0 L-PMT",bin_adc,min_adc,max_adc,range_para(16,0),range_para(16,1));

  //------- TDC vs ADC Prjection ----------//

 T->Project("tp[16][0]","R.s0.rt_c:R.s0.ra_c",cut_trig &&cut_time[i]);//Right PMT
    T->Project("tp[16][1]","R.s0.lt_c:R.s0.la_c",cut_trig &&cut_time[i]);//Left PMT
     //------------ Fit Function -------------------//

     ffit[16][0]=new TF1("ffit[16][0]","[0]*(1./sqrt(x)-1./sqrt([1]))",range_para(16,0),range_para(16,1));
     ffit[16][1]=new TF1("ffit[16][1]","[0]*(1./sqrt(x)-1./sqrt([1]))",range_para(16,0),range_para(16,1));








	

     //==== S0 ======//
    


       ffit[16][0]->SetParLimits(0,-5.0e-7,-1.0e-7);   
       ffit[16][0]->SetParLimits(1,1.0e0,1.0e1);
       ffit[16][1]->SetParLimits(0,-5.0e-7,-2.5e-7);   
       ffit[16][1]->SetParLimits(1,1.0e-3,1.0e4);

    
     tp[16][0]->Fit("ffit[16][0]","","");
     tp[16][1]->Fit("ffit[16][1]","","");
     p0[16][0]=ffit[16][0]->GetParameter(0);
     p0[16][1]=ffit[16][1]->GetParameter(0);
     p1[16][0]=ffit[16][0]->GetParameter(1);
     p1[16][1]=ffit[16][1]->GetParameter(1);
   

     p0[16][0]=p0[i][0]*1.0e9;
     p1[16][0]=p1[i][0]*1.0e5;
     p0[16][1]=p0[i][1]*1.0e9;
     p1[16][1]=p1[i][1]*1.0e5;
    



     


 //------ Project hist TDC vs ADC w/ Correction ----------------//
  
     T->Project("htdc_adc_c[16][0]",Form("R.s0.rt_c-(%Lf)*1.0e-9*(1./sqrt(R.s0.ra_c)-1./sqrt(%Lf*1.0e-5)):R.s0.ra_c",p0[16][0],p1[16][0]),cut_trig && cut_time[16]);//Right PMT
     T->Project("htdc_adc_c[16][1]",Form("R.s0.lt_c-(%Lf)*1.0e-9*(1./sqrt(R.s0.la_c)-1./sqrt(%Lf*1.0e-5)):R.s0.la_c",p0[16][1],p1[16][1]),cut_trig && cut_time[16]);//Right PMT
    



     //------- TDC hist w/o correction ---------//
     htdc[16][0]=new TH1F("htdc[16][0]","S0 R-PMT TDC hisht w/o correction",bin_tdc,range_para(16,0),range_para(16,1));
     htdc[16][1]=new TH1F("htdc[16][1]","S0 L-PMT TDC hisht w/o correction",bin_tdc,range_para(16,0),range_para(16,1));

 T->Project("htdc[16][0]","R.s0.rt_c",cut_trig && cut_time[16]);   	        
 T->Project("htdc[16][1]","R.s0.lt_c",cut_trig && cut_time[16]);   		    

		


 //------- TDC hist w/ correction ---------//

 htdc_c[16][0]=new TH1F("htdc_c[16][0]","S0 R-PMT TDC hisht w/ correction",bin_tdc,range_para(16,2),range_para(16,3));
 htdc_c[16][1]=new TH1F("htdc_c[16][1]","S0 L-PMT TDC hisht w/ correction",bin_tdc,range_para(16,2),range_para(16,3));

 T->Project("htdc_c[16][0]",Form("R.s0.rt_c-(%Lf)*1.0e-9*(1./sqrt(R.s0.ra_c)-1./sqrt(%Lf*1.0e-5))",p0[16][0],p1[16][0]),cut_trig && cut_time[16]);
 T->Project("htdc_c[16][1]",Form("R.s0.lt_c-(%Lf)*1.0e-9*(1./sqrt(R.s0.la_c)-1./sqrt(%Lf*1.0e-5))",p0[16][1],p1[16][1]),cut_trig && cut_time[16]);

 


 ftdc[16][0]=new TF1("ftdc[16][0]","gaus",range_para(16,0),range_para(16,1));
 ftdc_c[16][0]=new TF1("ftdc_c[16][0]","gaus",range_para(16,2),range_para(16,3));
  htdc[16][0]->Fit("ftdc[16][0]","","",range_para(16,0),range_para(16,1));
  htdc_c[16][0]->Fit("ftdc_c[16][0]","","",range_para(16,2),range_para(16,3));
  sig[16][0]=ftdc[16][0]->GetParameter(2);
 sig_c[16][0]=ftdc_c[16][0]->GetParameter(2);


 ftdc[16][1]=new TF1("ftdc[16][1]","gaus",range_para(16,0),range_para(16,1));
 ftdc_c[16][1]=new TF1("ftdc_c[16][1]","gaus",range_para(16,2),range_para(16,3));
  htdc[16][1]->Fit("ftdc[16][1]","","",range_para(16,0),range_para(16,1));
  htdc_c[16][1]->Fit("ftdc_c[16][1]","","",range_para(16,2),range_para(16,3));
  sig[16][1]=ftdc[16][1]->GetParameter(2);
 sig_c[16][1]=ftdc_c[16][1]->GetParameter(2);




 htof[i]=new TH1F(Form("htof[%d]",i),"S0 TOF w/o TW correction",200,0.19e-6,0.210e-6);

 T->Project(Form("htof[%d]",i),"(R.s0.rt_c+R.s0.lt_c)/2.0",cut_trig && cut_time[16]);

 htof_c[i]=new TH1F(Form("htof_c[%d]",i),"S0 TOF w/ TW correction",200,-1.0e-8,1.0e-8);

 T->Project(Form("htof_c[%d]",i),Form("((R.s0.rt_c-(%Lf)*1.0e-9*(1./sqrt(R.s0.ra_c)-1./sqrt(%Lf*1.0e-5)))+(R.s0.lt_c-(%Lf)*1.0e-9*(1./sqrt(R.s0.la_c)-1./sqrt(%Lf*1.0e-5))))/2.0",p0[16][0],p1[16][0],p0[16][1],p1[16][1]),cut_trig && cut_time[16]);


 //====== TDC mean : ADC hist ==============//
 htdcm_adc_c[i][0]=new TH2F(Form("htdc_adc_c[%d][0]",i),"TDC mean : ADC R-PMT w/o correciton",200,-1.0e-8,1.0e-8,4500,500,5000);

 T->Project(Form("htdcm_adc_c[%d][0]",i),Form("((R.s0.rt_c-(%Lf)*1.0e-9*(1./sqrt(R.s0.ra_c)-1./sqrt(%Lf*1.0e-5)))+(R.s0.lt_c-(%Lf)*1.0e-9*(1./sqrt(R.s0.la_c)-1./sqrt(%Lf*1.0e-5))))/2.0:R.s0.ra_c",p0[16][0],p1[16][0],p0[16][1],p1[16][1]),cut_trig && cut_time[16]);


 ctdcm_adc[i]=new TCanvas(Form("ctdcm_adc[%d]",i),Form("ctdcm_adc[%d]",i));
 ctdcm_adc[i]->cd(1);
 htdcm_adc_c[i][0]->Draw();








  //----------------- Comment out -----------------//
   cout<<"======== S0 R-PMT parameters ========"<<endl;
   cout<<"p0 :"<<p0[16][0]*1.0e-9<<"  p1 :"<<p1[16][0]*1.0e-5<<endl;    
   cout<<"======== S0 L-PMT parameters ========"<<endl;
   cout<<"p0 :"<<p0[16][1]*1.0e-9<<"  p1 :"<<p1[16][1]*1.0e-5<<endl;    


 //---- Comment out TDC correction get fit parameters --------//



   
 cout<<"========= S0 CH TDC correction get fit parameters ==========="<<endl;
 cout<<"R-PMT TDC w/o correction sigma "<<sig[16][0]<<endl;
 cout<<"L-PMT TDC w/o correction sigma "<<sig[16][1]<<endl;
 cout<<"R-PMT TDC w/ correction sigma "<<sig_c[16][0]<<endl;
 cout<<"L-PMT TDC w/ correction sigma "<<sig_c[16][1]<<endl;
 
   			 
     }
 }


  

}



double range_para(int i,int j){

  int npara=4;// number of parameters
  double par[chmax][npara];
  double param;
  //param 0 min_tdc, 1 max_tdc, 2 min_tdc_c, 3 max_tdc_c //

  // initial parameters //
  for(int k=0;k<chmax;k++){
par[k][0]=-2.00e-6, par[k][1]=2.0e-6, par[k][2]=-2e-6, par[k][3]=2.e-6;
  }


  par[8][0]=0.18e-6, par[8][1]=0.22e-6, par[8][2]=-0.02e-6, par[8][3]=0.02e-6;
  par[16][0]=0.180e-6, par[16][1]=0.23e-6, par[16][2]=-1e-8, par[16][3]=1.e-8;
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

  // initial parameters //
 

  for(int x=0;x<chmax;x++){
       fpar[k][0][0]= -6.0e-9, fpar[k][0][1]=1.0e-3,  fpar[k][1][0]= -1.25e-8, fpar[k][1][1]=3.67e-3;
  }


  fpar[8][0][0]= -6.0e-9, fpar[8][0][1]=1.0e-3,  fpar[8][1][0]= -1.25e-8, fpar[8][1][1]=3.67e-3;


  fpar[16][0][0]= -2.55e-7, fpar[16][0][1]=1.61,  fpar[16][1][0]= -7.40e-7, fpar[16][1][1]=1.0e2;

  for(int x=0;x<ch;x++){
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


double peak(int i,int j){

  char* cut_time[chmax];
  double xmin,xmax;
  xmin=0.0;
  xmax=10000.;
  int xbin;
  xbin=10000;
  TH1F* hist[chmax];
  // TProfile* TP[chmax];
 


  hist[i]=new TH1F(Form("hist[%d]",i),Form("hist[%d]",i),xbin,xmin,xmax);
  // TP[i]=new TProfile(Form("TP[%d]",i),"",);
 cut_time[i]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",i,i);
 TCut cut_trig="DR.evtypebits==112";
 //  T->Project(Form("hist[%d]",i), Form("R.s2.rt_c[%d]",i), cut_trig && cut_time[i]);




  return 0;
}
 // cut_time[ch]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",ch,ch);
