// 10/4/2018 Author Itabshi
// time walk correction R-HRS S2

//====== Define Function ===========//
double range_para(int i,int j);
double range_adc(int i,int j);
long double fit_ini(int i,int j,int k);

//=================================//



// Grobal Function //
 int chmax=17; // channel of S2 PMT 
 int ch=7;// defolt ch 8;
 double tdc=0.5e-9;//TDC converse ch->sec [ss/ch]


void hrs_S2_wt_correction(){
  gStyle->SetOptStat(111111111);
//-------- TTree data input ---------------//
 TChain *T = new TChain("T");
 for(Int_t i = 93500; i<93504;i++) /// these run goes 93495 to 93504
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));
 }
  
  
 //========== Defined Parameters =======================//
 //========= Time walk parameters =================//
 TH2F* htof_adc[chmax][2];
 TCanvas* cw[chmax];
 TCanvas* cwo[chmax];
 TCanvas* c0[chmax];
 TCanvas* ctdc[chmax];
 TCanvas* ct[chmax];
 TCanvas* ctof[chmax];
 TProfile* tp[chmax][2];
 TF1* ffit[chmax][2];
 TH2F* htdc_adc[chmax][2];
 TH2F* htdc_adc_c[chmax][2];
 TH1F* htdc[chmax][2];
 TH1F* htdc_c[chmax][2];
 TF1* ftdc[chmax][2];
 TF1* ftdc_c[chmax][2];
 double min_tof,max_tof,min_adc,max_adc;
 int bin_tof,bin_adc;
 
 
 double min_tof_c,max_tof_c,min_adc_c,max_adc_c;
 int bin_tof_c,bin_adc_c;

 
 double min_tdc,max_tdc;
 // min_tdc=-0.10e-6;
 //max_tdc=0.22e-6;
 int bin_tdc;
 long double p0[chmax][2],p1[chmax][2];
 long double sig[chmax][2],sig_c[chmax][2];


 char* cut_time[chmax];
 char *adc[chmax][2];
 char *tdc[chmax][2];

 int test=0;
 //================== TOF parameters ========================// 
 TCanvas* ctdc_m[chmax];
 TH1F* htof[chmax];
 TH1F* htof_c[chmax];
 TF1* ftof[chmax];
 TF1* ftof_c[chmax];
 TH1F* htdc_m[chmax];
 TF1* ftdc_m[chmax];
 TF1* ftdc_mc[chmax];
 double tmin,tmax,tcmin,tcmax;
 tmin=-1.0e-8;
 tmax=1.0e-8;
 tcmin=-1.0e-8;
 tcmax=1.0e-8;
 int tbin=400;
 int tcbin=400;

 double pm[chmax],pmc[chmax];

 //============================================================//
 //========== Time Walk correction ============================//
 //============================================================//


 for(int i=0;i<chmax;i++){
   
  
    if(i==ch){
       // Definition TCanvas //
     cwo[i]=new TCanvas(Form("cwo[%d]",i),Form("cwo[%d]",i));
                         cwo[i]->Divide(1,2);
     cw[i]=new TCanvas(Form("cw[%d]",i),Form("cw[%d]",i));
                         cw[i]->Divide(1,2);
     c0[i]=new TCanvas(Form("c0[%d]",i),Form("c0[%d]",i));
                         c0[i]->Divide(1,2);
     ctdc[i]=new TCanvas(Form("ctdc[%d]",i),Form("ctdc[%d]",i));
			 ctdc[i]->Divide(2,2);
     ctof[i]=new TCanvas(Form("ctof[%d]",i),Form("ctof[%d]",i));
			 ctof[i]->Divide(1,2);

   for(int j=0;j<2;j++){


     //======== Parameter Setting ==================//

 min_tof=range_para(i,0);
 max_tof=range_para(i,1);
 min_tof_c=range_para(i,2);
 max_tof_c=range_para(i,3);
 min_adc=range_para(i,4);
 max_adc=range_para(i,5); 
 min_adc_c=range_para(i,4);
 max_adc_c=range_para(i,5);
 bin_tof_c=50;
 bin_tof=50;
 bin_tdc=300;
 bin_adc=500;
 bin_adc_c=500;
       //Cut conditon //
 cut_time[i]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",i,i);
 cut_time[16]="R.s0.lt_c>0 && R.s0.rt_c>0";
 TCut cut_trig="DR.evtypebits==112";
 //TCut cut_trig="DR.evtypebits==16";
 cut_time[i]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",i,i);
 cut_time[16]="R.s0.lt_c>0 && R.s0.rt_c>0";

    // ADC & TDC inf//
 adc[i][0]=Form("R.s2.ra_c[%d]",i);
 adc[i][1]=Form("R.s2.la_c[%d]",i);
 adc[16][0]="R.s0.ra_c";
 adc[16][1]="R.s0.la_c";

 tdc[i][0]=Form("R.s2.rt_c[%d]",i);
 tdc[i][1]=Form("R.s2.lt_c[%d]",i);
 tdc[16][0]="R.s0.rt_c";
 tdc[16][1]="R.s0.lt_c";
 //======================================================//

     cwo[i]->cd(j+1);
     //--- TDC vs ADC 2D Hist --------//
     // htdc_adc[i][j]=new TH2F(Form("htdc_adc[%d][%d]",i,j),Form("TDC vs Scinti CH %d,[%d]-PMT adc w/o correlation ",i,j),bin_adc,min_adc,max_adc,bin_tof,min_tof,max_tof);

     htdc_adc[i][j]=new TH2F(Form("htdc_adc[%d][%d]",i,j),Form("TDC vs Scinti CH %d,[%d]-PMT adc w/o correlation ",i,j),bin_adc,range_para(i,4),range_para(i,5),bin_tof,range_para(i,0),range_para(i,1));
    
     //htdc_adc[i][j]=new TH2F(Form("htdc_adc[%d][%d]",i,j),Form("TDC vs Scinti CH %d,[%d]-PMT adc w/o correlation ",i,j),bin_adc,100.,500.,bin_tof,range_para(i,0),range_para(i,1));
   
 T->Project(Form("htdc_adc[%d][%d]",i,j),Form("%s:%s",tdc[i][j],adc[i][j]),cut_trig && cut_time[i]);//Right PMT
 htdc_adc[i][j]->GetXaxis()->SetRangeUser(range_para(i,4),range_para(i,5));

                       //----- ADC vs TDC w/o correction----------//
                    
			 htdc_adc[i][j]->Draw();


		      		 
 c0[i]->cd(j+1);

 // tp[i][j]=new TProfile(Form("tp[%d][%d]",i,j),Form("TProfile Scinti ch %d [%d]-PMT",i,j),bin_adc,range_para(i,4),range_para(i,5),bin_tof,range_para(i,0),range_para(i,1));
 tp[i][j]=new TProfile(Form("tp[%d][%d]",i,j),Form("TProfile Scinti ch %d [%d]-PMT",i,j),bin_adc,range_para(i,4),range_para(i,5),range_para(i,0),range_para(i,1));

		       //,bin_adc,min_adc,max_adc,min_tof,max_tof);
T->Project(Form("tp[%d][%d]",i,j),Form("%s:%s",tdc[i][j],adc[i][j]),cut_trig && cut_time[i]);


//ffit[i][j]=new TF1(Form("ffit[%d][%d]",i,j),"[0]*(1./sqrt(x)-1./sqrt([1]))",min_adc,max_adc);

 ffit[i][j]=new TF1(Form("ffit[%d][%d]",i,j),"[0]*(1./sqrt(x)-1./sqrt([1]))",range_para(i,4),range_para(i,5));


  ffit[i][j]->SetParameter(0,fit_ini(i,j,0));   
 ffit[i][j]->SetParameter(1,fit_ini(i,j,1));   

 //ffit[i][j]->SetParLimits(0,0.9*fit_ini(i,j,0),1.1*fit_ini(i,j,0));   
 //ffit[i][j]->SetParLimits(1,0.6*fit_ini(i,j,1),1.4*fit_ini(i,j,1));   
 tp[i][j]->Fit(Form("ffit[%d][%d]",i,j),"","");



		      //----- ADC vs TDC TProfile ----------//

			
			 htdc_adc[i][j]->Draw();
			 tp[i][j]->Draw("same");
		     	 ffit[i][j]->Draw("same");

     p0[i][j]=ffit[i][j]->GetParameter(0);
     p1[i][j]=ffit[i][j]->GetParameter(1);
   
    p0[i][j]=p0[i][j]*1.0e9;
     p1[i][j]=p1[i][j]*1.0e5;
    


     //------ Project hist TDC vs ADC w/ Correction ----------------//


     // htdc_adc_c[i][j]=new TH2F(Form("htdc_adc_c[%d][%d]",i,j),Form("TDC vs Scinti CH %d,[%d]-PMT adc w/ correlation ",i,j),bin_adc_c,min_adc_c,max_adc_c,bin_tof_c,min_tof_c,max_tof_c);
     htdc_adc_c[i][j]=new TH2F(Form("htdc_adc_c[%d][%d]",i,j),Form("TDC vs Scinti CH %d,[%d]-PMT adc w/ correlation ",i,j),bin_adc_c,range_para(i,4),range_para(i,5),bin_tof_c,range_para(i,2),range_para(i,3));
      

 T->Project(Form("htdc_adc_c[%d][%d]",i,j),Form("(%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5)):(%s)",tdc[i][j],p0[i][j],adc[i][j],p1[i][j],adc[i][j]),cut_trig && cut_time[i]);

                        
		    //----- ADC vs TDC w/ correction----------//

                         cw[i]->cd(j+1);
			 htdc_adc_c[i][j]->Draw();

                       

 //------- TDC hist w/o correction ---------//
  htdc[i][j]=new TH1F(Form("htdc[%d][%d]",i,j),Form("Scinti CH %d [%d]-PMT TDC hisht w/o correction",i,j),bin_tdc,min_tdc,max_tdc);
    
  T->Project(Form("htdc[%d][%d]",i,j),Form("%s",tdc[i][j]),cut_trig && cut_time[i]);

;


                    //--------- TDC hist --------------//

			 ctdc[i]->cd(j+1);
			 htdc[i][j]->SetFillColor(1);
		
			 //			 htdc[i][j]->GetXaxis()->SetRangeUser(min_tof,max_tof);
                         htdc[i][j]->Draw();




   		    
 //------- TDC hist w/ correction ---------//


                        //--- w/ correction --//
			 ctdc[i]->cd(j+3);


 htdc_c[i][j]=new TH1F(Form("htdc_c[%d][%d]",i,j),Form("Scinti CH %d [%d]-PMT TDC hisht w/ correction",i,j),bin_tdc,min_tdc,max_tdc);

  T->Project(Form("htdc_c[%d][%d]",i,j),Form("(%s)-(%Lf)*1.0e-9*(1./sqrt(%s)-1./sqrt(%Lf*1.0e-5))",tdc[i][j],p0[i][j],adc[i][j],p1[i][j]),cut_trig && cut_time[i]);



			 //	 htdc_c[i][j]->GetXaxis()->SetRangeUser(min_tof_c,max_tof_c);
                	 htdc_c[i][j]->SetFillColor(2);
                         htdc_c[i][j]->Draw();
                        

		
			 

 //--------TDC Fit function ------------------//    
 ftdc[i][j]=new TF1(Form("ftdc[%d][%d]",i,j),"gaus",min_tof,max_tof);
 ftdc_c[i][j]=new TF1(Form("ftdc_c[%d][%d]",i,j),"gaus",min_tdc,max_tdc);
 htdc[i][j]->Fit(Form("ftdc[%d][%d]",i,j),"","",min_tof,max_tof);
 htdc_c[i][j]->Fit(Form("ftdc_c[%d][%d]",i,j),"","",min_tdc,max_tdc);
 sig[i][j]=ftdc[i][j]->GetParameter(2);
 sig_c[i][j]=ftdc_c[i][j]->GetParameter(2);
   













   }


   //================== TDC mean Hist ===========================//

TCut cut_trig="DR.evtypebits==112";
 ctof[i]->cd(1);
 htof[i]=new TH1F(Form("htof[%d]",i),Form("S2 ch%d TDC hist  w/o TW correction",i),100,19e-8,21e-8);
 T->Project(Form("htof[%d]",i),Form("(R.s2.rt_c[%d]+R.s2.lt_c[%d])/2.0",i,i),cut_trig && cut_time[i]);
 htof[i]->Draw();

 htof_c[i]=new TH1F(Form("htof_c[%d]",i),Form("S2 ch%d TDC hist  w/ TW correction",i),100,-1.0e-8,1.0e-8);
 T->Project(Form("htof_c[%d]",i),Form("((R.s2.rt_c[%d]-(%Lf)*1.0e-9*(1./sqrt(R.s2.ra_c[%d])-1./sqrt(%Lf*1.0e-5)))+(R.s2.lt_c[%d]-(%Lf)*1.0e-9*(1./sqrt(R.s2.la_c[%d])-1./sqrt(%Lf*1.0e-5))))/2.0",i,p0[i][0],i,p1[i][0],i,p0[i][1],i,p1[i][1]),cut_trig && cut_time[i]);
 
ctof[i]->cd(2);
 htof_c[i]->Draw();

 
 //================== Fit TDC mean Hits ============================//
/*
 ftdc_m[i]=new TF1(Form("ftdc_m[%d]",i),"gaus",-1.0e-8,1.0e-8);
 ftdc_mc[i]=new TF1(Form("ftdc_mc[%d]",i),"gaus",-1.0e-8,1.0e-8);

 htof[i]->Fit(Form("ftdc_m[%d]",i),"","");
 htof_c[i]->Fit(Form("ftdc_mc[%d]",i),"","");
 pm[i]=ftdc_m[i]->GetParameter(2);
 pmc[i]=ftdc_mc[i]->GetParameter(2);
 */

//----------------- Comment out -----------------//
   cout<<"======== S2 R-PMT ch"<<i<<"parameters ========"<<endl;
   cout<<"p0 :"<<p0[i][0]*1.0e-9<<"  p1 :"<<p1[i][0]*1.0e-5<<endl;    
   cout<<"======== S2 L-PMT ch"<<i<<"parameters ========"<<endl;
   cout<<"p0 :"<<p0[i][1]*1.0e-9<<"  p1 :"<<p1[i][1]*1.0e-5<<endl;    


			
    }
    }
			 
}


 

double range_para(int i,int j){

  int npara=6;// number of parameters
  double par[chmax][npara];
  double param;
  //param 0 min_tdc, 1 max_tdc, 2 min_tdc_c, 3 max_tdc_c //
 //param 4 min_adc, 5 max_adc //
  //=== Inital parameters========//
  for(int k=0;k<chmax;k++){
    par[k][0]=0.18e-6, par[k][1]=0.22e-6, par[k][2]=-0.2e-6, par[k][3]=0.22e-6;//TDC
   if(k==16){//S0
  par[k][4]=500.0, par[k][5]=5000.;//ADC
    }else{ //S2
      par[k][4]=0.0, par[k][5]=500.;//ADC
    }
   }

  //===== Set Parameters ========//
  par[8][0]=0.18e-6, par[8][1]=-0.22e-6, par[8][2]=-0.18e-6,par[8][3]=0.18e-6;
  par[16][0]=0.180e-6, par[16][1]=0.23e-6, par[16][2]=-0.18e-6, par[16][3]=0.25e-6;

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
  fpar[k][0][0]= -6.0e-8, fpar[k][0][1]=1.1e5,  fpar[k][1][0]= -1.25e-8, fpar[k][1][1]=3.67e5;
   }



  fpar[8][0][0]= -5.70e-9, fpar[8][0][1]=0.81e5,  fpar[8][1][0]= -1.25e-8, fpar[8][1][1]=3.67e5;

  //fpar[16][0][0]= -2.55e-7, fpar[16][0][1]=1.61,  fpar[16][1][0]= -7.40e-7, fpar[16][1][1]=1.0e2;
fpar[16][0][0]= -2.55e-7, fpar[16][0][1]=3.61,  fpar[16][1][0]= -3.40e-7, fpar[16][1][1]=1.0e2;

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

