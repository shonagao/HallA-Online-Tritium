// 10/4/2018 Author Itabshi
// time walk correction R-HRS S2

//====== Define Function ===========//
double range_para(int i,int j);
long double fit_ini(int i,int j,int k);

//=================================//



// Grobal Function //
 int chmax=16; // channel of S2 PMT 
 int ch=8;// defolt ch 8;
 double tdc=0.5e-9;//TDC converse ch->sec [ss/ch]


void hrs_tof_tdcwt_correction(){
  //  gStyle->SetOptStat(111111);


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
 min_tof=range_para(ch,0);
 max_tof=range_para(ch,1);
 bin_tof=200;
 min_adc=0.;
 max_adc=500.;
 bin_adc=2500.;

 double min_tof_c,max_tof_c,min_adc_c,max_adc_c;
 int bin_tof_c,bin_adc_c;
 min_tof_c=range_para(ch,2);
 max_tof_c=range_para(ch,3);
 bin_tof_c=2000;
 min_adc_c=0.;
 max_adc_c=500.;
 bin_adc_c=2500.;
 
 double min_tdc,max_tdc;
 min_tdc=-0.10e-6;
 max_tdc=0.22e-6;
 int bin_tdc=500;
 long double p0[chmax][2],p1[chmax][2];
 long double sig[chmax][2],sig_c[chmax][2];


 char* cut_time[chmax];
 cut_time[ch]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",ch,ch);
 TCut cut_trig="DR.evtypebits==112";




 for(int i=0;i<chmax;i++){
     if(i==ch){

   for(int j=0;j<2;j++){
  htdc_adc[i][j]=new TH2F(Form("htdc_adc[%d][%d]",i,j),Form("TDC vs S2 CH %d,[%d] PMT adc w/o correlation ",i,j),bin_adc,min_adc,max_adc,bin_tof,min_tof,max_tof);

 htdc_adc_c[i][j]=new TH2F(Form("htdc_adc_c[%d][%d]",i,j),Form("TDC vs S2 CH %d,[%d] PMT adc w/ correlation ",i,j),bin_adc_c,min_adc,max_adc_c,bin_tof_c,min_tof_c,max_tof_c);
      }

     cwo[i]=new TCanvas(Form("cwo[%d]",i),Form("cwo[%d]",i));
     cw[i]=new TCanvas(Form("cw[%d]",i),Form("cw[%d]",i));
     c0[i]=new TCanvas(Form("c0[%d]",i),Form("c0[%d]",i));
     ctdc[i]=new TCanvas(Form("ctdc[%d]",i),Form("ctdc[%d]",i));

 cut_time[i]=Form("R.s2.lt_c[%d]>0 && R.s2.rt_c[%d]>0",i,i);

     //------ Project hist TDC vs ADC w/o Correction ----------------//
  
     T->Project(Form("htdc_adc[%d][0]",i),Form("R.s2.rt_c[%d]:R.s2.ra_c[%d]",i,i),cut_trig && cut_time[i]);//Right PMT
     T->Project(Form("htdc_adc[%d][1]",i),Form("R.s2.lt_c[%d]:R.s2.la_c[%d]",i,i),cut_trig && cut_time[i]);//Left PMT


     //-------------- TProfile Function-----------------//
     tp[i][0]=new TProfile(Form("tp[%d][0]",i),Form("TProfile S2 ch %d R-PMT",i),bin_adc,min_adc,max_adc,min_tof,max_tof);
   tp[i][1]=new TProfile(Form("tp[%d][1]",i),Form("TProfile S2 ch %d L-PMT",i),bin_adc,min_adc,max_adc,min_tof,max_tof);



  //------- TDC vs ADC Prjection ----------//
   T->Project(Form("tp[%d][0]",i),Form("R.s2.rt_c[%d]:R.s2.ra_c[%d]",i,i),cut_trig &&cut_time[i]);//Right PMT
   T->Project(Form("tp[%d][1]",i),Form("R.s2.lt_c[%d]:R.s2.la_c[%d]",i,i),cut_trig && cut_time[i]);//Left PMT


     //------------ Fit Function -------------------//
     ffit[i][0]=new TF1(Form("ffit[%d][0]",i),"[0]*(1./sqrt(x)-1./sqrt([1]))",min_adc,max_adc);
     ffit[i][1]=new TF1(Form("ffit[%d][1]",i),"[0]*(1./sqrt(x)-1./sqrt([1]))",min_adc,max_adc);
    
     /*   ffit[i][0]->SetParameter(0,-5.0e-6);   
     ffit[i][0]->SetParameter(1,1e-3);
     ffit[i][1]->SetParameter(0,-1.25e-8);   
     ffit[i][1]->SetParameter(1,3.5e-3);
     */
     ffit[i][0]->SetParameter(0,fit_ini(i,0,0));   
     ffit[i][0]->SetParameter(1,fit_ini(i,0,1));
     ffit[i][1]->SetParameter(0,fit_ini(i,1,0));   
     ffit[i][1]->SetParameter(1,fit_ini(i,1,1));
    
     tp[i][0]->Fit(Form("ffit[%d][0]",i),"","");
     tp[i][1]->Fit(Form("ffit[%d][1]",i),"","");
     p0[i][0]=ffit[i][0]->GetParameter(0);
     p0[i][1]=ffit[i][1]->GetParameter(0);
     p1[i][0]=ffit[i][0]->GetParameter(1);
     p1[i][1]=ffit[i][1]->GetParameter(1);
   

     p0[i][0]=p0[i][0]*1.0e9;
     p1[i][0]=p1[i][0]*1.0e5;
     p0[i][1]=p0[i][1]*1.0e9;
     p1[i][1]=p1[i][1]*1.0e5;
     //------ Project hist TDC vs ADC w/ Correction ----------------//
  
     T->Project(Form("htdc_adc_c[%d][0]",i),Form("R.s2.rt_c[%d]-(%Lf)*1.0e-9*(1./sqrt(R.s2.ra_c[%d])-1./sqrt(%Lf*1.0e-5)):R.s2.ra_c[%d]",i,p0[i][0],i,p1[i][0],i),cut_trig && cut_time[i]);//Right PMT
     T->Project(Form("htdc_adc_c[%d][1]",i),Form("R.s2.lt_c[%d]-(%Lf)*1.0e-9*(1./sqrt(R.s2.la_c[%d])-1./sqrt(%Lf*1.0e-5)):R.s2.la_c[%d]",i,p0[i][1],i,p1[i][1],i),cut_trig && cut_time[i]);//Left  PMT
     // T->Project(Form("htdc_adc_c[%d][1]",i),Form("R.s2.lt_c[%d]:R.s2.la_c[%d]",i,i),cut_trig);//Left PMT


     //------- TDC hist w/o correction ---------//
     htdc[i][0]=new TH1F(Form("htdc[%d][0]",i),Form("S2 CH %d R-PMT TDC hisht w/o correction",i),bin_tdc,min_tdc,max_tdc);
     htdc[i][1]=new TH1F(Form("htdc[%d][1]",i),Form("S2 CH %d L-PMT TDC hisht w/o correction",i),bin_tdc,min_tdc,max_tdc);

 T->Project(Form("htdc[%d][0]",i),Form("R.s2.rt_c[%d]",i),cut_trig && cut_time[i]);
 T->Project(Form("htdc[%d][1]",i),Form("R.s2.lt_c[%d]",i),cut_trig && cut_time[i]);   		    
 //------- TDC hist w/ correction ---------//
 htdc_c[i][0]=new TH1F(Form("htdc_c[%d][0]",i),Form("S2 CH%d R-PMT TDC hisht w/ correction",i),bin_tdc,min_tdc,max_tdc);
 htdc_c[i][1]=new TH1F(Form("htdc_c[%d][1]",i),Form("S2 CH%d L-PMT TDC hisht w/ correction",i),bin_tdc,min_tdc,max_tdc);

 T->Project(Form("htdc_c[%d][0]",i),Form("R.s2.rt_c[%d]-%Lf*1.0e-9*(1./sqrt(R.s2.ra_c[%d])-1./sqrt(%Lf*1.0e-5))",i,p0[i][0],i,p1[i][0]),cut_trig && cut_time[i]);
 T->Project(Form("htdc_c[%d][1]",i),Form("R.s2.lt_c[%d]-%Lf*1.0e-9*(1./sqrt(R.s2.la_c[%d])-1./sqrt(%Lf*1.0e-5))",i,p0[i][1],i,p1[i][1]),cut_trig && cut_time[i]);   		      
    


 for(int j=0;j<2;j++){
 ftdc[i][j]=new TF1(Form("ftdc[%d][%d]",i,j),"gaus",min_tof,max_tof);
 ftdc_c[i][j]=new TF1(Form("ftdc_c[%d][%d]",i,j),"gaus",min_tdc,max_tdc);
 htdc[i][j]->Fit(Form("ftdc[%d][%d]",i,j),"","",min_tof,max_tof);
 htdc_c[i][j]->Fit(Form("ftdc_c[%d][%d]",i,j),"","",min_tdc,max_tdc);
 sig[i][j]=ftdc[i][j]->GetParameter(2);
 sig_c[i][j]=ftdc_c[i][j]->GetParameter(2);


 }




  //----------------- Comment out -----------------//
   cout<<"======== S2 R-PMT ch"<<i<<"parameters ========"<<endl;
   cout<<"p0 :"<<p0[i][0]*1.0e-9<<"  p1 :"<<p1[i][0]*1.0e-5<<endl;    
   cout<<"======== S2 L-PMT ch"<<i<<"parameters ========"<<endl;
   cout<<"p0 :"<<p0[i][1]*1.0e-9<<"  p1 :"<<p1[i][1]*1.0e-5<<endl;    


 //---- Comment out TDC correction get fit parameters --------//

 cout<<"========= CH"<<i<<" TDC correction get fit parameters ==========="<<endl;
 cout<<"R-PMT TDC w/o correction sigma "<<sig[i][0]<<endl;
 cout<<"L-PMT TDC w/o correction sigma "<<sig[i][1]<<endl;
 cout<<"R-PMT TDC w/ correction sigma "<<sig_c[i][0]<<endl;
 cout<<"L-PMT TDC w/ correction sigma "<<sig_c[i][1]<<endl;



             //-----------TCanvas -----------------------//
		    
		      //----- ADC vs TDC TProfile ----------//
                         c0[i]->Divide(1,2);
			 c0[i]->cd(1);
			 htdc_adc[i][0]->Draw();
			 tp[i][0]->Draw("same");
			 ffit[i][0]->Draw("same");
                         c0[i]->cd(2);
			 htdc_adc[i][1]->Draw();
			 tp[i][1]->Draw("same");
			 ffit[i][1]->Draw("same");
		
                       //----- ADC vs TDC w/o correction----------//
                         cwo[i]->Divide(1,2);
                         cwo[i]->cd(1);
			 htdc_adc[i][0]->Draw();
                         cwo[i]->cd(2);
			 htdc_adc[i][1]->Draw();
		
     
		    //----- ADC vs TDC w/ correction----------//
                         cw[i]->Divide(1,2);
                         cw[i]->cd(1);
			 htdc_adc_c[i][0]->Draw();
                         cw[i]->cd(2);
			 htdc_adc_c[i][1]->Draw();
		      //--------- TDC hist --------------//
			 ctdc[i]->Divide(2,2);
			 ctdc[i]->cd(1);
			 htdc[i][0]->SetFillColor(1);
			 htdc_c[i][0]->SetFillColor(2);
			 htdc[i][0]->GetXaxis()->SetRangeUser(min_tof,max_tof);
                         htdc[i][0]->Draw();
			 //  htdc_c[i][0]->Draw("same");
                         ctdc[i]->cd(2);
			 htdc[i][1]->SetFillColor(1);
			 htdc_c[i][1]->SetFillColor(2);
			 htdc[i][1]->GetXaxis()->SetRangeUser(min_tof,max_tof);
			 htdc[i][1]->Draw();
		//	 htdc_c[i][1]->Draw("same");
                        //--- w/ correction --//
			 ctdc[i]->cd(3);
			 htdc_c[i][0]->GetXaxis()->SetRangeUser(min_tof_c,max_tof_c);
                         htdc_c[i][0]->Draw();
                         ctdc[i]->cd(4);
			 htdc_c[i][1]->GetXaxis()->SetRangeUser(min_tof_c,max_tof_c);
			 htdc_c[i][1]->Draw();
			 
			 }
 }





}



double range_para(int i,int j){

  int ch=16;//S2 Max Channel
  int npara=4;// number of parameters
  double par[chmax][npara];
  double param;
  //param 0 min_tdc, 1 max_tdc, 2 min_tdc_c, 3 max_tdc_c //
  par[8][0]=0.18e-6, par[8][1]=0.22e-6, par[8][2]=-0.02e-6, par[8][3]=0.02e-6;

  for(int k=0;k<ch;k++){
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
  int ch=16;
  int npara=2;
  double fpar[chmax][2][npara];
  double f_ini_param;

  fpar[8][0][0]= -6.0e-9, fpar[8][0][1]=1.0e-3,  fpar[8][1][0]= -1.25e-8, fpar[8][1][1]=3.67e-3;


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
