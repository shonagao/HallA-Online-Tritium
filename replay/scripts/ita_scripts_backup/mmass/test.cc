const double c=299792458e-9;// [m/ns]
const double mk=493.7e-3;// Kaon mass [GeV/c^2]
const double me=0.511e-3;// electron mass [GeV/c^2] 
const double ml=1115.7e-3;//Lambda mass [GeV/c^2]
const double mn=939.6e-3; // neutron mass [GeV/c^2]
const double mpi=139.6e-3;// pion mass [GeV/c^2]
const double tdc_time=56.23e-3;//[ns]
const double Ee=4.3;// [GeV]
const double mtr=938.27e-3;// proton mass [GeV/c^2]
const double PI=3.14;
void test(){
  double mh,mh_c;
  double pk=1.82;// GeV
  double pe_=2.1;
  double pe,theta_r,theta_l,rad_r,rad_l,theta_lc,theta_rc,rad_rc,rad_lc;

  double Ee_=sqrt(pow(pe_,2)+pow(me,2));
  double Ek=sqrt(pow(pk,2)+pow(mk,2));
  theta_r=0.0;
  theta_l=0.0; 
  rad_r=theta_r*PI/180.;
  rad_l=theta_l*PI/180.;
    
  theta_rc=13.2;
  theta_lc=-13.2;
  rad_rc=theta_rc*PI/180.;
  rad_lc=theta_lc*PI/180.;

  pe=sqrt(pow(Ee,2)-pow(me,2));
  mh=sqrt(pow(Ee+mtr-Ee_-Ek,2)-(pow(pe,2)+pow(pe_,2)+pow(pk,2)-2*pe*pe_*cos(rad_r)-2*pe*pk*cos(rad_l)+2*pk*pe_*cos(rad_r-rad_l)));

mh_c=sqrt(pow(Ee+mtr-Ee_-Ek,2)-(pow(pe,2)+pow(pe_,2)+pow(pk,2)-2*pe*pe_*cos(rad_lc)-2*pe*pk*cos(rad_rc)+2*pk*pe_*cos(rad_rc-rad_lc)));

 cout<<"Missing mass w/o deg Correction :"<<mh<<endl;
 cout<<"Missing mass w/  deg Correction :"<<mh_c<<endl;



}
