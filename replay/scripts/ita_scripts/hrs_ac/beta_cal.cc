
// beta calculation macro //

void beta_cal(){

  double mk,mp,mpi; //[MeV]
  double pk,beta_p,beta_pi,beta_k;//[MeV]
  mk=493.677;
  mpi=139.57;
  mp=938.272;
  
  pk=1800.;
  beta_p=pk/sqrt(pow(pk,2)+pow(mp,2));
  beta_pi=pk/sqrt(pow(pk,2)+pow(mpi,2));
  beta_k=pk/sqrt(pow(pk,2)+pow(mk,2));

  cout<<endl<<"==== beta inforamtion ======="<<endl;
  cout<<"Kaon mom [MeV] is "<<pk<<endl;
  cout<<"Kaon beta is "<<beta_k<<endl;
  cout<<"Pion beta is "<<beta_pi<<endl;
  cout<<"Proton beta is "<<beta_p<<endl;
}
