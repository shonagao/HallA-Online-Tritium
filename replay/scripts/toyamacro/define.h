/*
  define_hh
*/

#include <iostream>
#include <fstream>

#ifndef define_hh
#define define_hh 1
const int CID_F1S2=1;
const int CID_F1S0=2;
const int CID_RF=3;
const int CID_FbS2=4;
const int CID_FbS0=5;
const int nS2=16;
const int nS0=1;
const int nRF=1;
const int MaxChar= 100;
const double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
const double s2ns = 1.e+9;          // conv. factor sec to nano sec
const double ns2s = 1.e-9;          // conv. factor nano sec to sec

//const double PI = 4.*atan(1.);
//const double deg_to_rad = PI / 180.;
//const double rad_to_deg = 180. / PI;
//const double sigma_to_fwhm = 2.*sqrt(2.*log(2.));
//const double fwhm_to_sigma = 1./sigma_to_fwhm;
const double cm_to_barn = 1e+24;
const double alpha = 1./137.035999074; // fine structure constant
const double hbar = 6.58211928*1e-22;  // Planc constant (reduced) (MeV x s)
const double hbarc = 197.3269718;      // conversion constant (MeV x fm)
const double kb = 8.6173324*1e-5;      // Boltzmann constant
const double e = 1.602176565*1e-19;    // electron charge magnitude (C)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)
const double re = 2.817e-13;           // classical electron radius (cm)
const double Na = 6.02214129*1e+23;    // Avogadro constant
const double Me = 0.510998928;         // electron     mass (GeV/c2)
const double Mmu = 0.1056583715;       // muon         mass (GeV/c2)
const double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
const double Mpi0= 0.1349766;          // charged pion mass (GeV/c2)
const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
const double Mp = 0.938272046;         // proton       mass (GeV/c2)
const double Mn = 0.939565379;         // proton       mass (GeV/c2)
const double Mu = 0.931494061;         // proton       mass (GeV/c2)
const double ML = 1.115683;            // Lambda       mass (GeV/c2)
const double MS0 = 1.192642;           // Sigma Zero   mass (GeV/c2)
const double MSm = 1.197449;           // Sigma Minus  mass (GeV/c2)
const double MSp = 1.18937;            // Sigma Plus   mass (GeV/c2)
#endif
