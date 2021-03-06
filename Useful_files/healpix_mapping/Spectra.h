#ifndef SPECTRA_H_
#define SPECTRA_H_

#include "TSpline.h"
#include <string>
#include "TRandom3.h"


using namespace std;

class Spectra {

private:
  TRandom3 Rand3;
  double maxflux;   // max flux value
//  static const int NSPECTRA_MAX=300;  // why need this??
  static const int E_bin_max = 50;
  int E_bin;   // initialize # of energy bins (max = 50)
  
//  double energy[E_bin_max]; // energies that correspond to the fluxes in the previous array  
//  double EdNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E*dN/dE/dA/dt
//  double E2dNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E^2*dN/dE/dA/dt
  
  void GetFlux(string filename);    // read neutrino flux EdNdEdAdt (in GeV) from filename file

  TGraph *gEdNdEdAdt;   //graph for EdNdEdAdt flux
  TGraph *gE2dNdEdAdt;  //graph for E2dNdEdAdt flux

  TSpline3 *sEdNdEdAdt; //spline of EdNdEdAdt
  TSpline3 *sE2dNdEdAdt;    //spline of E2dNdEdAdt
  int EXPONENT; // set flux model

public:  

  double energy[E_bin_max]; // energies that correspond to the fluxes in the previous array  
  double EdNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E*dN/dE/dA/dt
  double E2dNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E^2*dN/dE/dA/dt
  
  Spectra(int EXPONENT); // constructor  
  
  double GetNuEnergy(); // get the neutrino energy which follows neutrino flux. 

  TGraph *GetGEdNdEdAdt();
  TGraph *GetGE2dNdEdAdt();

  TSpline3 *GetSEdNdEdAdt();
  TSpline3 *GetSE2dNdEdAdt();

  double *Getenergy();
  double *GetEdNdEdAdt();
  double *GetE2dNdEdAdt();
  double GetEdNdEdAdt(double E_val);    // return flux value from TSpline
  double GetE2dNdEdAdt(double E_val);   // return flux value from TSpline

  double Getmaxflux();
  
  int GetE_bin();   // return energy bin number


  int IsSpectrum(); // return 1 or 0 depend on EXPONENT value
  int IsMonoenergetic();    // return 1 or 0 depend of EXPONENT value

  // destructor

}; //class Spectra

#endif
