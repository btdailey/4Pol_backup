////////////////////////////////////////////////////////////////////////////////////////////////
//class Primaries:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PRIMARIES_H_
#define PRIMARIES_H_

#include "TRandom3.h" 
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TF3.h"
#include "TH2D.h"
#include "Settings.h"
#include "position.hh"


class Vector;
class Position;
class Interaction;
class Primaries;
class IceModel;
class Counting;

//using namespace std;
using std::string;
//
using std::cout;




class Y {
    //Code from Connolly Calc 2011.
private:
    TF1* ffrac; // fraction of events in the low y region.
    TF1* fC1_high[2][2]; // parameterization of C1 in high y region
    
    TF1 *fC1_low; // parameterization of C1 in the low y region.
    
    TF1 *fC2; // parameterization of C2 in the low y region.
    
    TF3* fy0_low; // for picking y0. 
    TF2* fy0_high;  // for picking y0.
    TRandom3 Rand3;
    
public:
    Y();
    double pickY(int NU,int CURRENT,double e); // pick an inelasticity
    // NU=0: nubar, NU=1: nu
    // CURRENT=0: CC, CURRENT-1: NC
    const static double miny_low=0.00002;//2.E-5
    const static double maxy_low=0.001;
    const static double miny_high=0.001;
    const static double maxy_high=1.;
};//Y

class Primaries {
    
private:
    TRandom3 Rand3;
    
    TH2D *m_hsigma;
    TCanvas *m_csigma;
    Y *m_myY;
    int run_old_code;
    
public:
    double ymin_low, ymax_low, ymin_high, ymax_high;
    double A_low[4];//same for any nu_nubar and current type.
    double A0_high[2][2];
    double A1_high[2][2];
    double A2_high[2][2];
    double A3_high[2][2];
    double b0, b1; 
    
    TF1* m_fy[2][2];
    TF1* m_fsigma[2][2];
    
    double c0[2][2];
    double c1[2][2];
    double c2[2][2];
    double c3[2][2];
    double c4[2][2];
    
    static const int NSIGMAS=2;// number of possible cross section models
    // 0=Gandhi et al.
    // 1=Connolly et al. 2011
    double mine[NSIGMAS];// minimum energy for cross section parametrizations, in eV
    double maxe[NSIGMAS]; //max
    
    Primaries();//constructor //default
    ~Primaries();//destructor //default
    //*primary1 must be manually deleted in icemc for deconstructor to actually be called.
    
    int GetSigma(double pnu,double& sigma,double &len_int_kgm2,Settings *settings1,int nu_nubar,int currentint);//not static
    double Gety(Settings *settings1,double pnu,int nu_nubar,int currentint);
    double Getyweight(double pnu,double y,int nu_nubar,int currentint);
    string GetCurrent();
    string GetNuFlavor();
protected:
};//Primaries

class Interaction  {
    
private:
    
    
    
    Vector tmp_banana; //Intermediate vector
    
    //For banana plot
    
    // static const double RADDEG_TMP=3.14159/180.;
    static const double nu_banana_theta_angle=-0.413 * 3.14159/180.;// don't let me use RADDEG which is annoying 
    
    
    static const double altitude_nu_banana=-400.;//Depth of interaction of banana neutrino
    
    
    static const double lat_nu_banana=0.; 
    static const double lon_nu_banana=0.;
    
    
    static const double banana_slopey=0.;//Turn slopyness off for banana plots (SLOPEY)
    static const double nu_banana_phi_angle=0. * 3.14159/180.; 
    
    
    
public:
    
    static const double phi_nu_banana=3.14159/4; //Location in phi
    
    static const double banana_observation_distance=600000.;//How far from the surface above the interaction are we when we measure the voltages? (meters) Note: Should be at least 100000 for best results.
    static const double theta_nu_banana=170.*3.14159/180.;//Location of banana neutrino in theta
    double banana_phi_obs;
    Vector banana_obs; //Vector from the neutrino interaction to the observation point
    Interaction(string inttype,Primaries *primary1,Settings *settings1,int whichray,Counting *count1);
    void PickAnyDirection();
    
    int noway;
    int wheredoesitleave_err;
    int neverseesice;
    int wheredoesitenterice_err;
    int toohigh;
    int toolow;
    
    double pathlength_inice;
    
    Vector nnu;  // direction of neutrino (+z in south pole direction)
    double costheta_nutraject; //theta of nnu with earth center to balloon as z axis 
    double phi_nutraject; //phi of nnu with earth center to balloon as z axis
    double weight_nu;//Weight for neutrino that survives to posnu
    double weight_nu_prob;//Weight for neutrino that survives to posnu and interacts in the ice
    double weight_tau;//Weight for a tau neutrino to interact, create a tau, and the tau to survive to the ice. NOT USED
    double weight_tau_prob;//Weight for tau neutrino to interact, create a tau, tau survives and decays in the ice.
  
    
    Position r_in; // position where neutrino enters the earth
    Position r_enterice; // position where neutrino enters the ice
    Position nuexit; // place where neutrino would have left the earth
    Position nuexitice; // place where neutrino would have left the ice
    double chord;  // chord in m from earth entrance to rock-ice boundary
    double logchord; // log_10 of chord length earth entrance to where it enters ice
    double weight_bestcase; // what weight1 would be if whole earth had density of crust - for quick and dirty calculation of best case scenario
    double chord_kgm2_bestcase; // the chord the neutrino would traverse if it all was crust density
    double chord_kgm2_ice; // from ice entrance to interaction point
    double d1;  //same as chord in m (earth entrance to rock-ice boundary)
    double d2;  // ice-rock boundary to interaction point in m
    
    
    static const double pnu_banana=2.00E19;
    static const double banana_y=0.2;//Elasticity.  0.2 is an average number.
    double banana_weight;//Weight measurement locations to account for phase space
    double banana_theta_obs;
    double banana_volts;//Total voltage measured at a spot on the sky
    static const double banana_signal_fluct=0.;//Turn off noise for banana plots (settings1->SIGNAL_FLUCT) (shouldn't matter)
    static const double banana_sigma=0.;//NSIGMA in the case of a banana plot
    
    
    void  setNuFlavor(Primaries *primary1,Settings *settings1,int whichray,Counting *count1);
    void setCurrent(Primaries *primary1);
    Position posnu;
    Position posnu_down;
    string  nuflavor;                   // neutrino flavor
    string  current;                    //  CC or NC?
    int nuflavorint;                // Added by Stephen for output purposes
    int currentint;                 // Ditto - Stephen
    
    
    double surface_over_banana_nu;
    
    string banana_flavor; //Force interaction to be a muon neutrino
    string banana_current;  //Force interaction to be a neutral current
    
    Vector nnu_banana; //Forced neutrino direction 
    
    Position nu_banana;  //The forced interaction point of the neutrino for the banana plots
    Position nu_banana_surface; //The location of the surface above the forced neutrino interaction point  
    
    // phase space weighting
    double dtryingdirection; //weighting factor: how many equivalent tries each neutrino counts for after having reduced angular phase space for possibly detectable events
    double dnutries; //product of dtryingdirection and dtryingposition
    
    double altitude_int;// depth of interaction
    double altitude_int_mirror;//depth of the mirror point of interaction.
    
    double r_fromballoon[2]; // distance from interaction to balloon for each ray
    
    double r_fromballoon_db; // same, for double bangs
    double r_exit2bn; // exit to balloon
    double r_exit2bn_measured; // exit to balloon deduced from measured theta
    int iceinteraction;// whether or not there is an interaction in the ice
    
    
protected:
};//Interaction
#endif
