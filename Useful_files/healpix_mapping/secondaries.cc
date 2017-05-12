#include "vector.hh"
#include "TRandom3.h"
#include "Settings.h"
#include "vector.hh"
#include "position.hh"
#include "signal.hh"
#include "Primaries.h"
#include "secondaries.hh"
#include "icemodel.hh"
#include "Tools.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "TH1F.h"
#include "Constants.h"
#include "Settings.h"
#include "TTreeIndex.h"
#include "TChain.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h" 
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"

using std::cout;
using std::stringstream;
using std::setprecision;
using std::accumulate;
using std::max_element;
using std::partial_sum;
using std::max;

 Secondaries::Secondaries() {
	//For Total Tau Survival probability equation
	//n.b. not in SI units.
	////from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
	//Equation 16  &  used in Equation 30. 
	/////////////////////|Units////|Description////////////////////////
	B0=1.2*pow(10,-7); //| m^2/kg  |
	B1=0.16*pow(10,-7);//| m^2/kg  | }parameterization using a logarithmic dependence on energy for B,
	E0=pow(10,19);     //| eV      | the tau elecromagnetic energy loss parameter.
	mT=1.777E9;	   //| eV      |Mass of Tau
	cT=0.00008693; 	   //| m       |Tau Decay length (86.93 microMeters)
	                   //|         |
	Mn=1.672622E-24;   //| g       |nucleon/ proton mass in grams,also equal to 0.938 GeV. 
	A=1.;              //| none    |constant that sets the total probability to unity
	//these last two constanst from Connolly Calc 2011, used in d_dzPsurvNu().
  
  flavors[0]="nue";
  flavors[1]="numu";
  flavors[2]="nutau"; // the gps path of the anita-lite flight

  SECONDARIES=1; // include secondary interactions
  TAUDECAY=1; // include secondary interactions
  // This is just the initialization, it is set in ReadInputs


    // reading in tauola data file for tau decays
   tauolainfile.open("data/tau_decay_tauola.dat",ifstream::in);

  InitTauola();
  
  TAUFRAC=.5; //fraction of tau neutrino-cc current events where the primare interaction point is the first bang   

 

  count_nfb=0;
  secondary_e_noncons=0;

  for (int i=0;i<7;i++) {
    Tools::Zero(dsdy_muon_brems[i],NPROB_MAX);
    Tools::Zero(dsdy_muon_epair[i],NPROB_MAX);
    Tools::Zero(dsdy_muon_pn[i],NPROB_MAX);
    
    Tools::Zero(y_muon_brems[i],NPROB_MAX);
    Tools::Zero(y_muon_epair[i],NPROB_MAX);
    Tools::Zero(y_muon_pn[i],NPROB_MAX);
    
    Tools::Zero(dsdy_tauon_brems[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_epair[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_pn[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_hadrdecay[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_edecay[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_mudecay[i],NPROB_MAX);
  


    Tools::Zero(y_tauon_brems[i],NPROB_MAX);
    Tools::Zero(y_tauon_epair[i],NPROB_MAX);
    Tools::Zero(y_tauon_pn[i],NPROB_MAX);
    Tools::Zero(y_tauon_hadrdecay[i],NPROB_MAX);
    Tools::Zero(y_tauon_edecay[i],NPROB_MAX);
    Tools::Zero(y_tauon_mudecay[i],NPROB_MAX);
  } //for (Tools::Zeroing)

  Tools::Zero(int_muon_brems,7);
  Tools::Zero(int_muon_epair,7);
  Tools::Zero(int_muon_pn,7);
  
  Tools::Zero(int_tauon_brems,7);
  Tools::Zero(int_tauon_epair,7);
  Tools::Zero(int_tauon_pn,7);
  Tools::Zero(int_tauon_hadrdecay,7);
  Tools::Zero(int_tauon_edecay,7);
  Tools::Zero(int_tauon_mudecay,7);

  // Read probability distributions for secondary interactions

  ReadSecondaries();

	
	
}//Secondaries Constructor

 void Secondaries::readData(string nuflavor,string secndryType, double (*y)[NPROB_MAX], double (*dsdy)[NPROB_MAX])
{
  
  stringstream senergy;
  
  ifstream ifile;
  string suffix=".vec";
  if(nuflavor=="tauon")
    suffix="_tau.vec";
  
  for(int index=0;index<7;index++)
    {senergy.str("");
      double energy=18+0.5*index;
      int precision=(index%2==0)?2:3;
      senergy << setprecision(precision) << energy;
      string path="secondary/"+nuflavor+"/dsdy_"+secndryType+"_1e"+senergy.str()+suffix;
      //cout << "openning file " << path.c_str() << endl;
      ifile.open(path.c_str());
      NPROB=0;
      while(!ifile.eof())
	{
	  ifile >> y[index][NPROB] >> dsdy[index][NPROB];
	  NPROB++;
	  if(NPROB>=NPROB_MAX)
	    {
	      // cerr << " ERROR in reading in y_muon_brem. \n";
	      break;
	    }
	 
	}
      ifile.close();
    }
  
}

 void Secondaries::ReadSecondaries() {
  // reading in data for secondary interactions
  
  cout<<"Reading in data on secondary interactions.\n";

  readData("muons","brems",y_muon_brems,dsdy_muon_brems);
  readData("muons","epair",y_muon_epair,dsdy_muon_epair);
  readData("muons","pn",y_muon_pn,dsdy_muon_pn);
  readData("tauon","brems",y_tauon_brems,dsdy_tauon_brems);
  readData("tauon","epair",y_tauon_epair,dsdy_tauon_epair);
  readData("tauon","pn",y_tauon_pn,dsdy_tauon_pn);
  readData("tauon","hadrdecay",y_tauon_hadrdecay,dsdy_tauon_hadrdecay);
  readData("tauon","edecay",y_tauon_edecay,dsdy_tauon_edecay);
  readData("tauon","mudecay",y_tauon_mudecay,dsdy_tauon_mudecay);
  //cout << "NPROB=" << NPROB << ",  NPROB_MAX=" << NPROB_MAX << endl;
 for(int j=0;j<7;j++) {
    // integrating prob. distributions.
    int_muon_brems[j]=accumulate(dsdy_muon_brems[j],dsdy_muon_brems[j]+NPROB_MAX,0.);//very important to keep the initial value the same type as the elements type
    int_muon_epair[j]=accumulate(dsdy_muon_epair[j],dsdy_muon_epair[j]+NPROB_MAX,0.);
    int_muon_pn[j]=accumulate(dsdy_muon_pn[j],dsdy_muon_pn[j]+NPROB_MAX,0.);
    int_tauon_brems[j]=accumulate(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+NPROB_MAX,0.);
    int_tauon_epair[j]=accumulate(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+NPROB_MAX,0.);
    int_tauon_pn[j]=accumulate(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+NPROB_MAX,0.);
    int_tauon_hadrdecay[j]=accumulate(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+NPROB_MAX,0.);
    int_tauon_edecay[j]=accumulate(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+NPROB_MAX,0.);
    int_tauon_mudecay[j]=accumulate(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+NPROB_MAX,0.);
    
    // maximum value of prob. dist.
    max_muon_brems=*max_element(dsdy_muon_brems[j],dsdy_muon_brems[j]+NPROB_MAX);
    //cout << "max_muon_brems=" << max_muon_brems << endl;//fenfang
    max_muon_epair=*max_element(dsdy_muon_epair[j],dsdy_muon_epair[j]+NPROB_MAX);
    max_muon_pn=*max_element(dsdy_muon_pn[j],dsdy_muon_pn[j]+NPROB_MAX);   
    max_tauon_brems=*max_element(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+NPROB_MAX);
    max_tauon_epair=*max_element(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+NPROB_MAX);
    max_tauon_pn=*max_element(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+NPROB_MAX);
    max_tauon_hadrdecay=*max_element(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+NPROB_MAX);
    max_tauon_edecay=*max_element(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+NPROB_MAX);
    max_tauon_mudecay=*max_element(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+NPROB_MAX);
     
    // minimum value of prob. dist.
    min_muon_brems=Tools::dMinNotZero(dsdy_muon_brems[j],NPROB_MAX);
    min_muon_epair=Tools::dMinNotZero(dsdy_muon_epair[j],NPROB_MAX);
    min_muon_pn=Tools::dMinNotZero(dsdy_muon_pn[j],NPROB_MAX);   
    min_tauon_brems=Tools::dMinNotZero(dsdy_tauon_brems[j],NPROB_MAX);
    min_tauon_epair=Tools::dMinNotZero(dsdy_tauon_epair[j],NPROB_MAX);
    min_tauon_pn=Tools::dMinNotZero(dsdy_tauon_pn[j],NPROB_MAX);
    min_tauon_hadrdecay=Tools::dMinNotZero(dsdy_tauon_hadrdecay[j],NPROB_MAX);
    min_tauon_edecay=Tools::dMinNotZero(dsdy_tauon_edecay[j],NPROB_MAX);
    min_tauon_mudecay=Tools::dMinNotZero(dsdy_tauon_mudecay[j],NPROB_MAX);
     
    if (min_muon_brems<=0)
      cout << "Minimum probability is <=0!\n";
    
    partial_sum(dsdy_muon_brems[j],dsdy_muon_brems[j]+NPROB_MAX,y_cumulative_muon_brems[j]);
    partial_sum(dsdy_muon_epair[j],dsdy_muon_epair[j]+NPROB_MAX,y_cumulative_muon_epair[j]);
    partial_sum(dsdy_muon_pn[j],dsdy_muon_pn[j]+NPROB_MAX,y_cumulative_muon_pn[j]);
    partial_sum(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+NPROB_MAX,y_cumulative_tauon_brems[j]);
    partial_sum(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+NPROB_MAX,y_cumulative_tauon_epair[j]);
    partial_sum(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+NPROB_MAX,y_cumulative_tauon_pn[j]);
    partial_sum(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+NPROB_MAX,y_cumulative_tauon_hadrdecay[j]);
    partial_sum(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+NPROB_MAX,y_cumulative_tauon_mudecay[j]);
    partial_sum(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+NPROB_MAX,y_cumulative_tauon_edecay[j]);
     
    for (int i=0;i<NPROB_MAX;i++) {
       y_cumulative_muon_brems[j][i]      /= y_cumulative_muon_brems[j][NPROB_MAX-1];
       y_cumulative_muon_epair[j][i]      /= y_cumulative_muon_epair[j][NPROB_MAX-1];
       y_cumulative_muon_pn[j][i]         /= y_cumulative_muon_pn[j][NPROB_MAX-1];
       y_cumulative_tauon_brems[j][i]     /= y_cumulative_tauon_brems[j][NPROB_MAX-1];
       y_cumulative_tauon_epair[j][i]     /= y_cumulative_tauon_epair[j][NPROB_MAX-1];
       y_cumulative_tauon_pn[j][i]        /= y_cumulative_tauon_pn[j][NPROB_MAX-1];
       y_cumulative_tauon_hadrdecay[j][i] /= y_cumulative_tauon_hadrdecay[j][NPROB_MAX-1];
       y_cumulative_tauon_mudecay[j][i]   /= y_cumulative_tauon_mudecay[j][NPROB_MAX-1];
       y_cumulative_tauon_edecay[j][i]    /= y_cumulative_tauon_edecay[j][NPROB_MAX-1];
    } //for

  }
  cout<<"Finished reading secondary interaction data.\n"; 
  
  /*string istring;
  char buffer[50];
  int n;  // counter
  

  ifstream ifile;
  int index;
  cout<<"Reading in data on secondary interactions.\n";
  for (int j=0;j<2;j++) {
    // for each energy
    for (int i=18;i<=21;i++) {
      if (!(i==21 && j==1)) {
      if (j==0)
	n=sprintf (buffer, "secondary/muons/dsdy_brems_1e%d.vec", i);
      if (j==1)
	n=sprintf (buffer, "secondary/muons/dsdy_brems_1e%d.5.vec", i);

      istring=buffer;

      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_muon_brems[index][NPROB] >> dsdy_muon_brems[index][NPROB];
	NPROB++;
      }
      ifile.close();   
      }
    }
  }

  for (int j=0;j<2;j++) {
    for (int i=18;i<=21;i++) {
      if (!(i==21 && j==1)) {
      if (j==0)
	n=sprintf (buffer, "secondary/muons/dsdy_epair_1e%d.vec", i);
      if (j==1)
	n=sprintf (buffer, "secondary/muons/dsdy_epair_1e%d.5.vec", i);

      istring=buffer;
      
      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_muon_epair[index][NPROB] >> dsdy_muon_epair[index][NPROB];
	NPROB++;
      }
      ifile.close();   
      }
    }
  }

  for (int j=0;j<2;j++) {
    for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {      
      if (j==0)
	n=sprintf (buffer, "secondary/muons/dsdy_pn_1e%d.vec", i);
      if (j==1)
	n=sprintf (buffer, "secondary/muons/dsdy_pn_1e%d.5.vec", i);

      istring=buffer;
      
      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_muon_pn[index][NPROB] >> dsdy_muon_pn[index][NPROB];
	NPROB++;
      }
      ifile.close();   
       }
    }
  }

   for (int j=0;j<2;j++) {
    for (int i=18;i<=21;i++) {
      if (!(i==21 && j==1)) {
      if (j==0)
	n=sprintf (buffer, "secondary/tauon/dsdy_brems_1e%d_tau.vec", i);
      if (j==1)
	n=sprintf (buffer, "secondary/tauon/dsdy_brems_1e%d.5_tau.vec", i);

      istring=buffer;
      
      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_tauon_brems[index][NPROB] >> dsdy_tauon_brems[index][NPROB];
	NPROB++;
      }
      ifile.close();   
      }
    }
   }

   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {       

       if (j==0)
	 n=sprintf (buffer, "secondary/tauon/dsdy_epair_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "secondary/tauon/dsdy_epair_1e%d.5_tau.vec", i);
       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_epair[index][NPROB] >> dsdy_tauon_epair[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {  
       if (j==0)
	 n=sprintf (buffer, "secondary/tauon/dsdy_pn_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "secondary/tauon/dsdy_pn_1e%d.5_tau.vec", i);
       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_pn[index][NPROB] >> dsdy_tauon_pn[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {

       if (j==0)
	 n=sprintf (buffer, "secondary/tauon/dsdy_hadrdecay_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "secondary/tauon/dsdy_hadrdecay_1e%d.5_tau.vec", i);       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_hadrdecay[index][NPROB] >> dsdy_tauon_hadrdecay[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       
       if (!(i==21 && j==1)) {
       if (j==0)
	 n=sprintf (buffer, "secondary/tauon/dsdy_edecay_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "secondary/tauon/dsdy_edecay_1e%d.5_tau.vec", i);       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_edecay[index][NPROB] >> dsdy_tauon_edecay[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {
       if (j==0)
	 n=sprintf (buffer, "secondary/tauon/dsdy_mudecay_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "secondary/tauon/dsdy_mudecay_1e%d.5_tau.vec", i);       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_mudecay[index][NPROB] >> dsdy_tauon_mudecay[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
 
   // for filling vectors with y values distributed so that they follow
   // dsdy distributions.

   for (int j=0;j<7;j++) {

     int_muon_brems[j]=0;
     int_muon_epair[j]=0;
     int_muon_pn[j]=0;   
     int_tauon_brems[j]=0;
     int_tauon_epair[j]=0;
     int_tauon_pn[j]=0;
     int_tauon_hadrdecay[j]=0;
     int_tauon_edecay[j]=0;
     int_tauon_mudecay[j]=0;
     
     // integrating prob. distributions.
     for (int i=0;i<NPROB_MAX;i++) {
       int_muon_brems[j]+=dsdy_muon_brems[j][i];
       //cout << "int_muon_brems is " << int_muon_brems[j] << "\n";
       int_muon_epair[j]+=dsdy_muon_epair[j][i];
       int_muon_pn[j]+=dsdy_muon_pn[j][i];
       int_tauon_brems[j]+=dsdy_tauon_brems[j][i];
       int_tauon_epair[j]+=dsdy_tauon_epair[j][i];
       int_tauon_pn[j]+=dsdy_tauon_pn[j][i];
       int_tauon_hadrdecay[j]+=dsdy_tauon_hadrdecay[j][i];
       int_tauon_edecay[j]+=dsdy_tauon_edecay[j][i];
       int_tauon_mudecay[j]+=dsdy_tauon_mudecay[j][i];
     }

     // maximum value of prob. dist. 
     max_muon_brems=dMax(dsdy_muon_brems[j],NPROB_MAX);
     max_muon_epair=dMax(dsdy_muon_epair[j],NPROB_MAX);
     max_muon_pn=dMax(dsdy_muon_pn[j],NPROB_MAX);   
     max_tauon_brems=dMax(dsdy_tauon_brems[j],NPROB_MAX);
     max_tauon_epair=dMax(dsdy_tauon_epair[j],NPROB_MAX);
     max_tauon_pn=dMax(dsdy_tauon_pn[j],NPROB_MAX);
     max_tauon_hadrdecay=dMax(dsdy_tauon_hadrdecay[j],NPROB_MAX);
     max_tauon_edecay=dMax(dsdy_tauon_edecay[j],NPROB_MAX);
     max_tauon_mudecay=dMax(dsdy_tauon_mudecay[j],NPROB_MAX);
     
     // minimum value of prob. dist.
     min_muon_brems=Tools::dMinNotZero(dsdy_muon_brems[j],NPROB_MAX);
     min_muon_epair=Tools::dMinNotZero(dsdy_muon_epair[j],NPROB_MAX);
     min_muon_pn=Tools::dMinNotZero(dsdy_muon_pn[j],NPROB_MAX);   
     min_tauon_brems=Tools::dMinNotZero(dsdy_tauon_brems[j],NPROB_MAX);
     min_tauon_epair=Tools::dMinNotZero(dsdy_tauon_epair[j],NPROB_MAX);
     min_tauon_pn=Tools::dMinNotZero(dsdy_tauon_pn[j],NPROB_MAX);
     min_tauon_hadrdecay=Tools::dMinNotZero(dsdy_tauon_hadrdecay[j],NPROB_MAX);
     min_tauon_edecay=Tools::dMinNotZero(dsdy_tauon_edecay[j],NPROB_MAX);
     min_tauon_mudecay=Tools::dMinNotZero(dsdy_tauon_mudecay[j],NPROB_MAX);
     
     if (min_muon_brems<=0)
       cout << "Minimum probability is <=0!\n";

     
     // for each y bin in dsdy curve, fill vector y_muon_brem with as
     // many of y's as you need to get the right distribution.
     for (int i=0;i<NPROB_MAX;i++) {
       y_cumulative_muon_brems[j][i]      = dSum(dsdy_muon_brems[j],i+1);
       y_cumulative_muon_epair[j][i]      = dSum(dsdy_muon_epair[j],i+1);
       y_cumulative_muon_pn[j][i]         = dSum(dsdy_muon_pn[j],i+1);
       y_cumulative_tauon_brems[j][i]     = dSum(dsdy_tauon_brems[j],i+1);
       y_cumulative_tauon_epair[j][i]     = dSum(dsdy_tauon_epair[j],i+1);
       y_cumulative_tauon_pn[j][i]        = dSum(dsdy_tauon_pn[j],i+1);
       y_cumulative_tauon_hadrdecay[j][i] = dSum(dsdy_tauon_hadrdecay[j],i+1);
       y_cumulative_tauon_mudecay[j][i]   = dSum(dsdy_tauon_mudecay[j],i+1);
       y_cumulative_tauon_edecay[j][i]    = dSum(dsdy_tauon_edecay[j],i+1);
     } //for

     // normalize the distributions
     for (int i=0;i<NPROB_MAX;i++) {
       y_cumulative_muon_brems[j][i]      /= y_cumulative_muon_brems[j][NPROB_MAX-1];
       y_cumulative_muon_epair[j][i]      /= y_cumulative_muon_epair[j][NPROB_MAX-1];
       y_cumulative_muon_pn[j][i]         /= y_cumulative_muon_pn[j][NPROB_MAX-1];
       y_cumulative_tauon_brems[j][i]     /= y_cumulative_tauon_brems[j][NPROB_MAX-1];
       y_cumulative_tauon_epair[j][i]     /= y_cumulative_tauon_epair[j][NPROB_MAX-1];
       y_cumulative_tauon_pn[j][i]        /= y_cumulative_tauon_pn[j][NPROB_MAX-1];
       y_cumulative_tauon_hadrdecay[j][i] /= y_cumulative_tauon_hadrdecay[j][NPROB_MAX-1];
       y_cumulative_tauon_mudecay[j][i]   /= y_cumulative_tauon_mudecay[j][NPROB_MAX-1];
       y_cumulative_tauon_edecay[j][i]    /= y_cumulative_tauon_edecay[j][NPROB_MAX-1];
     } //for
     
   }

   cout<<"Finished reading secondary interaction data.\n";
  */
} //end method ReadSecondaries


 void Secondaries::GetSecondaries(Settings *settings1,string nuflavor,double plepton,double &em_secondaries_max,double &had_secondaries_max,int &n_interactions,TH1F *hy) {


  em_secondaries_max=0.;
  had_secondaries_max=0.;

  int i=(int)((log10(plepton)-18.)*2.);
  if (i>6)
    i=6;
  if (i<0)
    i=0;

  int n_brems,n_epair,n_pn; // number of interactions of each type.
  int index_y; // index along the horizontal axis of ped's plots
  double rnd1=1000.;
  double rnd2=1000.;  // random numbers for throwing at dart board
  double y = 0; // inelasticity
 
  string whichtype; // which type of interaction corresponds to that index
  


  if (nuflavor=="numu") {   
    n_brems=gRandom->Poisson(int_muon_brems[i]); // pick number of brem interactions
    n_epair=gRandom->Poisson(int_muon_epair[i]); // # of pair production
    n_pn=gRandom->Poisson(int_muon_pn[i]); // # photonuclear interactions   
    
    n_interactions+=(n_brems+n_epair+n_pn);	


    for (int j=0;j<n_brems+n_epair+n_pn;j++) {
      rnd1=gRandom->Rndm();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn))
	whichtype="brems";
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn))
	whichtype="epair";
      else
	whichtype="pn";



      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;

      if (whichtype=="brems") {	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_brems[i],NPROB,rnd1,y);
      }
      else if (whichtype=="epair") {	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_epair[i],NPROB,rnd1,y);	
      }
      else if (whichtype=="pn") {
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_pn[i],NPROB,rnd1,y);
      }
     
      if (y*plepton>max(em_secondaries_max,had_secondaries_max)) {  // if this is the largest interaction for this event so far
	if (whichtype=="brems" || whichtype=="epair") {  // save it
	  em_secondaries_max=y*plepton;

	}
	if (whichtype=="pn") 
	  had_secondaries_max=y*plepton;
	 
		
      }
    } // loop over secondary interactions
  } // end if it was a muon neutrino
  if (nuflavor=="nutau") {
    n_brems=gRandom->Poisson(int_tauon_brems[i]);
    n_epair=gRandom->Poisson(int_tauon_epair[i]);
    n_pn=gRandom->Poisson(int_tauon_pn[i]);

    n_interactions+=(n_brems+n_epair+n_pn); // increment number of secondary interactions.

    for (int j=0;j<n_brems+n_epair+n_pn;j++) { // loop over secondary interactions. 
      
      rnd1=gRandom->Rndm();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn))
	whichtype="brems";
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn))
	whichtype="epair";
      else
	whichtype="pn";
  
      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;

      if (whichtype=="brems") {  // bremstrahlung interaction
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_brems[i],NPROB,rnd1,y);
      }
      if (whichtype=="epair") { // pair production
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_epair[i],NPROB,rnd1,y);
      }
      if (whichtype=="pn") {
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_pn[i],NPROB,rnd1,y);
      }

      if (settings1->HIST==1 && !settings1->ONLYFINAL && hy->GetEntries()<settings1->HIST_MAX_ENTRIES)
	hy->Fill(y);
      if (y*plepton>max(em_secondaries_max,had_secondaries_max)) { // if this is the biggest secondary signal yet,
	if (whichtype=="brems" || whichtype=="epair") // save it.
	  em_secondaries_max=y*plepton;
	if (whichtype=="pn")
	  had_secondaries_max=y*plepton;
      }
    }
   

    if (TAUDECAY) {
      n_interactions++; // increment number of interactions, for plotting.

      rnd1=gRandom->Rndm();
      if (rnd1<0.65011)  // pick which type of decay it is.
	whichtype="hadrdecay";
      if (rnd1>=0.65011 && rnd1<0.8219)
	whichtype="mudecay";
      if (rnd1>=0.8219)
	whichtype="edecay";
           
      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;     
      
      if (whichtype=="hadrdecay") { // hadronic decay
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_hadrdecay[i],NPROB,rnd1,y);	
      }
      else if (whichtype=="edecay") { // e decay	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_edecay[i],NPROB,rnd1,y);
      }
      else if (whichtype=="mudecay") { // mu decay
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_mudecay[i],NPROB,rnd1,y);
      }
      
     
      if (y*plepton>max(em_secondaries_max, had_secondaries_max)) {  // if this is the biggest interaction yet,    
	if (whichtype=="edecay") // save it.
	  em_secondaries_max=y*plepton;
	if (whichtype=="hadrdecay")
	  had_secondaries_max=y*plepton;
      } //if     
    } //if (TAUDECAY)
  } //if (nutau)

} //GetSecondaries



 int Secondaries::GetEMFrac(Settings *settings1,string nuflavor,
		     string current,
		     string taudecay,	      
		     double y,
		     TH1F *hy,
				  double pnu,				  
				  int inu,


		     double& emfrac,
		     double& hadfrac,
			    int& n_interactions, int taumodes1,double ptauf) {




  if (current=="cc")
    plepton=(1.-y)*pnu;
  else
    plepton=0.;
  
  if (nuflavor=="nue" && current=="cc") {
    emfrac=1.-y;
    hadfrac=y;
  }
  else if(nuflavor=="numu" && current=="cc") {
    emfrac=1.E-10;
    hadfrac=y;
  }
  else if(nuflavor=="nutau" && current=="cc") {
    // behaves like a muon
    if(taumodes1 ==1){//taumodes==1; tau created somewhere in rock and decays at posnu.
      
      this->GetEMFracDB(emfrac,hadfrac);
      
    }
    else if (taumodes1 == 0){
      emfrac=1.E-10;
      hadfrac=y;
    }
    


  }
  else if (current=="nc") {
    emfrac=1.E-10;
    hadfrac=y;
  }


  em_secondaries_max =emfrac; // initialize search for maximum signal among primary, secondary interactions.
  had_secondaries_max=hadfrac;

  
  
  if (SECONDARIES==1 && current=="cc" && settings1->FORSECKEL!=1) {

    while (1) {

      GetSecondaries(settings1,nuflavor,plepton,em_secondaries_max,had_secondaries_max,n_interactions,hy); // find how much em and hadronic energies comes from secondary interactions.  keep picking until you get a bunch of secondary interactions that conserve energy

      if (em_secondaries_max+had_secondaries_max<=plepton*(1.+1.E-5)) // if conserves energy, break.
	break;
      else {
	secondary_e_noncons++; //Record how many times we come up with something that doesn't conserve energy
	em_secondaries_max=emfrac;
	had_secondaries_max=hadfrac;
      } //else
    } //while(1)

    if ((em_secondaries_max+had_secondaries_max)>(emfrac+hadfrac)*pnu) { // if maximum signal from secondaries is larger than
                                                                         // signal from primary interaction
      emfrac=em_secondaries_max/pnu; // then use that one.
      hadfrac=had_secondaries_max/pnu;
      if (emfrac <= 1.E-10)
	emfrac=1.E-10;
      if (hadfrac<= 1.E-10)
	hadfrac=1.E-10;
    } //if
  } //if (charged current, secondaries on)

  if (nuflavor=="numu" && current=="cc" && n_interactions==0)
    cout << "Look at this one.  inu is " << inu << "\n";
  


  if ((y<0 || y>1) && y != -999.) 
    cout <<  "illegal y=" << y << "\n";
          
  if (emfrac+hadfrac>1.00001) {
    cout << "error emfrac,hadfrac=" << emfrac << " " << hadfrac << " " << emfrac+hadfrac << "\n";
    cout << "nuflavor,taudecay=" << nuflavor << " " << taudecay << "\n";
  } //if
  
  return 1;

} //GetEMFrac


//----------------------------------------------------------
//InitTauola()
//Initializes the tau decay information

 void Secondaries::InitTauola() {
   for(int k=0;k<5;k++)
    tauolainfile >> tauola[0][k];
   for(int i=1;i<N_TAUOLA;i++)
    for(int j=0;j<6;j++)
      tauolainfile >> tauola[i][j];
      
  //cout<<"first line of tauola is "<<tauola[0][0]<<","<<tauola[0][1]<<","<<tauola[0][2]<<","<<tauola[0][3]<<"\n";
 
  return;
}//InitTauola

void Secondaries::GetTauDecay(string nuflavor,string current,string& taudecay,double& emfrac_db,double& hadfrac_db) {

  if (!(nuflavor=="nutau" || current=="cc" || interestedintaus))
    return;

  // if nu_tau choose tau decay type
  
  double rnd = gRandom->Rndm();
  int decay = static_cast<int>(rnd*(N_TAUOLA-2)+1);
  
  hadfrac_db=tauola[decay][3];
  emfrac_db=tauola[decay][4];
  
  if(tauola[decay][1]!=0)
    taudecay="m";
  else
    taudecay="e";
  

  if(taudecay=="m")
    secondbang=false; //for all muon decays, the interaction point chosen is the neutrino interaction since we don't detect the decay if
  //the tau decays into a muon.
  else {
    double rnd=gRandom->Rndm();
    if(rnd>TAUFRAC) {
      secondbang=true;
      count_nfb++;
    } else
      secondbang=false;
  }
  

} //GetTauDecay

//-----------------------------------------------------
//GetEMFracDB()
//Gets the emfrac_db and hadfrac_db for a doublebang

 void Secondaries::GetEMFracDB(double& emfrac_db, double& hadfrac_db) {

   
  double rnd = gRandom->Rndm();
  int decay = static_cast<int>(rnd*(N_TAUOLA-2)+1);
  hadfrac_db=tauola[decay][3];
  emfrac_db=tauola[decay][4];

  return;
}//GetEMFracDB

//------------------------------------------------------
//GetDBViewAngle()
//Gets the viewangle of the second bang

 double Secondaries::GetDBViewAngle(const Vector &refr, const Vector &nnu) {

  return ((nnu.ChangeCoord(refr)).Angle(z_axis));

}//GetDBViewAngle

//------------------------------------------------------
//GetFirstBang()
//Gets the position of the first bang when the interaction point is the tau decay point

//  void Secondaries::GetFirstBang(const Position &r_in, const Vector &nnu, Position &posnu, double len_int_kgm2, double chord, double &nuentrancelength) {
  
//   double weightbang;
//   double junk1;
//   double junk2;
//   double junk3;
//   int junk4,junk5,junk6;
//   double myair=0;

//   Vector r_out = r_in + chord*nnu;

//   antarctica->Getchord(len_int_kgm2,r_in,r_out,
// 		  junk1,weightbang,junk2,myair,junk3,junk4,junk5,junk6);
//   double r1,r2;
//   if(weightbang>.999)
//     r2=gRandom->Rndm()*chord;
//   else {
//     do {
//       r1=gRandom->Rndm();
//       r2=gRandom->Rndm()*chord;
//       r_out = r_in + r2*nnu;
//       antarctica->Getchord(len_int_kgm2,r_in,r_out,
// 		      junk1,weightbang,junk2,myair,junk3,junk4,junk5,junk6);
//     }
//     while(r1>1-weightbang);
//   }
//   posnu = r_in + r2*nnu;
//   nuentrancelength=r2;

//   return;
// }//GetFirstBang

//---------------------------------------------------------
//NFBWeight()
//Gets the weight of the tau decay for second bang events
  double Secondaries::NFBWeight(double ptau, double taulength) {
  
  double gamma=ptau/MTAU;
  double D=TAUDECAY_TIME*CLIGHT*gamma;

  return exp(-taulength/D);

}
void Secondaries::Picky(double *y_cumulative,int NPROB,double rnd,double& y) {
  for (int i=0;i<NPROB;i++) {
    if (y_cumulative[i]<=rnd && y_cumulative[i+1]>rnd) {
      y=(double)i/(double)NPROB;
      continue; // once you found the right bin, stop looping.
    } //if
  } //for
} //Picky

//-------------------------------------------------------------
double Secondaries::GetTauWeight(Primaries *primary1, Settings *settings1,IceModel*antarctica1, double pnu, int nu_nubar, 
				 int currentint, double& ptauf, const Position posnu, const Position earth_in,
				 int& crust_entered, // 1 or 0
				 int& mantle_entered, // 1 or 0
				 int& core_entered,double  *prob_at_z_vector, double *distance_vector,
				 double *prob_at_z_vector1, double *distance_vector1, double *Energy_vector, double *Energy_vector1,				     double *nusurvarray, double *dEtauidEtaufarray,
				 double& tauchord, double *avgdensityarray, 
				 double *densityarray,int inu,double& TauWeight, double& weight_prob,int NNU,
				 const Position &r_enterice, const Position &nuexitice,int *arraytrigger,
				 int *looparray,double *etaufarray,double *TauWeightarray, double *PDFarray){
	
  // Settings *settings1=new Settings();
  //Primaries *primary1=new Primaries();
  Vector chord3;
  Vector nchord;
  
  int  N=1E3;
 

 //Find the chord, its length and its unit vector.
  chord3 = nuexitice - earth_in;//posnu-earth_in; 
  double Distance=chord3.Mag();
  nchord = chord3 / Distance;
  tauchord=Distance;
    
  
  double Etaui,dEtau,Etau_final;
   double y=0;
   double yweight;//d_dzPnu_surv;
  
   double  zdistance; //m
   double zdistance1=0;//m
   zdistance=0.;	
   TauWeight=0;//total tau survival probability.
   weight_prob=0;
   double TauWeight_tmp=0;
   double weight_prob_tmp=0;
   double prob_at_z=0; //this will be used to get weight1
  
   
   double Prob_Nu_Surv;
   double last_step_prob=0;
  
   double tau_surv;
   double tau_surv1;
   double sigma = 0;
   double len_int_kgm2 =0;
   
   primary1->GetSigma(pnu,sigma,len_int_kgm2,settings1,nu_nubar,currentint);
	
	
  
   double step=Tools::dMin(len_int_kgm2/densities[1]/10,25); //how big is the step size
   
   
   double avgdensity =0;//initilize average density.
   double density_total=0;//initilize running sum
   double density_now=0;//density at this step
   double density_tau_total=0;
   double avgtaudensity=0;
   Position posnunow;
   Position postaunow;
   Position postaustart;
   double lat;
   double lon;
   double tau_distance;
   double altitude_tau;
   
   double Etau_now=Etau_final;
   double Emin=1E15;
   double lat_tau=0;
   int i=0;
   int j=0;
   int i1=0;
   double dPdz_now=0;
   double dPdz_last=1;
   double taudensity=0;
   double dEtauidEtauf=0;
   double yintegralsum=0;
   double integral_sum=0;
   double integral_total=0;
   vector<double> mydensityvector;
  
   vector<double> myavgdensityvector;
   
   
   int etauint =0;
   double startingz=0;

   double totaltaudistance=0;
   int totalnusteps=0;
   Vector nchord1;
   double enter_ice_mag = r_enterice.Distance(earth_in);
   //first we will fill vectors for the density and avgdensity for neutrino
   for(double taudistance=0;taudistance<=Distance;taudistance+=step){
     nchord1=taudistance*nchord;
     postaunow=earth_in+nchord1;
     lat_tau=postaunow.Lat();
     altitude_tau = postaunow.Mag()-antarctica1->Geoid(lat_tau);
     density_now=antarctica1->GetDensity(altitude_tau,postaunow,posnu,crust_entered,mantle_entered,core_entered);
     mydensityvector.push_back(density_now);//filled with density at that step
     
     density_total+=density_now;
     avgdensity=density_total/(totalnusteps+1);//avgdensity
     myavgdensityvector.push_back(avgdensity);//avgdensity up to that step for neutrino
     totalnusteps++;//total number of steps neutrino will take through earth
     
   }
   //Begin looping over final energy states for the tau
   //cout<<"pnu is "<<pnu<<"\n";
   for(int clear=0; clear<2000;clear++){
     distance_vector[clear]=0;
     distance_vector1[clear]=0;
     prob_at_z_vector[clear]=0;
     prob_at_z_vector1[clear]=0;
     }
   for(double logEtau_final=log10(Emin);logEtau_final<=log10(pnu);logEtau_final+=.01){//integral over energy (in log space?)
     Etau_final= pow(10,logEtau_final);
     //cout<<"logEtau_final is "<<logEtau_final<<"\n";
     vector<double> myenergyvector;//vector to hold initial energies
     vector<double> myPsurvvector;//vector to hold the chance tau would survive from that step
     
     int totalsteps=0;
     i=0;
     TauWeight_tmp=0;
     myenergyvector.push_back(Etau_final);
     Etau_now=Etau_final;
     double gamma = Etau_final/mT;
     postaustart=nuexitice;
     //calculate the initial energy needed at the step so the tau will end at the correct final energy
     for(int energysteps=totalnusteps;energysteps>0;energysteps--){//start at nuexit, work backwards
       density_now=mydensityvector[energysteps];
       
       Etau_now=Etau_now + (B0+B1*log(Etau_now/E0))*density_now*Etau_now*step;//E_i=E_last +dE/dz*dz
       if(Etau_now <=pnu){
	 totaltaudistance=(totalnusteps-energysteps)*step;
	 totalsteps++;//number of steps tau can take
	 myenergyvector.push_back(Etau_now);
       }//Etau_now<=pnu
       else if(Etau_now >pnu){//Initial energy cannot go above pnu
	 break;
       }
       p++;
     } //energy steps
     //calculate the chance the tau would survive from one step to the end.
     for(int k1=totalsteps;k1>=0;k1--){//tau surv vector
       tau_surv=1;
       for(int k2=k1-1;k2>=0;k2--){
	 Etau_now=myenergyvector[k2];//energy vector starts at the endpoint and goes to earth_in;
	 if (Etau_now >0)
	   tau_surv=tau_surv*(1-step*mT/(cT*Etau_now));
	 else
	   tau_surv = 0;
       }
       myPsurvvector.push_back(tau_surv);
     }//tau surv
    
     //startingz=r_enterice.Distance(earth_in);
     
      startingz = Distance-totaltaudistance;

   //////////////////Integral over distance/////////////////////////////////////////
     for(zdistance = startingz; zdistance<=Distance; zdistance +=step){
        int nustep = (int)zdistance/step;
	int tauenergystep=totalsteps-i;//step number to get the correct initial energy at that point;
	int tausurvstep =i;
	
	nchord1 =  zdistance*nchord;
	posnunow = earth_in + nchord1; //vector pointing to the step we are currently on.
    


	avgdensity = myavgdensityvector[nustep];
	taudensity=mydensityvector[nustep];
	Etaui=myenergyvector[tauenergystep];//Energy vector is filled backwards (final to initial), pull out correct energy
	y=1.-Etaui/pnu;
	

	tau_surv = myPsurvvector[i];
	
	if(zdistance >=enter_ice_mag){
	  tau_surv=1;
	}

	if (Etaui>=pnu || Etaui<Etau_final || Etaui!=Etaui){//to catch anything that might make it through. Precaution
	prob_at_z =0;
     }
      else{

	yweight=primary1->Getyweight(pnu,y,nu_nubar,currentint);

	Prob_Nu_Surv = avgdensity/len_int_kgm2*exp(-zdistance*avgdensity*1./len_int_kgm2);
	
	dEtauidEtauf = exp(B1*taudensity*(Distance-zdistance))*Etaui/Etau_final;
	
	prob_at_z = Prob_Nu_Surv*yweight*1./pnu*dEtauidEtauf*step*tau_surv;
	
	if(zdistance<=enter_ice_mag){
	   prob_at_z*=1.-exp(-1.*(posnu.Distance(nuexitice)/(gamma*cT)));//chance to decay in ice
	}
	else if(zdistance>enter_ice_mag){
	  prob_at_z*=1.-exp(-1.*(posnunow.Distance(nuexitice)/(gamma*cT)));//chance to decay in ice if already in ice
	}
	TauWeight_tmp+=prob_at_z;
	
	if(logEtau_final>15.999 && logEtau_final < 16.001){
	  if( i<2000){
	    distance_vector[i]=Distance-zdistance;
	    // cout<<"Distance_vector["<<i<<"] is "<<distance_vector[i]<<"\n";
	    //cout<<"Energy_vector["<<i<<"] is "<<Energy_vector[i]<<"\n";
	    //cout<<"Etaui is "<<Etaui<<"\n";
	    prob_at_z_vector[i]=prob_at_z;
	    Energy_vector[i]=Etaui;
	  }//i<2000
	  if(i>1999 && i<4000){
	    Energy_vector1[i1]=Etaui;
	    distance_vector1[i1]=Distance-zdistance;
	    prob_at_z_vector1[i1]=prob_at_z;
	    i1++;
	  }
	}//logEtaufinal>
	i++;
      }//etaui<pnu
     }//zdist loop
     
     //TauWeight_tmp*=5.5E15;//dP/dE_final * dE_final
     TauWeight_tmp*=log(10)*Etau_final;//dP/dE ->dP/d(log10(E))
     TauWeight_tmp*=.01; //dP/dlog10(E) * dlog10(E)
   
     TauWeight+=TauWeight_tmp;
     
     //prob_at_z_vector.push_back(row);
     etaufarray[etauint]=Etau_final;
     PDFarray[etauint] = TauWeight;
     etauint++;
     
   }//Etau_final loop
  
   // }//xloop
   
    double xrandom= TauWeight*gRandom->Rndm();
    //cout<<"TauWeight is "<<TauWeight<<" xrandom is "<<xrandom<<"\n";
   for(int loopthrough =0; loopthrough<=etauint;loopthrough++){
     if(xrandom >PDFarray[loopthrough] && xrandom <PDFarray[loopthrough+1]){
	ptauf=etaufarray[loopthrough];
	break;
     }//if xrandom
     
   }//loopthrough
   //cout<<"ptauf is "<<ptauf<<"\n";
   return 1;
} //GetTauWeight

double Secondaries::probabilityTauSurv(double ptaui, double ptau_final,double density){
	//from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
	//Equation 16  & Equation 30.
 
/*	double B0,B1,E0;//parameterization using a logarithmic dependence on energy 
	//for B, the tau elecromagnetic energy loss parameter. 
	double mT;//mass of the Tau in Gev
	double cT;//Tau Decay length in cm
*/	double p = density;
        double probSurv1;
	double mT = 1.77684E9; //mass of Tau in eV
	double cT = 0.00008693; // m
	double A = mT*B1/(cT*p*pow(B0,2));//eV
	double B = mT/(cT*B0*p);//eV
	double C = (1/ptau_final)*(1+log(ptau_final/E0));//1/eV
	double D = (1/ptaui)*(1+log(ptaui/E0));//1/eV
	double F = (1/ptau_final)-(1/ptaui);//1/eV
	
	//	cout<< "A,B,C are "<<A<<","<<B<<","<<C<<"\n";
	//	cout<< "D,F are "<<D<<","<<F<<"\n";
	probSurv1 = exp(A*(C-D)-(B*F));

	return probSurv1;//Correct
}

double Secondaries::dTauSurvdx(double& probSurv1,double ptau_final,double step){//double ptaui, double ptau_final,double density,double distance){
	//from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
	//Equation 16  & Equation 30.
        double mT = 1.77684E9; //mass of Tau in eV
	double cT = 0.00008693; // m
        //double dTauSurvdx=-probSurv1*mT/(cT*ptau_final);
        probSurv1=probSurv1-(probSurv1*mT*step/(cT*ptau_final));
	return 0;//Correct
}

double Secondaries::TauEnergyInitial(double ptau_final, double Distance, double z_distance, double density){
	//equation found using Equation 13 Case (III),
	//from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
        double p = density;
        double ptaui;
	double zprime=(Distance-z_distance);

	ptaui = E0*exp((log(ptau_final/E0)+B0/B1*(1-exp(-B1*p*zprime)))*exp(B1*p*zprime));
	
	return ptaui;
}
double Secondaries::TauEnergyFinal(double ptaui,double step,double density){
        double ptau_final;
        double zprime=step;
	ptau_final=E0*exp(-(B0/B1)*(1-exp(-B1*density*zprime))+log(ptaui/E0)*exp(-B1*density*zprime));
	return ptau_final;
}
double Secondaries::FindBin(double ptau_final){
       
         double logptauf=0;
         double binnumber=0;
        
	 logptauf= log10(ptau_final);
	 binnumber= logptauf-14;//going to start binning at 10^14 eV
	 binnumber=binnumber/.1;//bin sizes are .1 (in log space);
	 binnumber=floor(binnumber);
	 return binnumber;
  
}
