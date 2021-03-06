#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include "TTreeIndex.h"
#include "TChain.h"
#include "TH1.h"
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
#include "TLatex.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
#include "Math/InterpolationTypes.h"
#include "Math/Interpolator.h"
#include "signal.h"
#include "TGraph.h"
#include "TGraphPolar.h"
#include "TGraph2D.h"
#include "TMarker.h"
#include "TArc.h"
#include "TPaletteAxis.h"
#include "TImage.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TPaveStats.h"
#include "TVirtualFitter.h"

#include "poly34.h"
#include "Constants.h"
#include "Settings.h"
#include "position.hh"
#include "earthmodel.hh"
#include "Tools.h"
#include "vector.hh"
#include "roughness.hh"
#include "anita.hh"
#include "balloon.hh"
#include "icemodel.hh"
#include "trigger.hh"
#include "Spectra.h"
#include "signal.hh"
#include "secondaries.hh"
#include "ray.hh"
#include "counting.hh"
#include "Primaries.h"

// HEALPix C++ module includes.
//#include "healpix_map.h"
#include "healpix_base.h"
#include "chealpix.h"
//#include "datatypes.h"
//#include "healpix_data_io.h"
//#include "alm.h"

#include "analysis_info_4pol.h"

analysis_info_4pol *pol4_Ptr = 0;

using namespace std;
TStyle* RootStyle();
TStyle *color=RootStyle();
class EarthModel;
class Position;


int whole_flag=0;//0 for fractional, 1 for while events
int k_value = 4;
long n_side =pow(2,k_value);
//Healpix_Base heal_pix;//(n_side,RING);
//pointing point_here;


int n_pix_int=12*n_side*n_side; //Total number of pixels
vector<double> BASE(n_pix_int,0.);
double BASE2;
vector<double> weight_frac(n_pix_int,0.);
vector<int> num_frac(n_pix_int,0);
vector<int> numevents(n_pix_int,0.);
//double PI = 3.14159264;

int NNU;        // number of neutrinos
int WHICHPATH=7;
 
double RANDOMISEPOL=0.;

double fit_10(double *x, double *par){
  double fitval = par[0]*pow(10,par[1]*1+par[2]);
  return fitval;

}
void FractionalArea(vector <vector <double> > phis_theta,vector<int> error_bin, vector<int> &area_pix, vector<double> &areas);
void FractionalArea_cone(vector<int> error_bin,vector<int> &area_pix, vector<double> &areas, double eta, double H,double surface,double heading, double lat,double lon);
int ErrorPoints(double *x,double *y,double *z,double eta, double theta, double phi,double H,double surface);
void GetLatLonfromCartesian(double *Lats, double *Lons, double heading, double lat, double lon, double *x, double *y, double *z);
void GetCartesianfromLatLon(double *Lats, double *Lons, double heading, double lat, double lon, double *x, double *y, double *z,double R);
void GetCartesianfromLatLon_vec(vector<double> Lats, vector<double> Lons, double heading, double lat, double lon, vector<double> &x,vector< double> &y, vector<double> &z,double R);
int HealPixBinFromCone(double x, double y,double z, double eta, double H,double surface,double heading, double lat,double lon);
void ToConeCoordinates(double *xfinal, double *yfinal, double *zfinal,double eta1, double *x,double *y,double *z,double H,double surface);
 void ToConeCoordinates_vec(vector<double> &xfinal, vector<double> &yfinal, vector<double> &zfinal,double eta1, vector<double> x,vector<double> y,vector<double> z,double H,double surface);
void LatLon2phitheta(double lat,double lon, double &phi, double &theta);
void PhiTheta2LatLon(double phi, double theta, double &lat, double &lon);
void SphericaltoCart(double phi, double theta, double &x, double &y);
void initializeBaselist(double latBases[1000], double lonBases[1000]);
void Cart2Sphere(double& phi, double& theta, double x, double y);
void GetTimePerBin(IceModel *antarctica);
double GetDistance(double Lat1, double Lon1, double Lat2, double Lon2);
int GetBrianBin(double lat, double lon);
void BintoLatLon(int binnumber, double& lat, double& lon);

double test_anita_x=0.;
double test_anita_y=0.;
double test_anita_z=0.;

double test_anita_lat=0.;
double test_anita_lon=0.;

vector< vector <double> > boundary_x;
vector< vector <double> > boundary_y;
vector< vector <double> > boundary_z;
vector< vector <double> > boundary_x_prime;
vector< vector <double> > boundary_y_prime;
vector< vector <double> > boundary_z_prime;

vector< vector< double> > ellipse_coordinates;


int main(int argc, char **argv) {
  // RootStyle()->cd();
  gStyle=color;
 
  double icethck;
 
  vector<double> lon_vector;
  vector<double> lat_vector;
  vector<double> icethick_vector;
  vector<double> phi_vector;
  vector<double> theta_vector;
  vector<double> bedmap_x;
  vector<double> bedmap_y;
  double highest_point=-1;
  int highest_point_bin=0;
  double theta_1;
  double lon_ice;
  double lat_ice;
  double min_lon=0;
  double max_lon=0;
  double min_lat=0;
  double max_lat=0;
  double min_ice=10000;
  double max_ice=0;
  double min_theta=7;
  double max_theta=0;
  double min_phi=7;
  double max_phi=0;
  double x_min=600;
  double x_max=-600;
  double y_min=600;
  double y_max=-600;
  int icemodel = 1;
  int earthmodel=0;
  int weightabsorption =1;

  // double theta_1;
  double phi_1;
  double PI = 3.14159264;
 
  

  ///////////healpix stuff
 double xholder_temp[100]={0.};
  double yholder_temp[100]={0.};
  double zholder_temp[100]={0.};
  double xholder_temp1[100]={0.};
  double yholder_temp1[100]={0.};
  double zholder_temp1[100]={0.};
  double xholder[100]={0.};
  double yholder[100]={0.};
  double zholder[100]={0.};
  double Lats[100]={0.};
  double Lons[100]={0.};
  double Lats_point[100]={0.};
  double Lons_point[100]={0.};
  double eta = 70.;
  double theta_error =1*PI/180;
  double phi_error =2*PI/180;
  double H = 32.;
  double lat_trial = -80;
  double lon_trial = 30;
  double heading = 180;
  
  double min_peakVal=1000;
  double min_Hilbert=1000;
  //loop through events
  //take sourceLat and Lon of events, get theta and phi put into pointing?   
  //then use ang2pix  to get pix, add event to Base[pixel]
     double Lat_point;
     double Lon_point;
     double anitaLat1;
     double anitaLon1;
     double phi_event;
     double theta_event;

     double anitaLat_event;
     double anitaLon_event;
     double Lat_event;
     double Lon_event;
     double x_map;
      double y_map;
      int EventEntry;
      int errorellipse_ctr=0;

      int errorpoints_return;
      double min_theta_error=100.;
      double max_theta_error=-100.;
      double min_phi_error=100.;
      double max_phi_error=-100.;
      double min_x_error=10000.;
      double max_x_error=-10000.;
      double min_y_error=10000.;
      double max_y_error=-10000.;
       double min_x_error_plot=10000.;
      double max_x_error_plot=-10000.;
      double min_y_error_plot=10000.;
      double max_y_error_plot=-10000.;
      vector<double> error_phi (100,0.);
      vector<double> error_theta (100,0.);
      vector<int> error_bin (100,0);
      //vector< vector <double> > peakVal_vector (n_pix_int+1, vector<double>(1));
      //vector< vector <double> > peakHilbert_vector (n_pix_int+1, vector<double>(1));
      //vector< vector <double> > SNR_vector (n_pix_int+1, vector<double>(1));
      
      //vector< vector <double> > healpix_bin_weights (n_pix_int+1, vector<double>(1));
      //vector< vector <double> > bedmap_flagged;

      //vector< vector <double> > ratio_vector (n_pix_int+1, vector<double>(1));
      //vector< vector <double> > polFrac_vector (n_pix_int+1, vector<double>(1));
      //vector< vector <int> > eventNumber_vector (n_pix_int+1, vector<int>(1));

      vector<double> peakVal_vector0;
      vector<double> peakHilbert_vector0; 
      vector<double> SNR_vector0;
      
      vector<double> healpix_bin_weights0; 
      vector<double> bedmap_flagged0;

      vector<double> ratio_vector0; 
      vector<double> polFrac_vector0;
      vector<int> eventNumber_vector0;
      double BASE0;
      // vector< vector <double> > areas_vector (n_pix_int_1, vector<double>(1));

      vector<double> SNR20;
      vector<double> peakVal20;
      vector<double> bedmap_holder(3,0);
      vector <int> pixelnum_vector;
      vector <int> area_pix;
      vector <double> areas;
      int n_old=0;
      int doublebinctr=0;
      int doublebinctr_cone=0;
      int mountainctr=0;
      int mountainctr2bin=0;
      int small_frac_ctr=0;
      double max_distance=0.;
      double min_distance=100;
      double distance=0.;
      double distance_perp=0.;
      double distance_hyp=0.;
      double Lon_furthest;
      double Lat_furthest;
      double x_furthest;
      double y_furthest;
      int ice_bin=0;
      double x_point;
      double y_point;
      double anita_pointx;
      double anita_pointy;
      int flag_ctr=0;
      int twobinFlag=0;
      int mountainFlag=0;
      int shadowctr=0;
      int shadowctr2=0;
      int twobinctr=0;
      int binctr=0;
      int pixel_far_ellipse;
      int pixel_bedmap_bin;
      int num_binned=0;

      double anitaLat_plot;
      double anitaLon_plot;
      double eventLat_plot;
      double eventLon_plot;
      double Lats_plot[100];
      double Lons_plot[100];

      int areas_frac=0;

      int pointctr=0;
      int pointctr1=0;
      TH1D *hdistance_from_end = new TH1D("distance_from_end",";Distance (km);Number of Bins",50,0,20);
      
      TH2F *hboundary = new TH2F("boundary",";X_coord;Y_coord",20,-10.,10.,20,-10.,10.);
      TH2F *hellipse = new TH2F("ellipse",";X_coord;Y_coord",20,-10.,10.,20,-10.,10.);

      double snrCoherent;

      hellipse->SetMarkerColor(kRed);
      //      cout<<"number of events is 10% sample is "<<nevents0<<"\n";

      ///////////////get boundary values//////////////
      double step_size = 2*PI/100.;
      double boundary_val;
      double boundary_val1;
      int offset;
      
     
      vector<double> Lats_boundary_tmp;
      vector<double> Lons_boundary_tmp;
      vector<double> Lats_boundary1_tmp;
      vector<double> Lons_boundary1_tmp;
      
      vector< vector<double> >Lats_boundary;
      vector< vector<double> >Lons_boundary;
      
      vector< vector<double> >Lats_boundary1;
      vector< vector<double> >Lons_boundary1;
           

  ////////////////////
  //read in root file
  string filter_name;

  double weight;
   UInt_t eventNumber_icemc;
   int n_icemc=0;
   double lat_icemc;
   double lon_icemc;
  
   double arrival_times[48];
   double peak2peak_max=0.;
  
   string temp_icemc = "/data/anita/btdailey/icemc_events/Kotera_max3/SimulatedAnitaEventFile.root";
   char *rootfile_icemc = Form(temp_icemc.c_str());
   TChain *event_Tree = new TChain("eventTree");
   event_Tree->Add(rootfile_icemc);
   event_Tree->SetBranchAddress("inu",&eventNumber_icemc);
   event_Tree->SetBranchAddress("weight",&weight);
   event_Tree->SetBranchAddress("Lat_icemc",&lat_icemc);
   event_Tree->SetBranchAddress("Lon_icemc",&lon_icemc);
   //event_Tree->SetBranchAddress("arrival_times",arrival_times);
   int nevents_icemc = event_Tree->GetEntries();

 

  //read in root file
   
  int  run=1;
   // string temp = "/data/anita/btdailey/final_filter/10sample/geom_4pol_partial_1213/output%d_0.root";
   string temp = "/data/anita/btdailey/final_filter/simulated/Kotera_4pol_partial_large_0301/output%d_0.root";
  
  
  char *rootfile = Form(temp.c_str(),run);
  
 
  TChain *pol4_Tree = new TChain("analysis_info_4pol");
  TChain *Pointed_Tree =  new TChain("tdataTracing");
  pol4_Tree->Add(rootfile);
  Pointed_Tree->Add(rootfile);
  //extras
 for(int extra=3000;extra<330000;extra+=3000){
    
       sprintf(rootfile,"/data/anita/btdailey/final_filter/simulated/Kotera_4pol_partial_large_0301/output%d_%d.root",run,extra);
	 TFile *fpTest = TFile::Open(rootfile);
	 if(!fpTest){ 
	   //cout<<"broke at extra "<<extra<<"\n";
	   
	   break;
	 }
	 else {
	   delete fpTest;
	   
	  
	  
	   pol4_Tree->Add(rootfile);
	   Pointed_Tree->Add(rootfile);
	 }
       }//extra
  
     
  int nevents0 = pol4_Tree->GetEntries();
  
  pol4_Tree->SetBranchAddress("pol4_Ptr",&pol4_Ptr);

  ///////////healpix stuff
  vector<double> *peakVal_vector2=0;
  vector<double> *peakHilbert_vector2=0;
  vector<double> *SNR_vector2=0;
  
  vector<double> *healpix_bin_weights2=0;
  vector<int> *eventNumber_vector2=0;
  vector<int> skipped_bins;

  vector< vector <double> > peakVal_vector (n_pix_int+1, vector<double>(1));
  vector< vector <double> > peakHilbert_vector (n_pix_int+1, vector<double>(1));
  vector< vector <double> > SNR_vector (n_pix_int+1, vector<double>(1));
  
  vector< vector <double> > healpix_bin_weights (n_pix_int+1, vector<double>(1));
  vector< vector <int> > eventNumber_vector (n_pix_int+1, vector<int>(1));

  //temp = "HealPix_partial_1215.root";
  temp = "HealPix_sim_partial_0301.root";
 
  rootfile = Form(temp.c_str(),filter_name.c_str());
    
  TChain *Binned_tree = new TChain("Binned_Tree");
 
       
 
  Binned_tree->Add(rootfile);
    
     
  int nevents1 = Binned_tree->GetEntries();
  Binned_tree->SetBranchAddress("eventNumber_vector",&eventNumber_vector2);
  Binned_tree->SetBranchAddress("peakVal_vector",&peakVal_vector2);
  Binned_tree->SetBranchAddress("peakHilbert_vector",&peakHilbert_vector2);
  Binned_tree->SetBranchAddress("SNR_vector",&SNR_vector2);
  Binned_tree->SetBranchAddress("healpix_bin_weights",&healpix_bin_weights2);
  Binned_tree->SetBranchAddress("BASE",&BASE2);
  cout<<"nevents1 is "<<nevents1<<"\n";


  for(int m=2990;m<n_pix_int;m++){
  // for(int m=3058;m<3059;m++){
     Binned_tree->GetEvent(m);
    
     BASE[m]=BASE2;
     for(int n=1;n<(*peakVal_vector2).size();n++)
       {
	 peakVal_vector[m].push_back((*peakVal_vector2)[n]);
	 peakHilbert_vector[m].push_back((*peakHilbert_vector2)[n]);
	 SNR_vector[m].push_back((*SNR_vector2)[n]);
	 healpix_bin_weights[m].push_back((*healpix_bin_weights2)[n]);
	 eventNumber_vector[m].push_back((*eventNumber_vector2)[n]);
	 // if(SNR_vector[m][n]>300) cout<<"eventnumber "<<eventNumber_vector[m][n]<<" has SNR = "<<SNR_vector[m][n]<<"\n";
	 //if((*SNR_vector2)[n] >300) cout<<"eventNumber is "<<(*eventNumber_vector2)[n]<<" has SNR ="<<(*SNR_vector2)[n]<<"\n";

	 if(eventNumber_vector[m][n]==8441658) cout<<"event is "<<eventNumber_vector[m][n]<<" pixel is "<<m<<" weight is "<<healpix_bin_weights[m][n]<<"\n";
       }
     if(BASE[m]>0) cout<<"pix is "<<m<<" num events is "<<BASE[m]<<"\n";
     //if(BASE[m] > max_base) max_base = BASE[m];
    /*
    for(int k=0;k<5;k++){
      cout<<"peakVal is "<<peakVal_vector[m][k]<<" snr is "<<SNR_vector[m][k]<<"\n";   
    }
    */
  }


  long pixel_num;
  long pixel_num_old;
  long pixel_num_event;
  double theta;
  double phi;
     //pix2ang_nest(n_side,pixel_num,theta,phi);
  int eventNumber_cut=0;
  double peakHilbertCoherent=0.;
  double peakVal=0.;
  double ratioFirstToSecondPeak=0.;
  double thetaMap=0.;
  double phiMap=0.;
  double polFractionCoherent=0.;
  double snrcoherent=0.;
  int hwTriggerFlag=0;
  int varnerFlag=0;
  int eventTracedFlag2=0;

   int varnerFlag2;
  
  double lon;
  double lat;

  double anitaLat;
  double anitaLon;
  double anitaAlt;
  double anitaHeading;
  double distance_from_source;
  
 
  double healpixctr=0;
      double payloadBlastctr=0;
   double hilbertctr=0;
   double ratiopeaksctr=0;
   double crosscorrctr=0;
   double polfractionctr=0;
   double rotatedctr=0;
   double elevationctr=0;
   double triggeredctr=0;
   double tracedctr=0;
   double badNoisectr=0;
   
   double varnerevents=0;
   double datactr=0.;
   double datactr2=0;

   double area_ctr=0.; 

   int badNoiseFlag;
     double ratio_last=0;
   double peakVal_last=0;
   double hilbert_last=0;
   double polfrac_last=0;
   double rotated_last=0;
   double elevation_last=0;
   double traced_last=0;
   double triggered_last=0;
   double varner_last=0;
   double badNoise_last=0;

   double snrPeak;
   
   float didIFilter;
   float didIFilterHoriz;
   float didIFilterAboveSatelliteHoriz;
   float didIFilterAboveSatellite;

   int mainrfcmflag;
     int bigenoughpeakflag;
     int dcoffsetflag;
     int shorttraceflag;
     int nadirrfcmflag;

    
     int numbins = 50;
   int index_mult = 1;
   int index=0;
   fstream myfile;
   int healpix_bin;
   double cut_val;
   double y_int;
   vector<int> healpix_bin_vector;
   vector<double> cut_val_vector;
   vector<double> num_pass_in_bin(3070,0.);

   int rotatedFlag=0;
   myfile.open("Cut_values_0501.txt");
   
   while (myfile >>healpix_bin >> cut_val){
     healpix_bin_vector.push_back(healpix_bin);
     cut_val_vector.push_back(cut_val);

   }
   cout<<"number of events in sample is "<<nevents0<<"\n";
      //////////////////
   int inserted_flag=0;
   int polarization=0;//Vpol=0, HPol=1;
   for(int m=0;m<nevents0;m++){
     
     //for(int m=849360;m<nevents0;m++){
     //for(int m=1369290;m<nevents0;m++){
       // for(int m=0;m<0;m++){
       	if (m % (nevents0 / 100) == 0){
	  cout << m << " events. " <<(double(m)/double(nevents0)) * 100 << "% complete.\n";
	}
	//get events
	//cout<<"m is "<<m<<"\n";
      pol4_Tree->GetEvent(m);
      Pointed_Tree->GetEvent(m);
       
       //cout<<"Anita is at "<<AnitaLat<<" "<<AnitaLon<<" "<<AnitaAlt<<"\n";
       eventNumber_cut = pol4_Ptr->eventNumber;
       for(int n=n_icemc;n<nevents_icemc;n++){
	 event_Tree->GetEvent(n);
	 //cout<<"n is "<<n<<"\n";
	 //cout<<"eventNumber_icemc is "<<eventNumber_icemc<<" weight is "<<weight<<"\n";
	 if(eventNumber_cut==eventNumber_icemc){
	   n_icemc = n;
	   break;
	 }
       }//n
       area_pix.clear();
       areas.clear();
       //cout<<"pixel_event_num was "<<pixel_num_event<<" ";
       	for(int j=3008;j<n_pix_int;j++){
	  //cout<<"j is "<<j<<"\n";
	  // cout<<"area_pix, areas size is "<<area_pix.size()<<" "<<areas.size()<<"\n";
	  for(int k=0;k<eventNumber_vector[j].size();k++){
	    //cout<<"eventnumber is "<<eventNumber_vector[j][k]<<"\n";
	    if(eventNumber_vector[j][k]==eventNumber_cut){
	      area_pix.push_back(j);
	      areas.push_back(healpix_bin_weights[j][k]);
	      // cout<<"eventnumber is "<<eventNumber_vector[j][k]<<" "<<eventNumber_cut<<" pix and area is "<<j<<" "<<healpix_bin_weights[j][k]<<"\n";
	    }
	  }
	}



	for(int j=0;j<area_pix.size();j++){
	  if(j==0){
	    pixel_num_event = area_pix[0];
	   
	  }

	  if(areas[j] > areas[0]){
	    pixel_num_event = area_pix[j];
	  }
	  //cout<<"changed to "<<pixel_num_event;
	}
	//	cout<<"\n";
	
	//cout<<"pixel_event_num is now "<<pixel_num_event<<"\n";

       inserted_flag=0;
       if(eventNumber_cut == 4586631 || eventNumber_cut == 8381355 || eventNumber_cut == 9362397 ||  eventNumber_cut == 10345208 || eventNumber_cut == 11207106 || eventNumber_cut == 12943715 || eventNumber_cut == 13662401 || eventNumber_cut == 15406954 || eventNumber_cut == 19480638 || eventNumber_cut == 20564174 || eventNumber_cut == 24699044 || eventNumber_cut == 25887362){
	 inserted_flag=1;
       }      
       peakHilbertCoherent=pol4_Ptr->peakHilbertCoherent[polarization];
       peakVal=pol4_Ptr->peakVal[polarization];
       ratioFirstToSecondPeak=pol4_Ptr->ratioFirstToSecondPeak[polarization];
       thetaMap=pol4_Ptr->thetaMap[polarization];
       phiMap = pol4_Ptr->phiMap[polarization];
       snrcoherent = pol4_Ptr->SNR_coherent[polarization];
       polFractionCoherent=pol4_Ptr->polFractionCoherent[polarization];
       anitaLat = pol4_Ptr->anitaLat;
       anitaLon = pol4_Ptr->anitaLon;
       anitaAlt = pol4_Ptr->anitaAlt;
       anitaHeading = pol4_Ptr->heading; 

       anitaLat = 90-anitaLat;
	 anitaLon = anitaLon+180;

       //if(peakHilbertCoherent != peakHilbertCoherent) cout<<"eventnumber is "<<eventNumber_cut<<" snr is "<<snrcoherent<<"\n";
       hwTriggerFlag=pol4_Ptr->hwTriggerFlag[polarization];
       snrPeak = pol4_Ptr->snrPeakAfterFilter[polarization];
       
       varnerFlag = pol4_Ptr->varnerFlag[polarization];
       didIFilter = pol4_Ptr->didIFilter[polarization];
       eventTracedFlag2 = pol4_Ptr->eventTracedFlag[polarization];
       varnerFlag2 = pol4_Ptr->varnerFlag2[polarization];
      
      
       mainrfcmflag=  pol4_Ptr->mainrfcmflag;
       bigenoughpeakflag=pol4_Ptr->bigenoughpeakflag;
       dcoffsetflag=pol4_Ptr->dcoffsetflag;
       shorttraceflag=pol4_Ptr->shorttraceflag;
       nadirrfcmflag=pol4_Ptr->nadirrfcmflag;
       badNoiseFlag = pol4_Ptr->noiseFlag[polarization];


       	lat = pol4_Ptr->sourceLat[polarization];
	lon = pol4_Ptr->sourceLon[polarization];
	
	//////////Get Pixel Num for Event//////
	lat = 90-lat;
	lon = lon+180;
	
	LatLon2phitheta(lat,lon,phi,theta);
	
	SphericaltoCart(phi,theta,x_map,y_map);

	pixel_num=ang2pix_ring(n_side,theta,phi);
	//pixel_num_event = pixel_num;
	
	y_int=100000;

	 for(int i =0;i<healpix_bin_vector.size();i++){
	     if(healpix_bin_vector[i] == pixel_num_event) y_int = cut_val_vector[i];
	   }
	 //cout<<"pixel_num, y_int is "<<pixel_num<<" "<<y_int<<"\n";

	 if(inserted_flag==1){
	  cout<<"event_number is "<<eventNumber_cut<<" values are "<<ratioFirstToSecondPeak<<" "<<peakVal<<" "<<peakHilbertCoherent<<" "<<thetaMap<<" "<<polFractionCoherent<<"\n";

	}


       if ( !mainrfcmflag || shorttraceflag !=0 || dcoffsetflag != 0 || bigenoughpeakflag !=1){
	 cout<<"quality cut event \n";
	 continue;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by quality \n";
       }
	

       //PAYLOAD BLASTS
       if(pol4_Ptr->payloadblastflag ==1){
	 payloadBlastctr+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by payload \n";
	 continue;
       }

       if(y_int> 500){
	 healpixctr+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by wrong bin. Healpix bin is "<<pixel_num_event<<" \n";
	 continue;
       }
       
       index = (int) snrPeak;

      
       if(index >= numbins){
	 index=numbins;
	 
       }

       if( snrcoherent < -38*peakVal + y_int) rotatedFlag=0;
       else rotatedFlag=1;
       

       float limit=0.9;

       
       if (didIFilterHoriz>1 || didIFilter>1){
	 limit=.85;
       }
       
      
       
       //do ctr for  last cuts////

       if(peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && rotatedFlag==1 && polFractionCoherent>0.3 && badNoiseFlag!=1){
	 if(ratioFirstToSecondPeak<=(1/limit)){
	     ratio_last+=weight;//ratio last
	   }
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && rotatedFlag==1 && polFractionCoherent>0.3  && badNoiseFlag!=1){
	  if(peakVal<0.075){
	    peakVal_last+=weight;//peakVal last
	   }
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075  && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && rotatedFlag==1 && polFractionCoherent>0.3  && badNoiseFlag!=1){
	  if(peakHilbertCoherent <=15){
	    hilbert_last+=weight;//hilbert last
	   }
       }
     if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && polFractionCoherent>0.3  && badNoiseFlag!=1){
	 if(rotatedFlag==0){
	   rotated_last+=weight; //rotated last
	 }
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && rotatedFlag==1  && badNoiseFlag!=1){
	 if(polFractionCoherent<=.3){
	   polfrac_last+=weight;//rotated last	
       }
      
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && rotatedFlag==1 && polFractionCoherent>0.3  && badNoiseFlag!=1){
	  if(thetaMap<=-35 || thetaMap>=0){
	    elevation_last+=weight;//elevation last
       }
      
       }

       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && hwTriggerFlag !=0 && varnerFlag !=1 && rotatedFlag==1 && polFractionCoherent>0.3  && badNoiseFlag!=1){
	 if(eventTracedFlag2!=1){
	   traced_last+=weight;//traced Last
	 }
       }

       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && varnerFlag !=1 && rotatedFlag==1 && polFractionCoherent>0.3  && badNoiseFlag!=1){
	 if(hwTriggerFlag==0){
	   //cout<<"eventnumber being cut by hwtrigger is "<<eventNumber<<" hwTriggerFlag is "<<hwTriggerFlag<<"\n";
	   triggered_last+=weight;//triggered_last
	 }
       }
        if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag!=0 && rotatedFlag==1 && polFractionCoherent>0.3  && badNoiseFlag!=1){
	 if(varnerFlag==1){
	   //cout<<"eventnumber being cut by hwtrigger is "<<eventNumber<<" hwTriggerFlag is "<<hwTriggerFlag<<"\n";
	   varner_last+=weight;//triggered_last
	 }
       }
	if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag!=0 && rotatedFlag==1 && polFractionCoherent>0.3  && varnerFlag!=1){
	 if(badNoiseFlag==1){
	   //cout<<"eventnumber being cut by hwtrigger is "<<eventNumber<<" hwTriggerFlag is "<<hwTriggerFlag<<"\n";
	   badNoise_last+=weight;//triggered_last
	 }
       }
      
       



       /////////////////////////////////////Stuff for plotting and in order cuts
      

       if(ratioFirstToSecondPeak <=(1/limit)){
	 ratiopeaksctr+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by ratio \n";
	 continue;
       }
       
       if(peakVal<.075){
	 //	 if(peakVal<.03){
	 crosscorrctr+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by peakVal \n";
	 continue;
       }
     
       if(peakHilbertCoherent <=15){
	 hilbertctr+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by hilbert \n";
	 continue;
       }
       	     
       
        if(polFractionCoherent<=.3){
	 polfractionctr+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by polfrac \n";
	 continue;
       }
       
	if(snrcoherent < -38*peakVal + y_int){
	  rotatedctr+=weight;
	  if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by rotated. pixel_num and y_int were "<<pixel_num_event<<" "<<y_int<<" \n";
	  continue;
	}
       
      
      

       if(thetaMap<=-35 || thetaMap>=0){
	 elevationctr+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by elevation \n";
	 continue;
       }
      
  
        if(eventTracedFlag2!=1){
	 tracedctr+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by traced \n";
	 continue;
       }
      
  
       if(hwTriggerFlag==0){
	 triggeredctr+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by triggered \n";
	 continue;
       }
      
   
       if(varnerFlag==1){// || varnerFlag2==1){
	 varnerevents+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by varner \n";
	 continue;
       }

        if(badNoiseFlag==1){
	 badNoisectr+=weight;
	 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by bad noise \n";
	 continue;
       }
	if(areas.size() ==0){
	  cout<<"event "<<eventNumber_cut<<" failed placing in bin (on horizon), remove for same reasons! \n";
	  healpixctr+=weight;
	  continue;
	}
	cout<<"passed cuts "<<eventNumber_cut<<" "<<pixel_num_event<<"\n";
	cout<<"ratio of peaks is "<<ratioFirstToSecondPeak<<"\n";
	cout<<"polfrac is "<<polFractionCoherent<<"\n";
	//area_pix.clear();
	//areas.clear();

	/*	for(int j=3008;j<n_pix_int;j++){
	  //cout<<"j is "<<j<<"\n";
	 
	  for(int k=0;k<eventNumber_vector[j].size();k++){
	    //cout<<"eventnumber is "<<eventNumber_vector[j][k]<<"\n";
	    if(eventNumber_vector[j][k]==eventNumber_cut){
	      area_pix.push_back(j);
	      areas.push_back(healpix_bin_weights[j][k]);
	      cout<<"eventnumber is "<<eventNumber_vector[j][k]<<" "<<eventNumber_cut<<" pix and area is "<<j<<" "<<healpix_bin_weights[j][k]<<"\n";
	    }
	  }
	}
	*/
	

	  
	datactr2+=weight;

	

	for(int k=0;k<areas.size();k++){
	   y_int=10000000;
	   if(area_pix[k]>0){
	     cout<<"pixel_num_event, area_pix, area is "<<pixel_num_event<<"("<<y_int<<") "<<area_pix[k]<<" "<<areas[k]<<"\n";
	     for(int i =0;i<healpix_bin_vector.size();i++){
	       if(healpix_bin_vector[i] == area_pix[k]) y_int = cut_val_vector[i];
	   }
	     cout<<"y_int now "<<y_int<<"\n";
	     if(snrcoherent >= -38*peakVal + y_int){
	       if(areas[k] <.5){
		 area_ctr+=areas[k];
		 cout<<"failed area cut! \n";
		 if(inserted_flag==1) cout<<"inserted event "<<eventNumber_cut<<" cut by area \n";
		 continue;
	       }//<.5
	       datactr+=weight*areas[k];
	        num_pass_in_bin[area_pix[k]]+=weight*areas[k];
	       
	       cout<<"event "<<eventNumber_cut<<" pixel_num is "<<area_pix[k]<<" weight is "<<weight*areas[k]<<" datactr is now "<<datactr<<"\n";
	     }
	     else{
	       cout<<"failed rotated cut! \n";
	     }
	     
	   }
	 }

   }//m

    cout<<"payload blast cut is: "<<payloadBlastctr<<"\n";
    cout<<"Number of events in 10% sample: "<<nevents0 - payloadBlastctr<<"\n";
    cout<<"not in correct healpix "<<healpixctr<<"\n";
   cout<<"ratio of peaks cut is: "<<ratiopeaksctr<<" :"<<ratio_last<<" \n";
   cout<<"cross corr cut is: "<< crosscorrctr<<" :"<<peakVal_last<<"\n";
   cout<<"number of events with hilbert <15 is: "<<hilbertctr<<" :"<<hilbert_last<<"\n";
   cout<<"polfraction cut is: "<< polfractionctr<<" :"<<polfrac_last<<"\n";
   cout<<"rotated cross cut is: "<<rotatedctr<<" :"<<rotated_last<<" \n";
   cout<<"elevation cut is: "<<elevationctr<<" :"<<elevation_last<<"\n";
   cout<<"traced cut is: "<<tracedctr<<" :"<<traced_last<<"\n";
   cout<<"triggered cut is: "<<triggeredctr<<" :"<<triggered_last<<"\n";
   cout<<"varner events is: "<<varnerevents<<" :"<<varner_last<<"\n";
   cout<<"badnoisectr is: "<<badNoisectr<<" :"<<badNoise_last<<"\n";
   cout<<"not enough area in event: "<<area_ctr<<"\n";
   cout<<"number that pass all cuts is: "<<datactr2<<" frac is "<<datactr<<"\n";

   for(int k=3008;k<3070;k++){
     cout<<"bin "<<k<<" had "<<num_pass_in_bin[k]<<" passing events \n";
   }
   
 
  return 0;
}
////////////////
void FractionalArea(vector <vector <double> > phis_theta,vector<int> error_bin,vector<int> &area_pix, vector<double> &areas){
  //sort(phis_theta.begin(),phis_theta.end());

  int pixel_num;

  int sizer = (int) phis_theta.size();
     
  for(int i =0;i<sizer;i++){
    //cout<<"phis and theta are "<<phis_theta[i][0]<<" "<<phis_theta[i][1]<<"\n";

  }
  double phi_start = phis_theta[0][0];
  double phi_end = phis_theta[sizer-1][0];
  double theta_start = phis_theta[0][1];
  double theta_end = phis_theta[sizer-1][1];
  
  vector< double > top_half_phi;
  vector< double > top_half_theta;
  
  vector< double > bottom_half_phi;
  vector< double > bottom_half_theta;
  double distance1;
  double distance2;
  int top_ctr=0;
  int bottom_ctr=0;
  
  top_half_phi.push_back(phis_theta[0][0]);
  top_half_theta.push_back(phis_theta[0][1]);
  
  bottom_half_phi.push_back(phis_theta[0][0]);
  bottom_half_theta.push_back(phis_theta[0][1]);
  
  
  if(phis_theta[1][1]>phis_theta[2][1]){//set next point to differentiate between top and bottom
    
    top_half_phi.push_back(phis_theta[1][0]);
    top_half_theta.push_back(phis_theta[1][1]);
    
    bottom_half_phi.push_back(phis_theta[2][0]);
    bottom_half_theta.push_back(phis_theta[2][1]);
    
  }
  else{
    
    top_half_phi.push_back(phis_theta[2][0]);
    top_half_theta.push_back(phis_theta[2][1]);
    
    bottom_half_phi.push_back(phis_theta[1][0]);
    bottom_half_theta.push_back(phis_theta[1][1]);
    
  }
  
  top_ctr=1;
  bottom_ctr=1;
  
  for(int k=3;k<100;k++){//separate the points
    distance1 = abs(phis_theta[k][1]-top_half_theta[top_ctr]);
    distance2 = abs(phis_theta[k][1]-bottom_half_theta[bottom_ctr]);
    
    if(distance2 > distance1){
      
      top_half_phi.push_back(phis_theta[k][0]);
      top_half_theta.push_back(phis_theta[k][1]);
      top_ctr++;
    }
    else if(distance1 > distance2){
      
      bottom_half_phi.push_back(phis_theta[k][0]);
      bottom_half_theta.push_back(phis_theta[k][1]);
      bottom_ctr++;
    }
    
    
  }//k==3
  
  //create TGraph that will be used for interpolation
  
  /*  cout<<"top half points are: \n";
  for(int i =0;i<top_half_phi.size();i++){
    cout<<top_half_phi[i]<<" "<<top_half_theta[i]<<"\n";

  }

   cout<<"bottom points are: \n";
  for(int i =0;i<bottom_half_phi.size();i++){
    cout<<bottom_half_phi[i]<<" "<<bottom_half_theta[i]<<"\n";

  }
  */
  TGraph *top = new TGraph(top_half_theta.size(),&(top_half_phi[0]),&(top_half_theta[0]));
  TGraph *bottom = new TGraph(bottom_half_phi.size(),&(bottom_half_phi[0]),&(bottom_half_theta[0]));
  
  
  double top_val;
  double bot_val;
  int n_steps=100;
  double step_sizer = (phi_end - phi_start)/n_steps;
  int offset;//needed because the equations of boundary had phi MOD(Pi/2)
  int offset_last=(int)(2*phi_start/PI);
  double boundary_val=0.;
  double boundary_val1=0.;
  double total_area=0.;
  double area1=0.;//area of first pixel we are located in 
  double area2=0.;//area of second pixel we are located in
  double area3=0.;
  double area4=0.;
  
  // vector<int> area_pix;//to hold bins numbers we are in
  area_pix.push_back(error_bin[0]);
  area_pix.push_back(0);//2
  area_pix.push_back(0);//3
  area_pix.push_back(0);//4
  
  //vector<double> areas;//to hold fractional areas
  areas.push_back(0);
  areas.push_back(0);
  areas.push_back(0);
  areas.push_back(0);
  
  for(int k=1;k<100;k++){//get all pixel numbers that error ellipse occurs in
    if(area_pix[1]==0){
      if(error_bin[k]!=area_pix[0]){
	area_pix[1]=error_bin[k];
      }
    }
    
    if(area_pix[2]==0){
      if(error_bin[k]!=area_pix[0] && error_bin[k] !=area_pix[1]){
	area_pix[2]=error_bin[k];
      }
    }
    
    if(area_pix[3]==0){
      if(error_bin[k]!=area_pix[0] && error_bin[k] !=area_pix[1] && error_bin[k] !=area_pix[2]){
	area_pix[3]=error_bin[k];
      }
    }
    
  }//k==1
  
  
  
  int pixel_num_first= ang2pix_ring(n_side,theta_start,phi_start);
  cout<<"phi_start and theta_start is "<<phi_start<<" "<<theta_start<<"\n";
  cout<<"pixel_num_first is "<<pixel_num_first<<"\n";
  double n_side_term = 3*pow(n_side,2.);
  double offset_next;
  
  vector <double> boundary_vals;//will hold theta values for boundaries inside error ellipse
  double area;
  double holder_theta=0.;
  //cout<<"phi_start,end are "<<phi_start<<" "<<phi_end<<"\n";
  //loop over all phis in error ellipse and find distances of ellipse theta from boundary theta, use those distances to calculate area in each healpix bin
  for(double nn=phi_start;nn<phi_end;nn+=step_sizer){
    double n = nn;
    
    top_val = top->Eval(n);//find theta value of top half
    bot_val = bottom->Eval(n);//find theta value on bottom half
    
    if(n>2*PI){//if we go above 2PI, wrap back around to 0
      n= n - 2*PI;//phi must be between 0 and 2Pi
      
    }
    boundary_vals.clear();//
    offset = (int)(2*n/PI);
    offset_next = (int)(2*(n+step_sizer)/PI);
    
  
    
    
    for(int k=1;k<=k_value;k++){
      boundary_val = 1-pow(k,2.)/(3*pow(n_side,2.))*pow(PI/(2*(n-offset*PI/2)),2.);
      boundary_val1 = 1-pow(k,2.)/(3*pow(n_side,2.))*pow(PI/(2*(n-offset*PI/2)-PI),2.);
      
      boundary_val = -1*boundary_val;//dealing with southern hemisphere
      boundary_val1 = -1*boundary_val1;
      
      
      boundary_val=acos(boundary_val);
      boundary_val1 = acos(boundary_val1);
      
      if(boundary_val>bot_val && boundary_val<top_val){//lies inside ellipse
	cout<<"boundary k is crossed "<<k<<"\n";
	boundary_vals.push_back(boundary_val);
      }
      if(boundary_val1>bot_val && boundary_val1<top_val){//lies inside ellipse
	cout<<"boundary1 k is crossed "<<k<<"\n";
	boundary_vals.push_back(boundary_val1);
      }
    }//k_values
    sort(boundary_vals.begin(),boundary_vals.end());
    
    if(boundary_vals.size()>0){//boundary cross through ellipse
      for(int k=0;k<boundary_vals.size();k++){
	
	if(k==0){
	  area = boundary_vals[k]-bot_val;
	  holder_theta = boundary_vals[k]+bot_val;
	  holder_theta = holder_theta/2;
	  
	}
	else{
	  area = boundary_vals[k]-boundary_vals[k-1];
	  holder_theta = boundary_vals[k]+boundary_vals[k-1];
	  holder_theta = holder_theta/2;
	}
	
	area = area*step_sizer;
	
	pixel_num=ang2pix_ring(n_side,holder_theta,n);

	for(int kk=0;kk<4;kk++){
	  if(pixel_num==area_pix[kk]){
	    
	    areas[kk]+=area;
	    // cout<<"phi is "<<n<<" pixel_num is "<<pixel_num<<" area being added is "<<area<<"\n";
	  }
	}
	
      }//k
      
      area = top_val - boundary_vals[boundary_vals.size()-1];
      area =area*step_sizer;
      pixel_num=ang2pix_ring(n_side,top_val,n);
      
      for(int kk=0;kk<4;kk++){
	if(pixel_num==area_pix[kk]){
	  
	  areas[kk]+=area;
	  //cout<<"phi is "<<n<<" pixel_num is "<<pixel_num<<" area is "<<area<<"\n";
	}
      }
    }//boundary_vals.size()>0
    else if(offset != offset_next){//crossing pi/2 boundary
      cout<<"offset problem \n";
      double partial_step;
      
      partial_step = offset_next*PI/2 - n;
      pixel_num=ang2pix_ring(n_side,top_val,n);
      area = top_val - bot_val;
      area =area*partial_step;
      for(int kk=0;kk<4;kk++){
	
	if(pixel_num==area_pix[kk]){
	  areas[kk]+=area;
	  //cout<<"phi is "<<n<<" pixel_num is "<<pixel_num<<" area is "<<area<<"\n";
	}
      }
      
      holder_theta =n+step_sizer;//find theta at next step to get correct bin
      holder_theta = top->Eval(holder_theta);//evaluate to get theta
      partial_step = n+step_sizer - offset_next*PI/2;
      pixel_num=ang2pix_ring(n_side,holder_theta,n+step_sizer);
      area = top_val - bot_val;
      area =area*partial_step;
      for(int kk=0;kk<4;kk++){
	
	if(pixel_num==area_pix[kk]){
	  areas[kk]+=area;
	}
      }
      
      
    }
    else{
      pixel_num=ang2pix_ring(n_side,top_val,n);
      area = top_val - bot_val;
      area =area*step_sizer;
      
      for(int kk=0;kk<4;kk++){
	
	if(pixel_num==area_pix[kk]){
	  areas[kk]+=area;
	  //cout<<"phi is "<<n<<" pixel_num is "<<pixel_num<<" area is "<<area<<"\n";
	}
      }
    }
    total_area+=(top_val-bot_val)*step_sizer;
    offset_last=offset;
  }//nn
  
   for(int i=0;i<4;i++){

    areas[i]=areas[i]/total_area;
    //cout<<"pix is "<<area_pix[i]<<" has area of "<<areas[i]<<"\n";
  }
  
  delete top;
  delete bottom;



}
////////////////////
////////////////
void FractionalArea_cone(vector<int> error_bin,vector<int> &area_pix, vector<double> &areas, double eta, double H,double surface,double heading, double lat,double lon){
  //sort(phis_theta.begin(),phis_theta.end());

  int pixel_num;
  //cout<<"in frac area cone \n";
  int sizer = (int) ellipse_coordinates.size();
  sort(ellipse_coordinates.begin(),ellipse_coordinates.end());
  double xmin=10.;
  double xmax=-10.;
  double ymin=10.;
  double ymax=-10.;
  for(int i =0;i<sizer;i++){
    
    //cout<<"phis and theta are "<<phis_theta[i][0]<<" "<<phis_theta[i][1]<<"\n";
    
    if(ellipse_coordinates[i][0] <xmin) xmin = ellipse_coordinates[i][0];
    if(ellipse_coordinates[i][1] <ymin) ymin = ellipse_coordinates[i][1];

    if(ellipse_coordinates[i][0] >xmax) xmax = ellipse_coordinates[i][0];
    if(ellipse_coordinates[i][1] >ymax) ymax = ellipse_coordinates[i][1];
    

  }
  
  double x_start = ellipse_coordinates[0][0];
  double x_end = ellipse_coordinates[sizer-1][0];
  double y_start = ellipse_coordinates[0][1];
  double y_end = ellipse_coordinates[sizer-1][1];
  double z_start = ellipse_coordinates[0][2];
  double z_end = ellipse_coordinates[sizer-1][2];

  // cout<<"x_start, end are "<<x_start<<" "<<x_end<<"\n";
  
  vector< double > top_half_x;
  vector< double > top_half_y;
  
  vector< double > bottom_half_x;
  vector< double > bottom_half_y;
  
  vector< double > top_half_z;
  vector< double > bottom_half_z;

  vector<double> boundary_x_tmp;
  vector<double> boundary_y_tmp;

  vector<double> boundary_x_prime_tmp;
  vector<double> boundary_y_prime_tmp;

  double distance1;
  double distance2;
  int top_ctr=0;
  int bottom_ctr=0;
  
  top_half_x.push_back(ellipse_coordinates[0][0]);
  top_half_y.push_back(ellipse_coordinates[0][1]);
  
  bottom_half_x.push_back(ellipse_coordinates[0][0]);
  bottom_half_y.push_back(ellipse_coordinates[0][1]);
  
  top_half_z.push_back(ellipse_coordinates[0][2]);
  bottom_half_z.push_back(ellipse_coordinates[0][2]);
  
  if(ellipse_coordinates[1][1]>ellipse_coordinates[2][1]){//set next point to differentiate between top and bottom
    
    top_half_x.push_back(ellipse_coordinates[1][0]);
    top_half_y.push_back(ellipse_coordinates[1][1]);
    top_half_z.push_back(ellipse_coordinates[1][2]);
    bottom_half_x.push_back(ellipse_coordinates[2][0]);
    bottom_half_y.push_back(ellipse_coordinates[2][1]);
    bottom_half_z.push_back(ellipse_coordinates[2][2]);
  }
  else{
    
    top_half_x.push_back(ellipse_coordinates[2][0]);
    top_half_y.push_back(ellipse_coordinates[2][1]);
    top_half_z.push_back(ellipse_coordinates[2][2]);
    bottom_half_x.push_back(ellipse_coordinates[1][0]);
    bottom_half_y.push_back(ellipse_coordinates[1][1]);
    bottom_half_z.push_back(ellipse_coordinates[2][2]);
  }
  
  top_ctr=1;
  bottom_ctr=1;
  
  for(int k=3;k<100;k++){//separate the points
    distance1 = abs(ellipse_coordinates[k][1]-top_half_y[top_ctr]);
    distance2 = abs(ellipse_coordinates[k][1]-bottom_half_y[bottom_ctr]);
    
    if(distance2 > distance1){
      
      top_half_x.push_back(ellipse_coordinates[k][0]);
      top_half_y.push_back(ellipse_coordinates[k][1]);
      top_half_z.push_back(ellipse_coordinates[k][2]);
      top_ctr++;
    }
    else if(distance1 > distance2){
      
      bottom_half_x.push_back(ellipse_coordinates[k][0]);
      bottom_half_y.push_back(ellipse_coordinates[k][1]);
      bottom_half_z.push_back(ellipse_coordinates[k][2]);
      bottom_ctr++;
    }
    
    
  }//k==3
  
  //create TGraph that will be used for interpolation
  /*
  cout<<"top half points are: \n";
  for(int i =0;i<top_half_x.size();i++){
    cout<<top_half_x[i]<<" "<<top_half_y[i]<<"\n";

  }

   cout<<"bottom points are: \n";
  for(int i =0;i<bottom_half_x.size();i++){
    cout<<bottom_half_x[i]<<" "<<bottom_half_y[i]<<"\n";

  }
  */
  TGraph *top = new TGraph(top_half_x.size(),&(top_half_x[0]),&(top_half_y[0]));
  TGraph *bottom = new TGraph(bottom_half_x.size(),&(bottom_half_x[0]),&(bottom_half_y[0]));
  TGraph *z_cone = new TGraph(top_half_z.size(),&(top_half_x[0]),&(top_half_z[0]));

  TGraph *boundary[k_value];
  TGraph *boundary_prime[k_value]; 
  /*
  for(int i=0;i<boundary_phi[0].size();i++){
    cout<<"cone phi, theta are "<<boundary_phi[0][i]<<" "<<boundary_theta[0][i]<<"\n";
  }
  cout<<"prime \n";
   for(int i=0;i<boundary_phi_prime[0].size();i++){
    cout<<"cone phi', theta' are "<<boundary_phi_prime[0][i]<<" "<<boundary_theta_prime[0][i]<<"\n";
  }
  */
   double x_val;
   double y_val;

   double x_val_next;
   double y_val_next;

   int boundary_flag[k_value];
   int boundary_prime_flag[k_value];
  for(int i=0;i<k_value;i++){
    boundary_x_tmp.clear();
    boundary_y_tmp.clear();
    boundary_flag[i]=0;
    boundary_prime_flag[i]=0;
    for(int k=0;k<boundary_x[i].size()-1;k++){
      x_val = boundary_x[i][k];
      y_val = boundary_y[i][k];

      x_val_next = boundary_x[i][k+1];
      y_val_next = boundary_y[i][k+1];
      
      //if boundary passes through . 
      if(x_val < xmax && x_val_next > xmin){//crosses in phi

	if(y_val <ymax && y_val_next > ymin){//crosses ellipse in theta
	  boundary_x_tmp.push_back(x_val);
	  boundary_x_tmp.push_back(x_val_next);

	  boundary_y_tmp.push_back(y_val);
	  boundary_y_tmp.push_back(y_val_next);
	}


      }

       
      if(x_val > xmin && x_val < xmax && y_val > ymin && y_val < ymax){//point lies inside
	boundary_x_tmp.push_back(x_val);
	boundary_y_tmp.push_back(y_val);
      }  
      
    }

    if(boundary_x_tmp.size()>0) boundary_flag[i]=1;
    boundary[i] = new TGraph(boundary_x_tmp.size(),&boundary_x_tmp[0],&boundary_y_tmp[0]);
    

    
    boundary_x_tmp.clear();
    boundary_y_tmp.clear();
    for(int k=0;k<boundary_x_prime[i].size()-1;k++){
      x_val = boundary_x_prime[i][k];
      y_val = boundary_y_prime[i][k];
      
      x_val_next = boundary_x_prime[i][k+1];
      y_val_next = boundary_y_prime[i][k+1];
      
      //if boundary passes through . 
      if(x_val < xmax && x_val_next > xmin){//crosses in phi
	
	if(y_val <ymax && y_val_next > ymin){//crosses ellipse in theta
	  boundary_x_tmp.push_back(x_val);
	  boundary_x_tmp.push_back(x_val_next);
	  
	  boundary_y_tmp.push_back(y_val);
	  boundary_y_tmp.push_back(y_val_next);
	}
      
	
      }
    }
    //SET FLAG FOR IF FILLED!
    if(boundary_x_tmp.size()>0) boundary_prime_flag[i]=1;
    boundary_prime[i] = new TGraph(boundary_x_tmp.size(),&boundary_x_tmp[0],&boundary_y_tmp[0]);
  }
  
  
  double top_val;
  double bot_val;
  int n_steps=100;
  double step_sizer = (x_end - x_start)/n_steps;
  int offset;//needed because the equations of boundary had phi MOD(Pi/2)
  int offset_last=0;//(int)(2*phi_start/PI);
  double boundary_val=0.;
  double boundary_val1=0.;
  double total_area=0.;
  double area1=0.;//area of first pixel we are located in 
  double area2=0.;//area of second pixel we are located in
  double area3=0.;
  double area4=0.;
  
  // vector<int> area_pix;//to hold bins numbers we are in
  area_pix.push_back(error_bin[0]);
  area_pix.push_back(0);//2
  area_pix.push_back(0);//3
  area_pix.push_back(0);//4
  
  //vector<double> areas;//to hold fractional areas
  areas.push_back(0);
  areas.push_back(0);
  areas.push_back(0);
  areas.push_back(0);
  
  for(int k=1;k<100;k++){//get all pixel numbers that error ellipse occurs in
    if(area_pix[1]==0){
      if(error_bin[k]!=area_pix[0]){
	area_pix[1]=error_bin[k];
      }
    }
    
    if(area_pix[2]==0){
      if(error_bin[k]!=area_pix[0] && error_bin[k] !=area_pix[1]){
	area_pix[2]=error_bin[k];
      }
    }
    
    if(area_pix[3]==0){
      if(error_bin[k]!=area_pix[0] && error_bin[k] !=area_pix[1] && error_bin[k] !=area_pix[2]){
	area_pix[3]=error_bin[k];
      }
    }
    
  }//k==1
  
  
  
 
  //int pixel_num_first= ang2pix_ring(n_side,theta_start,phi_start);
  int pixel_num_first = HealPixBinFromCone(x_start,y_start,z_start,eta,H,surface, heading,lat,lon);
  //cout<<"pixel_num_first is "<<pixel_num_first<<"\n";
  double n_side_term = 3*pow(n_side,2.);
  double offset_next=offset;
  
  vector <double> boundary_vals;//will hold theta values for boundaries inside error ellipse
  double area;
  double holder_y=0.;
  double zcone;
  //cout<<"phi_start,end are "<<phi_start<<" "<<phi_end<<"\n";
  //loop over all phis in error ellipse and find distances of ellipse theta from boundary theta, use those distances to calculate area in each healpix bin

 
  for(double nn=x_start;nn<x_end;nn+=step_sizer){
    double n = nn;
    
    top_val = top->Eval(n);//find theta value of top half
    bot_val = bottom->Eval(n);//find theta value on bottom half
    zcone = z_cone->Eval(n);
    
    boundary_vals.clear();//
   
    
    for(int k=0;k<k_value;k++){
      
      if(boundary_flag[k]==1){
	boundary_val = boundary[k]->Eval(n);
	if(boundary_val>bot_val && boundary_val<top_val){//lies inside ellipse
	   boundary_vals.push_back(boundary_val);
	}
      }
      boundary_val1 = boundary_prime[k]->Eval(n);
      //cout<<"boundary_val,1 top, bottom are "<<boundary_val<<" "<<boundary_val1<<" "<<top_val<<" "<<bot_val<<"\n";
      if(boundary_prime_flag[k]==1){
	if(boundary_val1>bot_val && boundary_val1<top_val){//lies inside ellipse
	   boundary_vals.push_back(boundary_val1);
	}
      }
    }//k_values
    sort(boundary_vals.begin(),boundary_vals.end());
    
    if(boundary_vals.size()>0){//boundary cross through ellipse
      for(int k=0;k<boundary_vals.size();k++){
	
	if(k==0){
	  area = boundary_vals[k]-bot_val;
	  holder_y = boundary_vals[k]+bot_val;
	  holder_y = holder_y/2;
	  
	}
	else{
	  area = boundary_vals[k]-boundary_vals[k-1];
	  holder_y = boundary_vals[k]+boundary_vals[k-1];
	  holder_y = holder_y/2;
	}
	
	area = area*step_sizer;
	
	//	pixel_num=ang2pix_ring(n_side,holder_y,n);

	pixel_num =  HealPixBinFromCone(n, holder_y,zcone,eta,H, surface, heading,lat,lon);
	
	for(int kk=0;kk<4;kk++){
	  if(pixel_num==area_pix[kk]){
	    
	    areas[kk]+=area;
	    // cout<<"phi is "<<n<<" pixel_num is "<<pixel_num<<" area being added is "<<area<<"\n";
	  }
	}
	
      }//k
      
      area = top_val - boundary_vals[boundary_vals.size()-1];
      area =area*step_sizer;
      //pixel_num=ang2pix_ring(n_side,top_val,n);
      
      pixel_num =  HealPixBinFromCone(n, top_val,zcone,eta,H, surface, heading,lat,lon);
      //cout<<"pixel_num (top_val) is "<<pixel_num<<"\n";
      for(int kk=0;kk<4;kk++){
	if(pixel_num==area_pix[kk]){
	  
	  areas[kk]+=area;
	  //cout<<"phi is "<<n<<" pixel_num is "<<pixel_num<<" area is "<<area<<"\n";
	}
      }
    }//boundary_vals.size()>0
    else if(offset != offset_next){//crossing pi/2 boundary
      
      double partial_step;
      
      partial_step = offset_next*PI/2 - n;
      pixel_num=ang2pix_ring(n_side,top_val,n);
      area = top_val - bot_val;
      area =area*partial_step;
      for(int kk=0;kk<4;kk++){
	
	if(pixel_num==area_pix[kk]){
	  areas[kk]+=area;
	  //cout<<"phi is "<<n<<" pixel_num is "<<pixel_num<<" area is "<<area<<"\n";
	}
      }
      
      holder_y =n+step_sizer;//find theta at next step to get correct bin
      holder_y = top->Eval(holder_y);//evaluate to get theta
      partial_step = n+step_sizer - offset_next*PI/2;
      pixel_num=ang2pix_ring(n_side,holder_y,n+step_sizer);
      area = top_val - bot_val;
      area =area*partial_step;
      for(int kk=0;kk<4;kk++){
	
	if(pixel_num==area_pix[kk]){
	  areas[kk]+=area;
	}
      }
      
      
    }//offset != offset_next
    else{
      //pixel_num=ang2pix_ring(n_side,top_val,n);
      
      pixel_num= HealPixBinFromCone(n, top_val,zcone,eta,H, surface, heading,lat,lon);
      area = top_val - bot_val;
      area = area*step_sizer;
      
      for(int kk=0;kk<4;kk++){
	
	if(pixel_num==area_pix[kk]){
	  areas[kk]+=area;
	  
	  //cout<<"phi is "<<n<<" pixel_num is "<<pixel_num<<" area is "<<area<<"\n";
	}
      }
    }
    
    total_area+=(top_val-bot_val)*step_sizer;
    offset_last=offset;
  }//nn
  
   for(int i=0;i<4;i++){
    
    areas[i]=areas[i]/total_area;
    
  }
  
  delete top;
  delete bottom;



}
////////////////////
int ErrorPoints(double *xfinal, double *yfinal, double *zfinal,double eta1, double theta, double phi,double H,double surface){
  double eta = eta1*PI/180;
  double R = 6378.1+surface; //km
  H = H-surface;
  double sec2phi = 1/cos(phi);
  double sec2theta = 1/cos(theta);
  double csc2theta = 1/sin(theta);
  double csc2phi = 1/sin(phi);
  sec2phi = pow(sec2phi,2);
  sec2theta = pow(sec2theta,2);
  csc2theta = pow(csc2theta,2);
  csc2phi = pow(csc2phi,2);
  double tan2phi = pow(tan(phi),2);
  double tan2theta = pow(tan(theta),2);
  
  double RplusH2 = pow(R+H,2);
  double R2 = pow(R,2);

  double sqrtterm = pow(R+H,2)-pow(R,2);
  sqrtterm = sqrtterm/(pow(R+H,2)*pow(cos(eta),2));
  sqrtterm*=sec2phi;
  sqrtterm = 1-sqrtterm;
  sqrtterm = sqrt(sqrtterm);
  double maxy = tan(phi)*(R+H)*cos(eta)*pow(cos(phi),2)*(1-sqrtterm);//(1-sqrt(1-sec2phi*(pow(R+H,2)-pow(R,2))/pow((R+H)*cos(eta),2))));
  double miny = -1* tan(phi)*(R+H)*cos(eta)*pow(cos(phi),2)*(1-sqrtterm);
  double stepsize = (maxy-miny)/52;
  //cout<<"maxy, miny is "<<maxy<<" "<<miny<<"stepsize is "<<stepsize<<"\n";
  
  //cout<<"z is "<<miny/tan(phi)<<" "<<maxy/tan(phi)<<"\n";
  double H_temp = H+surface;
  double R_temp = R-surface;
  double distance_to_horizon = sqrt(H*(2*R+H));
  double angle_ground_to_horizon = acos(distance_to_horizon/(R+H));
  
  //cout<<"distance_to_horizon, angle_below is "<<distance_to_horizon<<" "<<angle_ground_to_horizon<<" eta is "<<eta<<"\n";

  double angle_to_horizon_from_sight = angle_ground_to_horizon - eta;

  //cout<<"angle_to_horizon_from_sight is "<<angle_to_horizon_from_sight<<"\n";

  double x_horizon = distance_to_horizon*sin(angle_to_horizon_from_sight);
  double z_horizon = distance_to_horizon*cos(angle_to_horizon_from_sight);
 


  //  cout<<"x_hor is "<<x_horizon<<"\n";


  if (maxy !=maxy || miny !=miny){
    
  }
  double a =pow(csc2theta,2);
  double b;
  double c;
  double d;
  double e;

  // cout<<"co-efficients are "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<"\n";


  double z[4];
  double z_b = 2*(R+H)*cos(eta);
  double z_a = sec2phi;
  double z_c0;
  double z_c1;
  double x[4];
  double zval0;
  double zval1;
  double xval0;
  double xval1;

  double xval_temp0;
  double xval_temp1;
  int output;
  

  double xprime0;
  double yprime0;
  double zprime0;

  double xdoubleprime0;
  double ydoubleprime0;
  double zdoubleprime0;

  double xprime1;
  double yprime1;
  double zprime1;

  double xdoubleprime1;
  double ydoubleprime1;
  double zdoubleprime1;

  int ctr=0;
  int output_return=10;
  TH2D *hxymap = new TH2D("xymap",";x;y",1000,0,300,1000,-300,300);
  double y_val;
  
  for(int ctr1=0;ctr1<50;ctr1++){
    y_val = (ctr1+1)*stepsize + miny;
    
    b = 4*(R+H)*csc2theta*sin(eta);
    c =2*pow(y_val,2)*csc2theta + 2*csc2theta*(RplusH2-R2)+4*RplusH2*pow(sin(eta),2)-4*RplusH2*pow(cos(eta),2)/pow(tan(theta),2);
    d =4*pow(y_val,2)*csc2phi*(R+H)*sin(eta)+4*(R+H)*sin(eta)*(RplusH2-R2);
    e = pow(y_val,4)*csc2phi+2*pow(y_val,2)*(RplusH2-R2)*csc2phi-4*pow(y_val,2)*RplusH2*pow(cos(eta),2)/pow(tan(phi),2)+pow(RplusH2-R2,2);
    
    
   
    //cout<<"co-efficients are "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" ";

   output= SolveP4(x, b/a,c/a,d/a,e/a);
   // cout<<"output is "<<output<<"\n";
   //cout<<"x solns are "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<"\n";

   //if(output<4) cout<<"eta is "<<eta1<<" output is "<<output<<"\n";
    if(output<output_return) output_return = output;
    if(output==2){
      xval0 = 1E10;
      xval1 = 1E10;
      if(abs(x[0])<abs(xval0)) xval0 = x[0];
      if(abs(x[1])<abs(xval0)) xval0 = x[1];
     
      if(abs(x[0])<abs(xval1) && abs(x[0])>abs(xval0)) xval1 = x[0];
      if(abs(x[1])<abs(xval1) && abs(x[1])>abs(xval0)) xval1 = x[1];
      
      
      z_c0 = R2-RplusH2 - 2*(R+H)*sin(eta)*xval0-(1-tan2phi/tan2theta)*pow(xval0,2);
      z_c1 = R2-RplusH2 - 2*(R+H)*sin(eta)*xval1-(1-tan2phi/tan2theta)*pow(xval1,2);
      /*
      z[0] = -z_b + sqrt(pow(z_b,2)-4*z_a*z_c0);
      z[0] = z[0]/(2*z_a);
      
      z[1] = -z_b - sqrt(pow(z_b,2)-4*z_a*z_c0);
      z[1] = z[1]/(2*z_a);
      
      z[2] = -z_b + sqrt(pow(z_b,2)-4*z_a*z_c1);
      z[2] = z[2]/(2*z_a);
      
      z[3] = -z_b - sqrt(pow(z_b,2)-4*z_a*z_c1);
      z[3] = z[3]/(2*z_a);

     
      
      zval0=1E10;
      zval1=1E10;
      if(abs(z[0])<abs(zval0)){ 
	zval0 = z[0];
	xval_temp0 = xval0;
      }
      if(abs(z[1])<abs(zval0)){
	zval0 = z[1];
	xval_temp0 = xval0;
      }
      if(abs(z[2])<abs(zval0)){
	zval0 = z[2];
	xval_temp0 = xval1;
	
      }
      if(abs(z[3])<abs(zval0)){
	zval0 = z[3];
	xval_temp0=xval1;
      }

     
      //cout<<"xval0 si "<<xval0<<"\n";
      
      if(abs(z[0])<abs(zval1) && abs(z[0])>abs(zval0)){
	zval1 = z[0];
	xval_temp1 = xval0;
      }
      if(abs(z[1])<abs(zval1) && abs(z[1])>abs(zval0)){
	zval1 = z[1];
	xval_temp1=xval0;
      }
      if(abs(z[2])<abs(zval1) && abs(z[2])>abs(zval0)){
	zval1 = z[2];
	xval_temp1=xval1;
      }
      if(abs(z[3])<abs(zval1) && abs(z[3])>abs(zval0)){
	zval1 = z[3];
	xval_temp1 = xval1;
      }
      
      xval0=xval_temp0;
      xval1=xval_temp1;
      */
      //cout<<" error is 2, x0,y0,z0 are "<<xval0<<" "<<y_val<<" "<<zval0<<"\n";
      //cout<<" error is 2, x1,y1,z1 are "<<xval1<<" "<<y_val<<" "<<zval1<<"\n";
      zval0 = sqrt(pow(xval0/tan(theta),2) + pow(y_val/tan(phi),2));
      zval1 = sqrt(pow(xval1/tan(theta),2) + pow(y_val/tan(phi),2));
      // cout<<"point is "<<xval0<<" "<<y_val<<" "<<zval0<<"\n";
      //cout<<"secondpoint is "<<xval1<<" "<<y_val<<" "<<zval1<<"\n";

     
      
      xval1 = x_horizon;
      zval1 = z_horizon;

    }
    else{
      xval0 = 1E10;
      xval1 = 1E10;
      if(abs(x[0])<abs(xval0)) xval0 = x[0];
      if(abs(x[1])<abs(xval0)) xval0 = x[1];
      if(abs(x[2])<abs(xval0)) xval0 = x[2];
      if(abs(x[3])<abs(xval0)) xval0 = x[3];
      
      if(abs(x[0])<abs(xval1) && abs(x[0])>abs(xval0)) xval1 = x[0];
      if(abs(x[1])<abs(xval1) && abs(x[1])>abs(xval0)) xval1 = x[1];
      if(abs(x[2])<abs(xval1) && abs(x[2])>abs(xval0)) xval1 = x[2];
      if(abs(x[3])<abs(xval1) && abs(x[3])>abs(xval0)) xval1 = x[3];
      
      // cout<<"xval is "<<xval0<<" "<<xval1<<"\n";
      
      z_c0 = R2-RplusH2 - 2*(R+H)*sin(eta)*xval0-(1-tan2phi/tan2theta)*pow(xval0,2);
      z_c1 = R2-RplusH2 - 2*(R+H)*sin(eta)*xval1-(1-tan2phi/tan2theta)*pow(xval1,2);
      
      z[0] = -z_b + sqrt(pow(z_b,2)-4*z_a*z_c0);
      z[0] = z[0]/(2*z_a);
      
      z[1] = -z_b - sqrt(pow(z_b,2)-4*z_a*z_c0);
      z[1] = z[1]/(2*z_a);
      
      z[2] = -z_b + sqrt(pow(z_b,2)-4*z_a*z_c1);
      z[2] = z[2]/(2*z_a);
      
      z[3] = -z_b - sqrt(pow(z_b,2)-4*z_a*z_c1);
      z[3] = z[3]/(2*z_a);
      

      // cout<<"x,y z0 is "<<xval0<<" "<<y_val<<" "<<sqrt(pow(xval0/tan(theta),2) + pow(y_val/tan(phi),2))<<"\n";
      //cout<<"x,y, z1 is "<<xval1<<" "<<y_val<<" "<<sqrt(pow(xval1/tan(theta),2) + pow(y_val/tan(phi),2))<<"\n";
      
      zval0 = sqrt(pow(xval0/tan(theta),2) + pow(y_val/tan(phi),2));
      zval1 = sqrt(pow(xval1/tan(theta),2) + pow(y_val/tan(phi),2));
      // cout<<"z values are "<<z[0]<<" "<<z[1]<<" "<<z[2]<<" "<<z[3]<<"\n"; 
	/* zval0=1E10;
      zval1=1E10;
      if(abs(z[0])<abs(zval0)){ 
	zval0 = z[0];
	xval_temp0 = xval0;
      }
      if(abs(z[1])<abs(zval0)){
	zval0 = z[1];
	xval_temp0 = xval0;
      }
      if(abs(z[2])<abs(zval0)){
	zval0 = z[2];
	xval_temp0 = xval1;
	
      }
      if(abs(z[3])<abs(zval0)){
	zval0 = z[3];
	xval_temp0=xval1;
      }
      
      //cout<<"xval0 si "<<xval0<<"\n";
      
      if(abs(z[0])<abs(zval1) && abs(z[0])>abs(zval0)){
	zval1 = z[0];
	xval_temp1 = xval0;
      }
      if(abs(z[1])<abs(zval1) && abs(z[1])>abs(zval0)){
	zval1 = z[1];
	xval_temp1=xval0;
      }
      if(abs(z[2])<abs(zval1) && abs(z[2])>abs(zval0)){
	zval1 = z[2];
	xval_temp1=xval1;
      }
      if(abs(z[3])<abs(zval1) && abs(z[3])>abs(zval0)){
	zval1 = z[3];
	xval_temp1 = xval1;
      }
      
      xval0=xval_temp0;
      xval1=xval_temp1;
	*/
      // cout<<"point is "<<xval0<<" "<<y_val<<" "<<zval0<<"\n";
      //cout<<"secondpoint is "<<xval1<<" "<<y_val<<" "<<zval1<<"\n";
    }
    // cout<<"xval0,y,z are "<<xval0<<" "<<y_val<<" "<<zval0<<"\n";
    //cout<<"xval1,y,z are "<<xval1<<" "<<y_val<<" "<<zval1<<"\n";
    //Transform x,y,z from cone coordinates to geo-centric coordinates
    
    xprime0 = xval0+(R+H)*sin(eta);
    yprime0 = y_val;
    zprime0 = -zval0 + (R+H)*cos(eta);

    xdoubleprime0 = xprime0*cos(eta)-zprime0*sin(eta);
    ydoubleprime0 = yprime0;
    zdoubleprime0 = xprime0*sin(eta)+zprime0*cos(eta);
    
    //cout<<"error ellipse x,y,z, are "<<xdoubleprime0<<" "<<ydoubleprime0<<" "<<zdoubleprime0<<"\n";
    double test_anita_x0;
    double test_anita_y0;
    double test_anita_z0;
    test_anita_x=0.;
    test_anita_y=0.;
    test_anita_z=0.;
    test_anita_x0 = test_anita_x+(R+H)*sin(eta);
    test_anita_y0 = test_anita_y;
    test_anita_z0 = -test_anita_z + (R+H)*cos(eta);

    test_anita_x = test_anita_x0*cos(eta)-test_anita_z0*sin(eta);
    test_anita_y = test_anita_y0;
    test_anita_z = test_anita_x0*sin(eta)+test_anita_z0*cos(eta);


    // if (output==4){
      xprime1 = xval1+(R+H)*sin(eta);
      yprime1 = y_val;
      zprime1 = -zval1 + (R+H)*cos(eta);
    
      xdoubleprime1 = xprime1*cos(eta)-zprime1*sin(eta);
      ydoubleprime1 = yprime1;
      zdoubleprime1 = xprime1*sin(eta)+zprime1*cos(eta);
      /* }
    else{
      xdoubleprime1=0;
      ydoubleprime1=0;
      zdoubleprime1=0;
    }
      */

    //cout<<"pts are "<<xdoubleprime0<<" "<<ydoubleprime0<<" "<<zdoubleprime0<<"\n";
    //cout<<"pts are "<<xdoubleprime1<<" "<<ydoubleprime1<<" "<<zdoubleprime1<<"\n";
    xfinal[ctr] = xdoubleprime0;
    yfinal[ctr] = ydoubleprime0;
    zfinal[ctr] = zdoubleprime0;
    
    xfinal[ctr+1] = xdoubleprime1;
    yfinal[ctr+1] = ydoubleprime1;
    zfinal[ctr+1] = zdoubleprime1;


    

    ctr+=2;
    /////////////////////////////
    //convert to lat and lon?

   
    
    
    hxymap->Fill(xdoubleprime0,ydoubleprime0);
    hxymap->Fill(xdoubleprime1,ydoubleprime1);



  }
  /*cout<<"anita x,y,z, are "<<test_anita_x<<" "<<test_anita_y<<" "<<test_anita_z<<"\n";
  cout<<"x, y, z, is : \n";
  for(int i=0;i<100;i++){
    cout<<xfinal[i]<<" "<<yfinal[i]<<" "<<zfinal[i]<<"\n";
    
  }
*/
  //cout<<"Ctr is "<<ctr<<"\n";
  char printer[250];
  // cout<<"eta is "<<eta1<<"\n";
  /* TCanvas *c0 = new TCanvas("c0","c0",800,800);
  hxymap->Draw();
  sprintf(printer,"errorellipse_%f.png",eta1);
  c0->Print(printer);
  */
  delete hxymap;
  
  return output_return;
}

int HealPixBinFromCone(double x, double y, double z,double eta,double H,double surface,double heading, double lat,double lon){
  eta = eta*PI/180;
  
  double R = 6378.1+surface; //km
  H = H-surface;
 
  
  double xprime;
  double yprime;
  double zprime;

  double xdoubleprime;
  double ydoubleprime;
  double zdoubleprime;
  

  xprime = x+(R+H)*sin(eta);
  yprime = y;
  zprime = -z + (R+H)*cos(eta);
  
  xdoubleprime = xprime*cos(eta)-zprime*sin(eta);
  ydoubleprime = yprime;
  zdoubleprime = xprime*sin(eta)+zprime*cos(eta);

 
  lat = 90-lat;
  
  lat = lat*PI/180.;
  lon = lon*PI/180.;
 
  heading =heading *PI/180.;
  
  double u;
  double uprime;
  double v;
  double vprime;
  double w;
  double wprime;
 
  double Lat_holder;
  double Lon_holder;
  double R1;
  u = ydoubleprime*sin(heading) + xdoubleprime*cos(heading);
  v = ydoubleprime*cos(heading) - xdoubleprime*sin(heading);
  w = zdoubleprime;
 
  uprime = u*sin(lat) - w*cos(lat);
  vprime = v;
  wprime = u*cos(lat) + w*sin(lat);
  R1 = sqrt(pow(uprime,2)+pow(vprime,2)+pow(wprime,2));
 
  Lat_holder = asin(wprime/R)*180/PI;
  Lat_holder = 90-Lat_holder;
  
  Lon_holder = asin(vprime/sqrt(pow(uprime,2)+pow(vprime,2)));
   if(uprime <0){
     Lon_holder=PI-Lon_holder;


   }
  
  Lon_holder = (Lon_holder+lon)*180/PI;
  
  if(Lon_holder>360) Lon_holder-=360;
  
  double phi1;
  double theta1;
  
  LatLon2phitheta(Lat_holder,Lon_holder-180,phi1,theta1);//90-Lats,180+Lons
  
 int  pixel_num=ang2pix_ring(n_side,theta1,phi1);
 return pixel_num;
}



void ToConeCoordinates(double *xfinal, double *yfinal, double *zfinal,double eta1, double *x,double *y,double *z,double H,double surface){
  double eta = eta1*PI/180;
  double R = 6378.1+surface; //km
  H = H-surface;
  
  double xprime0;
  double yprime0;
  double zprime0;

  double xdoubleprime0;
  double ydoubleprime0;
  double zdoubleprime0;

  for(int n=0;n<100;n++){
    xprime0=x[n]*cos(eta)+z[n]*sin(eta);
    yprime0=y[n];
    zprime0=-x[n]*sin(eta)+z[n]*cos(eta);

    xfinal[n] = xprime0 -(R+H)*sin(eta);
    yfinal[n] = yprime0;
    zfinal[n] = -zprime0 +(R+H)*cos(eta);


  }
  

}
 void ToConeCoordinates_vec(vector<double> &xfinal, vector<double> &yfinal, vector<double> &zfinal,double eta1, vector<double> x,vector<double> y,vector<double> z,double H,double surface){
  double eta = eta1*PI/180;
  double R = 6378.1+surface; //km
  H = H-surface;
  
  double xprime0;
  double yprime0;
  double zprime0;

  double xdoubleprime0;
  double ydoubleprime0;
  double zdoubleprime0;

  for(int n=0;n<100;n++){
    xprime0=x[n]*cos(eta)+z[n]*sin(eta);
    yprime0=y[n];
    zprime0=-x[n]*sin(eta)+z[n]*cos(eta);

    xfinal[n] = xprime0 -(R+H)*sin(eta);
    yfinal[n] = yprime0;
    zfinal[n] = -zprime0 +(R+H)*cos(eta);


  }
  

}

void GetLatLonfromCartesian(double *Lats, double *Lons, double heading, double lat, double lon, double *x, double *y, double *z){
  //cout<<"lat is "<<lat<<"\n";
  //cout<<"90-lat is "<<90-lat<<"\n";
  lat = 90-lat;
  //lon = lon-180;
  lat = lat*PI/180.;
  lon = lon*PI/180.;
  //cout<<"heading is "<<heading<<"\n";
  heading =heading *PI/180.;
  
  double u;
  double uprime;
  double v;
  double vprime;
  double w;
  double wprime;
  double R;
  //cout<<" anita x,y,z are "<<test_anita_x<<" "<<test_anita_y<<" "<<test_anita_z<<"\n";
  R = sqrt(pow(test_anita_x,2)+pow(test_anita_y,2)+pow(test_anita_z,2));
   
  u = test_anita_y*sin(heading) + test_anita_x*cos(heading);
  v = test_anita_y*cos(heading) - test_anita_x*sin(heading);
  w = test_anita_z;
  //cout<<"anita u v w are "<<u<<" "<<v<<" "<<w;
  
  
  uprime = u*sin(lat) - w*cos(lat);
  vprime = v;
  wprime = u*cos(lat) + w*sin(lat);
  //cout<<" anita u,v,w' are "<<uprime<<" "<<vprime<<" "<<wprime<<"\n";
  test_anita_lat = asin(wprime/R)*180/PI;
  test_anita_lat = 90-test_anita_lat;
  //cout<<" u' v' w' are "<<uprime<<" "<<vprime<<" "<<wprime<<"\n";
  
  test_anita_lon = asin(vprime/sqrt(pow(uprime,2)+pow(vprime,2)));
  
  if(uprime <0){
    //test_anita_lon=PI-test_anita_lon;


      }
    test_anita_lon = (test_anita_lon+lon)*180/PI;
    if(test_anita_lon > 360) test_anita_lon-=360;

    // cout<<"heading is "<<heading<<"\n";
    //cout<<"sin heading, cos heading are "<<sin(heading)<<" "<<cos(heading)<<"\n";
  for(int n=0;n<100;n++){
    
    //cout<<"x,y,z are "<<x[n]<<" "<<y[n]<<" "<<z[n];
    R = sqrt(pow(x[n],2)+pow(y[n],2)+pow(z[n],2));
   
    u = y[n]*sin(heading) + x[n]*cos(heading);
    v = y[n] *cos(heading) - x[n]*sin(heading);
    w = z[n];
    //cout<<" u v w are "<<u<<" "<<v<<" "<<w;
    
    
    uprime = u*sin(lat) - w*cos(lat);
    vprime = v;
    wprime = u*cos(lat) + w*sin(lat);
    
    
    Lats[n] = asin(wprime/R)*180/PI;
    //cout<<" u' v' w' are "<<uprime<<" "<<vprime<<" "<<wprime<<"\n";
   
    Lons[n] = asin(vprime/sqrt(pow(uprime,2)+pow(vprime,2)));
    /*if(uprime >0){
      Lons[n] = PI-Lons[n];
      }*/
    
     if(uprime <0){
       Lons[n]=PI-Lons[n];


      }
    Lons[n] = (Lons[n]+lon)*180/PI;
    if(Lons[n] > 360) Lons[n]-=360;
   
    //cout<<"x,y,z are "<<x[n]<<" "<<y[n]<<" "<<z[n];
    //cout<<" primes are "<<uprime<<" "<<vprime<<" "<<wprime;
    // cout<<" Lat and Lon are  "<<Lats[n]<<" "<<Lons[n]<<"\n";
    //cout<<" Lat and Lon is "<<Lats[n]<<" "<<Lons[n]<<"\n";
   
  }
  

}
void GetCartesianfromLatLon(double *Lats, double *Lons, double heading, double lat, double lon, double *x, double *y, double *z,double R){
  //cout<<"lat is "<<lat<<"\n";
  //cout<<"90-lat is "<<90-lat<<"\n";
  lat = 90-lat;
  //lon = lon-180;
  lat = lat*PI/180.;
  lon = lon*PI/180.;
  //cout<<"heading is "<<heading<<"\n";
  heading =heading *PI/180.;
  R = 6378.1+R; //km
  double u;
  double uprime;
  double v;
  double vprime;
  double w;
  double wprime;
  
  double lat_holder;
  double lon_holder;
  for(int n=0;n<100;n++){
 
    lat_holder = Lats[n];
    lon_holder = Lons[n];

    lat_holder *=PI/180.;
    lon_holder *=PI/180.;
    lon_holder -=lon;
    
    
    wprime = R*sin(lat_holder);
    double horizontal = sqrt(pow(R,2)-pow(wprime,2));
    vprime = sin(lon_holder)*horizontal;
    uprime = cos(lon_holder)*horizontal;
    
    v = vprime;
    u = uprime*sin(lat) + wprime*cos(lat);
    w = -uprime*cos(lat) + wprime*sin(lat);
    
    x[n] = u*cos(heading) - v*sin(heading);
    y[n] = u*sin(heading) + v*cos(heading);
    z[n]= w;
  }
  
 
}

void GetCartesianfromLatLon_vec(vector<double> Lats, vector<double> Lons, double heading, double lat, double lon, vector<double> &x,vector< double> &y, vector<double> &z,double R){
  //cout<<"lat is "<<lat<<"\n";
  //cout<<"90-lat is "<<90-lat<<"\n";
  lat = 90-lat;
  //lon = lon-180;
  lat = lat*PI/180.;
  lon = lon*PI/180.;
  //cout<<"heading is "<<heading<<"\n";
  heading =heading *PI/180.;
  R = 6378.1+R; //km
  
  double u;
  double uprime;
  double v;
  double vprime;
  double w;
  double wprime;
  
  double lat_holder;
  double lon_holder;
  int sizer = Lats.size();
  


  for(int n=0;n<sizer;n++){
 
    lat_holder = Lats[n];
    lon_holder = Lons[n];

    lat_holder *=PI/180.;
    lon_holder *=PI/180.;
    lon_holder -=lon;
    
    
    wprime = R*sin(lat_holder);
    double horizontal = sqrt(pow(R,2)-pow(wprime,2));
    vprime = sin(lon_holder)*horizontal;
    uprime = cos(lon_holder)*horizontal;
    
    v = vprime;
    u = uprime*sin(lat) + wprime*cos(lat);
    w = -uprime*cos(lat) + wprime*sin(lat);
    
    x[n] = u*cos(heading) - v*sin(heading);
    y[n] = u*sin(heading) + v*cos(heading);
    z[n]= w;
  }
  
 
}

void LatLon2phitheta(double lat,double lon, double &phi, double &theta){
 
  double NLAT = 90;
  double NLON=180;
  double MAXTHETA=180;
  double PI = 3.14159264;
  double RADDEG = 180./PI;
  double DEGRAD = PI/180.;
  // theta = ((lat+0.5)/(NLAT*MAXTHETA))*RADDEG;
   theta=lat*DEGRAD; 
 
   
  // phi = ((-1*(lon+0.5)+NLON)*2*PI/NLON)-PI/2;
   /* phi = (lon+0.5)*DEGRAD;
  // phi = phi;
  phi-=(2*PI);
  phi= -phi;
   */
   if(lon>=360){
     lon = lon-360;
   }
   phi=lon*DEGRAD;
  
  /* if(lat>20 && lat <21 &&((lon>85 && lon<86)||(lon>94 &&  lon<95)))
      cout<<"lat and lon are "<<lat<<","<<lon<<" and theta phi are "<<theta<<","<<phi<<"\n";
  */
}//latlon2rtheta

void PhiTheta2LatLon(double phi, double theta, double &lat, double &lon){
  double NLAT = 90;
  double NLON=180;
  double MAXTHETA=180;
  double PI = 3.14159264;
  double RADDEG = 180./PI;
  double DEGRAD = PI/180.;
  lat = theta*RADDEG;

  lon = phi*RADDEG;


}


void SphericaltoCart(double phi, double theta, double &x, double &y){
  double  r= 6378.1;
  y = -r*sin(theta)*cos(phi);
  x = -r*sin(theta)*sin(phi);


 

  /*
  if(theta>=.0726193 && theta<=.0726195){
    if(phi >=1.40991 && phi<=1.40993){
      cout<<"for theta, phi "<<theta<<","<<phi<<" x and y are "<<x<<","<<y<<"\n\n";
    }
    if(phi>=1.69676 && phi<=1.69678){
      cout<<"for theta, phi "<<theta<<","<<phi<<" x and y are "<<x<<","<<y<<"\n\n";
    }
  }
  */
}//spherical to cart
void Cart2Sphere(double& phi, double& theta, double x, double y){
  double  r= 6378.1;
  double PI = 3.14159264;
  double RADDEG = 180./PI;
  double z = -1*sqrt(pow(r,2)-pow(x,2)-pow(y,2));
  phi = atan2(x,y);
  phi = phi+PI;

  theta=acos(z/r);
  
  
  
}//cart to sphere
void initializeBaselist(double latBases[1000], double lonBases[1000]){  
  cout<<"initializing base list..."<<endl;
  std::string name;
  double input_lat;
  double input_lon;
  char lat_sign; //A character: this one should always be 'S'.
  char lon_sign; //A character: either 'E' or 'W'
  int misc_field;
  int ctr=0;
  std::string filename="/home/btdailey/analysis/baseLocations/all_base_locations_new.txt";
  std::ifstream base_file(filename.c_str());
  if (base_file.fail())
    {
      std::cerr<<"Error!  Could not open "
	       <<("all_base_locations.txt")<<" to read in base locations!\n";
    } //end if
  while (base_file >> name >> input_lat >> lat_sign >> input_lon >> lon_sign >> misc_field )
    {
      if (lat_sign == 'S')
	input_lat *= -1.;
      if (lon_sign == 'W')
	input_lon = (input_lon*-1.);//+360;
      latBases[ctr]=90-input_lat;
      lonBases[ctr]=input_lon+180;
      ctr++;
      //cout<<"basename: "<<name<<", ctr: "<<ctr<<endl;
    }
  
  //now do pseudo bases
  std::string filenamePseudo="/home/btdailey/analysis/baseLocations/pseudoBases.txt";
  std::ifstream base_filePseudo(filenamePseudo.c_str());
  if (base_filePseudo.fail())
    {
      std::cerr<<"Error!  Could not open "
	       <<("pseudoBases.txt")<<" to read in base locations!\n";
    } //end if
  while (base_filePseudo >> name >> input_lat >> lat_sign >> input_lon >> lon_sign >> misc_field )
    {
      
      if (lat_sign == 'S')
	input_lat *= -1.;
      if (lon_sign == 'W')
	input_lon = (input_lon*-1.);//+360;
      latBases[ctr]=input_lat+90;
      lonBases[ctr]=input_lon+180;
      ctr++;
    }
 
  cout<<"done initializing base list"<<endl;
}//base list


void GetTimePerBin(IceModel *antarctica){
  Balloon *bn1 = new Balloon();
  
  bn1->InitializeBalloon();
  double t_0;
  double t_1;
  double total_time;
  double distance;
  double last_distance=680000;//m
  double min_distance=0;
  
  double horizon = 674.2; //horion ~675km 
  double anita_lon;
  double anita_lat;
  unsigned int realtime;
  
  ofstream MyNewFile;
  MyNewFile.open("timeperbin.txt");
  MyNewFile<<"Lat \t Lon \t Total_time \t minimum distance\n";


  ofstream MyNewFile1;
  MyNewFile1.open("AnitaPosandTime.txt");
  MyNewFile1<<"Lat \t Lon \t Realtime \n";

  double max_lat = 35;
  double max_lon = 360;
  double anita_lat_max=-90;
  double anita_lat_min=90;
  double anita_lon_max = -370;
  double anita_lon_min = 370;
  //cout<<"number of entries in chain is "<<bn1->flightdatachain->GetEntries()<<"\n";
  int nentries = bn1->flightdatachain->GetEntries();
  nentries=22483;//data seems to have more then 1 flight path. cut it short
  cout<<"nentries is "<<nentries<<"\n";
  for(double lat=0;lat<max_lat;lat++){
    cout<<"lat is "<<lat<<". "<<(lat/max_lat)*100<<"% done\n";    
    	
    for(double lon=0;lon<max_lon;lon++){
      //cout<<"lon is "<<lon<<"\n";
      total_time=0;
      min_distance=1E9;
      last_distance=horizon+1;
      bn1->PickBalloonPosition(antarctica,0);
      for(int i=0;i <nentries;i++){

        bn1->PickBalloonPosition(antarctica,i);//go around the flight path
	
	anita_lat=bn1->latitude;//get lat
	anita_lon=bn1->longitude;//get lon
	realtime=bn1->realTime_flightdata;//get time


	

	//cout<<"realtime is "<<realtime<<"\n";
	
	anita_lat +=90;
	anita_lon +=180;
	//	cout<<"anita lat is "<<anita_lat<<"\n";
	//cout<<"anita lon is "<<anita_lon<<"\n";
	
	 if(anita_lon<anita_lon_min){
	    anita_lon_min=anita_lon;
	  }
	  if(anita_lon>anita_lon_max){
	    anita_lon_max=anita_lon;
	  }


	  if(anita_lat<anita_lat_min){
	    anita_lat_min=anita_lat;
	  }
	  if(anita_lat>anita_lat_max){
	    anita_lat_max=anita_lat;
	  }


	distance = GetDistance(anita_lat, anita_lon, lat, lon);
	//cout<<"distance is "<<distance<<"\n";
	if(distance<min_distance){
	  min_distance=distance;
	}
	if(lon==0 && lat==0){
	  MyNewFile1<<anita_lat<<" \t "<<anita_lon<<" \t "<<realtime<<"\n";
	  
	 

	}

	if(distance <=horizon && last_distance > horizon){
	  t_0= (double) realtime;
	  
	 
	 
	  
	}
	
	if(distance >= horizon && last_distance <horizon){
	  t_1 = (double) realtime;
	  total_time += (t_1 - t_0);
	}
	else if(i==nentries-1){//reached end of flight
	  if(distance <=horizon){//if still in horizon, find final time
	    total_time += realtime-t_0;
	  }
	  else if(min_distance>horizon){//still not in horizon, so no time
	    total_time = 0; 
	  }
	  
	}//last anita flight point

	last_distance=distance;
	

      }//i
      MyNewFile1.close();
      MyNewFile<<lat<<"\t"<<lon<<"\t"<<total_time<<"\t"<<min_distance<<"\n";

    }//lon
  }//lat
  
  
  MyNewFile.close();
  
  cout<<"anita lat goes from "<<anita_lat_min<<" to "<<anita_lat_max<<"\n";
  cout<<"anita lon goes from "<<anita_lon_min<<" to "<<anita_lon_max<<"\n";


}//gettimeperbin


double GetDistance(double Lat1, double Lon1, double Lat2, double Lon2){
  double r;
  double theta1;
  double phi1;
  double theta2;
  double phi2;

  double x1;
  double y1;
  double x2;
  double y2;

  double deltax;
  double deltay;

  LatLon2phitheta(Lat1,Lon1,phi1,theta1);
  SphericaltoCart(phi1,theta1,x1,y1);

  LatLon2phitheta(Lat2,Lon2,phi2,theta2);
  SphericaltoCart(phi2,theta2,x2,y2);

  cout<<"x1,y1 is "<<x1<<" "<<y1<<"\n";
  cout<<"x2,y2 is "<<x2<<" "<<y2<<"\n";
  deltax= x2-x1;
  deltay= y2-y1;

  deltax = pow(deltax,2);
  deltay = pow(deltay,2);
  
  r = sqrt(deltax+deltay);

  return r;

}


int GetBrianBin(double lat, double lon){
  int binnumber;

  binnumber = (360*lat)+lon;
  return binnumber;
}

void BintoLatLon(int binnumber, double& lat, double& lon){
  lat = (double)(binnumber%360);
  
  lon = binnumber-(lat*360);
 

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TStyle* RootStyle() {

  // const char* modified = "Borrowed and adapted from paus et al";

  TStyle *RootStyle = new TStyle("Root-Style","The Perfect Style for Plots ;-)");

#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                          // save the global style reference

  gStyle = RootStyle;
#endif
  // otherwise you need to call TROOT::SetStyle("Root-Style")

  // Paper size

  RootStyle->SetPaperSize(TStyle::kUSLetter);

  // Canvas

  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderSize(10);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH      (600);
  RootStyle->SetCanvasDefW      (600);
  RootStyle->SetCanvasDefX      (10);
  RootStyle->SetCanvasDefY      (10);

  // Pads

  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderSize  (10);
  RootStyle->SetPadBorderMode  (0);
  //  RootStyle->SetPadBottomMargin(0.13);
  RootStyle->SetPadBottomMargin(0.16);
  RootStyle->SetPadTopMargin   (0.08);
  RootStyle->SetPadLeftMargin  (0.15);
  RootStyle->SetPadRightMargin (.05);
  RootStyle->SetPadGridX       (0);
  RootStyle->SetPadGridY       (0);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);

  // Frames

  RootStyle->SetFrameFillStyle ( 0);
  RootStyle->SetFrameFillColor ( 0);
  RootStyle->SetFrameLineColor ( 1);
  RootStyle->SetFrameLineStyle ( 0);
  RootStyle->SetFrameLineWidth ( 2);
  RootStyle->SetFrameBorderSize(10);
  RootStyle->SetFrameBorderMode( 0);


  // Histograms

  RootStyle->SetHistFillColor(0);
  RootStyle->SetHistFillStyle(1);
  RootStyle->SetHistLineColor(1);
  RootStyle->SetHistLineStyle(0);
  RootStyle->SetHistLineWidth(2);

  // Functions

  RootStyle->SetFuncColor(1);
  RootStyle->SetFuncStyle(0);
  RootStyle->SetFuncWidth(1);

  //Legends 

  RootStyle->SetStatBorderSize(0);
  RootStyle->SetStatFont      (42);
  RootStyle->SetOptStat       (111111);
  //RootStyle->SetOptStat       ("ne");
  RootStyle->SetStatColor     (0);
  RootStyle->SetStatX         (1);
  RootStyle->SetStatY         (1);
  RootStyle->SetStatFontSize  (0.06);
  RootStyle->SetStatW         (0.2);
  RootStyle->SetStatH         (0.15);
  
  // Labels, Ticks, and Titles

  RootStyle->SetTickLength ( 0.015,"X");
  RootStyle->SetTitleSize  ( 0.055,"X");
  RootStyle->SetTitleOffset( 1.000,"X");
  RootStyle->SetTitleBorderSize(0);
  //  RootStyle->SetTitleFontSize((float)3.);
  RootStyle->SetLabelOffset( 0.015,"X");
  RootStyle->SetLabelSize  ( 0.040,"X");
  RootStyle->SetLabelFont  ( 42   ,"X");

  RootStyle->SetTickLength ( 0.015,"Y");
  RootStyle->SetTitleSize  ( 0.045,"Y");
  RootStyle->SetTitleOffset( 1.65,"Y");
  RootStyle->SetLabelOffset( 0.015,"Y");
  RootStyle->SetLabelSize  ( 0.040,"Y");
  RootStyle->SetLabelSize  ( 0.040,"Z");
  RootStyle->SetLabelFont  ( 42   ,"Y");
  RootStyle->SetLabelFont  ( 42   ,"Z");
  RootStyle->SetTitleOffset(1.00  ,"Z");
  RootStyle->SetTitleSize  ( 0.055,"Z");
  RootStyle->SetTitleFont  (42,"XYZ");
  RootStyle->SetTitleColor  (1);

  // Options

  RootStyle->SetOptFit     (1);

  RootStyle->SetMarkerStyle(20);
  RootStyle->SetMarkerSize(.8);

  //  cout << ">> Style initialized with the Root Style!" << endl;
  //  cout << ">> " << modified << endl << endl;
  return RootStyle;
}
