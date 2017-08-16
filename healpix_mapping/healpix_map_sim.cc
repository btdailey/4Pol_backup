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
vector<double> BASE_allcuts(n_pix_int,0.);
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
void pix2ang(int nside, int pix, double &theta, double &phi);
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

  double x;
  double y;
  int r_limit =2000;

  int x_coord=0;
  int y_coord=0;

  int x_limit_min=x_coord-r_limit;
  int x_limit_max=x_coord+r_limit;
  int y_limit_min=y_coord-r_limit;
  int y_limit_max=y_coord+r_limit;
 
  TMarker *point_sources[100000];
  TMarker *point_boundaries[100000];
  vector< vector<int> > baseNumbers;
  vector <int> baseNumber_holder;
  vector <int> cluster_number;
  int cluster_info[1000][3];
  vector< vector <double> > phis_theta(100, vector<double>(2));
  vector<int> bedmap_bin;
  vector<int> bedmap_bin_flag(1200000,0);
  int bin;
  int bin1=0;
  double surface_above_geoid;
  vector<double> surface_height;
  //TH2F *iceweight= new TH2F("iceweight",";Lat;Lon",40,0,40,180,0,360);//lat on x, lon on y
   TH2F *iceweight= new TH2F("iceweight",";x(km);y(km)",500,-r_limit,r_limit,500,-r_limit,r_limit);//lat on x, lon on y
  // TH2F *iceweight= new TH2F("iceweight",";x(km);y(km)",50,x_limit_min,x_limit_max,50,y_limit_min,y_limit_max);//lat on x, lon on y
   cout<<"here before antarctica! \n";
  IceModel *antarctica = new IceModel(icemodel, earthmodel, weightabsorption);//instance of icemodel->rename to antarctica (follow icemc?)
  //GetTimePerBin(antarctica);
  cout<<"here \n";
  double surface;
  double geoid;



  for (int i=0;i<antarctica->nRows_ice;i++) {
	for (int j=0;j<antarctica->nCols_ice;j++) {
	  //cout<<"lon, lat are "<<lon_ice<<" "<<lat_ice<<"\n";
	    antarctica->IceENtoLonLat(j,i,lon_ice,lat_ice);
	    icethck=antarctica->IceThickness(lon_ice,lat_ice);
	    surface_above_geoid =antarctica->SurfaceAboveGeoid(lon_ice,lat_ice);
	    surface = antarctica->Surface(lon_ice,lat_ice);
	    lat_ice=180-lat_ice;
	    //lon_ice = 180+lon_ice;
	    lon_vector.push_back(lon_ice);
	    lat_vector.push_back(lat_ice);
	    icethick_vector.push_back(icethck);

	   
	    surface_height.push_back(surface_above_geoid/1000);
	   
	    if(surface_above_geoid/1000. > highest_point){
	      highest_point=surface_above_geoid/1000.;
	      highest_point_bin = bin1;

	    }
	    
	    
	    geoid = antarctica->Geoid(lat_ice);
	   
	    bedmap_bin.push_back(bin1);
	    bin1++;
	    // cout<<"lon_ice is "<<lon_ice<<" and lat is "<<lat_ice<<"\n";

	    if(lon_ice<min_lon){
	      min_lon=lon_ice;
	    }
	    if(lon_ice>max_lon){
	      max_lon=lon_ice;
	    }


	    if(lat_ice<min_lat){
	      min_lat=lat_ice;
	    }
	    if(lat_ice>max_lat){
	      max_lat=lat_ice;
	    }


	    if(icethck<min_ice){
	      min_ice=icethck;
	    }
	    if(icethck>max_ice){
	      max_ice=icethck;
	    }
	    
	    LatLon2phitheta(lat_ice,lon_ice,phi_1,theta_1);
	    phi_vector.push_back(phi_1);
	    theta_vector.push_back(theta_1);
	    if(phi_1<min_phi){
	      min_phi=phi_1;
	    }
	    if(phi_1>max_phi){
	      max_phi=phi_1;
	    }
	    

	    if(theta_1<min_theta){
	      min_theta=theta_1;
	    }
	    if(theta_1>max_theta){
	      max_theta=theta_1;
	    }


	    SphericaltoCart(phi_1,theta_1,x,y);
	    bedmap_x.push_back(x);
	    bedmap_y.push_back(y);
	    
	    if(x<x_min){
	      x_min=x;
	    }
	    if(x>x_max){
	      x_max=x;
	    }
	    

	    if(y<y_min){
	      y_min=y;
	    }
	    if(y>y_max){
	      y_max=y;
	    }


	    //iceweight->Fill(x,y,icethck/1000);
	    bin = iceweight->FindBin(x,y);
	    //cout<<"icethck is "<<icethck<<" max and min are "<<max_ice<<" "<<min_ice<<"\n";
	    
	    iceweight->SetBinContent(bin,icethck/1000.);
	}
    }

  TCanvas *c1 = new TCanvas("c1","c1",880,800);
  iceweight->Draw("zcol");
  //gStyle -> SetOptStat("eou");
  gStyle -> SetOptStat(0000000000);
  gStyle -> SetOptFit(1110);
  //gStyle->SetOptFit(1);
  /*Int_t palette[6];
  palette[0]=17;
  palette[1]=16;
  palette[2]=15;
  palette[3]=14;
  palette[4]=13;
  palette[5]=12;
  // palette[6]=1;
  
 
  
  //gStyle->SetPalette(6,palette);
  gStyle->SetNdivisions(150,"z");
  gPad->Update();*/
  TPaletteAxis *paletteaxis = (TPaletteAxis*)iceweight->GetListOfFunctions()->FindObject("palette");
   paletteaxis->SetY2NDC(0.38);
   paletteaxis->SetY1NDC(0.16);
   paletteaxis->SetX1NDC(0.15);
   paletteaxis->SetX2NDC(0.25);
  //c1->SetGrayscale();
  
  c1->Update();

  char printer[256];
  c1->Update();
  sprintf(printer,"antarctica.png");
  c1->Print(printer);



  ///////////healpix stuff
  
  
  /////////////////////
 
  //read in root file
  string filter_name;
 
  string temp = "/data/anita/btdailey/passingCuts/sim_passedcuts_0301.root";
  char *rootfile;
  rootfile = Form(temp.c_str(),filter_name.c_str());
    
  TChain *pol4_Tree = new TChain("analysis_info_4pol");
  TChain *extraTree = new TChain("weightTree");
       
 
  pol4_Tree->Add(rootfile);
  extraTree->Add(rootfile);
    
     
  int nevents0 = pol4_Tree->GetEntries();
  
  pol4_Tree->SetBranchAddress("pol4_Ptr",&pol4_Ptr);

 


  long pixel_num;
  long pixel_num_old;
  long pixel_num_event;
  double theta;
  double phi;
     //pix2ang_nest(n_side,pixel_num,theta,phi);
     
  double peakHilbertCoherent;
  double peakVal;
  double ratioFirstToSecondPeak;
  double polFractionCoherent;
  
     
  double lon;
  double lat;
  double thetaMap;
  double phiMap;
  double anitaLat;
  double anitaLon;
  double anitaAlt;
  double anitaHeading;
  double distance_from_source;
  double weight;
  extraTree->SetBranchAddress("anitaLat",&anitaLat);
  extraTree->SetBranchAddress("anitaLon",&anitaLon);
  extraTree->SetBranchAddress("anitaAlt",&anitaAlt);
  extraTree->SetBranchAddress("heading",&anitaHeading);
  extraTree->SetBranchAddress("weight",&weight);
  
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
      vector< vector <double> > peakVal_vector (n_pix_int+1, vector<double>(1));
      vector< vector <double> > peakHilbert_vector (n_pix_int+1, vector<double>(1));
      vector< vector <double> > SNR_vector (n_pix_int+1, vector<double>(1));
      vector< vector <double> > weight_vector (n_pix_int+1, vector<double>(1));
      vector< vector <double> > healpix_bin_weights (n_pix_int+1, vector<double>(1));
      vector< vector <double> > bedmap_flagged;

      vector< vector <double> > ratio_vector (n_pix_int+1, vector<double>(1));
      vector< vector <double> > polFrac_vector (n_pix_int+1, vector<double>(1));
      vector< vector <int> > eventNumber_vector (n_pix_int+1, vector<int>(1));


       vector<double> peakVal_vector0;
      vector<double> peakHilbert_vector0; 
      vector<double> SNR_vector0;
      
      vector<double> healpix_bin_weights0; 
      vector<double> bedmap_flagged0;

      vector<double> ratio_vector0; 
      vector<double> polFrac_vector0;
       vector<int> eventNumber_vector0;
       vector<double> weight_vector0;
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
      cout<<"number of events is 10% sample is "<<nevents0<<"\n";

      ///////////////get boundary values//////////////
      double step_size = 2*PI/100.;
      double boundary_val;
      double boundary_val1;
      int offset;
      
      int ratiopeaksctr=0;
  
      int polfractionctr=0;
      
      
      vector<double> Lats_boundary_tmp;
      vector<double> Lons_boundary_tmp;
      vector<double> Lats_boundary1_tmp;
      vector<double> Lons_boundary1_tmp;
      
      vector< vector<double> >Lats_boundary;
      vector< vector<double> >Lons_boundary;
      
      vector< vector<double> >Lats_boundary1;
      vector< vector<double> >Lons_boundary1;
      
      double Lat_temp=0.;
      double Lon_temp=0.;
      int ctr_k;
      int ctr_k1;
      double x_tmp;
      double y_tmp;
      for(int k=1;k<=k_value;k++){
	ctr_k=0;
	ctr_k1=0;
	Lats_boundary_tmp.clear();
	Lons_boundary_tmp.clear();
	
	Lats_boundary1_tmp.clear();
	Lons_boundary1_tmp.clear();
	for(double n = 0;n<2*PI;n+=step_size){
	  offset = (int)(2*n/PI);
	  boundary_val = 1-pow(k,2.)/(3*pow(n_side,2.))*pow(PI/(2*(n-offset*PI/2)),2.);
	  boundary_val1 = 1-pow(k,2.)/(3*pow(n_side,2.))*pow(PI/(2*(n-offset*PI/2)-PI),2.);
	  
	  boundary_val = -1*boundary_val;//dealing with southern hemisphere
	  boundary_val1 = -1*boundary_val1;
	  
	  if(boundary_val <2./3.){
	    
	    boundary_val=acos(boundary_val);
	    
	    SphericaltoCart(n,boundary_val,x_tmp,y_tmp);
	    
	    PhiTheta2LatLon(n,boundary_val,Lat_temp,Lon_temp);
	    Lats_boundary_tmp.push_back(90-Lat_temp);
	    Lons_boundary_tmp.push_back(Lon_temp-180);
	    
	  }
	  
	  
	  if(boundary_val1 <2./3.){
	    boundary_val1 = acos(boundary_val1);
	    SphericaltoCart(n,boundary_val1,x_tmp,y_tmp);
	    
	    
	    PhiTheta2LatLon(n,boundary_val1,Lat_temp,Lon_temp);
	    Lats_boundary1_tmp.push_back(90-Lat_temp);
	    Lons_boundary1_tmp.push_back(Lon_temp-180);
	    
	  }
	  
	  ctr_k++;
	  
	}
	Lats_boundary.push_back(Lats_boundary_tmp);
	Lons_boundary.push_back(Lons_boundary_tmp);
	
	Lats_boundary1.push_back(Lats_boundary1_tmp);
	Lons_boundary1.push_back(Lons_boundary1_tmp);
      }//k==1




      ////////////ROOT OUTPUT
      char filename[150];//simCWdata/largesample/CW2/abby/
      sprintf(filename,"HealPix_sim_partial_0301.root");
      cout<<"outputting to file: "<<filename<<endl;
      TFile *rootfile_out = new TFile(filename,"RECREATE");
  
      TTree *Binned_tree = new TTree("Binned_Tree","variables used for optimization"); 
      Binned_tree->Branch("eventNumber_vector",&eventNumber_vector0);
      Binned_tree->Branch("peakHilbert_vector",&peakHilbert_vector0);
      Binned_tree->Branch("SNR_vector",&SNR_vector0);
      Binned_tree->Branch("peakVal_vector",&peakVal_vector0);
      Binned_tree->Branch("healpix_bin_weights",&healpix_bin_weights0);
      Binned_tree->Branch("ratio_vector",&ratio_vector0);
      Binned_tree->Branch("polFrac_vector",&polFrac_vector0);
      Binned_tree->Branch("weight_vector",&weight_vector0);
      //Binned_tree->Branch("areas_vector",&areas_vector);
      Binned_tree->Branch("BASE",&BASE0);

      //////////////////

      for(int m=0;m<nevents0;m+=1){
	if (m % (nevents0 / 100) == 0){
	  cout << m << " events. " <<(double(m)/double(nevents0)) * 100 << "% complete.\n";
	}
	area_pix.clear();
	areas.clear();
	bedmap_flagged.clear();
	areas_frac=0;

	pol4_Tree->GetEvent(m);
	extraTree->GetEvent(m);
	//cout<<"weight is "<<weight<<"\n";
	//cout<<"eventnumber is "<<pol4_Ptr->eventNumber<<"\n";
	//cout<<"anita is at "<<anitaLat<<" "<<anitaLon<<" "<<anitaAlt<<"\n";
	lat = pol4_Ptr->sourceLat[0];
	lon = pol4_Ptr->sourceLon[0];
	
	

	peakHilbertCoherent = pol4_Ptr->peakHilbertCoherent[0];
	snrCoherent = pol4_Ptr->SNR_coherent[0];
	peakVal = pol4_Ptr->peakVal[0];
	phiMap = pol4_Ptr->phiMap[0];
	thetaMap = pol4_Ptr->thetaMap[0];
	ratioFirstToSecondPeak=pol4_Ptr->ratioFirstToSecondPeak[0];
	thetaMap=pol4_Ptr->thetaMap[0];
      
	polFractionCoherent=pol4_Ptr->polFractionCoherent[0];
       
	float limit =0.85;
	
	if(ratioFirstToSecondPeak <=(1/limit)){
	  ratiopeaksctr++;
	 continue;
       }
       
	if(polFractionCoherent<=.3){
	  polfractionctr++;
	 continue;
       }


	//////////Get Pixel Num for Event//////
	 lat = 90-lat;
	 lon = lon+180;
	
	 LatLon2phitheta(lat,lon,phi,theta);
	 
	 SphericaltoCart(phi,theta,x_map,y_map);
	 phi_event = phi;
	 theta_event = theta;
	
	 pixel_num=ang2pix_ring(n_side,theta,phi);
	 pixel_num_event = pixel_num;
	
	 ///////////Get XY point for ANITA;
	 anitaLat = 90-anitaLat;
	 anitaLon = anitaLon+180;
	 
	 LatLon2phitheta(anitaLat,anitaLon,phi,theta);
	 
	 SphericaltoCart(phi,theta,anita_pointx,anita_pointy);
	


	 if(peakVal<min_peakVal) min_peakVal=peakVal;
	 if(peakHilbertCoherent < min_Hilbert) min_Hilbert = peakHilbertCoherent;
	 
	  
	 Lat_point = lat;
	 Lon_point = lon;
	 eta = 90+thetaMap;
	 

	 ///////////////Get Error Ellipse/////////
	 int event_bin;
	 double min_distance_event=10000;
	 
	 for(int k=0;k<bedmap_x.size();k++){//get Bedmap bin that event lies in to get surface height
	   distance = sqrt(pow(x_map -bedmap_x[k],2)+pow(y_map-bedmap_y[k],2));
	   if(distance < min_distance_event){
	     event_bin=k;
	     min_distance_event=distance;
	   }
	 }
	 
	 
	 double event_bin_height = surface_height[event_bin]; //get that surface height
	 
	 anitaLat_event = anitaLat;
	 anitaLon_event = anitaLon;
	 Lat_event = Lat_point;
	 Lon_event = Lon_point;
	 pixelnum_vector.clear();
	 for(int k=0;k<100;k++){
	   Lats[k]=0.;
	   Lons[k]=0.;
	   xholder_temp[k]=0.;
	   yholder_temp[k]=0.;
	   zholder_temp[k]=0.;
	 }
	 //cout<<"\n\n";
	 
	 errorpoints_return = ErrorPoints(xholder_temp,yholder_temp,zholder_temp,eta,theta_error, phi_error,anitaAlt/1000,event_bin_height);//get error ellipse
	 
	 if(errorpoints_return <4) errorellipse_ctr++;//keep track of number of events that only intersect earth 2 times

	   
	   if(errorpoints_return >2){//use only points that intersect the earth the correct number of times. 
	     
	      //set Lat and Lon for error ellipse
	     
	     GetLatLonfromCartesian(Lats,Lons,anitaHeading - phiMap, anitaLat, anitaLon, xholder_temp, yholder_temp, zholder_temp);
	     
	     GetCartesianfromLatLon(Lats,Lons,anitaHeading-phiMap, anitaLat, anitaLon, xholder_temp, yholder_temp, zholder_temp,event_bin_height);
	      
	     ToConeCoordinates(xholder,yholder,zholder, eta, xholder_temp,yholder_temp,zholder_temp, anitaAlt/1000.,0);

	     
	    
	     vector <double> holder_coordinates;
	     vector <double> holder_x;
	     vector <double> holder_y;
	     vector <double> holder_z;
	    
	     
	     ellipse_coordinates.clear();
	     for(int i=0;i<100;i++){
 
	       holder_coordinates.push_back(xholder[i]);
	       holder_coordinates.push_back(yholder[i]);
	       holder_coordinates.push_back(zholder[i]);
	      
	       ellipse_coordinates.push_back(holder_coordinates);
	       holder_coordinates.clear();
	       
	     }

	    
	       //////////////get boundaries in payload coordinates and then cone coordinates (based on position and heading of payload)/////////////
	       int sizer;// = Lats_boundary[k].size();
	       vector<double> x_holder_tmp;
	       vector<double> y_holder_tmp;
	       vector<double> z_holder_tmp;
	       
	       vector<double> x_boundary;
	       vector<double> y_boundary;
	       vector<double> z_boundary;
		 
	       vector<double> x_boundary1;
	       vector<double> y_boundary1;
	       vector<double> z_boundary1;
	       for(int k=0;k<k_value;k++){
		 x_holder_tmp.clear();
		 y_holder_tmp.clear();
		 z_holder_tmp.clear();
		 
		 x_boundary.clear();
		 y_boundary.clear();
		 z_boundary.clear();
		
		 //cout<<"k is "<<k+1<<"\n";
		 sizer = Lats_boundary[k].size();
		 //cout<<"sizer is "<<sizer<<"\n";
		 for(int i=0;i<sizer;i++){
		   x_holder_tmp.push_back(0.);
		   y_holder_tmp.push_back(0.);
		   z_holder_tmp.push_back(0.);
		   
		   x_boundary.push_back(0.);
		   y_boundary.push_back(0.);
		   z_boundary.push_back(0.);
		   
		  
		 }
		 
		 GetCartesianfromLatLon_vec(Lats_boundary[k],Lons_boundary[k],anitaHeading-phiMap, anitaLat, anitaLon, x_holder_tmp, y_holder_tmp, z_holder_tmp,event_bin_height);
		
		 ToConeCoordinates_vec(x_boundary, y_boundary,z_boundary, eta, x_holder_tmp,y_holder_tmp,z_holder_tmp, anitaAlt/1000.,event_bin_height);

		 
		 for(int i=0;i<x_boundary.size();i++){
		   //cout<<"k is "<<k<<" x,y,z are "<<x_boundary[i]<<" "<<y_boundary[i]<<" "<<z_boundary[i]<<"\n";
		   
		   holder_x.push_back(x_boundary[i]);
		   holder_y.push_back(y_boundary[i]);
		   holder_z.push_back(z_boundary[i]);
		 }
		 boundary_x.push_back(holder_x);
		 boundary_y.push_back(holder_y);
		 boundary_z.push_back(holder_z);
		 holder_x.clear();
		 holder_y.clear();
		 holder_z.clear();
		 
	       }//k==0 < k_value
	       for(int k=0;k<k_value;k++){
		 x_holder_tmp.clear();
		 y_holder_tmp.clear();
		 z_holder_tmp.clear();
		 
		
		 x_boundary1.clear();
		 y_boundary1.clear();
		 z_boundary1.clear();
		 //cout<<"k is "<<k+1<<"\n";
		 sizer = Lats_boundary1[k].size();
		 //cout<<"sizer is "<<sizer<<"\n";
		 for(int i=0;i<sizer;i++){
		   x_holder_tmp.push_back(0.);
		   y_holder_tmp.push_back(0.);
		   z_holder_tmp.push_back(0.);
		   
		  
		   
		   x_boundary1.push_back(0.);
		   y_boundary1.push_back(0.);
		   z_boundary1.push_back(0.);
		 }
		 
		 GetCartesianfromLatLon_vec(Lats_boundary1[k],Lons_boundary1[k],anitaHeading-phiMap, anitaLat, anitaLon, x_holder_tmp, y_holder_tmp, z_holder_tmp,event_bin_height);
		 
		 ToConeCoordinates_vec(x_boundary1, y_boundary1,z_boundary1, eta, x_holder_tmp,y_holder_tmp,z_holder_tmp, anitaAlt/1000.,event_bin_height);

		 for(int i=0;i<x_boundary1.size();i++){
		   //cout<<"k is "<<k<<" x,y,z are "<<x_boundary1[i]<<" "<<y_boundary1[i]<<" "<<z_boundary1[i]<<"\n";
		   
		   holder_x.push_back(x_boundary1[i]);
		   holder_y.push_back(y_boundary1[i]);
		   holder_z.push_back(z_boundary1[i]);
		 }
		 boundary_x_prime.push_back(holder_x);
		 boundary_y_prime.push_back(holder_y);
		 boundary_z_prime.push_back(holder_z);
		 holder_x.clear();
		 holder_y.clear();
		 holder_z.clear();
	     
		 
		
	       }//k==0


		 
	       //////////done with getting boundary values	 
	

	       /////////xholder, yholder are for plotting///////
	     for(int n=0;n<100;n++){
	       if(errorpoints_return>2){
		 //cout<<"Lats and Lons of ellipse are  "<<90-Lats[n]<<" "<<180+Lons[n]<<" Lat_event, Lon_event "<<Lat_event<<" "<<Lon_event<<"\n";
	       }
	       LatLon2phitheta(90-Lats[n],180+Lons[n],phi,theta);
	       
	       SphericaltoCart(phi,theta,xholder[n],yholder[n]);
	       
	     }

	     
	     //cout<<"m is "<<m<<"\n";
	     min_distance = 1000000.;
	     max_distance=0.;
	     min_x_error=10000.;
	     max_x_error=-10000.;
	     min_y_error = 10000.;
	     max_y_error = -10000.;
	   
	     for(int k=0;k<100;k++){
	       
	       if(xholder[k]<min_x_error) min_x_error = xholder[k];
	       if(xholder[k] > max_x_error) max_x_error = xholder[k];

	       if(yholder[k]<  min_y_error) min_y_error = yholder[k];
	       if(yholder[k] > max_y_error) max_y_error = yholder[k];
	      
	       distance = sqrt(pow(xholder[k]-anita_pointx,2)+pow(yholder[k]-anita_pointy,2));
	       
	       if(distance > max_distance && errorpoints_return>2) {
		
		 max_distance = distance;
		 Lat_furthest = Lats[k];//furthest point from balloon
		 Lon_furthest = Lons[k];
		 x_furthest = xholder[k];
		 y_furthest = yholder[k];
		
		  LatLon2phitheta(90-Lats[k],180+Lons[k],phi_1,theta_1);
		 
		 pixel_far_ellipse=ang2pix_ring(n_side,theta_1,phi_1);
		
	       }
	     }//k==0
	     // cout<<"anitaLat,Lon,alt phi, theta, x,y are "<<anitaLat<<" "<<anitaLon<<" "<<anitaAlt<<" "<<phi<<" "<<theta<<" "<<anita_pointx<<" "<<anita_pointy<<"\n";
	     

	     //////////////////////////////////////////////////////////////////////////////////////////

	     /////use error points to find which healpix bins the error ellipse lies in
	     for(int k=0;k<100;k++){
	       if(Lats[k]!=Lats[k] || Lons[k] !=Lons[k]) continue;
	       //if(k==0) cout<<" ellipse 0 point is "<<Lats[k]<<" "<<Lons[k]<<"\n";
	       LatLon2phitheta(90-Lats[k],180+Lons[k],phi,theta);//90-Lats,180+Lons
	      
	       pixel_num=ang2pix_ring(n_side,theta,phi);
	       
	       error_phi[k]=phi;
	       error_theta[k]=cos(theta);
	       error_bin[k]=pixel_num;
	       if(k==0) pixelnum_vector.push_back(pixel_num);
	       for(int l=0;l< pixelnum_vector.size();l++){
		 if(pixel_num == pixelnum_vector[l]){
		     break;
		 }
		 if(l==pixelnum_vector.size()-1){
		   pixelnum_vector.push_back(pixel_num);
		 }
	       }
		  

	       
 
	       if(cos(theta)>max_theta_error) max_theta_error = cos(theta);
	       if(cos(theta)<min_theta_error) min_theta_error = cos(theta);
	       
	       if(phi > max_phi_error) max_phi_error = phi;
	       if(phi < min_phi_error) min_phi_error = phi;

	       
	       

	       //cout<<"phi and theta are "<<phi<<" "<<theta<<"\n";
	       if(k>0 && abs(phi-phis_theta[k-1][0])>6){ 
		 // cout<<"phi is "<<phi<<" last phi was "<<phis_theta[k-1][0]<<" ";

		 if(phi < phis_theta[k-1][0]){
		     phi= phi+2*PI;
		   }
		   else{
		      phi= phi-2*PI;
		   }
		 //cout<<"phi_now is "<<phi<<"\n";
	       }
	       phis_theta[k][0]=phi;
	       phis_theta[k][1]=theta;
	       
	       //cout<<"k is "<<k<<" pixel_num is "<<pixel_num<<" theta, phi are "<<theta<<" "<<phi<<"\n";
	       if(errorpoints_return >2){
		 if(k==0) pixel_num_old = pixel_num;
		 if(pixel_num != pixel_num_old){
		   
		   // cout<<"m is "<<m<<" point is "<<k<<" Lat, Lon, theta, phi are "<<Lats[k]<<" "<<Lons[k]<<" "<<theta<<" "<<phi<<" error point "<<errorpoints_return<<" is located in pixel "<<pixel_num<<" original pixel is "<<pixel_num_old<<"\n";
		 }//pixel_num !=*/
	       }//errorpoints_return
	       
	     }//k


	     /////////Now find amount of area in each bin
	     area_pix.clear();
	     areas.clear();
	     if(errorpoints_return >2){
	       //cout<<"eta, anitaAlt,anitaHeading,phiMao, anitaLat,Lon are "<<eta<<" "<<anitaAlt<<" "<<anitaHeading<<" "<<phiMap<<" "<<anitaLat<<" "<<anitaLon<<"\n";
	       FractionalArea_cone(error_bin,area_pix, areas,eta, anitaAlt/1000.,event_bin_height,anitaHeading - phiMap, anitaLat, anitaLon);
	     }
	     // cout<<"cone area_pix, areas are :";
	     if(area_pix[1] >0){
	       //cout<<"m in different bins! m is "<<m<<"\n";
	       doublebinctr_cone++;
	       areas_frac=1;
	       
	       
	     }
	     if(areas_frac==1){
	        for(int i=0;i<4;i++){
		 
		 if(area_pix[i]>0){
		   //cout<<"area is "<<area_pix[i]<<" "<<areas[i]<<"\n";
		   weight_frac[area_pix[i]]+=areas[i];
		   num_frac[area_pix[i]]++;
		 }
		}
	     }
	    
	     if(errorpoints_return >2){
	       for(int i=0;i<4;i++){
		 
		 if(area_pix[i]>0){
		  
		    /////FOR PLOTTING EVENTS///////
		     anitaLat_plot= anitaLat;
		     anitaLon_plot= anitaLon;
		     eventLat_plot= Lat_event; 
		     eventLon_plot= Lon_event;
		     min_x_error_plot = min_x_error;
		     max_x_error_plot = max_x_error;
		     min_y_error_plot = min_y_error;
		     max_y_error_plot = max_y_error;
		     
		     
		     
		     LatLon2phitheta(eventLat_plot,eventLon_plot,phi,theta);
		     pixel_num=ang2pix_ring(n_side,theta,phi);
		     SphericaltoCart(phi,theta,x_point,y_point);
		    
		     point_sources[pointctr] = new TMarker(x_point,y_point,40);
		     point_sources[pointctr]->SetMarkerStyle(21);
		     point_sources[pointctr]->SetMarkerSize(1);
		     point_sources[pointctr]->SetMarkerColor(kBlack);
		     pointctr++;
		     if(areas_frac ==1){
		       point_boundaries[pointctr1] = new TMarker(x_point,y_point,40);
		       point_boundaries[pointctr1]->SetMarkerStyle(21);
		       point_boundaries[pointctr1]->SetMarkerSize(1);
		       point_boundaries[pointctr1]->SetMarkerColor(kRed);
		       pointctr1++;
		     }
		   
		     //cout<<area_pix[i]<<" "<<areas[i]<<"\n";
		   
		   healpix_bin_weights[area_pix[i]].push_back(areas[i]);
		   numevents[area_pix[i]]=numevents[area_pix[i]]+1;
		   
		   num_binned++;
		  
		   if(whole_flag ==0){
		     BASE[area_pix[i]]+=(weight*areas[i]);
		     peakVal_vector[area_pix[i]].push_back(peakVal);
		     if(peakVal >= 0.075){
		       BASE_allcuts[area_pix[i]]+=areas[i];
		     }
		     peakHilbert_vector[area_pix[i]].push_back(peakHilbertCoherent*distance_from_source/1.E6);
		     SNR_vector[area_pix[i]].push_back(snrCoherent);
		     ratio_vector[area_pix[i]].push_back(ratioFirstToSecondPeak);
		     polFrac_vector[area_pix[i]].push_back(polFractionCoherent);
		     eventNumber_vector[area_pix[i]].push_back(pol4_Ptr->eventNumber);
		     weight_vector[area_pix[i]].push_back(weight);
		   }
		 }
	       }//i==4
	     }//>2
	     


	     //if boundary is in 1 bin//
	     if(errorpoints_return >2){
	       if(whole_flag ==1){
		 BASE[pixel_num_event]+=weight;
		 peakVal_vector[pixel_num_event].push_back(peakVal);
		 peakHilbert_vector[pixel_num_event].push_back(peakHilbertCoherent*distance_from_source/1.E6);
		 healpix_bin_weights[pixel_num_event].push_back(1);
		 SNR_vector[pixel_num_event].push_back(snrCoherent);
		 eventNumber_vector[pixel_num_event].push_back(pol4_Ptr->eventNumber);
		 weight_vector[pixel_num_event].push_back(weight);
	       }
	     }
	   }

	  
      }//for m

      cout<<"ratio of peaks cut is: "<<ratiopeaksctr<<"\n";
  
      cout<<"polfraction cut is: "<< polfractionctr<<"\n";
       for(int m=0;m<n_pix_int;m++){
	healpix_bin_weights0= healpix_bin_weights[m];
	
	BASE0=BASE[m];
	peakVal_vector0=peakVal_vector[m];
	peakHilbert_vector0=peakHilbert_vector[m];
	SNR_vector0=SNR_vector[m];
	ratio_vector0=ratio_vector[m];
	polFrac_vector0=polFrac_vector[m];
	eventNumber_vector0=eventNumber_vector[m];
	weight_vector0 = weight_vector[m];
	Binned_tree->Fill();

      }


       //Binned_tree->Fill();
      cout<<"total events binned is "<<num_binned<<"\n"; 

      rootfile_out->Write();
      rootfile_out->Close();
      //}
    
        ////////GOT ALL INFO FROM EVENTS. NOW CAN USE FOR OPTIMIZATION and PLOTTING//////////
      double max_base=1.;
      for(int m=0;m<n_pix_int;m++){
	if(BASE_allcuts[m]>0) cout<<"pix is "<<m<<" num events is "<<BASE_allcuts[m]<<"\n";
	if(BASE_allcuts[m] > max_base) max_base = BASE_allcuts[m];
      }
      for(int m=0;m<n_pix_int;m++){
	if(num_frac[m]>0) cout<<"pix is "<<m<<" num frac is "<<num_frac[m]<<" weight_frac is "<<weight_frac[m]<<"\n";
	//if(BASE[m] > max_base) max_base = BASE[m];
      }
      
      max_base = 1.1*max_base;
      cout<<"max_base is "<<max_base<<"\n";
      double LatMC = 90+77.85;
      double LonMC = 166.67+180;
      //double phi_1,theta_1;
      double x_MC,y_MC;
      LatLon2phitheta(LatMC,LonMC,phi_1,theta_1);
      pixel_num=ang2pix_ring(n_side,theta_1,phi_1);
      SphericaltoCart(phi_1,theta_1,x_MC,y_MC);

      cout<<"McMurdo is at "<<x_MC<<" "<<y_MC<<"("<<phi_1<<","<<theta_1<<") pix is "<<pixel_num<<"\n";
      int max_healpixmap=2000;
      int min_healpixmap=-2000;
      int numbin_healpixmap=1500;
      TH2D *hhealpix_map = new TH2D("healpix_map",";x(km);y(km);Number of Events",numbin_healpixmap,min_healpixmap,max_healpixmap,numbin_healpixmap,min_healpixmap-300,max_healpixmap-300);
      //TH2D *hhealpix_map = new TH2D("healpix_map",";x(km);y(km);Number of Events",numbin_healpixmap,min_healpixmap,max_healpixmap,numbin_healpixmap,2300,1700);
      hhealpix_map->SetMarkerSize(20);
      double min_x_map = 350;
      double min_y_map = -700;
      double max_x_map = 650;
      double max_y_map = -400;
      TH2F *hhealpix_error = new TH2F("healpix_map error ",";x(km);y(km);Number of Events",200,500,800,200,300,700);
      hhealpix_error->SetMarkerSize(0.5);
     

      LatLon2phitheta(eventLat_plot,eventLon_plot,phi,theta);
      pixel_num=ang2pix_ring(n_side,theta,phi);
      SphericaltoCart(phi,theta,x_point,y_point);
      cout<<"eventLat, lon, phi, theta are "<<eventLat_plot<<" "<<eventLon_plot<<" "<<phi<<" "<<theta<<"\n";
      cout<<"x_point, y_point is "<<x_point<<" "<<y_point<<"\n";
      TMarker *point;
      cout<<"Lat_event, Lon, are "<<Lat_event<<" "<<Lon_event<<"\n";
      cout<<"x_point,y_point is "<<x_point<<" "<<y_point<<"\n";
      point = new TMarker(x_point,y_point,40);
      point->SetMarkerStyle(21);
      point->SetMarkerSize(1);
      point->SetMarkerColor(kRed);



      LatLon2phitheta(anitaLat_plot,anitaLon_plot,phi,theta);
      pixel_num=ang2pix_ring(n_side,theta,phi);
      SphericaltoCart(phi,theta,anita_pointx,anita_pointy);
      
      TMarker *anitapoint;
      cout<<"anitaLat, Lon are "<<anitaLat_event<<" "<<anitaLon_event<<"\n";
      anitapoint = new TMarker(anita_pointx,anita_pointy,40);
      anitapoint->SetMarkerStyle(21);
      anitapoint->SetMarkerColor(kMagenta);
      anitapoint->SetMarkerSize(2);
      cout<<"anita is located at "<<anita_pointx<<" "<<anita_pointy<<"\n";
      
      /////////fill healpix_map
      double map_step_size = max_healpixmap - min_healpixmap;
      map_step_size = map_step_size/numbin_healpixmap;
      cout<<"min_val is "<<min_healpixmap<<" max_Val is "<<max_healpixmap<<" step size is "<<map_step_size<<"\n";
      for(double x_val=min_healpixmap;x_val<max_healpixmap;x_val+=map_step_size){//lat
	for(double y_val=min_healpixmap-300;y_val<max_healpixmap;y_val+=map_step_size){//lat
	  
	  Cart2Sphere(phi,theta,x_val,y_val);
	  //cout<<"x_val, y_val are "<<x_val<<" "<<y_val<<" phi, theta are "<<phi<<" "<<theta<<"\n";
	  bin = hhealpix_map->GetBin(x_val,y_val);
	  pixel_num=ang2pix_ring(n_side,theta,phi);
	  SphericaltoCart(phi,theta,x_map,y_map);
	  
	  // cout<<"bin is "<<bin<<"\n";
	  // if(pixel_num==3008 || pixel_num==3010 || pixel_num==3012 || pixel_num==3028 || pixel_num==3029 || pixel_num==3030 || pixel_num==3031 || pixel_num==3032 || pixel_num==3033 ||pixel_num==3037 || pixel_num==3045 ||pixel_num==3046 || pixel_num==3048 || pixel_num==3051 ||pixel_num==3051|| pixel_num==3053 || pixel_num==3057 || pixel_num==3061 || pixel_num==3062 || pixel_num==3063 || pixel_num==3066 || pixel_num==3069){
	  hhealpix_map->Fill(x_map,y_map,BASE_allcuts[pixel_num]);
	  
	  //hhealpix_map->SetBinContent(bin,BASE[pixel_num]);
	  
	  //cout<<"x_map,y_map,pixel,BASE is "<<x_map<<" "<<y_map<<" "<<pixel_num<<" "<<BASE[pixel_num]<<"\n";
	  //  }
	}
       }
      
    
       cout<<"passed healpix map stuff \n";
       
      
      hhealpix_map->SetMinimum(1E-5);
      hhealpix_map->SetMaximum(ceil(max_base));
      gStyle->SetPalette(55);
      //gStyle->SetNdivisions(10,"z");
      gStyle=color;
      //gStyle->SetPalette(55);
      gStyle->SetPadRightMargin (.20);
      gStyle->SetOptTitle(0); //this will disable the title for all coming histograms
      TH2F *haxes = new TH2F("",";x(km);y(km)",200,-1500,-500,200,-200,700);
      TCanvas *c2 = new TCanvas("c2","c2",880,800);
      c2->SetLogz();
      //haxes->Draw();

      double pixel_theta;
      double pixel_phi;
      double x_pos;
      double y_pos;
      TPaveText *pt;
      hhealpix_map->Draw("zcol");
      

        TFile *flightpath_file = new TFile("anita2flightpath.root"); 
      flightpath_file -> ls(); 
      TGraph *flightpath_graph = (TGraph*)flightpath_file -> Get("Graph");
      flightpath_graph->SetLineColor(kViolet-1);
      flightpath_graph->SetLineWidth(4);    
      flightpath_graph->Draw("L same");  


        for(int m=2988;m<3072;m++){
	  // if(m==3008 || m==3010 || m==3012 || m==3028 || m==3029 || m==3030 || m==3031 || m==3032 || m==3033 ||m==3037 || m==3045 ||m==3046 || m==3048 || m==3051 ||m==3051|| m==3053 || m==3057 || m==3061 || m==3062 || m==3063 || m==3066 || m==3069){
	  if(BASE_allcuts[m]>1E-2){
	    pix2ang(n_side,m,pixel_theta,pixel_phi);
	    //cout<<"pixel_theta, phi are "<<pixel_theta<<" "<<pixel_phi<<"\n ";
	    SphericaltoCart(pixel_phi, pixel_theta, x_pos, y_pos);
	    //cout<<"x_pos, y_pos is "<<x_pos<<" "<<y_pos<<"\n";
	    
	    pt = new TPaveText(x_pos-100,y_pos-100,x_pos+100,y_pos+100);
	    // pt->SetTextColor(kRed);
	    //pt->SetFillColor(kRed);
	    pt->SetFillStyle(4000);
	    sprintf(printer,"%i",m);
	    pt->AddText(printer);
	    pt->Draw("same");
	    //delete pt;
	     }
	  // }//pixel_num
      }

    

      for(int j=0;j<pointctr;j++){
	//point_sources[j]->Draw("same");
      }
      for(int j=0;j<pointctr1;j++){
	//point_boundaries[j]->Draw("same");
      }
    
      c2->Print("healpix_map_sim_0501.png");
      c2->Print("healpix_map_sim_0501.eps");
      

     

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

void pix2ang(int nside, int pix, double &theta, double &phi){


  static const double halfpi=1.570796326794896619231321691639751442099;
  long npix=12*nside*nside;
  double fact2  = 4./npix;
  double z;
    
  int ip = npix - pix;
  int iring = (int)(0.5*(1+sqrt(2*ip-1))); /* counted from South pole */
  int iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

  z = -1.0 + (iring*iring)*fact2;
  phi = (iphi-0.5) * halfpi/iring;
  theta=acos(z);

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
