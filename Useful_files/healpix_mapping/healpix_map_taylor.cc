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

using namespace std;
TStyle* RootStyle();
TStyle *color=RootStyle();
class EarthModel;
class Position;

int k_value = 4;
long n_side =pow(2,k_value);
//Healpix_Base heal_pix;//(n_side,RING);
//pointing point_here;


int n_pix_int=12*n_side*n_side; //Total number of pixels
vector<double> BASE(n_pix_int,0.);
//double PI = 3.14159264;

int NNU;        // number of neutrinos
int WHICHPATH=7;
 
double RANDOMISEPOL=0.;

double fit_10(double *x, double *par){
  double fitval = par[0]*pow(10,par[1]*1+par[2]);
  return fitval;

}
void FractionalArea(vector <vector <double> > phis_theta,vector<int> error_bin, vector<int> &area_pix, vector<double> &areas);
int ErrorPoints(double *x,double *y,double *z,double eta, double theta, double phi,double H,double surface);
void GetLatLonfromCartesian(double *Lats, double *Lons, double heading, double lat, double lon, double *x, double *y, double *z);

void LatLon2phitheta(double lat,double lon, double &phi, double &theta);
void SphericaltoCart(double phi, double theta, double &x, double &y);
void initializeBaselist(double latBases[1000], double lonBases[1000]);
void Cart2Sphere(double& phi, double& theta, double x, double y);
void GetTimePerBin(IceModel *antarctica);
double GetDistance(double Lat1, double Lon1, double Lat2, double Lon2);
int GetBrianBin(double lat, double lon);
void BintoLatLon(int binnumber, double& lat, double& lon);

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
  int r_limit =3000;

  int x_coord=0;
  int y_coord=0;

  int x_limit_min=x_coord-r_limit;
  int x_limit_max=x_coord+r_limit;
  int y_limit_min=y_coord-r_limit;
  int y_limit_max=y_coord+r_limit;
 

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
   TH2F *surfaceweight = new TH2F("surfaceweight",";x(km);y(km)",500,-r_limit,r_limit,500,-r_limit,r_limit);//lat on x, lon on y
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
	    lon_ice =lon_ice;
	    lon_vector.push_back(lon_ice);
	    lat_vector.push_back(lat_ice);
	    icethick_vector.push_back(icethck);

	   
	    surface_height.push_back(surface_above_geoid/1000.);
	    
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
	    surfaceweight->SetBinContent(bin,surface_above_geoid/1000.);
	    //cout<<"icethck, surface_above_geoid is "<<icethck/1000<<" "<<surface_above_geoid/1000.<<"\n";
	}
    }

  


  TCanvas *c1 = new TCanvas("c1","c1",880,800);
  iceweight->Draw("zcol");
  gStyle -> SetOptStat("n");
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

  TCanvas *c1a = new TCanvas("c1a","c1a",880,800);
  surfaceweight->Draw("zcol");
  gStyle -> SetOptStat("n");
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
  TPaletteAxis *paletteaxis1 = (TPaletteAxis*)surfaceweight->GetListOfFunctions()->FindObject("palette");
   paletteaxis1->SetY2NDC(0.38);
   paletteaxis1->SetY1NDC(0.16);
   paletteaxis1->SetX1NDC(0.15);
   paletteaxis1->SetX2NDC(0.25);
  //c1->SetGrayscale();
  
  c1a->Update();

  
  c1a->Update();
  sprintf(printer,"antarcticaheight.png");
  c1a->Print(printer);



  ///////////healpix stuff
  
  //read in root file
  string filter_name;
  int run=13;
  string temp = "/data/anita/btdailey/filter_outputs/pulser_final_geom/%s/output%d_0_0.root";
  char *rootfile;
  filter_name="interpolated";
   rootfile = Form(temp.c_str(),filter_name.c_str(),run);
    
     TChain *PointedTree = new TChain("tdataPointed");
     TChain *TracingTree  = new TChain("tdataTracing");
     TChain *dataTree3  = new TChain("ndata3");
     TChain *dataTree4  = new TChain("ndata4");
     TChain *dataTree2  = new TChain("ndata2");
     TChain *dataTree  = new TChain("ndata");
     PointedTree->Add(rootfile);
     dataTree->Add(rootfile);
     dataTree3->Add(rootfile);
     dataTree2->Add(rootfile);
     dataTree4->Add(rootfile);
     TracingTree->Add(rootfile);
     for(int extra=14;extra<261;extra+=1){
    
       sprintf(rootfile,"/data/anita/btdailey/filter_outputs/pulser_final_geom/%s/output%d_0_0.root",filter_name.c_str(),extra);
       TFile *fpTest = TFile::Open(rootfile);
       if(!fpTest){ 
	 // cout<<"broke at extra "<<extra<<"\n";
	 continue;
	 //break;
       }
       else {
	 delete fpTest;

	 PointedTree->Add(rootfile);
	 dataTree->Add(rootfile);
	 dataTree3->Add(rootfile);
	 dataTree4->Add(rootfile);
	 dataTree2->Add(rootfile);
	 TracingTree->Add(rootfile);
	 }
       }//extra
      
     int nevents0 = PointedTree->GetEntries();
     cout<<"nevents0 is "<<nevents0<<"\n";
     int nevents1 = TracingTree->GetEntries();
     int nevents2 = dataTree->GetEntries();
     cout<<"nevents1 and 2 are "<<nevents1<<"\n";
     double lon;
     double lat;
     double thetaMap;
     double phiMap;
     double anitaLat;
     double anitaLon;
     double anitaAlt;
     double anitaHeading;
     double distance_from_source;


     float ratioFirstToSecondPeak=0;
     //float thetaMap=0;
     float polFractionCoherent=0;
     float hwTriggerFlag=0;
     
     float varnerFlag2;
     float varnerFlag;
     float nadirFlag;
     
     float payloadBlastFlag;
     
     float didIFilter;
     float didIFilterHoriz;
     float didIFilterAboveSatelliteHoriz;
     float didIFilterAboveSatellite;
     float peakHilbertCoherent;
     float peakVal;
     int eventNumber;
     int eventNumber2;
     double sourceLon_tracing;
     double sourceLat_tracing;

     int hilbertctr=0;
     int ratiopeaksctr=0;
     int crosscorrctr=0;
     int polfractionctr=0;
     int rotatedctr=0;
     int elevationctr=0;
     int triggeredctr=0;
     int tracedctr=0;
     int eventTracedFlag2=0;
     int dataeventsctr=0;
     int varnerctr=0;
     int datactr=0;

     dataTree->SetBranchAddress("peakHilbertCoherent",&peakHilbertCoherent);
     dataTree->SetBranchAddress("peakVal",&peakVal);
     dataTree->SetBranchAddress("ratioFirstToSecondPeak",&ratioFirstToSecondPeak);
     dataTree->SetBranchAddress("thetaMap",&thetaMap);
     
     dataTree3->SetBranchAddress("hwTriggerFlag",&hwTriggerFlag);
     dataTree3->SetBranchAddress("varnerFlag2",&varnerFlag2);
     dataTree3->SetBranchAddress("varnerFlag",&varnerFlag);
     dataTree3->SetBranchAddress("nadirFlag",&nadirFlag);
     
     dataTree4->SetBranchAddress("polFractionCoherent",&polFractionCoherent);
     dataTree4->SetBranchAddress("payloadBlastFlag",&payloadBlastFlag);
     dataTree4->SetBranchAddress("didIFilterHoriz",&didIFilterHoriz);
     dataTree4->SetBranchAddress("didIFilterAboveSatelliteHoriz",&didIFilterAboveSatelliteHoriz);
     dataTree4->SetBranchAddress("didIFilterAboveSatellite",&didIFilterAboveSatellite);
      
     dataTree2->SetBranchAddress("didIFilter",&didIFilter);

     TracingTree->SetBranchAddress("eventNumber",&eventNumber);
     TracingTree->SetBranchAddress("sourceLon",&sourceLon_tracing);
     TracingTree->SetBranchAddress("sourceLat",&sourceLat_tracing);
     TracingTree->SetBranchAddress("eventTracedFlag",&eventTracedFlag2);

     PointedTree->SetBranchAddress("eventNumber",&eventNumber2);
     PointedTree->SetBranchAddress("sourceLon",&lon);
     PointedTree->SetBranchAddress("sourceLat",&lat);
     PointedTree->SetBranchAddress("thetaMap",&thetaMap);
     PointedTree->SetBranchAddress("phiMap",&phiMap);
     PointedTree->SetBranchAddress("anitaLatitude",&anitaLat);
     PointedTree->SetBranchAddress("anitaLongitude",&anitaLon);
     PointedTree->SetBranchAddress("anitaAltitude",&anitaAlt);
     PointedTree->SetBranchAddress("heading",&anitaHeading);
     PointedTree->SetBranchAddress("distance_from_source",&distance_from_source);
     long pixel_num;
     long pixel_num_old;
     long pixel_num_event;
     double theta;
     double phi;
     //pix2ang_nest(n_side,pixel_num,theta,phi);
     
     
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
      vector< vector <double> > peakVal_vector (n_pix_int, vector<double>(1));
      vector< vector <double> > peakHilbert_vector (n_pix_int, vector<double>(1));
      vector< vector <double> > bedmap_flagged;
      vector<double> bedmap_holder(3,0);
      vector <int> pixelnum_vector;
      vector <int> area_pix;
      vector <double> areas;
      int n_old=0;
      int doublebinctr=0;
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


      double anitaLat_plot;
      double anitaLon_plot;
      double eventLat_plot;
      double eventLon_plot;
      double Lats_plot[100];
      double Lons_plot[100];

      TH1D *hdistance_from_end = new TH1D("distance_from_end",";Distance (km);Number of Bins",50,0,20);
      
      double LatTD = 90+77.8803;
      double LonTD = 158.459+180;
      //double phi_1,theta_1;
      double x_TD,y_TD;
      int bin_TD;
      vector <int> bin_shadowing;
      
      
      LatLon2phitheta(LatTD,LonTD,phi_1,theta_1);
      pixel_num=ang2pix_ring(n_side,theta_1,phi_1);
      SphericaltoCart(phi_1,theta_1,x_TD,y_TD);

      cout<<"McMurdo is at "<<x_TD<<" "<<y_TD<<"("<<phi_1<<","<<theta_1<<") pix is "<<pixel_num<<"\n";

       min_distance=10000;
       for(int k=0;k<bedmap_x.size();k++){
	
	 distance = sqrt(pow(x_TD -bedmap_x[k],2)+pow(y_TD-bedmap_y[k],2));
	 if(distance < min_distance){
	   bin_TD=k;
	   min_distance=distance;
	 }
	 
	 
       }

       cout<<"TD located in bedmap bin "<<bin_TD<<" with height of "<<surface_height[bin_TD]<<"\n";
      
      for(int m=0;m<nevents0;m++){
	bin_shadowing.clear();
	//m=306026;
	//m=370578;
	//m= 774317;
	//m = 30781;
	//m = 32778;
	if(m < 105301) m = 105301;
	if(m > 105301 && m < 105306) m=105306;
	if(m > 105306 && m < 105322) m=105322;
	if(m > 105322 && m < 105324) m= 105324;
	if(m > 105324 && m < 105326) m= 105326;
	if(m > 105326 && m < 105407) m=105407;
	
	area_pix.clear();
	areas.clear();
	bedmap_flagged.clear();
	PointedTree->GetEvent(m);

	//if(anitaAlt > 35000){
	 //cout<<"lat, lon are "<<lat<<" "<<lon<<" tracing are "<<sourceLat_tracing<<" "<<sourceLon_tracing<<"\n";
	 lat = 90-lat;
	 lon = lon+180;
	 //cout<<"lat, lon are "<<lat<<" "<<lon<<"\n";
	 
	 LatLon2phitheta(lat,lon,phi,theta);
	 SphericaltoCart(phi,theta,x_map,y_map);
	 phi_event = phi;
	 theta_event = theta;
	 //cout<<"phi and theta for event are "<<phi_event<<" "<<theta_event<<"\n";
	 /*
	 phi = PI/4.;
	 theta = 1-pow(1,2.)/(3*pow(n_side,2.))*pow(PI/(2*phi),2.);
	 theta = acos(-1*theta);
	       //cout<<"phi and theta are "<<phi<<" "<<theta<<"\n";
	 */
	 // cout<<"pixel num for event \n";
	 pixel_num=ang2pix_ring(n_side,theta,phi);
	 pixel_num_event = pixel_num;
	 //cout<<"lat, lon, theta, phi "<<lat<<" "<<lon<<" "<<theta<<" "<<phi<<" "<<x_map<<" "<<y_map<<" "<<pixel_num<<"\n";
	 //cout<<"pixel_num  = "<<pixel_num<<"\n";
	 //cout<<"phi and theta are "<<phi<<" "<<theta<<" pixel is "<<pixel_num<<"\n";
	
	 
	 //if(pixel_num==3051){
	
	 anitaLat = 90-anitaLat;
	 anitaLon = anitaLon+180;
	 //cout<<"anitaLat, lon are "<<anitaLat<<" "<<anitaLon<<"\n";
	 
	 LatLon2phitheta(anitaLat,anitaLon,phi,theta);
	 
	 SphericaltoCart(phi,theta,anita_pointx,anita_pointy);
	 //if(anita_pointx > 500) cout<<"other side of TD, m is "<<m<<"\n";
	 //cout<<"anita is at "<<anita_pointx<<" "<<anita_pointy<<"\n";
	   for(int n=n_old;n<nevents1;n++){
	     //cout<<"n is "<<n<<"\n";
	     dataTree->GetEvent(n);
	     TracingTree->GetEvent(n);
	     dataTree3->GetEvent(n);
	     dataTree2->GetEvent(n);
	     dataTree4->GetEvent(n);
	     //cout<<"n is "<<n<<" eventNumber is "<<eventNumber<<" eventNumber2 is "<<eventNumber2<<"\n";
	    
	     if(eventNumber==eventNumber2){
	       n_old = n;
	       break;
	     }
	   }//n
	   
	   float limit=0.9;
	  
	   if (didIFilterHoriz>1 || didIFilter>1){
	     limit=.85;
	   }
	   else if (didIFilterAboveSatellite>0 || didIFilterAboveSatelliteHoriz>0){
	     limit=.85;
	   }
	   
	   if(ratioFirstToSecondPeak <=(1/limit)){
	     ratiopeaksctr++;
	     continue;
	   }
	   
	   if(peakVal<.075){
	     crosscorrctr++;
	     continue;
	   }
	   
	   
	   if(peakHilbertCoherent <=15){
	     hilbertctr++;
	    continue;
	   }
	   
	   if(polFractionCoherent<=.3){
	     polfractionctr++;
	     continue;
	   }
	   
	  // if(peakHilbertCoherent<-350*peakVal+57.14){
	  //   rotatedctr++;
	  //   continue;
	  // }
	   
	   if(thetaMap<=-35 || thetaMap>=0){
	     elevationctr++;
	     continue;
	   }
	   
	   if(eventTracedFlag2!=1){
	     tracedctr++;
	     continue;
	   }
	   
	   if(hwTriggerFlag==0){
	     triggeredctr++;
	     continue;
	   }
	   
	   if(varnerFlag==1){// || varnerFlag2==1){
	     varnerctr++;
	     continue;
	   }
	   
	   //cout<<"m is "<<m<<"\n";
	   //BASE[pixel_num]++;
	   datactr++;

	   if(peakVal<min_peakVal) min_peakVal=peakVal;
	   if(peakHilbertCoherent < min_Hilbert) min_Hilbert = peakHilbertCoherent;
	   // pixelnum_vector[pixel_num].push_back(pixel_num);
	   peakVal_vector[pixel_num].push_back(peakVal);
	   //cout<<"distance_from_source is "<<distance_from_source<<"\n";
	   peakHilbert_vector[pixel_num].push_back(peakHilbertCoherent*distance_from_source/1.E6);
	   // }
	 //cout<<"pixel_num is "<<pixel_num<<"\n";
	   //get error ellipse of one point
	   // cout<<"m is "<<m<<"\n";
	  
	  


	  
	   Lat_point = lat;
	   Lon_point = lon;
	   eta = 90+thetaMap;
	   //cout<<"distance from source is "<<distance_from_source<<"\n";
	   //cout<<"Lat_point, Lon_point is "<<Lat_point<<" "<<Lon_point<<"\n";
	   //cout<<"thetaMap,phiMap, anitaLat,anitaLon,anitaAlt,heading are "<<thetaMap<<" "<<phiMap<<" "<<anitaLat<<" "<<anitaLon<<" "<<anitaAlt<<" "<<anitaHeading<<"\n";
	   //cout<<"m is "<<m<<"\n";
	   
	     //BASE[pixel_num_event]++;
	   
	   //if(m==102131){
	   //if(m==125297){
	   // if(m==203693){//3 bins
	   //if(m==7643){

	   int event_bin;
	   double min_distance_event=10000;
	   
	   for(int k=0;k<bedmap_x.size();k++){
	     distance = sqrt(pow(x_map -bedmap_x[k],2)+pow(y_map-bedmap_y[k],2));
	     if(distance < min_distance_event){
	       event_bin=k;
	       min_distance_event=distance;
	     }
	   }
	   
	    
	   double event_bin_height = surface_height[event_bin];
	   //cout<<"event_bin_height is "<<event_bin_height<<"\n";
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

	   errorpoints_return = ErrorPoints(xholder_temp,yholder_temp,zholder_temp,eta,theta_error, phi_error,anitaAlt/1000,event_bin_height);
	     if(errorpoints_return >2){
	     if(errorpoints_return <4) errorellipse_ctr++;
	     if(errorpoints_return <4){
	       //cout<<"eta is "<<eta<<" theta was "<<thetaMap<<"\n";
	     }
	     //cout<<"anitaLat and Lon are "<<anitaLat<<" "<<anitaLon<<"\n";
	     //cout<<"phi_map, anitaHeading is "<<phiMap<<" "<<anitaHeading<<"\n";
	    
	     GetLatLonfromCartesian(Lats,Lons,anitaHeading-phiMap, anitaLat, anitaLon, xholder_temp, yholder_temp, zholder_temp);//phiMap-anitaHeading

	     
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
	     //cout<<"m is "<<m<<"\n";
	     //cout<<"anita x,y is "<<anita_pointx<<" "<<anita_pointy<<"\n";
	     //cout<<"x,y location are "<<x_map<<" "<<y_map<<"\n";
	     for(int k=0;k<100;k++){
	       //cout<<"x_holder, y_holder are "<<xholder[k]<<" "<<yholder[k]<<"\n";
	       if(xholder[k]<min_x_error) min_x_error = xholder[k];
	       if(xholder[k] > max_x_error) max_x_error = xholder[k];

	       if(yholder[k]<  min_y_error) min_y_error = yholder[k];
	       if(yholder[k] > max_y_error) max_y_error = yholder[k];
	       /*if(m==306026){
		 cout<<"x_holder, y_holder are "<<xholder[k]<<" "<<yholder[k]<<"\n";

		 }*/
	       //cout<<"Lats and Lons of ellipse are  "<<Lats[k]<<" "<<Lons[k]<<"\n";
	       //cout<<"x and y of ellipse are "<<xholder[k]<<" "<<yholder[k]<<"\n";
	       // distance = GetDistance(anitaLat,anitaLon,Lats[k],Lons[k]);
	       distance = sqrt(pow(xholder[k]-anita_pointx,2)+pow(yholder[k]-anita_pointy,2));
	       
	       if(distance > max_distance && errorpoints_return>2) {
		
		 max_distance = distance;
		 Lat_furthest = Lats[k];//furthest point from balloon
		 Lon_furthest = Lons[k];
		 x_furthest = xholder[k];
		 y_furthest = yholder[k];
		 // cout<<"Lats, lons is "<<Lats[k]<<" "<<Lons[k]<<"\n";
		 // LatLon2phitheta(90-Lats[k],180+Lons[k],phi_1,theta_1);
		  LatLon2phitheta(90-Lats[k],180+Lons[k],phi_1,theta_1);
		  //cout<<"doing error ellipse pixel \n";
		  //cout<<"Lats, lons is "<<Lats[k]<<" "<<Lons[k]<<"\n";
		  //cout<<"theta_1,phi_1 are "<<theta_1<<" "<<phi_1<<"\n";
		 pixel_far_ellipse=ang2pix_ring(n_side,theta_1,phi_1);
		 // cout<<"theta_1, phi_1 are "<<theta_1<<" "<<phi_1<<"\n";
		 //cout<<"pixel_far_ellipse is "<<pixel_far_ellipse<<"\n";
		 //Cart2Sphere(phi_1, theta_1, xholder[k], yholder[k]);
		 //pixel_far_ellipse = ang2pix_ring(n_side,theta_1,phi_1);
		 //cout<<"theta_1, phi_1 are "<<theta_1<<" "<<phi_1<<"\n";
		 //cout<<"pixel_far_ellipse is "<<pixel_far_ellipse<<"\n";
     
		 //cout<<"distance is "<<distance<<" lat, lon, x,y are "<<Lat_furthest<<" "<<Lon_furthest<<" "<<x_furthest<<" "<<y_furthest<<"\n";
	       }
	     }//k==0
	     //if(anita_pointx > 500 && max_distance>100) cout<<"other side of TD, m is "<<m<<" max_distance is "<<max_distance<<"\n";
	     // cout<<"anitaLat,Lon,alt phi, theta, x,y are "<<anitaLat<<" "<<anitaLon<<" "<<anitaAlt<<" "<<phi<<" "<<theta<<" "<<anita_pointx<<" "<<anita_pointy<<"\n";
	     
	     double sin_rot = sin(phi_error);
	     double cos_rot = cos(phi_error);

	     double x1;
	     double x2;
	     double y1;
	     double y2;
	     

	     
	     x1 = cos_rot*(x_furthest-anita_pointx) - sin_rot*(y_furthest-anita_pointy);//rotate by 2
	     y1 = sin_rot*(x_furthest-anita_pointx) + cos_rot*(y_furthest-anita_pointy);

	     x2 = cos_rot*(x_furthest-anita_pointx) + sin_rot*(y_furthest-anita_pointy);//rotate by -2
	     y2 = -1*sin_rot*(x_furthest-anita_pointx) + cos_rot*(y_furthest-anita_pointy);

	     x1+=anita_pointx;
	     x2+=anita_pointx;

	     y1+=anita_pointy;
	     y2+=anita_pointy;
	      
	     /*cout<<"x_furthest, y_furthest are "<<x_furthest<<" "<<y_furthest<<"\n";
	     cout<<"x1, y1 are "<<x1<<" "<<y1<<"\n";
	     cout<<"x2, y2 are "<<x2<<" "<<y2<<"\n";
	     */
	     if(x_furthest==anita_pointx) cout<<"vertical line for stepping! blargh \n";
	     double slope0 = (y_furthest-anita_pointy)/(x_furthest-anita_pointx);
	     double slope1 = (y1-anita_pointy)/(x1-anita_pointx);
	     double slope2 = (y2-anita_pointy)/(x2-anita_pointx);

	     int lower_flag;
	     int fill_flag=0;
	     double y_int1;
	     double y_int2;
	     //cout<<"x_furthest, y_furthest are "<<x_furthest<<" "<<y_furthest<<"\n";
	     double delta_x = x_furthest - anita_pointx;
	     double delta_y = y_furthest - anita_pointy;
	     
	     distance = sqrt(pow(delta_x,2)+pow(delta_y,2));
	     delta_x = delta_x/distance;
	     delta_y = delta_y/distance;
	     // cout<<"delta_x, delta_y is "<<delta_x<<" "<<delta_y<<"\n";

	     if(y2 < y1) lower_flag =2;
	     else lower_flag =1;
	    
	     // cout<<"lower_flag is "<<lower_flag<<"\n";
	     //cout<<"bedmap_x.size() is "<<bedmap_x.size()<<"\n";

	     // cout<<"distance to furthest point is "<<(x_furthest-anita_pointx)*delta_x + (y_furthest-anita_pointy)*delta_y<<"\n";
	     int bin_furthest;
	     int contains_highest=0;
	     //cout<<"highest point is "<<highest_point<<" occurs in bin "<<highest_point_bin<<"\n";
	     min_distance=10000;
	     for(int k=0;k<bedmap_x.size();k++){
	       bedmap_bin_flag[k]=0;
	       distance = sqrt(pow(x_furthest -bedmap_x[k],2)+pow(y_furthest-bedmap_y[k],2));
	       if(distance < min_distance){
		 bin_furthest=k;
		 min_distance=distance;
	       }

	      
	     }
	    
	     double bedmap_bin_height = surface_height[bin_furthest];
	    
	     double B_x=0.;
	     double B_y=0.;
	     double C_x=0.;
	     double C_y=0.;
	     double holder;
	     for(int k=0;k<bedmap_x.size();k++){
	       fill_flag=0;
	       distance_hyp = sqrt(pow(bedmap_x[k]-anita_pointx,2)+pow(bedmap_y[k]-anita_pointy,2));
	       distance = (bedmap_x[k]-anita_pointx)*delta_x + (bedmap_y[k]-anita_pointy)*delta_y;
	       distance_perp = -1*(bedmap_x[k]-anita_pointx)*delta_y + (bedmap_y[k]-anita_pointy)*delta_x;
	       // cout<<"distance_perp is "<<distance_perp<<"\n";
	       //distance_perp = sqrt(pow(distance_hyp,2)-pow(distance,2));
	       //cout<<"distance_perp is "<<distance_perp<<"\n";
	       /*  B_x= bedmap_x[k]-anita_pointx -(bedmap_x[k]-anita_pointx)*delta_x;
	       B_y= bedmap_y[k]-anita_pointy -(bedmap_y[k]-anita_pointy)*delta_y;
	       C_x = bedmap_x[k]-anita_pointx;
	       C_y = bedmap_y[k]-anita_pointy;

	       holder = B_x*C_x + B_y*C_y;
	       if(holder<0) distance_perp *=-1;*/
	       
	       
	       if(bedmap_x[k] > anita_pointx && bedmap_x[k] < x_furthest){
		  
		 y_int1 = slope1*(bedmap_x[k]-anita_pointx)+anita_pointy;
		 y_int2 = slope2*(bedmap_x[k]-anita_pointx)+anita_pointy;
		 
		 if(lower_flag==1){
		   if(bedmap_y[k] > y_int1 && bedmap_y[k]<y_int2){
		     
		     bedmap_bin_flag[k]=1;
		     bedmap_holder[0]=k;
		     bedmap_holder[1]=distance;
		    
		     bedmap_holder[2]=distance_perp;
		     fill_flag=1;
		     
		   }
		 }
		 if(lower_flag==2){
		   if(bedmap_y[k] < y_int1 && bedmap_y[k]>y_int2){
		    
		     bedmap_bin_flag[k]=1;
		     bedmap_holder[0]=k;
		     bedmap_holder[1]=distance;
		     bedmap_holder[2]=distance_perp;
		     fill_flag=1;
		   }
		 }
	       }
	       else if(bedmap_x[k] < anita_pointx && bedmap_x[k] > x_furthest){
		
		
		 y_int1 = slope1*(bedmap_x[k]-anita_pointx)+anita_pointy;
		 y_int2 = slope2*(bedmap_x[k]-anita_pointx)+anita_pointy;
		
		 if(lower_flag==1){
		   if(bedmap_y[k] > y_int1 && bedmap_y[k]<y_int2){
		     
		     bedmap_bin_flag[k]=1;
		     bedmap_holder[0]=k;
		     bedmap_holder[1]=distance;
		     bedmap_holder[2]=distance_perp;
		     fill_flag=1;
		   }
		 }
		 if(lower_flag==2){
		   if(bedmap_y[k] < y_int1 && bedmap_y[k]>y_int2){
		     
		     bedmap_bin_flag[k]=1;
		     bedmap_holder[0]=k;
		     bedmap_holder[1]=distance;
		     bedmap_holder[2]=distance_perp;
		     fill_flag=1;
		   }
		 }
	       }
	       //cout<<"distance/ perp is "<<distance<<" "<<distance_perp<<"\n";
	       if(fill_flag==1 && distance < max_distance){
		 //if(m==306026){
		 // cout<<"x,y,height are "<<bedmap_x[bedmap_holder[0]]<<" "<<bedmap_y[bedmap_holder[0]]<<" "<<surface_height[bedmap_holder[0]]<<"\n";
		 // }
		 if(bedmap_holder[0]==highest_point_bin){
		   contains_highest=1;
		 }
		 bedmap_flagged.push_back(bedmap_holder);
	       }
	     }
	     //cout<<"bedmap_flagged size is "<<bedmap_flagged.size()<<"\n";
	     sort(bedmap_flagged.begin(),bedmap_flagged.end());
	     //cout<<"distance from anita are ";
	     
	    
	     slope0 = -1*((anitaAlt/1000)-bedmap_bin_height)/(max_distance);
	    
	     /*
	     for(int k=0;k<bedmap_flagged.size();k++){
	       // cout<<bedmap_flagged[k][0]<<" for bin "<<bedmap_flagged[k][1]<<"\n"; 
	       y_int1 = slope0*(bedmap_flagged[k][0])+anitaAlt/1000;
	       //cout<<"k is "<<k<<" y_int1 is "<<y_int1<<" distance was "<<bedmap_flagged[k][0]<<"\n";
	       if(y_int1+bedmap_bin_height < surface_height[(int) bedmap_flagged[k][1]]){
		 //cout<<"k is "<<bedmap_flagged[k][1]<<" is to tall \n";
		   flag_ctr++;
		   break;
		   
		 }
	     }
	     */
	     
	    

	     /*int flag_ctr=0;
	     for(int k=0;k<bedmap_x.size();k++){
	       if(bedmap_bin_flag[k]>0){
		 cout<<"surface height is "<<surface_height[k]<<"\n";
		 
		 flag_ctr++;
	       }
	     }
	     
	     cout<<"number of bedmap bins inside region is "<<flag_ctr<<"\n";
	     */
	     ///////////////////
	     for(int k=0;k<100;k++){
	       if(Lats[k]!=Lats[k] || Lons[k] !=Lons[k]) continue;
	       
	       LatLon2phitheta(90-Lats[k],180+Lons[k],phi,theta);//90-Lats,180+Lons
	       //cout<<"error ellipse points \n";
	       pixel_num=ang2pix_ring(n_side,theta,phi);
	       //cout<<"theta, phi are "<<theta<<" "<<phi<<"\n";
	       //Cart2Sphere(phi,theta,xholder[k],yholder[k]);
	       //pixel_num=ang2pix_ring(n_side,theta,phi);
	      
	       //cout<<"pixel_nums are "<<pixel_num<<"\n";
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
	     
	     if(errorpoints_return >2){
	       BASE[pixel_num_event]++;
	       sort(phis_theta.begin(),phis_theta.end());
	       FractionalArea(phis_theta,error_bin,area_pix,areas);
	     }
	     //}//m==#
	     
	     twobinFlag=0;
	     mountainFlag=0;
	if(area_pix[1] >0){
	   doublebinctr++;
	   twobinFlag=1;
	}
	if((areas[0] <0.1 || areas[1]<0.1) && areas[1]>1E-3 && areas[2]<1E-10 && areas[3]<1E-10){
	  
	  small_frac_ctr++;
	  }
	
	/*	if(phi_event < min_phi_error || phi_event >max_phi_error){
		 if(cos(theta_event) < min_theta_error || cos(theta_event) > max_theta_error){
		   //cout<<"mountain! m : "<<m<<" theta,phi "<<theta_event<<" "<<phi_event<<" mins maxes "<<min_theta_error<<" "<<max_theta_error<<" "<<min_phi_error<<" "<<max_phi_error<<"\n";
		   mountainctr++;
		 }

		 }*/

	int event_diff_bin=0;
	for(int i=0;i<areas.size();i++){
	       if(area_pix[i] == pixel_num_event) break;
	       if(i==areas.size()-1){
		 event_diff_bin=1;
	       }
	}
	
	if(x_map < min_x_error || x_map >max_x_error){
	  if( y_map< min_y_error || y_map > max_y_error){
	    mountainFlag=1;
	    /*cout<<"mountain event: "<<m<<"\n";
	    cout<<"error points return is "<<errorpoints_return<<"\n";
	    cout<<"anita is at "<<anita_pointx<<" "<<anita_pointy<<"\n";
	    cout<<"source is at "<<x_map<<" "<<y_map<<"\n";
	    for(int k=0;k<100;k++){
	      cout<<"x y and for error is "<<xholder[k]<<" "<<yholder[k]<<"\n";
	    }
	    */
	    //cout<<"mountain! m : "<<m<<" theta,phi "<<theta_event<<" "<<phi_event<<" mins maxes "<<min_theta_error<<" "<<max_theta_error<<" "<<min_phi_error<<" "<<max_phi_error<<"\n";
	    //cout<<"mountain? m  is "<<m<<"\n";
	    //cout<<"x_map,y_map, x error, y eror are "<<x_map<<" "<<y_map<<" "<<min_x_error<<" "<<max_x_error<<" "<<min_y_error<<" "<<max_y_error<<"\n";
	    mountainctr++;
	    if(mountainFlag==1 && event_diff_bin==1){
	      mountainctr2bin++;
	      cout<<" mountain ctr event diff bin. m is "<<m<<"\n";
	    }
	  }
	  
	}
	/*	cout<<"m is "<<m<<" errorpoints_return is "<<errorpoints_return<<"\n";
		cout<<"anita_pointx,y are "<<anita_pointx<<" "<<anita_pointy<<"\n";
		cout<<"x_map, y_map are "<<x_map<<" "<<y_map<<"\n";
		cout<<"x_furthest,y "<<x_furthest<<" "<<y_furthest<<"\n";
		cout<<"minx_error, max is "<<min_x_error<<" "<<max_x_error<<"\n";
		cout<<"min/max y is "<<min_y_error<<" "<<max_y_error<<"\n";
		cout<<"pixel_num_ellipse is "<<pixel_far_ellipse<<"\n";
		cout<<"pixel_num_event is "<<pixel_num_event<<"\n";*/
	for(int i=0;i<areas.size();i++){
	  
	  if(area_pix[i]>0){
	    // cout<<"area_pix is "<<area_pix[i]<<"\n";
	    //cout<<"filling BASE with area "<<areas[i]<<"\n";
	    //BASE[area_pix[i]]+=areas[i];
	    
	  }
	}

	
	double distance_to_event;
	double slope_from_anita;
	double slope_from_bin;
	int flip=0;
	int pixel_num_bedmap;
	int pixel_num_bedmap_temp;
	int bin_bedmap;
	int bin_bedmap_temp;
	int closer_flag=0;
	
	if((mountainFlag==1 && event_diff_bin==1) || twobinFlag==1){
	//if(contains_highest==1){
	   distance_to_event = (x_map-anita_pointx)*delta_x + (y_map-anita_pointy)*delta_y;
	   //if(m<1 && mountainFlag==2){
	     flip=0;
	     shadowctr2++;
	     /*   
	     for(int k3=0;k3<100;k3++){
	     cout<<"x_holder, y_holder are "<<xholder[k3]<<" "<<yholder[k3]<<"\n";
	     
	   }
	   for(int k3=0;k3<bedmap_flagged.size();k3++){
	     cout<<"x,y,height are "<<bedmap_x[bedmap_flagged[k3][0]]<<" "<<bedmap_y[bedmap_flagged[k3][0]]<<" "<<surface_height[bedmap_flagged[k3][0]]<<"\n";
	   }
	   cout<<"m is "<<m<<"\n";
	     */
	   /* cout<<"anita_pointx,y are "<<anita_pointx<<" "<<anita_pointy<<"\n";
	   cout<<"source x an y are "<<x_map<<" "<<y_map<<"\n";
	   cout<<"AnitaAlt is "<<anitaAlt<<"\n";
	   cout<<"xfurthest, y are "<<x_furthest<<" "<<y_furthest<<"\n";
	   */   
	   for(int k=0;k<bedmap_flagged.size();k++){
	     if(bedmap_flagged[k][1] > distance_to_event && bedmap_flagged[k][1]<max_distance){
	       bin_bedmap = (int) bedmap_flagged[k][0];
	       pixel_num_bedmap = ang2pix_ring(n_side,theta_vector[bin_bedmap],phi_vector[bin_bedmap]);
	       slope_from_anita = surface_height[bin_bedmap]-anitaAlt/1000;
	       slope_from_anita = slope_from_anita/(bedmap_flagged[k][1]);
	       /*  if(m==306026){
		 cout<<"anita_pointx,y are "<<anita_pointx<<" "<<anita_pointy<<"\n";
		 cout<<"cource x an y are "<<x_map<<" "<<y_map<<"\n";
		 cout<<"source Location is "<<lat<<" "<<lon<<"\n";
		 cout<<"bedmap_x[k] is "<<bedmap_x[bin_bedmap]<<" "<<bedmap_y[bin_bedmap]<<"\n";
		 cout<<"pixel_num event is "<<pixel_num_event<<"\n";
		 cout<<"pixel nums are "<<pixel_num_bedmap<<"\n";
		 cout<<"slope from anita is "<<slope_from_anita<<"\n";
		 cout<<"surface height 1 is "<< surface_height[(int) bedmap_flagged[k][0]]<<"\n";
		 }*/
	       for(int k2=0;k2<bedmap_flagged.size();k2++){
		 if(bedmap_flagged[k2][1] > 0 && bedmap_flagged[k2][1]<bedmap_flagged[k][1]-2.5){
		   
		   closer_flag=0;
		     if(bedmap_flagged[k][2]>0 && (bedmap_flagged[k2][2] < bedmap_flagged[k][2]-2.5) && bedmap_flagged[k2][2]>0){
		     closer_flag=1;
		   }
		   else if(bedmap_flagged[k][2]<0 && (bedmap_flagged[k2][2] > bedmap_flagged[k][2]+2.5) && bedmap_flagged[k2][2]<0){
		     closer_flag=1;
		     }
		   //cout<<"closer flag is "<<closer_flag<<" perp distance original, temp is "<< bedmap_flagged[k][2]<<" "<<bedmap_flagged[k2][2]<<"\n";
		   if(closer_flag==1 ){
		     bin_bedmap_temp = (int) bedmap_flagged[k2][0];
		     pixel_num_bedmap_temp = ang2pix_ring(n_side,theta_vector[bin_bedmap_temp],phi_vector[bin_bedmap_temp]);
		     slope_from_bin = surface_height[(int) bedmap_flagged[k][0]]-surface_height[(int) bedmap_flagged[k2][0]];
		     slope_from_bin = slope_from_bin / (bedmap_flagged[k][1]-bedmap_flagged[k2][1]);
		     /* cout<<"m is "<<m<<"\n";
		     cout<<"anita_pointx,y are "<<anita_pointx<<" "<<anita_pointy<<"\n";
		     cout<<"source x an y are "<<x_map<<" "<<y_map<<"\n";
		     cout<<"slope_from_bin, slope_from_anita is "<<slope_from_bin<<" "<<slope_from_anita<<"\n";
		     cout<<"pixel_num_bedmap is "<<pixel_num_bedmap<<" "<<pixel_num_bedmap_temp<<"\n";
		     cout<<"AnitaAlt is "<<anitaAlt<<"\n";
		     cout<<"surface_height of bin is "<<surface_height[bin_bedmap]<<"\n";
		     cout<<"surface_height of shadower is "<<surface_height[bin_bedmap_temp]<<"\n";
		     cout<<"distance to bin is "<<bedmap_flagged[k][1]<<"\n";
		     cout<<"distance to shadow is "<<bedmap_flagged[k2][1]<<"\n";
		     cout<<"perp distance is "<<bedmap_flagged[k][2]<<" "<<bedmap_flagged[k2][2]<<"\n";
		      
		        cout<<"bedmap_x[k] is "<<bedmap_x[bin_bedmap]<<" "<<bedmap_y[bin_bedmap]<<"\n";
		     cout<<"bedmap_x[k1] is "<<bedmap_x[bin_bedmap_temp]<<" "<<bedmap_y[bin_bedmap_temp]<<"\n";
		     */
		   if(slope_from_bin < slope_from_anita){
		     bin_bedmap_temp = (int) bedmap_flagged[k2][0];
		     pixel_num_bedmap_temp = ang2pix_ring(n_side,theta_vector[bin_bedmap_temp],phi_vector[bin_bedmap_temp]);
		     if(pixel_num_bedmap != pixel_num_bedmap_temp){

		       bin_shadowing.push_back(pixel_num_bedmap_temp);
		       cout<<"\n";
		       cout<<"m is "<<m<<"\n";
		       for(int k3=0;k3<100;k3++){
			 cout<<"x_holder, y_holder are "<<xholder[k3]<<" "<<yholder[k3]<<"\n";

		       }
		       for(int k3=0;k3<bedmap_flagged.size();k3++){
			 cout<<"x,y,height are "<<bedmap_x[bedmap_flagged[k3][0]]<<" "<<bedmap_y[bedmap_flagged[k3][0]]<<" "<<surface_height[bedmap_flagged[k3][0]]<<"\n";
		       }
		       cout<<"x_furthest, y_furthest are "<<x_furthest<<" "<<y_furthest<<"\n";
		       cout<<"bedmap_distance is "<<bedmap_flagged[k][1]<<" "<<bedmap_flagged[k2][1]<<"\n";
		     cout<<"bedmap_perp is "<<bedmap_flagged[k][2]<<" "<<bedmap_flagged[k2][2]<<"\n";
		     cout<<"anita_pointx,y are "<<anita_pointx<<" "<<anita_pointy<<"\n";
		      cout<<"source x an y are "<<x_map<<" "<<y_map<<"\n";
		     cout<<"bedmap_x[k] is "<<bedmap_x[bin_bedmap]<<" "<<bedmap_y[bin_bedmap]<<"\n";
		     cout<<"bedmap_x[k1] is "<<bedmap_x[bin_bedmap_temp]<<" "<<bedmap_y[bin_bedmap_temp]<<"\n";
		     cout<<"ellipse exists in healpix bins: "<<area_pix[0]<<" "<<area_pix[1]<<" "<<area_pix[2]<<" "<<area_pix[3]<<"\n";
		     cout<<"pixel nums are "<<pixel_num_bedmap<<" "<<pixel_num_bedmap_temp<<"\n";
		     cout<<"bedmap_bin is "<<bin_bedmap<<" "<<bin_bedmap_temp<<"\n";
		     cout<<"theta phi bedmap are "<<theta_vector[bin_bedmap]<<" "<<phi_vector[bin_bedmap]<<"\n";
		     cout<<"theta phi bedmap temp are "<<theta_vector[bin_bedmap_temp]<<" "<<phi_vector[bin_bedmap_temp]<<"\n";
		     cout<<"slope from anita is "<<slope_from_anita<<"\n";
		     cout<<"surface height 1 is "<< surface_height[(int) bedmap_flagged[k][1]]<<"\n";
		     cout<<"distance from anita is "<<bedmap_flagged[k][1]<<"\n";
		     cout<<"slope from bin is "<<slope_from_bin<<"\n";
		     cout<<"height is "<<surface_height[(int) bedmap_flagged[k2][1]]<<"\n";
		     cout<<"distance from other bin is "<<bedmap_flagged[k][1]-bedmap_flagged[k2][1]<<"\n";
		       
		     if(flip==0){
		       cout<<"shadow occurs in event m = "<<m<<" mountain flag is "<<mountainFlag<<"\n";
		     }
		       flip=1;
		     }//different pixels?
		   }//slope?
		   }//perp_distance
		 }//between source location and far point on ellipse
	       }//k2
	     
	     }
	     
	   }
	   if(flip==1){
	     shadowctr++;
	   }
	   }//distance>300
	   // cout<<"\n\n";

	       // }//mountainFlag==1
	   
	     //cout<<"event hits "<<pixelnum_vector.size()<<" number of bins, return was "<<errorpoints_return<<" \n";
	     for(int i=0;i<areas.size();i++){
	       if(area_pix[i] == pixel_num_event) break;
	       
	       if(i== areas.size()-1){
		 // cout<<"oh no!  m is "<<m<<" pixel_num is "<<pixel_num_event<<" area_pix are "<<area_pix[0]<<" "<<area_pix[1]<<" "<<area_pix[2]<<" "<<area_pix[3]<<" errorpoint return was "<<errorpoints_return<<"\n";
		 if(m==306026){
		   anitaLat_plot= anitaLat;
		   anitaLon_plot= anitaLon;
		   eventLat_plot= Lat_event; 
		   eventLon_plot= Lon_event;
		   min_x_error_plot = min_x_error;
		   max_x_error_plot = max_x_error;
		   min_y_error_plot = min_y_error;
		   max_y_error_plot = max_y_error;
		   for(int n=0;n<100;n++){
		     //cout<<"x and y is "<<xholder[n]<<" "<<yholder[n]<<"\n";
		     Lats_plot[n]=Lats[n];
		     Lons_plot[n]=Lons[n];
		   }

		 }
		 /* cout<<"m is "<<m<<" errorpoints_return is "<<errorpoints_return<<"\n";
		 cout<<"anita_pointx,y are "<<anita_pointx<<" "<<anita_pointy<<"\n";
		 cout<<"x_map, y_map are "<<x_map<<" "<<y_map<<"\n";
		 cout<<"x_furthest,y "<<x_furthest<<" "<<y_furthest<<"\n";
		 cout<<"minx_error, max is "<<min_x_error<<" "<<max_x_error<<"\n";
		 cout<<"min/max y is "<<min_y_error<<" "<<max_y_error<<"\n";
		 cout<<"pixel_num_ellipse is "<<pixel_far_ellipse<<"\n";
		 cout<<"pixel_num_event is "<<pixel_num_event<<"\n";
		 */
		 for(int n=0;n<100;n++){
		 
		   if(errorpoints_return>2){
		     //cout<<"Lats and Lons of ellipse are  "<<Lats[n]<<" "<<Lons[n]<<" Lat_event, Lon_event "<<90-Lat_event<<" "<<Lon_event-180<<"\n";
		   }
		 }
	       }
	     }
	     
	     /* int two_bin_flagged=0;
	    
	     if(errorpoints_return > 2){
	       for(int k=0;k<bedmap_flagged.size();k++){
		  y_int1 = slope0*(bedmap_flagged[k][0])+anitaAlt/1000;
		 //cout<<"k is "<<k<<" y_int1 is "<<y_int1<<" distance was "<<bedmap_flagged[k][0]<<"\n";
		 if(y_int1 < surface_height[(int) bedmap_flagged[k][1]]){
		  
		   hdistance_from_end->Fill(max_distance - bedmap_flagged[k][0]);
		   
		   LatLon2phitheta(180-lat_vector[k],lon_vector[k],phi_1,theta_1);
		   //cout<<"lat_vector, lon_vector,phi,theta are "<<180-lat_vector[k]<<" "<<lon_vector[k]<<" "<<phi_1<<" "<<theta_1<<"\n";
		   pixel_bedmap_bin=ang2pix_ring(n_side,theta_1,phi_1);
		   
		    //cout<<"pixel_bedmap_bin is "<<pixel_bedmap_bin<<" "<<pixel_far_ellipse<<"\n";
		    if(pixel_bedmap_bin != pixel_far_ellipse && two_bin_flagged==0){
		      //cout<<"pixel_bedmap_bin is "<<pixel_bedmap_bin<<" "<<pixel_far_ellipse<<"\n";
		      twobinctr++;
		      two_bin_flagged=1;
		    }
		  

		   
		 }
	       }
	     }
	     */
	     // }
	     //}//contains highest
	     }//errorpoints >2
      }//for m
      
	 
      //}
     
     

      cout<<"ratio cut is :"<<ratiopeaksctr<<"\n";
      cout<<"cross corr cut is :"<<crosscorrctr<<"\n";
      cout<<"hilbert cut is :"<<hilbertctr<<"\n";
      cout<<"pol frac cut is :"<<polfractionctr<<"\n";
      cout<<"rotated cut is :"<<rotatedctr<<"\n";
      cout<<"elevation cut is :"<<elevationctr<<"\n";
      cout<<"traced cut is :"<<tracedctr<<"\n";
      cout<<"triggered cut is :"<<triggeredctr<<"\n";
      cout<<"varner cut is :"<<varnerctr<<"\n";
      
      cout<<"numevents passed :"<<datactr<<"\n";
      cout<<"\n errorellipse is "<<errorellipse_ctr<<"\n";
      cout<<"events in more than 1 bin: "<<doublebinctr<<"\n";
      cout<<"events with small frac are "<<small_frac_ctr<<"\n";
      cout<<"events with somethign in the way "<<flag_ctr<<"\n";
      cout<<"events with somethign in the way all in 1 bin: "<<binctr<<" over more than 1 bin "<<twobinctr<<"\n";
      cout<<"events where mountain is in way are "<<mountainctr<<" (of those, "<<mountainctr2bin<<" are in more tha1 bins \n";
      cout<<"events that have shadow are "<<shadowctr<<" out of "<<shadowctr2<<"\n";
      cout<<"min peakVal and hilbert are "<<min_peakVal<<" "<<min_Hilbert<<"\n";
      cout<<"peakVal_vector size is "<<peakVal_vector.size()<<"\n";
      
      double max_base=1.;
      for(int m=0;m<n_pix_int;m++){
	if(BASE[m]>0) cout<<"pix is "<<m<<" num events is "<<BASE[m]<<"\n";
	if(BASE[m] > max_base) max_base = BASE[m];
      }
      
      max_base = 1.1*max_base;

      double LatMC = 90+77.85;
      double LonMC = 166.67+180;
      //double phi_1,theta_1;
      double x_MC,y_MC;
      LatLon2phitheta(LatMC,LonMC,phi_1,theta_1);
      pixel_num=ang2pix_ring(n_side,theta_1,phi_1);
      SphericaltoCart(phi_1,theta_1,x_MC,y_MC);

      cout<<"McMurdo is at "<<x_MC<<" "<<y_MC<<"("<<phi_1<<","<<theta_1<<") pix is "<<pixel_num<<"\n";
      int max_healpixmap=3000;
      int min_healpixmap=-3000;
      int numbin_healpixmap=600;
      TH2D *hhealpix_map = new TH2D("healpix_map",";x(km);y(km);Number of Events",numbin_healpixmap,min_healpixmap,max_healpixmap,numbin_healpixmap,min_healpixmap,max_healpixmap);
      hhealpix_map->SetMarkerSize(20);
      double min_x_map = 350;
      double min_y_map = -700;
      double max_x_map = 650;
      double max_y_map = -400;
      TH2F *hhealpix_error = new TH2F("healpix_map error ",";x(km);y(km);Number of Events",200,500,800,200,300,700);
      hhealpix_error->SetMarkerSize(0.5);
      // for(int x=x_limit_min;x<x_limit_max;x+=1){
      //	for(int y=y_limit_min;y<y_limit_max;y+=1){
      
      /*for(int x=-3000;x<3000;x+=10){//lat
	for(int y=-3000;y<3000;y+=10){//lon
	  x_map=x;
	  y_map=y;
	  Cart2Sphere(phi, theta, x_map, y_map);
	  //cout<<"phi, theta is "<<phi<<" "<<theta<<"\n";
	  pixel_num=ang2pix_ring(n_side,theta,phi);
	  hhealpix_map->Fill(x_map,y_map,BASE[pixel_num]);
	  if(x_map >300 && x_map <350) cout<<"theta,phi /x,y are "<<theta<<" "<<phi<<" "<<x_map<<" "<<y_map<<"\n";
	  //cout<<"x_map, y_map is "<<x_map<<" "<<y_map<<" Base is "<<BASE[pixel_num]<<"\n";
	  //hhealpix_map->Fill(x_map,y_map,1);
	 
	}
	}*/
      

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
      
      //int bin;
      /* for(double Lat=-90;Lat<-60;Lat+=.5){//lat
	for(double Lon=0;Lon<360;Lon+=.5){//lon
	  
	  lat = 90-Lat;
	  lon = Lon+180;
	  // cout<<"lat and lon are "<<lat<<" "<<lon<<"\n";
	   LatLon2phitheta(lat,lon,phi,theta);
	   pixel_num=ang2pix_ring(n_side,theta,phi);
	   SphericaltoCart(phi,theta,x_map,y_map);
	   bin = hhealpix_map->GetBin(x_map,y_map);
	   // cout<<"bin is "<<bin<<"\n";
	   hhealpix_map->Fill(x_map,y_map,BASE[pixel_num]);
	   //hhealpix_map->SetBinContent(bin,BASE[pixel_num]);
	  
	  //cout<<"x_map,y_map,pixel,BASE is "<<x_map<<" "<<y_map<<" "<<pixel_num<<" "<<BASE[pixel_num]<<"\n";
	 
	  }
	  }*/
     
	double map_step_size = max_healpixmap - min_healpixmap;
	map_step_size = map_step_size/numbin_healpixmap;
	cout<<"min_val is "<<min_healpixmap<<" max_Val is "<<max_healpixmap<<" step size is "<<map_step_size<<"\n";
       for(double x_val=min_healpixmap;x_val<max_healpixmap;x_val+=map_step_size){//lat
	 for(double y_val=min_healpixmap;y_val<max_healpixmap;y_val+=map_step_size){//lat
	   
	   Cart2Sphere(phi,theta,x_val,y_val);
	   //cout<<"x_val, y_val are "<<x_val<<" "<<y_val<<" phi, theta are "<<phi<<" "<<theta<<"\n";
	   bin = hhealpix_map->GetBin(x_val,y_val);
	   pixel_num=ang2pix_ring(n_side,theta,phi);
	   SphericaltoCart(phi,theta,x_map,y_map);
	   
	   // cout<<"bin is "<<bin<<"\n";
	   // if(pixel_num==3051 || pixel_num==3036){
	   hhealpix_map->Fill(x_map,y_map,BASE[pixel_num]);
	     // }
	   //hhealpix_map->SetBinContent(bin,BASE[pixel_num]);
	  
	  //cout<<"x_map,y_map,pixel,BASE is "<<x_map<<" "<<y_map<<" "<<pixel_num<<" "<<BASE[pixel_num]<<"\n";
	 
	  }
       }
      
    
       cout<<"passed healpix map stuff \n";
       
       /*  cout<<"x_map,y_map for error are: \n";
       for(int n=0;n<100;n++){
	 lat = 90-Lats_plot[n];
	 lon = 180+Lons_plot[n];
	 //cout<<"lat and lon are "<<lat<<" "<<lon<<" Lats and Lons are "<<Lats[n]<<" "<<Lons[n]<<"\n";
	 LatLon2phitheta(lat,lon,phi,theta);
	 pixel_num=ang2pix_ring(n_side,theta,phi);
	 SphericaltoCart(phi,theta,x_map,y_map);
	 cout<<x_map<<" "<<y_map<<"\n";
	 //cout<<"xmap, y map are "<<x_map<<" "<<y_map<<"\n";
	 
	 hhealpix_error->Fill(x_map,y_map);
	 
	 
       }
       */
      hhealpix_map->SetMinimum(1E-2);
      hhealpix_map->SetMaximum(ceil(max_base));
      gStyle->SetPalette(55);
      //gStyle->SetNdivisions(10,"z");
      gStyle=color;
      //gStyle->SetPalette(55);
      gStyle->SetPadRightMargin (.20);
      gStyle->SetOptTitle(0); //this will disable the title for all coming histograms
      TH2F *haxes = new TH2F("",";x(km);y(km)",200,500,800,200,300,700);
      TCanvas *c2 = new TCanvas("c2","c2",880,800);
      c2->SetLogz();
      // haxes->Draw();
      hhealpix_map->Draw("zcol");
      //hhealpix_error->Draw("same");
      //point->Draw("same");
      //anitapoint->Draw("same");
      
      c2->Print("healpix_map_taylor.png");
      c2->Print("healpix_map_taylor.eps");

       TCanvas *cdistance = new TCanvas("cdistance","cdistance",880,800);
       hdistance_from_end->Draw();
       cdistance->Print("Distance_to_bad_bins.png");


      int bin_number = 3042;
      int numberevents = peakVal_vector[bin_number].size();
      TLine *cutLine = new TLine(0,30,.5,0);
      cutLine->SetLineWidth(2);
      cutLine->SetLineColor(kRed);
      /* TH2F *hRotatedCut = new TH2F("rotatedcut",";Peak Cross Corr Val;Corrected Peak Hilbert Value/1E6 (mV m);",100,0,1,100,0,100);
      for(int n=0;n<numberevents;n++){
	//	cout<<"peakVal, peakHilbert is "<<peakVal_vector[n]<<" "<<peakHilbert_vector[n]<<"\n";
	hRotatedCut->Fill(peakVal_vector[bin_number][n], peakHilbert_vector[bin_number][n]);
      }
      TCanvas *c3 = new TCanvas("c3","c3",800,800);
      hRotatedCut->Draw("zcol");
      cutLine->Draw("same");
      c3->Print("rotatedCut.png");
      */
      double slope=-10.;
      // double y_int;
      int num_cut=0;
      int num_cut_old=0;
      double otherside;
      int fit_start=0;
      int fit_end=0;
      int num_cut_max=0;
      int bad_fit_ctr=0;
      for(int sloper=10;sloper<25;sloper++){
	slope = -1*sloper;
	cout<<" slope is "<<slope<<"\n";
	bad_fit_ctr=0;
	TH1D *hChiSquare = new TH1D("ChiSquare",";ChiSquare/NDF;Number of Events",100,0.,200.);
	TH1D *hpVal = new TH1D("pVal",";p Val;Number of Events",50,0.,1.);
     
	//peakHilbertCoherent<-350*peakVal+57.14
	for(int m=0;m<n_pix_int;m++){
	  //cout<<"m is "<<m<<"\n";
	  fit_start=0;
	  fit_end=0;
	  num_cut_max=0;
	  num_cut=0;
	  num_cut_old=0;
	  //TH1D *hdiffplot = new TH1D("diffplot",";y_intercept;Number of Events Cut;",100,0.,100.);
	  
	  if(BASE[m]>1) {
	    numberevents = BASE[m];
	    // cout<<"Base num and size of vector is "<<BASE[m]<<" "<<peakVal_vector[m].size()<<"\n";
	    TH1D *hdiffplot = new TH1D("diffplot",";y_intercept;Number of Events Cut;",1000,0.,1000.);
	    TH2F *hRotatedCut = new TH2F("rotatedcut",";Peak Cross Corr Val;Corrected Peak Hilbert Value/1E6 (mV m);",100,0,1,100,0,100);
	    
	    for(double y_int=0;y_int<1000;y_int+=1){
	      num_cut=0;
	      for(int n=1;n<=numberevents;n++){
		hRotatedCut->Fill(peakVal_vector[m][n], peakHilbert_vector[m][n]);
		//cout<<"peakVal, hilbert is "<<peakVal_vector[m][n]<<" "<<peakHilbert_vector[m][n]<<"\n";
		//if(m==3042) cout<<"peakHilbert, peakVal is "<<peakHilbert_vector[m][n]<<" "<<peakVal_vector[m][n]<<"\n";
		//if(m==3066) cout<<"peakHilbert, otherside is "<<peakHilbert_vector[m][n]<<" "<<slope*peakVal_vector[m][n]+y_int<<"\n";
		if(peakHilbert_vector[m][n] < slope*peakVal_vector[m][n]+y_int){ 
		  //if(m==3066) cout<<"cut! \n";
		  num_cut++;
		}
	      }
	      
	      //cout<<"bin is "<<m<<" y_int is "<<y_int<<" num_cut is "<<num_cut<<" last_num_cut is "<<num_cut_old<<"\n";
	      
	      num_cut = num_cut-num_cut_old;
	      //cout<<"num_cut is "<<num_cut<<" num_cut_old is "<<num_cut_old<<" y_int is "<<y_int<<"\n";
	      if(num_cut==0){
		//continue;
	      }
	      if(num_cut > num_cut_max){
		num_cut_max=num_cut;
		fit_start = y_int+3;
		//cout<<"m is "<<m<<" num_cut is "<<num_cut<<" "<<fit_start<<"\n";
	      }
	      
	      hdiffplot->Fill(y_int,num_cut);
	      
	      num_cut_old = num_cut+num_cut_old;
	      //cout<<"num_cut_old is "<<num_cut_old<<" numberevents is "<<numberevents<<" y_int is "<<y_int<<"\n";
	      // if( (num_cut_old >0 && num_cut ==0 && num_cut_old > numberevents/2) || num_cut_old==numberevents){
	      if(num_cut_old==numberevents){
		//cout<<"num_cut_old is "<<num_cut_old<<" "<<numberevents<<" "<<y_int<<"\n";
		fit_end = y_int+1;
		break;
	      }
	      
	      
	      
	    }
	    
	    //if(m==3017) cout<<"m is "<<m<<" numberevents is "<<numberevents<<" fit_start and end is "<<fit_start<<" "<<fit_end<<"\n";
	    TF1 *fit;
	    
	    TH2F *haxes_diff;
	    int num_events = hdiffplot->GetEntries();
	    if(fit_end >fit_start+1){
	      // cout<<"here 1 \n";
	      haxes_diff = new TH2F("axes_diff",";y_int;Number of Events Cut",200,0,fit_end+5,200,.9,num_cut_max*2);
	      if(numberevents <100){
		hdiffplot->Fit("expo","QR LL","",fit_start,fit_end);
	      }
	      else{
		hdiffplot->Fit("expo","QR L","",fit_start,fit_end);
	      }
	      
	      fit = hdiffplot->GetFunction("expo");
	      fit->SetLineColor(kRed);
	      fit->SetLineWidth(2);
	      double chisquare = fit->GetChisquare();
	      double NDF = fit->GetNDF();
	      double pval = fit->GetProb();
	      hChiSquare->Fill(chisquare/NDF);
	      hpVal->Fill(pval);
	      if(pval >.95 || pval<0.05){
		bad_fit_ctr++;
		cout<<"m is "<<m<<" sloper is "<<sloper<<" numevents is "<<numberevents<<" pval is "<<pval<<"\n";
	      }
	      double param0 = fit->GetParameter(0);
	      double param1 = fit->GetParameter(1);
	      
	      //cout<<"chisquare is "<<chisquare<<" pVal is "<<pval<<"\n";
	      //delete fit;
	    }
	    
	    /*
	      TCanvas *c3 = new TCanvas("c3","c3",800,800);
	      hRotatedCut->Draw("zcol");
	      cutLine->Draw("same");
	      sprintf(printer,"rotatedCut_%i.png",m);
	      c3->Print(printer);
	    */
	    if(fit_end > fit_start+1 && sloper==22){ 
	      
	      //cout<<"here 2 \n";
	      TCanvas *c4 = new TCanvas("c4","c4",800,800);
	      c4->SetLogy();
	      
	      haxes_diff->Draw();
	      hdiffplot->Draw("sames");
	      c4->Update();
	      //hdiffplot->SetStats(1);
	      //TPaveStats *p2 = (TPaveStats*)fit->FindObject("stats");
	      
	      //sprintf(printer,"Entries:\t %i",num_events);
	      //p2->AddText(printer);
	      //p2->AddText("TEST! \n");
	      //p2->AddText("TEST");
	      //fit->SetStats(0);
	      //c4->Update();
	      sprintf(printer,"diffplot_%i_%i.png",m,sloper);
	      c4->Print(printer);
	      delete c4;
	      delete haxes_diff;
	    }
	    //cout<<"about to delete \n";
	    delete hdiffplot;
	    delete hRotatedCut;
	    
	    
	    //cout<<"deleted \n";
	    //delete fit;
	    //delete c3;
	  }//num_events >0

	
	}//m==0
     cout<<"bad_fit_ctr "<<bad_fit_ctr<<"\n";
     TCanvas *c5 = new TCanvas("c5","c5",800,800);
     //c5->SetLogy();
     hChiSquare->Draw();
     sprintf(printer,"ChiSquare_%i.png",(int)abs(slope));
     c5->Print(printer);

     TCanvas *c6 = new TCanvas("c6","c6",800,800);
     //c6->SetLogy();
     hpVal->Draw();
     sprintf(printer,"pVal_%i.png",(int)abs(slope));
     c6->Print(printer);
     delete c5;
     delete c6;
     delete hpVal;
     delete hChiSquare;
      }//sloper


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
	boundary_vals.push_back(boundary_val);
      }
      if(boundary_val1>bot_val && boundary_val1<top_val){//lies inside ellips
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
  double stepsize = (maxy-miny)/50;
  //cout<<"maxy, miny is "<<maxy<<" "<<miny<<"stepsize is "<<stepsize<<"\n";
  
  //cout<<"z is "<<miny/tan(phi)<<" "<<maxy/tan(phi)<<"\n";
 
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
    y_val = ctr1*stepsize + miny;
    
    b = 4*(R+H)*csc2theta*sin(eta);
    c =2*pow(y_val,2)*csc2theta + 2*csc2theta*(RplusH2-R2)+4*RplusH2*pow(sin(eta),2)-4*RplusH2*pow(cos(eta),2)/pow(tan(theta),2);
    d =4*pow(y_val,2)*csc2phi*(R+H)*sin(eta)+4*(R+H)*sin(eta)*(RplusH2-R2);
    e = pow(y_val,4)*csc2phi+2*pow(y_val,2)*(RplusH2-R2)*csc2phi-4*pow(y_val,2)*RplusH2*pow(cos(eta),2)/pow(tan(phi),2)+pow(RplusH2-R2,2);
    
    
   
    // cout<<"co-efficients are "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" ";

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

      zval0 = sqrt(pow(xval0/tan(theta),2) + pow(y_val/tan(phi),2));
      zval1 = sqrt(pow(xval1/tan(theta),2) + pow(y_val/tan(phi),2));
      // cout<<"point is "<<xval0<<" "<<y_val<<" "<<zval0<<"\n";
      //cout<<"secondpoint is "<<xval1<<" "<<y_val<<" "<<zval1<<"\n";

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
    //cout<<"output is "<<output<<"\n";
    //cout<<"xval0,y,z, are "<<xval0<<" "<<y_val<<" "<<zval0<<"\n";
    //cout<<"xval1,y,z, are "<<xval1<<" "<<y_val<<" "<<zval1<<"\n";

    xprime0 = xval0+(R+H)*sin(eta);
    yprime0 = y_val;
    zprime0 = -zval0 + (R+H)*cos(eta);

    xdoubleprime0 = xprime0*cos(eta)-zprime0*sin(eta);
    ydoubleprime0 = yprime0;
    zdoubleprime0 = xprime0*sin(eta)+zprime0*cos(eta);
    if (output==4){
      xprime1 = xval1+(R+H)*sin(eta);
      yprime1 = y_val;
      zprime1 = -zval1 + (R+H)*cos(eta);
    
      xdoubleprime1 = xprime1*cos(eta)-zprime1*sin(eta);
      ydoubleprime1 = yprime1;
      zdoubleprime1 = xprime1*sin(eta)+zprime1*cos(eta);
    }
    else{
      xdoubleprime1=0;
      ydoubleprime1=0;
      zdoubleprime1=0;
    }


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


void GetLatLonfromCartesian(double *Lats, double *Lons, double heading, double lat, double lon, double *x, double *y, double *z){
  //cout<<"lat is "<<lat<<"\n";
  //cout<<"90-lat is "<<90-lat<<"\n";
  lat = 90-lat;
  //lon = lon-180;
  lat = lat*PI/180.;
  lon = lon*PI/180.;
  heading =heading *PI/180.;
  
  double u;
  double uprime;
  double v;
  double vprime;
  double w;
  double wprime;
  double R;

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
  RootStyle->SetTitleOffset( 1.30,"Y");
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
