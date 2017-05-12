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
#include "Math/GSLIntegrator.h"
#include "Math/WrappedTF1.h"

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
//vector<double> *BASE=0;
vector<double> BASE(n_pix_int,0.);;
double BASE2;
vector<double> BASE_sim(n_pix_int,0.);;
double BASE2_sim;
vector<double> weight_frac(n_pix_int,0.);
vector<int> num_frac(n_pix_int,0);
vector<int> numevents(n_pix_int,0.);
//double PI = 3.14159264;

int NNU;        // number of neutrinos
int WHICHPATH=7;
 
double RANDOMISEPOL=0.;


int main(int argc, char **argv) {
  // RootStyle()->cd();
  gStyle=color;
 
  	gStyle->SetPalette(55);
      //gStyle->SetNdivisions(10,"z");
	gStyle=color;
	gStyle->SetOptStat(0000000000);
	gStyle->SetOptFit(0);

  ///////////healpix stuff
  vector<double> *peakVal_vector2=0;
  vector<double> *peakHilbert_vector2=0;
  vector<double> *SNR_vector2=0;
  
  vector<double> *healpix_bin_weights2=0;
  vector<int> *eventNumber_vector2=0;

  vector< vector <double> > peakVal_vector (n_pix_int+1, vector<double>(1));
  vector< vector <double> > peakHilbert_vector (n_pix_int+1, vector<double>(1));
  vector< vector <double> > SNR_vector (n_pix_int+1, vector<double>(1));
  
  vector< vector <double> > healpix_bin_weights (n_pix_int+1, vector<double>(1));
  vector< vector <int> > eventNumber_vector (n_pix_int+1, vector<int>(1));

  //read in root file
  string filter_name;
 
  //string temp = "HealPix_partial_1215.root";
  string temp = "HealPix_partial_0301.root";
  char *rootfile;
  rootfile = Form(temp.c_str(),filter_name.c_str());
    
  TChain *Binned_tree = new TChain("Binned_Tree");
 
  ofstream myfile;
  myfile.open("Cut_values_0301_test2.txt");

  ofstream myfile2;
  myfile2.open("background_values.txt");
  
 
  Binned_tree->Add(rootfile);
    
     
  int nevents0 = Binned_tree->GetEntries();
  Binned_tree->SetBranchAddress("eventNumber_vector",&eventNumber_vector2);
  Binned_tree->SetBranchAddress("peakVal_vector",&peakVal_vector2);
  Binned_tree->SetBranchAddress("peakHilbert_vector",&peakHilbert_vector2);
  Binned_tree->SetBranchAddress("SNR_vector",&SNR_vector2);
  Binned_tree->SetBranchAddress("healpix_bin_weights",&healpix_bin_weights2);
  Binned_tree->SetBranchAddress("BASE",&BASE2);
  cout<<"nevents0 is "<<nevents0<<"\n";

  //////Get Simulated events
  vector<double> *peakVal_vector_sim2=0;
  vector<double> *peakHilbert_vector_sim2=0;
  vector<double> *SNR_vector_sim2=0;
  
  vector<double> *healpix_bin_weights_sim2=0;
  vector<int> *eventNumber_vector_sim2=0;
  vector<double> *weight_vector=0;

  vector< vector <double> > peakVal_vector_sim (n_pix_int+1, vector<double>(1));
  vector< vector <double> > peakHilbert_vector_sim (n_pix_int+1, vector<double>(1));
  vector< vector <double> > SNR_vector_sim (n_pix_int+1, vector<double>(1));
  
  vector< vector <double> > healpix_bin_weights_sim (n_pix_int+1, vector<double>(1));
  vector< vector <int> > eventNumber_vector_sim (n_pix_int+1, vector<int>(1));
  // vector< vector <double> > weight_vector(n_pix_int+1, vector<double>(1));

  //read in root file
  string filter_name_sim;
 
  //string temp_sim = "HealPix_sim_partial_1215.root";
  string temp_sim = "HealPix_sim_partial_0301.root";
  char *rootfile_sim;
  rootfile_sim = Form(temp_sim.c_str(),filter_name_sim.c_str());
    
  TChain *Binned_tree_sim = new TChain("Binned_Tree");
 
       
 
  Binned_tree_sim->Add(rootfile_sim);
    
     
  int nevents_sim = Binned_tree_sim->GetEntries();
  Binned_tree_sim->SetBranchAddress("eventNumber_vector",&eventNumber_vector_sim2);
  Binned_tree_sim->SetBranchAddress("peakVal_vector",&peakVal_vector_sim2);
  Binned_tree_sim->SetBranchAddress("peakHilbert_vector",&peakHilbert_vector_sim2);
  Binned_tree_sim->SetBranchAddress("SNR_vector",&SNR_vector_sim2);
  Binned_tree_sim->SetBranchAddress("healpix_bin_weights",&healpix_bin_weights_sim2);
  Binned_tree_sim->SetBranchAddress("weight_vector",&weight_vector);
  Binned_tree_sim->SetBranchAddress("BASE",&BASE2_sim);
  cout<<"nevents_sim is "<<nevents_sim<<"\n";


  ////////GOT ALL INFO FROM EVENTS. NOW CAN USE FOR OPTIMIZATION and PLOTTING//////////
  double max_base=1.;
  double pval=0.;
 
  
  int hilbertFlag=0;
  int cut_flag=0;
  int cut_flag2=0;

  double num_cut_peak_val;
  int int_cut_rotated=0;
  double num_cut_rotated=0.;
  double num_cut_rotated_sim=0.;
  double num_cut_peak_val_sim;
  for(int m=2980;m<n_pix_int;m++){
  // for(int m=3063m<3064;m++){
    // if( m==3029 || m==3030 || m==3032 || m==3036 || m==3037 || m==3038  || m==3042 || m==3053 || m==3055 || m==3057 || m==3061 || m==3062  || m==3063  || m==3066  || m==3058  || m==3070) {
     Binned_tree->GetEvent(m);
     num_cut_peak_val=0.;
     num_cut_peak_val_sim=0.;
     int_cut_rotated=0;
     num_cut_rotated=0.;
     num_cut_rotated_sim;
     BASE[m]=BASE2;
     
     for(int n=1;n<(*peakVal_vector2).size();n++)
       {
	 peakVal_vector[m].push_back((*peakVal_vector2)[n]);
	 // cout<<"peakVal is "<<(*peakVal_vector2)[n]<<"\n";
	 peakHilbert_vector[m].push_back((*peakHilbert_vector2)[n]);
	 //cout<<"peakHilbert push_back is "<<(*peakHilbert_vector2)[n]<<"\n";
	 SNR_vector[m].push_back((*SNR_vector2)[n]);
	 healpix_bin_weights[m].push_back((*healpix_bin_weights2)[n]);
	 eventNumber_vector[m].push_back((*eventNumber_vector2)[n]);
	 // if(SNR_vector[m][n]>300) cout<<"eventnumber "<<eventNumber_vector[m][n]<<" has SNR = "<<SNR_vector[m][n]<<"\n";
	 // if((*SNR_vector2)[n] >300) cout<<"eventNumber is "<<(*eventNumber_vector2)[n]<<" has SNR ="<<(*SNR_vector2)[n]<<"\n";

	 if(peakVal_vector[m][n] <0.075) num_cut_peak_val+=healpix_bin_weights[m][n];
       }
     // if(BASE[m]>0) cout<<"pix is "<<m<<" num events is "<<BASE[m]<<"\n";
    if(BASE[m] > max_base) max_base = BASE[m];
 

    Binned_tree_sim->GetEvent(m);
    BASE_sim[m]=BASE2_sim;
    cout<<"peakVal_vector_sim2).size() is "<<(*peakVal_vector_sim2).size()<<"\n";
    double weighted_area=0.;
     for(int n=1;n<(*peakVal_vector_sim2).size();n++)
       {
	 
	 peakVal_vector_sim[m].push_back((*peakVal_vector_sim2)[n]);
	 peakHilbert_vector_sim[m].push_back((*peakHilbert_vector_sim2)[n]);
	 SNR_vector_sim[m].push_back((*SNR_vector_sim2)[n]);
	 weighted_area = (*healpix_bin_weights_sim2)[n];
	 weighted_area *=(*weight_vector)[n];
	 //healpix_bin_weights_sim[m].push_back((*healpix_bin_weights_sim2)[n]);
	 healpix_bin_weights_sim[m].push_back(weighted_area);
	 if(peakVal_vector_sim[m][n] <0.075) num_cut_peak_val_sim+=healpix_bin_weights_sim[m][n];

	 //eventNumber_vector_sim[m].push_back((*eventNumber_vector_sim2)[n]);
	 // if(SNR_vector[m][n]>300) cout<<"eventnumber "<<eventNumber_vector[m][n]<<" has SNR = "<<SNR_vector[m][n]<<"\n";
	 // if((*SNR_vector2)[n] >300) cout<<"eventNumber is "<<(*eventNumber_vector2)[n]<<" has SNR ="<<(*SNR_vector2)[n]<<"\n";
       }
     
     cout<<"pix is "<<m<<" num events is "<<BASE[m]<<" "<<BASE_sim[m]<<" "<<BASE[m]-num_cut_peak_val<<" sim is "<<BASE_sim[m]-num_cut_peak_val_sim<<" numcut peak val is "<<num_cut_peak_val<<" "<<num_cut_peak_val_sim<<"\n";


     // }
  }
  
  max_base = 1.1*max_base;
  
      //////OPTIMIZATION! 
      double slope=-38.;
      double max_val=0.;
      double cut_max=0.;
      double cut_Sup=0.;
      double cut_signal=0;
      double max_signal=0;
      double max_Sup=0;
      // double y_int;
      double num_cut=0.;
      double num_cut_old=0.;
      double min_cut=0.;
      double otherside;
      int fit_start=0;
      int fit_end=0;
      double num_cut_max=0;
      int bad_fit_ctr=0;
      int fullcut_flag=0;
      double min_weight=10.;
      double max_weight=0.;

      int multbinFlag=0;
      int numberevents=0;
      char printer[256];
      double max_snr=0.;
      double max_peakval=0.;
      double slope_fit=0.;
      double yint=0.;
      double num_background=0.;
      double num_signal =0.;
      double num_cut_peakVal=0.;
      double y_int_tmp=0.;
      double max_y_int=0.;
      int good_bin=1;
      double step_size=.1;
	//peakHilbertCoherent<-350*peakVal+57.14
      for(int m=2980;m<n_pix_int+1;m++){
      //for(int m=3013;m<3014;m++){
	  // if( m==3029 || m==3030 || m==3032 || m==3036 || m==3037 || m==3038  || m==3042 || m==3053 || m==3055 || m==3057 || m==3061 || m==3062  || m==3063  || m==3066  || m==3058) {
	 
	  //slope40
	//  if(m!=2990 && m!=3008 && m!=3009 && m!=3013 && m!=3014 && m!=3015 && m!=3016 && m!=3017 && m!=3018 && m!=3019 && m!=3027 && m!=3028 && m!=3034 && m!=3035 && m!=3036 && m!=3037 && m!=3038 && m!=3039 && m!=3040 && m!=3041 && m!=3042 && m!=3044 && m!=3046 && m!=3047 && m!=3050 && m!=3052 && m!=3055 && m!=3056 &&  m!= 3059 && m!=3060 && m!=3064 && m!=3065 && m!=3067 && m!=3068 && m!=3071){
	//slope38
	if(m==3008 || m==3010 || m==3012 || m==3029 || m==3030 || m==3031 || m==3032 || m==3033 ||m==3037 ||  m==3045 || m==3048 || m==3051 ||m==3052|| m==3053 || m==3057 || m==3061 || m==3062 || m==3063 || m==3066 || m==3069){  
	
	    // if(m==3009 || m==3014 || m==3015 || m==3027 || m==3029 || m==3030 || m==3032 || m==3033 || m==3036 || m==3037 || m==3038 || m==3039 || m==3040 || m==3041 || m==3042 || m==3044 || m==3045  || m==3053 || m==3055 || m==3057 || m==3061 || m==3062  || m==3063  || m==3066  || m==3058 || m==3069  || m==3070) {
	    /*
	    bad_fit_ctr=0;
	    SNR_vector[n_pix_int].clear();
	    SNR_vector[n_pix_int].push_back(1);
	    
	    peakVal_vector[n_pix_int].clear();
	    peakVal_vector[n_pix_int].push_back(1);
	    
	    healpix_bin_weights[n_pix_int].clear();
	    healpix_bin_weights[n_pix_int].push_back(1);
	    */
	   

	    //if(SNR_vector[n_pix_int].size() <=1) continue;
	  cout<<"DOING OPTIMIZATION ON BIN "<<m<<"\n";
	  
	  //	for(int m=3055;m<3026;m++){
	  cout<<"m is "<<m<<"\n";
	  max_snr=0;
	  max_peakval=0.;
	  max_signal=0;
	  max_Sup=0;
	  Binned_tree->GetEvent(m);
	  fit_start=0;
	  fit_end=0;
	  num_cut_max=0;
	  num_cut=0;
	  num_cut_old=0;
	  fullcut_flag=0;
	  max_y_int=0.;
	  y_int_tmp=0.;
	  //TH1D *hdiffplot = new TH1D("diffplot",";y_intercept;Number of Events Cut;",100,0.,100.);
	  

	  
	  numberevents = peakHilbert_vector[m].size();
	    // cout<<"Base num and size of vector is "<<BASE[m]<<" "<<peakVal_vector[m].size()<<"\n";
	    TH1D *hdiffplot = new TH1D("diffplot",";y_intercept;Number of Weighted Events Cut;",10000,0.,10000.);
	    
	   
	    hdiffplot->Sumw2();
	    TH2F *hRotatedCut;
	    TH2F *hRotatedCut_sim;
	    double max_cutval=500;
	    
	    int hbin_num=(int) (max_cutval/step_size);
	    TH1F *hopt = new TH1F("optimization",";Cut_Value;Signal/S_{up}",hbin_num,0,max_cutval);
	  TH1F *hnum = new TH1F("num",";Cut_Value;Signal",hbin_num,0,max_cutval);
	  TH1F *hdenom = new TH1F("denom",";Cut_Value;S_{up}",hbin_num,0,max_cutval);
	  TH1F *hbackground = new TH1F("background",";Cut_Value;Number of Background Expected",hbin_num,0,500);
	  
	    if(hilbertFlag==0){
	     hRotatedCut = new TH2F("rotatedcut",";Peak Cross Corr Val;SNR of Coherent Waveform;",200,0,1,2000,0,1000);
	     hRotatedCut_sim = new TH2F("rotatedcut_sim",";Peak Cross Corr Val;SNR of Coherent Waveform;",200,0,1,2000,0,1000);
	    }
	    else{
	     hRotatedCut = new TH2F("rotatedcut",";Peak Cross Corr Val;Peak Hilbert of Coherent Waveform;",1000,0,1,5000,0,500);
	    }
	     max_weight=0;
	     //cout<<"m is "<<m<<"\n";
	     if(numberevents<=1){
	       /* for(int n=0;n<numberevents;n++){
		 SNR_vector[n_pix_int].push_back(SNR_vector[m][n]);
		 peakVal_vector[n_pix_int].push_back(peakVal_vector[m][n]);
		 healpix_bin_weights[n_pix_int].push_back(healpix_bin_weights[m][n]);
		 peakHilbert_vector[n_pix_int].push_back(peakHilbert_vector[m][n]);

		 SNR_vector_sim[n_pix_int].push_back(SNR_vector_sim[m][n]);
		 peakVal_vector_sim[n_pix_int].push_back(peakVal_vector_sim[m][n]);
		 healpix_bin_weights_sim[n_pix_int].push_back(healpix_bin_weights_sim[m][n]);
		 peakHilbert_vector_sim[n_pix_int].push_back(peakHilbert_vector_sim[m][n]);

	       }
	       */
	       cout<<"m is "<<m<<" continuing! \n";
	       delete hdiffplot;
	       
	       
	       delete hRotatedCut;
	       delete hopt;
	       continue;
	     }
	     num_cut_peakVal=0.;
	     min_cut=100.;
	     for(int n=1;n<numberevents;n++){
		  if(peakVal_vector[m][n] < 0.075){
		    num_cut_peakVal+=healpix_bin_weights[m][n]; //Apply peak val cut now so these events appear in rotated plot
		 }
		 if(hilbertFlag==0){
		   if(SNR_vector[m][n]>max_snr){
		     //cout<<"SNR_vector is "<<SNR_vector[m][n]<<"\n";
		     //cout<<"max_snr is "<<max_snr<<"\n";
		     max_snr = SNR_vector[m][n];
		     max_peakval = peakVal_vector[m][n];
		     //cout<<"SNR_vector is "<<SNR_vector[m][n]<<"\n";
		     //cout<<"max_snr is "<<max_snr<<"\n";
		     y_int_tmp = SNR_vector[m][n]-slope*peakVal_vector[m][n];
		     if(y_int_tmp > max_y_int) max_y_int = y_int_tmp;
		   }
		 }
		 else{
		   //if(peakHilbert_vector[m][n]>max_snr) max_snr = peakHilbert_vector[m][n];
		 }

		 if(hilbertFlag==0){
		   hRotatedCut->Fill(peakVal_vector[m][n], SNR_vector[m][n],healpix_bin_weights[m][n]);
		   //if(m==3029) cout<<"n is "<<n<<" peakVal is "<<peakVal_vector[m][n]<<" "<<SNR_vector[m][n]<<"\n";
		 }
		 else{
		   hRotatedCut->Fill(peakVal_vector[m][n],peakHilbert_vector[m][n],healpix_bin_weights[m][n]);
		 }
	     }
	     
	     //double max_y_int = max_snr-slope*max_peakval;
	     //cout<<"max_y_int is "<<max_y_int<<" max_pakval, snr are "<<max_peakval<<" "<<max_snr<<"\n";
	     for(int n=1;n<SNR_vector_sim[m].size();n++){
	       

	       hRotatedCut_sim->Fill(peakVal_vector_sim[m][n], SNR_vector_sim[m][n],healpix_bin_weights_sim[m][n]);
	       
	     }
	     
	     for(double y_int=0;y_int<10000;y_int+=1){
	       num_cut=0.;
	       
	       for(int n=1;n<numberevents;n++){
		 cut_flag=0;
		 cut_flag2=0;
		 //cout<<"hilbert is "<<peakHilbert_vector[m][n]<<"\n";
		 if(peakVal_vector[m][n] < 0.075){
		   continue; //Apply peak val cut now so these events appear in rotated plot
		 }
		//cout<<"y_int is "<<y_int<<" n is "<<n<<" peakVal, SNR are "<<peakVal_vector[m][n]<<" "<<SNR_vector[m][n]<<"\n";
		 if(hilbertFlag==0){
		   if(SNR_vector[m][n] < slope*peakVal_vector[m][n]+y_int && SNR_vector[m][n] > slope*peakVal_vector[m][n]+y_int-1) cut_flag=1; 
		 }
		 else{
		   if(peakHilbert_vector[m][n] < slope*peakVal_vector[m][n]+y_int && peakHilbert_vector[m][n] > slope*peakVal_vector[m][n]+y_int-1) cut_flag=1; 
		 }


		 if(cut_flag==1){ 
		   // cout<<"filling! at y_int = "<<y_int<<"\n";
		   hdiffplot->Fill(y_int,healpix_bin_weights[m][n]);
		   if(healpix_bin_weights[m][n]>max_weight) max_weight = healpix_bin_weights[m][n];
		   
		 }
		 if(hilbertFlag==0){
		   if(SNR_vector[m][n] < slope*peakVal_vector[m][n]+y_int){
		     cut_flag2=1;
		   }
		 }
		 else{
		   if(peakHilbert_vector[m][n] < slope*peakVal_vector[m][n]+y_int){
		     cut_flag2=1;
		   }
		 }

		 if(cut_flag2==1){
		   if(whole_flag==1){
		     num_cut+=1;
		   }
		   else{
		     num_cut+=healpix_bin_weights[m][n];
		     if(healpix_bin_weights[m][n] < min_weight) min_weight = healpix_bin_weights[m][n];
		     
		   }
		 }//num_cut
	       }//n=1
	       
	       //cout<<"bin is "<<m<<" y_int is "<<y_int<<" num_cut is "<<num_cut<<" last_num_cut is "<<num_cut_old<<"\n";
	       // if(m==3051) cout<<"numcut is "<<num_cut<<" old is "<<num_cut_old<<"\n";
	       num_cut = num_cut-num_cut_old;
	       //cout<<"num_cut is "<<num_cut<<" num_cut_old is "<<num_cut_old<<" y_int is "<<y_int<<"\n";
	       if(num_cut==0){
		 //continue;
	       }
	       if(num_cut > num_cut_max){
		 num_cut_max=num_cut;
		 fit_start = y_int+1;
		 //cout<<"m is "<<m<<" num_cut is "<<num_cut<<" "<<fit_start<<"\n";
	       }
	       if(num_cut < min_cut && num_cut>0){
		 min_cut=num_cut;
	       }
	       
	       //hdiffplot->Fill(y_int,num_cut);
	       
	       num_cut_old = num_cut+num_cut_old;
	       
	       //cout<<"num_cut_old is "<<num_cut_old<<" numberevents is "<<numberevents<<" y_int is "<<y_int<<"\n";
	       // if( (num_cut_old >0 && num_cut ==0 && num_cut_old > numberevents/2) || num_cut_old==numberevents){
	       
	       if(num_cut_old > BASE[m]-num_cut_peakVal-.0001){
		 //cout<<"num_cut_old, BASE is "<<num_cut_old<<" "<<BASE[m]-num_cut_peakVal<<"\n";
		 //cout<<"num_cut_old is "<<num_cut_old<<" "<<numberevents<<" "<<y_int<<"\n";
		 fit_end = y_int+1;
		 //y_int = 1E5;
		 fullcut_flag=1;
		 break;
	       }
	      
	       
	     }//y_int
	     // cout<<"doing fit \n";
	    if(fit_end<fit_start) fit_end = fit_start;
	    //if(fullcut_flag==0) fit_end = 3*fit_start;
	    //if(m==3017) cout<<"m is "<<m<<" numberevents is "<<numberevents<<" fit_start and end is "<<fit_start<<" "<<fit_end<<"\n";
	    TF1 *fit;
	   
	    double likelihood_val;
	   
	    double edm, errdef;
	    int nvpar,nparx;
	    int num_events = hdiffplot->GetEntries();
	    //cout<<"m is "<<m<<" sloper is "<<sloper<<" fit start and end are "<<fit_start<<" "<<fit_end<<" num_cut_max is "<<num_cut_max<<" fullcut_flag is "<<fullcut_flag<<" num_cut_old is "<<num_cut_old<<" min_weight is "<<min_weight<<"\n";
	    cout<<"fit_start, fit_end is "<<fit_start<<" "<<fit_end<<"\n";
	    /*if(fit_end < fit_start+3){
	      cout<<"doing fit \n";
	       for(int n=0;n<numberevents;n++){
		 SNR_vector[n_pix_int].push_back(SNR_vector[m][n]);
		 peakVal_vector[n_pix_int].push_back(peakVal_vector[m][n]);
		 healpix_bin_weights[n_pix_int].push_back(healpix_bin_weights[m][n]);
		 peakHilbert_vector[n_pix_int].push_back(peakHilbert_vector[m][n]);
	       }
	       cout<<"m is "<<m<<" continuing! \n";
	       delete hdiffplot;
	       
	       
	       delete hRotatedCut;
	       continue;
	    }
	    */
	    if(fit_end <= fit_start+3) continue;


	    if(fit_end >fit_start+3 && numberevents >1){
	      cout<<"doing fit \n";
	      // haxes_diff = new TH2F("axes_diff",";y_int;Number of Events Cut",200,0,fit_end+5,200,.9,num_cut_max*2);
	      // if(numberevents <100){
	      hdiffplot->Fit("expo","QR WL","",fit_start,fit_end);
		//}
		//else{
		//	hdiffplot->Fit("expo","QR L","",fit_start,fit_end);
		//}
	      
	      fit = hdiffplot->GetFunction("expo");
	      fit->SetLineColor(kRed);
	      fit->SetLineWidth(2);
	      double chisquare = fit->GetChisquare();
	      double NDF = fit->GetNDF();
	      
	      
	      slope_fit = fit->GetParameter(1);
	      yint = fit->GetParameter(0);
	      //cout<<"slope_fit and yint is "<<slope_fit<<" "<<yint<<"\n";
	      TVirtualFitter *fitter = TVirtualFitter::Fitter(hdiffplot);
	      fitter->GetStats(likelihood_val, edm, errdef,nvpar,nparx);
	    }//fit_end >f_start+3

	    //GOT FIT
	    max_val=0.;
	    max_signal=0.;
	    max_Sup=0.;
	    double max_background=0.;
	    int bin_number=0;
	    double cut_val_end=5000;
	    for(double cutval=0;cutval<cut_val_end;cutval+=step_size){//loop over cutvals
	      
	      //number of background: integrate fit line from cutval to inf. fitline y=exp(yin+slopefit*x). Integral is: -1/slope*exp(yint+slopefit*cutval)
	      //cout<<"yint and slope_fit are "<<yint<<" "<<slope_fit<<"\n";
	      num_background = -1*exp(yint+slope_fit*cutval)/slope_fit;
	      
	      if(slope_fit >0){
		cout<<"slope isnt decay! m = "<<m<<"\n";
		good_bin=0;
		break;
		
	      }
	      num_signal=0.;
	      
	      bin_number = (int)( round(cutval/step_size)+1.);
	      // cout<<"cutval/stepsize is "<<cutval/step_size<<" binnumber is "<<bin_number<<"\n";
	      //cout<<"cut_val is "<<cutval<<" hopt->GetXaxis()->GetBinCenter(bin_number) is "<<hopt->GetXaxis()->GetBinCenter(bin_number)<<"\n";
	      //cout<<"num_background is "<<num_background<<"\n";
	     
	      
	      if(num_background <0) continue;
	      //cout<<"num_background originally is "<<num_background<<"\n";
	      num_background*=9;//SET TO 90% background to get correct poisson
	      //cout<<"num_backgorund now is "<<num_background<<"\n";
	      hbackground->SetBinContent(bin_number,num_background);
	      if(num_background >100) continue;

	      if(num_background > max_background) max_background = num_background;

	      for(int i=1;i< SNR_vector_sim[m].size();i++){
		if(peakVal_vector_sim[m][i] > 0.075) {
		  if( SNR_vector_sim[m][i] >= slope*peakVal_vector_sim[m][i]+cutval){
		    num_signal+=healpix_bin_weights_sim[m][i];
		  }
		}
	
	      }
	      //cout<<"cutval is "<<cutval<<" num_signal is "<<num_signal<<"\n";
	      
	      //get value of mean for 90% confidence
	      double end_int = 100;
	      TF1 *poisson = new TF1("mypoisson","(pow(x+[0],[0])/TMath::Gamma([0]+1))*exp(-1*([0]+x))",0.,end_int);//0 is lower limit of poisson, 100 is higher limit?
	      poisson->SetParameter(0,num_background);
	      //  cout<<"poisson at 1 is "<<poisson->Eval(0)<<"\n";
	      
	      //poisson->SetParameter(0,num_background);
	      //cout<<"num_background is "<<num_background<<"\n";
	      double hand_int=0.;
	      double delta_int= 0.1;//(num_background/10);
	      //double end_int=10000;
	      //cout<<"delta_int is "<<delta_int<<"\n";
	      /* for(double i=0.;i<end_int;i+=delta_int){
		
	      hand_int+=(poisson->Eval(i)*delta_int);
		
	      }
	      */
	      //cout<<"integral by machine is "<<poisson->Integral(0.,end_int)<<"\n";
	      double normalization=poisson->Integral(0.,end_int);
	      //cout<<"normalization is "<<normalization<<"\n";
	     
	      double Sup=0.;
	      double CL=0.9;
	      double start_int = num_background-20;
	      
	      //if(normalization <1E-10) continue;
	      
	      // if(num_background <1) end_int = 1;
	      double val=0.;
	      if(start_int <0) start_int =0;
	      for(double s=0;s<end_int;s+=delta_int){
		//cout<<"s is "<<s<<"\n";
		hand_int=0.;
		val = poisson->Integral(0.,s);
		/*	for(double i=0.;i<s;i+=delta_int){
		  // cout<<"i is "<<i<<" eval is "<<poisson->Eval(i)<<"\n";
		  hand_int+=(poisson->Eval(i)*delta_int);
		}
		*/
	      
		
		if(val > CL*normalization){	
		    Sup=s;
		    break;
		  
		  //cout<<"breaking sup = <<s<<"\n";
		  
		}
	      }
	     
	      
	      //cout<<"cut_val, num_signal,num_background, Sup,normalization, num_signal/Sup "<<cutval<<" "<<num_signal<<" "<<num_background<<" "<<Sup<<" "<<normalization<<" "<<num_signal/Sup<<"\n";
	      if(num_signal <.01){
		cut_val_end=cutval;
		cout<<"num_signal is "<<num_signal<<" cutval is "<<cutval<<"\n";
		cut_max=cutval/2;
	      }
	      if (num_signal/Sup >= max_val){
		max_val = num_signal/Sup;
		cut_max = cutval;
	      }
	      // cout<<"num_signal is "<<num_signal<<" max_signal is "<<max_signal<<"\n";
	      if(num_signal > max_signal){
		max_signal = num_signal;
		cut_signal = cutval;
	      }

	      if(Sup > max_Sup){
		max_Sup = Sup;
		cut_Sup = cutval;
	      }
	      //cout<<"cutval is "<<cutval<<"\n";
	      //cout<<" normalization is "<<normalization<<" Sup is "<<Sup<<" num_signal/Sup "<<num_signal/Sup<<"\n";
	      //cout<<"cutval+1,, signal/sup "<<cutval+1<<" "<<num_signal/Sup<<"\n";
	      // cout<<"cut_val is "<<cutval<<" "<<bin_number<<" "<<num_signal/Sup<<"\n";
	      hopt->SetBinContent(bin_number,num_signal/Sup);
	      hnum->SetBinContent(bin_number,num_signal);
	      hdenom->SetBinContent(bin_number,Sup);
	      
	      delete poisson;
	      good_bin=1;
	    }//cutval
	    if(good_bin==0) continue;
	    int yint_bin =hopt->GetMaximumBin();
	    //cout<<"yint_bin is "<<yint_bin<<"\n";
	    double y_int = hopt->GetXaxis()->GetBinCenter(yint_bin);
	    //cout<<"value of max is "<<hopt->GetMaximum()<<" y_int is "<<y_int<<"\n";
	    //cout<<"slope and y_int are "<<slope<<" "<<y_int<<" max_val is "<<max_val<<"\n";
	    myfile<<m<<"\t"<<y_int<<"\n";

	    
	    int_cut_rotated=0;
	    num_cut_rotated=0.;
	    num_cut_rotated_sim=0.;

	    
	    for(int n=1;n<numberevents;n++){
	      if(peakVal_vector[m][n] < 0.075) continue;

	      if(SNR_vector[m][n] >= slope*peakVal_vector[m][n]+y_int){
		    num_cut_rotated+=healpix_bin_weights[m][n]; //Apply peak val cut now so these events appear in rotated plot
		    cout<<"eventNumber is "<<eventNumber_vector[m][n]<<"\n";
		    int_cut_rotated++;
		}
	    }
	    
	     for(int n=1;n< SNR_vector_sim[m].size();n++){
		if(peakVal_vector_sim[m][n] > 0.075) {
		  if( SNR_vector_sim[m][n] >= slope*peakVal_vector_sim[m][n]+y_int){
		    num_cut_rotated_sim+=healpix_bin_weights_sim[m][n];
		  }
		}
	
	      }
	     
	     cout<<"number of surviving events in bin = "<<m<<" is "<<num_cut_rotated<<" integer events is "<<int_cut_rotated<<" 90% sample is "<< -9*exp(yint+slope_fit*y_int)/slope_fit<<" simulated set surviving is "<<num_cut_rotated_sim<<"\n";

	     myfile2<<"\t"<<-9*exp(yint+slope_fit*y_int)/slope_fit<<"\t"<<num_cut_rotated_sim<<"\n";
	    //if(y_int > max_val) max_val=y_int;

	    double x2 = -y_int/slope;
	    double y2=0.;
	    if(x2 >1) {
	      x2=1;
	      y2 = slope+y_int;
	    }
	    char printer[256];
	    TCanvas *c0 = new TCanvas("c0","c0",800,800);
	    TH2F *haxes = new TH2F("axes","Optimization Distribution;Cut Value;Signal/S_{up}",10,0,cut_max*2, 10,0, max_val*1.1);
	    haxes->Draw();
	    hopt->Draw("same C*");
	    sprintf(printer,"Optimization_%i.png",m);
	    c0->Print(printer);

	    hnum->SetMarkerColor(kBlue);
	    hdenom->SetMarkerColor(kRed);

	    TCanvas *c1 = new TCanvas("c1","c1",800,800);
	    TH2F *haxes1 = new TH2F("axes1","Signal Distribution;Cut Value;Signal",10,0,cut_max*2, 10,1E-1, max_signal*1.1);
	    c1->SetLogy();
	    haxes1->Draw();
	    hnum->Draw("same C*");
	    sprintf(printer,"Signal_%i.png",m);
	    c1->Print(printer);
	    //cout<<"1 \n";
	    TCanvas *c2 = new TCanvas("c2","c2",800,800);
	    TH2F *haxes2 = new TH2F("axes2","S_{up} Distribution;Cut Value;S_{up}",10,0.,cut_max*2, 10,1E-1, max_Sup*1.1);
	    c2->SetLogy();
	    haxes2->Draw();
	    hdenom->Draw("same C*");
	    sprintf(printer,"Sup_%i.png",m);
	    c2->Print(printer);
	    //cout<<"2 \n";
	    double maxer=0.;
	    if(max_Sup > max_signal) maxer=max_Sup;
	    else maxer = max_signal;
	     TCanvas *c4 = new TCanvas("c4","c4",800,800);
	    TH2F *haxes4 = new TH2F("axes","Signal, S_{up} Distributions;Cut Value;Counts",10,0,cut_max*2, 10,5E-1, maxer*3);
	    //c4->SetLogy();
	    haxes4->Draw();
	    hnum->SetLineColor(kBlue);
	    hdenom->SetLineColor(kRed);
	    hnum->Draw("same C*");
	    hdenom->Draw("same C*");

	    TLegend leg(.7,.8,1,1,"l");
	    leg.SetFillColor(0);
	    leg.AddEntry(hnum,"Signal");
	    leg.AddEntry(hdenom,"S_{up}");
	    leg.Draw("same");
	    sprintf(printer,"Signal_Sup_%i.png",m);
	    c4->Print(printer);


	    	TH2D *haxes_rotated; 
		if(y_int > max_snr) max_snr=y_int;
		haxes_rotated = new TH2D("rotated axes","Rotated Cross Correlation Distribution;Peak Val of Cross Correlation;SNR of Coherent Waveform",10,0,1,10,0,max_snr+5);
	

	    TLine *line1 = new TLine(0,y_int,x2,y2);
	    line1->SetLineColor(kRed);
	    line1->SetLineWidth(2);

	    TLine *line2 = new TLine(0.075,0,0.075,max_snr+5);
	    line2->SetLineColor(kRed);
	    line2->SetLineWidth(2);

	    
	    TCanvas *c3 = new TCanvas("c3","c3",1600,800);
	    TPaveText *pt = new TPaveText(.5,.6,.8,.9,"NDC");
	    pt->SetFillColor(0);
	    pt->SetTextSize(.02);
	    sprintf(printer,"# of Background Passing: %.3E",num_cut_rotated);
	    pt->AddText(printer);
	    
	    
	    sprintf(printer,"intercept: %g",y_int);
	    pt->AddText(printer);
	
	     TPaveText *pt2 = new TPaveText(.5,.6,.8,.9,"NDC");
	    pt2->SetFillColor(0);
	    pt2->SetTextSize(.02);
	    sprintf(printer,"# of Signal Passing: %.3E",num_cut_rotated_sim);
	    pt2->AddText(printer);


	    c3->Divide(2,1);
	    c3->cd(1);
	    haxes_rotated->Draw();
	    hRotatedCut->Draw("same zcol");
	    line1->Draw("same");
	    line2->Draw("same");
	    pt->Draw("same");

	    c3->cd(2);
	    haxes_rotated->Draw();
	    hRotatedCut_sim->Draw("same zcol");
	    line1->Draw("same");
	    line2->Draw("same");
	    pt2->Draw("same");
	    
	    
	    sprintf(printer,"rotatedCut_%i_optimized.png",m);
	    c3->Print(printer);


	    TCanvas *c5 = new TCanvas("c5","c5",800,800);
	    double max_exponent = log10(max_background);
	    max_exponent = floor(max_exponent);
	    max_exponent = pow(10,max_exponent);
	    cout<<"Cut_max is "<<cut_max<<"\n";
	    TH2F *haxes5 = new TH2F("axes5",";Cut Value;Number of Backgroundu Expected",10,0.,cut_max*2, 10,1E-6, max_background*max_exponent);
	    c5->SetLogy();
	    haxes5->Draw();
	    hbackground->Draw("same C*");
	    sprintf(printer,"background_%i.png",m);
	    c5->Print(printer);



	    TH2F *haxes_diff = new TH2F("axes_diff","Differential Plot;y_int;Number of Weighted Events Cut",200,0,fit_end+10,200,0.5*min_cut,num_cut_max*2);

	    TCanvas *c6 = new TCanvas("c6","c6",880,800);
	    c6->Divide(2,2);
	    c6->cd(1);
	    TPad *p1 = (TPad *)(c6->cd(1));
	    p1->SetTitle("Differential Plot");
	    p1->SetLogy();
	    haxes_diff->Draw();
	    hdiffplot->Draw("same");

	    c6->cd(2);
	     TPad *p2 = (TPad *)(c6->cd(2));
	     p2->SetTitle("Optimization Distribution");
	     haxes->Draw();
	    hopt->Draw("same C*");
	   
	    c6->cd(3);
	    TPad *p3 = (TPad *)(c6->cd(3)); 
	    p3->SetTitle("Signal Distribution");
	    p3->SetLogy();
	    haxes4->Draw();
	    hnum->SetLineColor(kBlue);
	    hdenom->SetLineColor(kRed);
	    hnum->Draw("same C*");
	    hdenom->Draw("same C*");

	    TLegend leg2(.7,.8,1,1,"");
	    leg2.SetFillColor(0);
	    leg2.AddEntry(hnum,"Signal");
	    leg2.AddEntry(hdenom,"S_{up}");
	    leg2.Draw("same");
	    
	    c6->cd(4);
	    TPad *p4 = (TPad *)(c6->cd(4)); 
	    haxes_rotated->Draw();
	    hRotatedCut->Draw("same zcol");
	    line1->Draw("same");
	    line2->Draw("same");
	    pt->Draw("same");

	    sprintf(printer,"Bin_%i.png",m);
	    c6->Print(printer);

	    delete hdiffplot;
	    delete hRotatedCut;
	    
	    delete hopt;
	    delete hnum;
	    delete hdenom;
	    delete hbackground;

	       }//only in good bins
	}//m
	myfile.close();
	myfile2.close();
	return 0;
}
////////////////


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
  RootStyle->SetPadRightMargin (.15);
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
