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
//vector<double> *BASE=0;
vector<double> BASE(n_pix_int,0.);;
double BASE2;
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
  vector<int> skipped_bins;
  vector<int> used_bins;

  vector< vector <double> > peakVal_vector (n_pix_int+1, vector<double>(1));
  vector< vector <double> > peakHilbert_vector (n_pix_int+1, vector<double>(1));
  vector< vector <double> > SNR_vector (n_pix_int+1, vector<double>(1));
  
  vector< vector <double> > healpix_bin_weights (n_pix_int+1, vector<double>(1));
  vector< vector <int> > eventNumber_vector (n_pix_int+1, vector<int>(1));
   ofstream myfile;
  myfile.open("Diff_fit_values_test_615_1.txt");

  ofstream myfile2;
  myfile2.open("Pval_info_626_1.txt");

  ifstream mycutvals;
  mycutvals.open("Cut_values_0501.txt");


  //read in root file
  string filter_name;
  // string temp = "HealPix_partial_1215.root";
  //string temp = "HealPix_partial_111.root";//removed events <50% weight
  //string temp = "HealPix_partial_HPol_120.root";//HPol
  // string temp = "HealPix_partial_0207.root";//HPol
   string temp = "HealPix_partial_0301.root";//HPol
  char *rootfile;
  rootfile = Form(temp.c_str(),filter_name.c_str());
    
  TChain *Binned_tree = new TChain("Binned_Tree");
 
       
 
  Binned_tree->Add(rootfile);
    
     
  int nevents0 = Binned_tree->GetEntries();
  Binned_tree->SetBranchAddress("eventNumber_vector",&eventNumber_vector2);
  Binned_tree->SetBranchAddress("peakVal_vector",&peakVal_vector2);
  Binned_tree->SetBranchAddress("peakHilbert_vector",&peakHilbert_vector2);
  Binned_tree->SetBranchAddress("SNR_vector",&SNR_vector2);
  Binned_tree->SetBranchAddress("healpix_bin_weights",&healpix_bin_weights2);
  Binned_tree->SetBranchAddress("BASE",&BASE2);
  cout<<"nevents0 is "<<nevents0<<"\n";

  ////////GOT ALL INFO FROM EVENTS. NOW CAN USE FOR OPTIMIZATION and PLOTTING//////////
  double max_base=1.;
  double pval=0.;
 
  
  int hilbertFlag=0;
  int cut_flag=0;
  int cut_flag2=0;

  string line;
  vector<int> healpixbin_cuts_vector;
  vector<double> cut_val_vector;
  while(getline(mycutvals,line)){

    istringstream ss(line);
    int number;
    double cut_val_temp;
    ss >> number >> cut_val_temp;
    healpixbin_cuts_vector.push_back(number);
    cut_val_vector.push_back(cut_val_temp);

  }



  for(int m=2990;m<n_pix_int;m++){
    //for(int m=3030;m<3035;m++){
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
    if(BASE[m] > max_base) max_base = BASE[m];
    /*
    for(int k=0;k<5;k++){
      cout<<"peakVal is "<<peakVal_vector[m][k]<<" snr is "<<SNR_vector[m][k]<<"\n";   
    }
    */
  }
  
  
  max_base = 1.1*max_base;
  
      //////OPTIMIZATION! 
      double slope=-20.;
      // double y_int;
      double num_cut=0.;
      double num_cut_old=0.;
      double otherside;
      int fit_start=0;
      int fit_end=0;
      double num_cut_max=0;
      int bad_fit_ctr=0;
      int fullcut_flag=0;
      double min_weight=10.;
      double max_weight=0.;
      double min_cut=0.;
      int num_cut_int=0;
      int num_cut_int_old=0;

      int multbinFlag=0;
      int numberevents=0;
      char printer[256];
      double max_snr=0.;
       double num_cut_peakVal=0.;
      double num_bins_plots=10000;
      double bin_size = 1;
      double fit_events=0.;
      double fit_error=0.;
      
       for(int sloper=38;sloper<39;sloper+=1){
	slope = -1*sloper;
	skipped_bins.clear();
	used_bins.clear();
	cout<<" slope is "<<slope<<"\n";
	bad_fit_ctr=0;
	SNR_vector[n_pix_int].clear();
	SNR_vector[n_pix_int].push_back(1);

	peakVal_vector[n_pix_int].clear();
	peakVal_vector[n_pix_int].push_back(1);
	
	healpix_bin_weights[n_pix_int].clear();
	healpix_bin_weights[n_pix_int].push_back(1);
	
	cout<<"SNR_vector[add] size is "<<SNR_vector[n_pix_int].size()<<"\n";
	TH1D *hChiSquare = new TH1D("ChiSquare",";ChiSquare/NDF;Number of Events",100,0.,200.);
	TH1D *hpVal = new TH1D("pVal",";p Val;Number of HealPix Bins",1000,0.,1.2);
	TH2D *hpvalDist = new TH2D("pvalDist",";PVal;Log(Number of entries in Bin;Number of Bins",1000,0,1.2,1000,1E-5,5);
	//peakHilbertCoherent<-350*peakVal+57.14
	//for(int m=2990;m<n_pix_int;m++){
	  for(int m=3031;m<3032;m++){
	  //if(m>3033) m=3055;
	  //cout<<"m is "<<m<<"\n";
	  max_snr=0;
	  Binned_tree->GetEvent(m);
	  fit_start=0;
	  fit_end=0;
	  num_cut_max=0.;
	  num_cut=0.;
	  num_cut_old=0.;
	  fullcut_flag=0;
	  num_cut_peakVal=0.;
	  min_cut=100.;
	  fit_events=0.;
	  fit_error=0.;
	  //TH1D *hdiffplot = new TH1D("diffplot",";y_intercept;Number of Events Cut;",100,0.,100.);
	  

	  
	  numberevents = peakHilbert_vector[m].size();
	  if(numberevents ==1) continue;
	    // cout<<"Base num and size of vector is "<<BASE[m]<<" "<<peakVal_vector[m].size()<<"\n";
	    TH1D *hdiffplot = new TH1D("diffplot",";y_intercept;Number of Events Cut;",num_bins_plots,0.,num_bins_plots*bin_size);
	    TH1D *hdiffplotWhole = new TH1D("diffplotWhole",";y_intercept;Number of Events Cut;",num_bins_plots,0.,num_bins_plots*bin_size);
	    TH1D *hdiffplotTotal = new TH1D("diffplotTotal",";y_intercept;Number of Events Cut;",num_bins_plots,0.,num_bins_plots*bin_size);

	   
	    hdiffplot->Sumw2();
	    TH2F *hRotatedCut;
	    if(hilbertFlag==0){
	      // hRotatedCut = new TH2F("rotatedcut",";Peak Cross Corr Val;SNR of Coherent Waveform;",200,0,1,500,0,500);
	      hRotatedCut = new TH2F("rotatedcut",";Peak Cross Corr Val;SNR of Coherent Waveform;",200,0,1,2000,0,1000);
	    }
	    else{
	     hRotatedCut = new TH2F("rotatedcut",";Peak Cross Corr Val;Peak Hilbert of Coherent Waveform;",100,0,1,1000,0,1000);
	    }
	     max_weight=0;
	     min_weight=100.;
	     cout<<"m is "<<m<<"\n";
	     if(numberevents<=0){
	       for(int n=0;n<numberevents;n++){
		 SNR_vector[n_pix_int].push_back(SNR_vector[m][n]);
		 peakVal_vector[n_pix_int].push_back(peakVal_vector[m][n]);
		 healpix_bin_weights[n_pix_int].push_back(healpix_bin_weights[m][n]);
		 peakHilbert_vector[n_pix_int].push_back(peakHilbert_vector[m][n]);
	       }
	       cout<<"m is "<<m<<" continuing! \n";
	       delete hdiffplot;
	       
	       delete hdiffplotWhole;
	       delete hdiffplotTotal;
	       delete hRotatedCut;
	       continue;
	     }
	     
	     for(int n=1;n<numberevents;n++){
		if(peakVal_vector[m][n] < 0.075){
		    num_cut_peakVal+=healpix_bin_weights[m][n]; //Apply peak val cut now so these events appear in rotated plot
		 }
		 //cout<<"hilbert is "<<peakHilbert_vector[m][n]<<"\n";
		 if(hilbertFlag==0){
		   if(SNR_vector[m][n]>max_snr) {
		     //cout<<"max_snr is from "<<eventNumber_vector[m][n]<<"\n";
		     max_snr = SNR_vector[m][n];
		   }
		 }
		 else{
		   if(peakHilbert_vector[m][n]>max_snr) max_snr = peakHilbert_vector[m][n];
		 }

		 if(hilbertFlag==0){
		   hRotatedCut->Fill(peakVal_vector[m][n], SNR_vector[m][n],healpix_bin_weights[m][n]);//for plotting
		 }
		 else{
		   hRotatedCut->Fill(peakVal_vector[m][n],peakHilbert_vector[m][n]);
		 }
	     }


	     for(double y_int=0;y_int<10000;y_int+=bin_size){
	       num_cut=0.;
	       num_cut_int=0;
	       
	       for(int n=1;n<numberevents;n++){
		 cut_flag=0;
		 cut_flag2=0;
		 if(peakVal_vector[m][n] < 0.075) continue; //Apply peak val cut now so these events appear in rotated plot
		 if(hilbertFlag==0){
		    if(SNR_vector[m][n] < slope*peakVal_vector[m][n]+y_int && SNR_vector[m][n] > slope*peakVal_vector[m][n]+y_int-bin_size) cut_flag=1;
		    //if(SNR_vector[m][n] < slope*peakVal_vector[m][n]+y_int ) cut_flag=1;
		 }
		 else{
		   if(peakHilbert_vector[m][n] < slope*peakVal_vector[m][n]+y_int && peakHilbert_vector[m][n] > slope*peakVal_vector[m][n]+y_int-bin_size) cut_flag=1; 
		 }

		 //hdiffplot->Fill(y_int,0);
		 if(cut_flag==1){ 
		   cout<<"y_int is "<<y_int<<" event is "<< eventNumber_vector[m][n]<<" weight is "<<healpix_bin_weights[m][n]<<"\n";
		   hdiffplot->Fill(y_int,healpix_bin_weights[m][n]);
		   if(healpix_bin_weights[m][n]>max_weight) max_weight = healpix_bin_weights[m][n];
		   hdiffplotTotal->Fill(y_int);
		   if(healpix_bin_weights[m][n]==1){
		     hdiffplotWhole->Fill(y_int);
		   }
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
		   num_cut_int++;
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
	       //num_cut_int-=num_cut_int_old;
	       //cout<<"num_cut is "<<num_cut<<" num_cut_old is "<<num_cut_old<<" y_int is "<<y_int<<"\n";
	       if(num_cut==0){
		 //continue;
	       }
	       if(num_cut > num_cut_max){
		 num_cut_max=num_cut;
		 fit_start = y_int+1;
		 fit_events+=num_cut;
		 //cout<<"m is "<<m<<" num_cut is "<<num_cut<<" "<<fit_start<<"\n";
	       }
	       if(num_cut < min_cut && num_cut>0){
		 min_cut=num_cut;
	       }
	       
	       //hdiffplot->Fill(y_int,num_cut);
	       cout<<"y_int is "<<y_int<<" num_cut is "<<num_cut<<"\n";
	       num_cut_old = num_cut+num_cut_old;
	       //num_cut_int_old +=num_cut_int;
	       //cout<<"num_cut_old is "<<num_cut_old<<" numberevents is "<<numberevents<<" y_int is "<<y_int<<"\n";
	       // if( (num_cut_old >0 && num_cut ==0 && num_cut_old > numberevents/2) || num_cut_old==numberevents){
	       // cout<<"num_cut,old,Base, are "<<num_cut<<" "<<num_cut_old<<" "<<BASE[m]<<"\n";
	       if(num_cut_old > BASE[m]-num_cut_peakVal-1E-10){
	       //if(num_cut_int_old == numberevents){
		 //cout<<"num_cut_old is "<<num_cut_old<<" "<<numberevents<<" "<<y_int<<"\n";
		 //cout<<"num_cut is "<<num_cut_old<<" BASE[m] is "<<BASE[m]<<"\n";
		 fit_end = y_int+1;
		 y_int = 1E5;
		 fullcut_flag=1;
		 cout<<"fit_events is "<<fit_events<<" num_cut_old is "<<num_cut_old<<"\n";
		 fit_events = num_cut_old - fit_events; 
		 cout<<"fit_events is now "<<fit_events<<"\n";
		 break;
	       }
	      
	       
	     }//y_int

	     TH2F *haxes_diff_0 = new TH2F("axes_diff",";y_int;Number of Events Cut",200,0,fit_end+10,200,0.1*min_cut,num_cut_max*2);
	       TCanvas *c4_0 = new TCanvas("c4_0","c4_0",800,800);
		c4_0->SetLogy();
		
		haxes_diff_0->Draw();
		hdiffplot->Draw("sames");
		sprintf(printer,"diffplot_%i_all_test.png",m,sloper);
		c4_0->Print(printer);

		delete haxes_diff_0;
		delete c4_0;
	

		 cout<<"NUMBER OF BINS USED IS "<<fit_end - fit_start<<"\n";
		cout<<"NUMBER OF EVENTS FIT OVER IS "<<fit_events<<"\n";

	     // cout<<"doing fit \n";
	    if(fit_end<fit_start) fit_end = fit_start;
	    cout<<"fit start, end are "<<fit_start<<" "<<fit_end<<"\n";
	    //if(fullcut_flag==0) fit_end = 3*fit_start;
	    //if(m==3017) cout<<"m is "<<m<<" numberevents is "<<numberevents<<" fit_start and end is "<<fit_start<<" "<<fit_end<<"\n";
	    TF1 *fit;
	    TF1 *fit_pseudo;
	    
	    TH2F *haxes_diff;
	    
	    TH1D *hpseudoVal = new TH1D("pseudoVal","likelihood value;Counts",num_bins_plots*10,0,10000);
	    
	    double likelihood_val;
	    double likelihood_val_pseudo;
	    double edm, errdef;
	    int nvpar,nparx;
	    int num_events = hdiffplot->GetEntries();
	   
	    //cout<<"m is "<<m<<" sloper is "<<sloper<<" fit start and end are "<<fit_start<<" "<<fit_end<<" num_cut_max is "<<num_cut_max<<" fullcut_flag is "<<fullcut_flag<<" num_cut_old is "<<num_cut_old<<" min_weight is "<<min_weight<<"\n";
	    cout<<"min_cut is "<<min_cut<<"\n";
	    if(min_weight >90) min_weight=1E-2;
	    haxes_diff = new TH2F("axes_diff",";y_int;Number of Events Cut",200,0,fit_end+10,200,0.1*min_cut,num_cut_max*2);

	    TH2D *haxes_rotated; 
	    if(hilbertFlag==0){
	      haxes_rotated = new TH2D("rotated axes",";Peak Val of Cross Correlation;SNR of Coherent Waveform",10,0,1,10,0,max_snr+5);
	    }
	    else{
	      haxes_rotated = new TH2D("rotated axes",";Peak Val of Cross Correlation;Peak hilbert of Coherent Waveform",10,0,1,10,0,max_snr+5);
	    }
	    
	    TLine *line1 = new TLine(0.075,0,0.075,max_snr+5);
	    line1->SetLineColor(kRed);
	    line1->SetLineWidth(2);

	    
	    TCanvas *c3 = new TCanvas("c3","c3",800,800);
	    haxes_rotated->Draw();
	    hRotatedCut->Draw("same zcol");
	    line1->Draw("same");
	    
	    //cutLine->Draw("same");
	    sprintf(printer,"rotatedCut_%i.png",m);
	    c3->Print(printer);
	    
	    

	    
	    if(fit_end < fit_start+3){
	      //cout<<"doing fit \n";
	       for(int n=0;n<numberevents;n++){
		 SNR_vector[n_pix_int].push_back(SNR_vector[m][n]);
		 peakVal_vector[n_pix_int].push_back(peakVal_vector[m][n]);
		 healpix_bin_weights[n_pix_int].push_back(healpix_bin_weights[m][n]);
		 peakHilbert_vector[n_pix_int].push_back(peakHilbert_vector[m][n]);
	       }
	       cout<<"m is "<<m<<" continuing! had "<<numberevents<<" background "<<" fit_start, end were "<<fit_start<<" "<<fit_end<<"\n";
	       skipped_bins.push_back(m);
	       delete hdiffplot;
	       
	       delete hdiffplotWhole;
	       delete hdiffplotTotal;
	       delete hRotatedCut;
	       myfile2<<"\n";
	       continue;
	    }

	    
	    if(fit_end <= fit_start+3 ){
	      skipped_bins.push_back(m);
	      cout<<"bad fits on m="<<m<<"\n";
	      myfile2<<"\n";
	      continue;
	    }
	    
	    
	    
	    //if(1){
	    double stats[4];
	    hdiffplot->GetStats(stats);
	    for(int k=0;k<fit_end*2+3;k++){
	      cout<<"bin: "<<k<<" value: "<<hdiffplot->GetBinCenter(k)<<" events: "<<hdiffplot->GetBinContent(k)<<" error: "<<hdiffplot->GetBinError(k)<<" sumw: "<<stats[0]<<" sumw2: "<<stats[1]<<"\n";
	    }
	    
	    if(fit_end > fit_start +3 && numberevents >1){
	      cout<<"doing fit \n";
	      // haxes_diff = new TH2F("axes_diff",";y_int;Number of Events Cut",200,0,fit_end+5,200,.9,num_cut_max*2);
	      // if(numberevents <100){
	      hdiffplot->Fit("expo","QR WL","",fit_start,fit_end);
	      //hdiffplotTotal->Fit("expo","QR","",fit_start,fit_end);
	      //hdiffplot->Fit("expo","QR WL","",fit_start,fit_start+3);
		//}
		//else{
		//	hdiffplot->Fit("expo","QR L","",fit_start,fit_end);
		//}
	      //cout<<"CHANGED TO TOTAL PLOT!!!!!!!!!!\n\n\n";
	       fit = hdiffplot->GetFunction("expo");
	       //fit = hdiffplotTotal->GetFunction("expo");
	      fit->SetLineColor(kRed);
	      fit->SetLineWidth(2);
	      double chisquare = fit->GetChisquare();
	      cout<<"chisquare is "<<chisquare<<"\n";
	      double NDF = fit->GetNDF();
	      
	      
	      double slope_fit = fit->GetParameter(1);
	      double yint = fit->GetParameter(0);
	      
	      myfile<<m<<"\t"<<slope_fit<<"\t"<<yint<<"\t"<<fit_start<<"\t"<<fit_end<<"\n";
	      double slope_error = fit->GetParError(1);
	      double yint_error = fit->GetParError(0);
	      myfile2<<m<<"\t"<<slope_fit<<"\t"<<slope_error<<"\t"<<yint<<"\t"<<yint_error<<"\t";
	      cout<<"slope_fit and yint is "<<slope_fit<<" "<<yint<<"\n";
	      cout<<"errors are "<<slope_error<<" "<<yint_error<<"\n";
	      TVirtualFitter *fitter = TVirtualFitter::Fitter(hdiffplot);
	      fitter->GetStats(likelihood_val, edm, errdef,nvpar,nparx);
	      //double pval1 = fitter->GetProb();
	      //cout<<"bin is "<<m<<"\n";
	      //cout<<"numberevents is "<<numberevents<<"\n";
	      //cout<<"likelihood, edm, errdef, nvpar, nparx are "<<likelihood_val<<" "<<edm<<" "<<errdef<<" "<<nvpar<<" "<<nparx<<"\n";
	      //cout<<"chisquare for fit is "<<chisquare<<" degrees of freedom "<<NDF<<" pval is "<<pval<<"\n";
	      cout<<"likelihood_val is "<<likelihood_val<<"\n";
	      fit_error = edm;
	     
		cout<<"Error in fitter is "<<fit_error<<"\n";
	      if(slope_fit >0) {
		cout<<"bin "<<m<<" is not a decay! \n";
		skipped_bins.push_back(m);
		myfile2<<"\n";
		continue;
	      }
	      double x_val;
	      double y_val;
	      double y_val_whole;
	      double yval_total;
	      double yval_fit;
	      double log_val;
	      double num_events=0.;
	      double num_events_whole=0;
	      double num_events_total=0;
	      double max_likelihood=likelihood_val+5;
	      int min_bin_x=0;
	      int max_bin_x;
	      for(int n=1;n<10001;n++){
		x_val = hdiffplot->GetBinCenter(n);
		//x_val = hdiffplotTotal->GetBinCenter(n);
		//cout<<"xval is "<<x_val<<" fit_start is "<<fit_start<<"\n";
		if(x_val >= fit_start && x_val <= fit_end){
		  if(min_bin_x==0) min_bin_x=n;
		  max_bin_x=n;
		  y_val = hdiffplot->GetBinContent(n);
		  // y_val = hdiffplotTotal->GetBinContent(n);
		  yval_fit = exp(yint+slope_fit*x_val);
		  yval_total =hdiffplotTotal->GetBinContent(n);
		  y_val_whole = hdiffplotWhole->GetBinContent(n);
		  num_events+=y_val;
		  num_events_whole+=y_val_whole;
		  num_events_total+=yval_total;
		  /* if(y_val>0){
		    log_val = yval_fit - y_val +y_val*log(y_val/yval_fit);
		  }
		  else{
		    log_val = yval_fit - y_val;
		  }

		  cout<<"y_val, yval_fit are "<<y_val<<" "<<yval_fit<<" log_val is "<<log_val<<"\n";
		  */
		}//x_val>
		  
		

	      }//n

	      double percentwhole =num_events_whole/num_events_total;
	      cout<<"m is "<<m<<" percentwhole is "<<percentwhole<<"\n"; 
	      double num_events_pseudo;
	      double y_rndm;
	      double x_rndm;
	      double y_min = exp(yint+slope_fit*fit_start);
	      double y_max = exp(yint+slope_fit*fit_end);

	      double rndm;
	      double weight;
	      cout<<"y_max, min are "<<y_max<<" "<<y_min<<"\n";
	      cout<<"yint, slope are "<<yint<<" "<<slope_fit<<"\n";
	      //cout<<"yint, slope_fit, fit_start, fit_end are "<<yint<<" "<<slope_fit<<" "<<fit_start<<" "<<fit_end<<"\n";
	      // cout<<"num_events is "<<num_events<<"\n";
	      //cout<<"num_events_whole is "<<num_events_whole<<"\n";
	      //cout<<"num_events_total is "<<num_events_total<<"\n";
	      num_events_pseudo = num_events_total;
	      cout<<"NUM_EVENTS_PSEUDO are "<<num_events_pseudo<<"\n";
	      cout<<"max_weight is "<<max_weight<<"\n";
	      for(int exp_num=0;exp_num<10000;exp_num++){
		//num_events_pseudo= gRandom->Poisson(num_events_total);
		
		TH1D *hpseudo = new TH1D("psuedo",";xval;yval",num_bins_plots,0.,10000.);
					
		hpseudo->Sumw2();
		//cout<<"num_events_pseudo is "<<num_events_pseudo<<"\n";
		
		for(int j=0;j<num_events_pseudo;j++){
		  y_rndm = gRandom->Rndm();
		  x_rndm = (y_max-y_min)*y_rndm + y_min;
		  //cout<<"y_rndm is "<<y_rndm<<" x_rdnm is "<<x_rndm;
		  x_rndm = log(x_rndm);
		 
		  x_rndm = (x_rndm - yint)/slope_fit;
		  //x_rndm = floor(x_rndm);
		  //cout<<" final x is "<<x_rndm<<"\n";
		  //cout<<"y_rndm, y_max,y_min, x_rndm is "<<y_rndm<<" "<<y_max<<" "<<y_min<<" "<<" "<<x_rndm<<"\n";
		 
		  //////whole or frac area?
		  rndm = gRandom->Rndm();
		  
		  if(rndm < (1-percentwhole) ) {
		    if(max_weight <1){
		      weight =  (max_weight)*gRandom->Rndm();
		      // if(weight>=1)  weight = (max_weight/num_events_total + max_weight)*gRandom->Rndm();
		    }
		    else {
		      weight = gRandom->Rndm();
		    }
		  }
		  
		  else weight=1;
		  
		  if(weight >1) weight=1;
		  //cout<<"j is "<<j<<" x_rndm is "<<x_rndm<<" weight is "<<weight<<"\n";
		  hpseudo->Fill(x_rndm,weight);

		  
		}//j

	
		//cout<<"fit_start, end are "<<fit_start<<" "<<fit_end<<"\n";
		TF1 *myexpo = new TF1("myexpo","exp([0]+[1]*x)");
		myexpo->FixParameter(0,yint);
		myexpo->FixParameter(1,slope_fit);
		//cout<<"yint and slope_fit are "<<yint<<" "<<slope_fit<<"\n";
		hpseudo->Fit("myexpo","QR WL","",fit_start,fit_end);
		
		  //myfit->SetParameter(0,yint);
		  //myfit->SetParameter(1,slope_fit);
		 //hpseudo->Fit("expo");
		//hpseudo->Fit("expo","QR WL","",fit_start,fit_end);
		   fit_pseudo = hpseudo->GetFunction("myexpo");

		   
		   //fit_pseudo->FixParameter(0,yint);
		   //hpseudo->Fit("expo","QR WL B","",fit_start,fit_end);
		   //cout<<"parametrs are "<<fit_pseudo->GetParameter(0)<<" "<<fit_pseudo->GetParameter(1)<<"\n";
		  /*TF1 *fit2 = hpseudo->GetFunction("expo");
		   fit2->SetLineColor(kRed);
		  fit2->SetLineWidth(2);
		  double y_max = hpseudo->GetMaximum();
		if(exp_num==0){
		  cout<<"num_events_psuedo is "<<num_events_pseudo<<"\n";
		  TH2F *axes_fit = new TH2F("axes_fit",";xval;yval",10,0,fit_end+2,10,1,1.1*y_max);
		  TCanvas *tester = new TCanvas("tester","tester",800,800);
		  tester->SetLogy();
		  axes_fit->Draw();
		  hpseudo->Draw("same");
		  sprintf(printer,"pseudo_%i_%i.png",m, sloper);
		  tester->Print(printer);

		  delete tester;
		  delete axes_fit;
		}
		  */ 
		  //cout<<"slope_fit = "<< fit_pseudo->GetParameter(1)<<" y_int is "<<fit_pseudo->GetParameter(0)<<"\n";
	     
		TVirtualFitter *fitter = TVirtualFitter::Fitter(hpseudo);
		fitter->GetStats(likelihood_val_pseudo, edm, errdef,nvpar,nparx);
		//cout<<"num_events_psuedo, likelihood is "<<num_events_pseudo<<" "<<likelihood_val_pseudo<<"\n";
		hpseudoVal->Fill(likelihood_val_pseudo);
		
		if(likelihood_val_pseudo > max_likelihood) max_likelihood = likelihood_val_pseudo;
		/*
		if(likelihood_val_pseudo < likelihood_val+10){
		  TH2F *haxes_diff_pseudo = new TH2F("axes_diff",";y_int;Number of Events Cut",200,0,fit_end+10,200,0.1*min_cut,num_cut_max*2);
		 TCanvas *c4_pseudo = new TCanvas("c4_pseudo","c4_pseudo",800,800);
		c4_pseudo->SetLogy();
		
		haxes_diff_pseudo->Draw();
		hpseudo->Draw("sames");
		sprintf(printer,"diffplot_%i_%i_%i.png",m,likelihood_val_pseudo,likelihood_val);
		c4_pseudo->Print(printer);
		}
		*/
		//cout<<"likelihood_val_pseudo is "<<likelihood_val_pseudo<<"\n";

		/*
			TF1 *myexpo2 = new TF1("myexpo2","exp([0]+[1]*x)",fit_start,fit_end);
		  myexpo2->FixParameter(0,yint);
		  myexpo2->FixParameter(1,slope_fit);
		  myexpo2->SetLineColor(kRed);

		  

		  TCanvas *c_pseudo = new TCanvas("c_pseudo","c_pseudo",800,800);
		  c_pseudo->SetLogy();
		  haxes_diff->Draw();
		  hpseudo->Draw("sames");
		  myexpo2->Draw("same");

		  TPaveText *pt_like = new TPaveText(.7,.7,.9,.9,"NDC");
		  pt_like->SetFillColor(0);
		  sprintf(printer,"likelihood: %g",likelihood_val_pseudo);
		  pt_like->AddText(printer);
		  pt_like->Draw("same");
		  sprintf(printer,"pseudo_experiment_%i_%i.png",exp_num,m);
		  c_pseudo->Print(printer);

		  delete myexpo2;
		  delete c_pseudo;
		  delete pt_like;
		*/

		delete hpseudo;
	      }//pseudo experiment


	      
	      //double bin_converter = 10000/(num_bins_plots*10);
	      //hpseudoVal->GetBinWithContent(likelihood_val,min_bin_x,1,1000);
	      // min_bin_x = (int)ceil(likelihood_val/bin_converter);// /0.01
	      min_bin_x=hpseudoVal->GetXaxis()->FindBin(likelihood_val);
		// cout<<"min_bin_x should be "<<hpseudoVal->GetXaxis()->FindBin(likelihood_val)<<"\n";
	      cout<<"min_bin_x is "<<min_bin_x<<" likelihood_val is "<<likelihood_val<<"\n";
	      double scale = hpseudoVal->Integral();
	      pval = hpseudoVal->Integral(min_bin_x+1,100000);
	      cout<<"pval, scale is "<<pval<<" "<<scale<<"\n";
	      
	      pval = pval/scale;
	      //scale = hpseudoVal->GetMaximum();
	      hpseudoVal->Scale(1/scale);
	      double scale_max = hpseudoVal->GetMaximum();
	      cout<<"pval is now "<<pval<<"\n";
	      cout<<"log10 is "<<log10(num_events_total)<<"\n";
	      //if(m==3012 || m==3030 || m==3031 || m==3032 || m==3033 || m==3045 || m==3048 || m==3051 || m==3053 || m==3057 || m==3061 || m==3062 || m==3063 || m==3066 || m==3069 || m==3054){//fill pval with bins we will use! 
	      if(pval >=0){
		hpvalDist->Fill(pval,log10(num_events_total));
	      }
	      //hChiSquare->Fill(chisquare/NDF);
	      hpVal->Fill(pval);
	      //	}
	      
	      //double param0 = fit->GetParameter(0);
	      //double param1 = fit->GetParameter(1);
	      
	      //cout<<"chisquare is "<<chisquare<<" pVal is "<<pval<<"\n";
	      //delete fit;
	      // if(pval>0.999 || pval <0.001){
		TCanvas *clikelihood = new TCanvas("clikelihood","clikelihood",880,800);
		TH2D *haxes_like = new TH2D("Axes_like",";Likelihood;Normalized Number of Counts",10,0,max_likelihood,10,0,scale_max*1.5);
		TLine *line = new TLine(likelihood_val,0,likelihood_val,scale_max*1.5);
		line->SetLineColor(kRed);
		haxes_like->Draw();
		line->Draw("same");
		hpseudoVal->Draw("same");
		sprintf(printer,"likelihood_%i_%i_test.png",m,sloper);
		clikelihood->Print(printer);

		delete clikelihood;
		delete line;
		//delete line1;
		delete haxes_like;
		delete hRotatedCut;
		
		myfile2<<pval<<"\t"<<fit_end-fit_start<<"\t"<<fit_events<<"\n";
		if(pval<0.05){
		bad_fit_ctr++;
		cout<<"m is "<<m<<" sloper is "<<sloper<<" numevents is "<<numberevents<<" pval is "<<pval<<"\n";
		skipped_bins.push_back(m);
		//myfile2<<"\n";
		continue;
	      }
		

		//cout<<"max snr is "<<max_snr<<"\n";
		/*	TH2D *haxes_rotated; 
		if(hilbertFlag==0){
		 haxes_rotated = new TH2D("rotated axes",";Peak Val of Cross Correlation;SNR of Coherent Waveform",10,0,1,10,0,max_snr+5);
		}
		else{
		  haxes_rotated = new TH2D("rotated axes",";Peak Val of Cross Correlation;Peak hilbert of Coherent Waveform",10,0,1,10,0,max_snr+5);
		}
		TLine *line1 = new TLine(0,50,0.2,0);
		line1->SetLineColor(kRed);
		line1->SetLineWidth(2);
		TCanvas *c3 = new TCanvas("c3","c3",800,800);
		haxes_rotated->Draw();
		hRotatedCut->Draw("same zcol");
		//line1->Draw("same");
		//cutLine->Draw("same");
		sprintf(printer,"rotatedCut_%i.png",m);
		c3->Print(printer);
		
		delete clikelihood;
		delete line;
		delete line1;
		delete haxes_like;
		delete hRotatedCut; 
		*/
		//}
	      
	      
	      // }//Fit_end > fit_Start
	   


	    /*
	      TCanvas *c3 = new TCanvas("c3","c3",800,800);
	      hRotatedCut->Draw("zcol");
	      cutLine->Draw("same");
	      sprintf(printer,"rotatedCut_%i.png",m);
	      c3->Print(printer);
	    */
		//myfile2<<pval<<"\t"<<fit_end-fit_start<<"\t"<<fit_events<<"\n";
	      cout<<"m is "<<m<<" pval is "<<pval<<"\n";
	      //if(1){
	      if(fit_end > fit_start+3 && numberevents >1){
		used_bins.push_back(m);
		//if(1){
	    //if(sloper==22){ 
	      //cout<<"here 2 \n";
	      //if(pval<0.001 || pval>.999){
	      	cout<<"PVAL IS "<<pval<<"\n";
		cout<<"NUMBER OF BINS USED IS "<<fit_end - fit_start<<"\n";
		cout<<"NUMBER OF EVENTS FIT OVER IS "<<fit_events<<"\n";
		cout<<"Error in fitter is "<<fit_error<<"\n";
		  

		double cut_val_bin=-1.;

		for(int kk=0;kk<healpixbin_cuts_vector.size();kk++){
		  if(healpixbin_cuts_vector[kk]==m){
		    cut_val_bin = cut_val_vector[kk];
		    //cout<<"m is "<<m<<" cut_val is "<<cut_val_bin<<"\n";
		  }
		}
		
		TLine *cut_val_line; 

		if(cut_val_bin >0){
		  cut_val_line = new TLine(cut_val_bin,0,cut_val_bin,num_cut_max*2);
		  cut_val_line->SetLineColor(kRed+2);
		}
		


		TCanvas *c4 = new TCanvas("c4","c4",800,800);
		c4->SetLogy();
		
		haxes_diff->Draw();
		hdiffplot->Draw("sames");
		if(cut_val_line > 0){
		  cut_val_line->Draw("sames");
		}
		c4->Update();
		//hdiffplot->SetStats(1);
		
		TPaveText *pt = new TPaveText(.7,.7,.9,.9,"NDC");
		pt->SetFillColor(0);
		sprintf(printer,"pval: %g",pval);
		pt->AddText(printer);
	     
		sprintf(printer,"slope: %g",slope_fit);
		pt->AddText(printer);


		sprintf(printer,"intercept: %g",yint);
		pt->AddText(printer);
	
		pt->Draw("same");
		// TLatex *myt = new TLatex(0.,0.,printer);
	      
		//  TLatex *myt1 = new TLatex(1.,1.,printer);
	      
		//TLatex *myt2 = new TLatex(.8,.5,printer);
	      //myt->AddText(printer);
	      
	      // myt->Draw("same");
	      //myt1->Draw("same");
	      //myt2->Draw("same");

	      
	      //TList *list = p2->GetListOfLines();
	      //list->Add(myt);
	      //p2->AddText(printer);
	      //p2->AddText("TEST! \n");
	      //p2->AddText("TEST");
	      //hdiffplot->SetStats(0);
	      c4->Modified();
	      //c4->Update();
		sprintf(printer,"diffplot_%i_%i_%i_test.png",m,sloper,whole_flag);
		c4->Print(printer);
		delete c4;
		delete pt;
		delete haxes_diff;
		// }
	      }
	    
	    }//fit_end > fit_start+3
	    //cout<<"about to delete \n";
	    delete hdiffplot;
	    //delete hRotatedCut;
	    delete hpseudoVal;
	     
	     
	    delete hdiffplotWhole;
	    delete hdiffplotTotal;
	    
	    //cout<<"deleted \n";
	    //delete fit;
	    //delete c3;
	 

	    
	}//m==0
	
	
	 
	gStyle->SetPalette(55);
      //gStyle->SetNdivisions(10,"z");
	gStyle=color;
	gStyle->SetOptStat(0000000000);
	TCanvas *cpvalDist = new TCanvas("cpvalDist","cpvalDist",880,800);
	hpvalDist->Draw();
	sprintf(printer,"pvalDist_%i_test.png",sloper);
	cpvalDist->Print(printer);
	delete cpvalDist;

     cout<<"bad_fit_ctr "<<bad_fit_ctr<<"\n";
     TCanvas *c5 = new TCanvas("c5","c5",800,800);
     //c5->SetLogy();
     hChiSquare->Draw();
     sprintf(printer,"ChiSquare_%i.png",(int)abs(slope));
     //c5->Print(printer);
     //TH2F *hpval_axes = new TH2F("pval_axes",";pval;Counts",10,0,1.1);
     TCanvas *c6 = new TCanvas("c6","c6",800,800);
     //c6->SetLogy();
     // hpval_axes->Draw();
     hpVal->Draw();
     sprintf(printer,"pVal_%i_test.png",(int)abs(slope));
     c6->Print(printer);
     delete c5;
     delete c6;
     delete hpVal;
     delete hpvalDist;
     delete hChiSquare;

     for(int k=0;k<skipped_bins.size();k++){
       cout<<"skipped "<<skipped_bins[k]<<" due to bad fit! \n";
     }
     cout<<"using "<<used_bins.size()<<" bins for analysis! \n";
     for(int k=0;k<used_bins.size();k++){
       cout<<"using "<<used_bins[k]<<"! \n";
     }

       }//sloper

     
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
