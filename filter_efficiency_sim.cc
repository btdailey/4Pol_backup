#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cmath>
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
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
#include <iomanip>
//#include <fftw3.h>

//#ifdef ANITA_UTIL_EXISTS
#include "UsefulAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "RawAnitaHeader.h"
#include "analysis_info_4pol.h"
RawAnitaHeader *fHeadPtr = 0;
UsefulAnitaEvent *fUsefulEventPtr;
UsefulAnitaEvent *fUsefulEventPtr_CW;
analysis_info_4pol *pol4_Ptr = 0;
//#endif

using namespace std;
TStyle* RootStyle();
TStyle *color=RootStyle();
/*
#include "Tools.h"
#include "Constants.h"
#include "Vector.h"
#include "Position.h"
#include "EarthModel.h"
#include "IceModel.h"
#include "Efficiencies.h"
#include "Spectra.h"
#include "Event.h"
#include "Trigger.h"
#include "Detector.h"
#include "Settings.h"
#include "counting.hh"
#include "Primaries.h"
#include "Report.h"

#include "Ray.h"

class EarthModel; //class
*/



int main() { 

   gStyle=color;
   gStyle->SetPalette(1);
   gStyle -> SetOptStat(0000000000000);

   
   int eventNumber_cut=0;
   double peakHilbertCoherent=0.;
   double peakVal=0.;
   double ratioFirstToSecondPeak=0.;
   double thetaMap=0.;
   double phiMap=0.;
   double polFractionCoherent=0.;
   int hwTriggerFlag=0;
   int varnerFlag=0;
   int eventTracedFlag2=0;

   int varnerFlag2;
  
   float nadirFlag;
   
   float payloadBlastFlag;
   
   float didIFilter;
   float didIFilterHoriz;
   float didIFilterAboveSatelliteHoriz;
   float didIFilterAboveSatellite;
    

   int eventctr=0;
   int eventctr_total=0;
   int eventctr_total2=0;
   int datactr=0;
   
   int eventpointed_ctr=0;
   
   int payloadBlastctr=0;
   int hilbertctr=0;
   int ratiopeaksctr=0;
   int crosscorrctr=0;
   int polfractionctr=0;
   int rotatedctr=0;
   int elevationctr=0;
   int triggeredctr=0;
   int tracedctr=0;
   
   int dataeventsctr=0;
   int varnerevents=0;

   int polarization=0;

   int eventNumber=0;

   double max_correlation[4][40][40];
   double RMS_ants[4][40];
   double time_delay[4][40][40];
   
   
   double secondThetaMap; // 4 
  double secondPhiMap; // 4 

   double peak_corr;
   int peak_corr_bin1;
   int peak_corr_bin2;
   double peak_corr_time;

   int realtime;

   int eventnumber_old=0;
   int realtime_old=0;
   double noiseBeforeFilter=0.;
   double snrPeak;
   double peak2peak=0.;

   double sourceLon;
   double sourceLat;
   double SNRs_after_filter[40];
   double peakVoltage[40];
   fHeadPtr = new RawAnitaHeader();
   
   
   int numbins = 200;
   int index_mult = 5;
   numbins=50;
   index_mult=1;
   TH1D *hnum_peak_vol;
   TH1D *hnum_peak_vol2;
   TH1D *hnum_peak_vol3;

   hnum_peak_vol = new TH1D("num_peak_vol","1st Peak;Vpeak/RMS;Number",numbins/index_mult,0,numbins);
  
    int ratio_last=0;
   int peakVal_last=0;
   int hilbert_last=0;
   int polfrac_last=0;
   int rotated_last=0;
   int elevation_last=0;
   int traced_last=0;
   int triggered_last=0;
   


   TH1D *hpeakVal;
   TH1D *hratiopeaks;
   TH1D *hhilbert;
   TH1D *hpolfrac;
   TH1D *hprojected;
   TH1D *htriggered;
   TH1D *helevation;
   TH1D *hrotated;
   TH1D *hfinal;

   int index_ctr=0;
   int index=0;
   TH1D *hdistance = new TH1D("distance",";Distance from source;Number of events",200,0,80E4);
   TH1D *hcross_corr = new TH1D("cross_Cor",";Correlation Value;Number of Events",100,0,1);
   TH2D *hcross_corr_SNR = new TH2D("cross_Cor snr",";Correlation Value;SNR",100,0,1,100,0,numbins);
   TH2D *hcross_corr_cut = new TH2D("cross_Cor cut",";Correlation Value;Peak Hilbert Value (mv)",200,0,1,250,0,500);
   TH2D *hcross_corr_cut_corrected = new TH2D("cross_Cor corrected",";Correlation Value;Corrected Peak Hilbert Value (mv m)",250,0,1,500,0,200E6);

    TH2D *hcross_corr_cut_SNR = new TH2D("cross_Cor cut SNR",";Correlation Value;SNR",200,0,1,numbins/index_mult*10,0,numbins);
    TH2D *hcross_corr_LR_corrected = new TH2D("cross_Cor LR corrected",";Correlation Value;Peak Hilbert Value (mv)",200,0,1,500,0,200E6);
   TH2D *hcross_corr_LR_SNR = new TH2D("cross_Cor Lr SNR",";Correlation Value;SNR",200,0,1,numbins/index_mult*10,0,numbins);
  
   TH2D *hcross_corr_LR = new TH2D("cross_Cor LR",";Correlation Value;Peak Hilbert Value (mv)",200,0,1,250,0,500);
   TH2D *hdeltaTheta = new TH2D("deltaTheta",";SNR;cos(#theta_{L})-cos(#theta_{R})",numbins/index_mult,0,numbins,100,0,1);
   TH2D *hdeltaPhi = new TH2D("deltaPhi",";SNR;#Delta #Phi_{L-R}",numbins/index_mult,0,numbins,100,0,20);
 

   TH2D *hpeakLvsR = new TH2D("hpeakLvsR",";Right Peak Corr;Left Peak Corr",100,0,1,100,0,1);
   TH2D *hpeaksratio = new TH2D("peaksratio",";SNR;#frac{Peak H}{Peak V}",numbins/index_mult,0,numbins,100,0,1);
    TH2D *hcorrratio = new TH2D("corrratio",";SNR;#frac{Peak H}{Peak V}",numbins/index_mult,0,numbins,100,0,1);


    TH2D *hintegrandvsSNR = new TH2D("integrand SNR",";SNR;Log(Peak Correlation between Antennas)",numbins/index_mult,0,numbins,100,0,7.2);

    TH2D *hdenomvsSNR = new TH2D("integrand SNR",";SNR;Log(Denominator of Correlation)",numbins/index_mult,0,numbins,100,0,7.2);

    TH2D *hCorrvstime = new TH2D("Corr time ",";Time Delay(ns);Log(Peak Correlation)",100,-20,20,100,0,7.2);



   TH1D *helevation_dist = new TH1D("elevation_dist",";Elevation Angle;Number of Events",180,-90,90);
   TH1D *helevation_traced = new TH1D("elevation_traced",";Elevation Angle;Number of Events",180,-90,90);
   helevation_traced->SetLineColor(kRed);
   //create histograms that are filled at the end of the code
   hratiopeaks = new TH1D("ratiopeaks","Signal peak After Filter; Peak-to-Peak; #frac{# passing cuts}{# passing trigger}",numbins/index_mult,0,numbins);
   hpeakVal = new TH1D("peakVal","SNR peak After Filter; SNR; Fraction passing the cut",numbins/index_mult,0,numbins);
   hpeakVal->SetLineColor(2);
   hhilbert = new TH1D("hilbert","SNR peak After Filter; SNR; Fraction passing the cut",numbins/index_mult,0,numbins);
   hhilbert->SetLineColor(3);
   hpolfrac = new TH1D("polfrac","SNR peak After Filter; SNR; Fraction passing the cut",numbins/index_mult,0,numbins);
   hpolfrac->SetLineColor(4);
   hrotated = new TH1D("rotated","SNR peak After Filter; SNR; Fraction passing the cut",numbins/index_mult,0,numbins);
   hrotated->SetLineColor(kOrange+2);
   helevation = new TH1D("elevation","SNR peak After Filter; SNR; Fraction passing the cut",numbins/index_mult,0,numbins);
   helevation->SetLineColor(6);
   hprojected = new TH1D("projected","SNR peak After Filter; SNR; Fraction passing the cut",numbins/index_mult,0,numbins);
   hprojected->SetLineColor(7);
   htriggered = new TH1D("triggered","SNR peak After Filter; SNR; Fraction passing the cut",numbins/index_mult,0,numbins);
   htriggered->SetLineColor(8);
    hfinal = new TH1D("final","SNR peak After Filter; SNR; Fraction passing the cut",numbins/index_mult,0,numbins);
   hfinal->SetLineColor(kViolet-1);



   //vectors that will be filled with number of V/RMS per bin
   vector<double> ratiopeaks_vec(numbins,0.);
   
   vector<double> peakval_vec(numbins,0.);
   
   vector<double> hilbert_vec(numbins,0.);
  
   vector<double> polfrac_vec(numbins,0.);
  
   vector<double> rotated_vec(numbins,0.);
  
   vector<double> elevation_vec(numbins,0.);
  
   vector<double> projected_vec(numbins,0.);
   
   vector<double> triggered_vec(numbins,0.);
   
   vector<double> final_vec(numbins,0.);
   //vectors used to normalization of the weight
   vector<double> volt_tracker(numbins,0.);
  
   


   //vector of sums
   vector<double> ratiopeaks_vec_sum(numbins/index_mult,0.);
  
   vector<double> peakval_vec_sum(numbins/index_mult,0.);
  
   vector<double> hilbert_vec_sum(numbins/index_mult,0.);
  
   vector<double> polfrac_vec_sum(numbins/index_mult,0.);
  
   vector<double> rotated_vec_sum(numbins/index_mult,0.);
  
   vector<double> elevation_vec_sum(numbins/index_mult,0.);
  
   vector<double> projected_vec_sum(numbins/index_mult,0.);
  
   vector<double> triggered_vec_sum(numbins/index_mult,0.);
  
   vector<double> final_vec_sum(numbins/index_mult,0.);

   vector<double> volt_tracker_sum(numbins/index_mult,0.);
   
   vector<double> weight_tracker(numbins/index_mult,0.);
   
  
    hilbertctr=0;
    ratiopeaksctr=0;
    crosscorrctr=0;
    polfractionctr=0;
    rotatedctr=0;
    elevationctr=0;
    triggeredctr=0;
    tracedctr=0;
    eventTracedFlag2=0;
    dataeventsctr=0;
    varnerevents=0;
    datactr=0;
     ratio_last=0;
     peakVal_last=0;
     hilbert_last=0;
     polfrac_last=0;
     rotated_last=0;
     elevation_last=0;
     traced_last=0;
     triggered_last=0;
   
   double weight;
  UInt_t eventNumber_icemc;
  int n_icemc=0;
  double lat_icemc;
  double lon_icemc;
  double distance_from_source;
  double arrival_times[48];
  double peak2peak_max=0.;
  double peakSNR_2;
  string temp_icemc = "/data/anita/btdailey/icemc_events/Kotera_max/SimulatedAnitaEventFile.root";
  char *rootfile_icemc = Form(temp_icemc.c_str());
  TChain *event_Tree = new TChain("eventTree");
  event_Tree->Add(rootfile_icemc);
  event_Tree->SetBranchAddress("inu",&eventNumber_icemc);
  event_Tree->SetBranchAddress("weight",&weight);
  event_Tree->SetBranchAddress("Lat_icemc",&lat_icemc);
  event_Tree->SetBranchAddress("Lon_icemc",&lon_icemc);
  event_Tree->SetBranchAddress("arrival_times",arrival_times);
  int nevents_icemc = event_Tree->GetEntries();
  ////////////////

  //read in root file
  string filter_name;
  int run=1;
  string temp = "/data/anita/btdailey/final_filter/simulated/Kotera_4pol/output%d_0.root";
  char *rootfile;
  filter_name="no_noise";
  rootfile = Form(temp.c_str(),run);
  
  TChain *FlagTree = new TChain("ndataflags");
  TChain *TracingTree = new TChain("tdataTracing");
  TChain *PointedTree = new TChain("tdataPointed");
  TChain *pol4_Tree = new TChain("analysis_info_4pol");
  
       
  FlagTree->Add(rootfile);
  TracingTree->Add(rootfile);
  PointedTree->Add(rootfile);
  pol4_Tree->Add(rootfile);
  //extras
     
     for(int extra=1000;extra<25000;extra+=1000){
    
       sprintf(rootfile,"/data/anita/btdailey/final_filter/simulated/Kotera_4pol/output%d_%d.root",run,extra);
	 TFile *fpTest = TFile::Open(rootfile);
	 if(!fpTest){ 
	   cout<<"broke at extra "<<extra<<"\n";
	   //continue;
	   break;
	 }
	 else {
	   delete fpTest;
	   
	  
	   FlagTree->Add(rootfile);
	   TracingTree->Add(rootfile);
	   PointedTree->Add(rootfile);
	   pol4_Tree->Add(rootfile);
	 }
       }//extra

     
     int nevents0 = pol4_Tree->GetEntries();//dataTree->GetEntries();
     int nevents1 = TracingTree->GetEntries();
     int nevents2 = PointedTree->GetEntries();
     
     cout<<"nevents0,1,2 are "<<nevents0<<" "<<nevents1<<" "<<nevents2<<"\n";
   
     dataeventsctr+=nevents0;
     eventctr_total2+=nevents2;
     
     //TracingTree->SetBranchAddress("eventTracedFlag",&eventTracedFlag2);
     TracingTree->SetBranchAddress("eventNumber",&eventNumber);
     TracingTree->SetBranchAddress("realTime",&realtime);
     TracingTree->SetBranchAddress("snrPeakAfterFilter",&snrPeak);
     TracingTree->SetBranchAddress("noiseBeforeFilter",&noiseBeforeFilter);
     TracingTree->SetBranchAddress("sourceLon",&sourceLon);
     TracingTree->SetBranchAddress("sourceLat",&sourceLat);
     TracingTree->SetBranchAddress("peak2peak",&peak2peak);
     TracingTree->SetBranchAddress("phiMap",&phiMap);
     TracingTree->SetBranchAddress("distance_from_source",&distance_from_source);
     
     pol4_Tree->SetBranchAddress("pol4_Ptr",&pol4_Ptr);

    double timer;
    
     cout<<"nevents0 is : "<<nevents0<<"\n";
     cout<<"nevents_icemc is "<<nevents_icemc<<"\n";
     vector< vector<double> > index_array(40, vector<double>(2.));
      double max_volts=-1000;
      double min_volts = 1000;
       
      //double peak2peak=0.;
      double SNRforindex=0.;

     for(int m=0;m<nevents0;m++){
      
	//get events
      
       TracingTree->GetEvent(m);
       pol4_Tree->GetEvent(m);


       phiMap = pol4_Ptr->phiMap[0];
       secondPhiMap = pol4_Ptr->secondPhiMap[0];
       secondThetaMap = pol4_Ptr->secondThetaMap[0];


       eventNumber_cut = eventNumber;
       peakHilbertCoherent=pol4_Ptr->peakHilbertCoherent[0];
       peakVal=pol4_Ptr->peakVal[0];
       ratioFirstToSecondPeak=pol4_Ptr->ratioFirstToSecondPeak[0];
       thetaMap=pol4_Ptr->thetaMap[0];
       // thetaMap = pol4_Ptr->peakThetaFinal[0];
       polFractionCoherent=pol4_Ptr->polFractionCoherent[0];
       
       hwTriggerFlag=pol4_Ptr->hwTriggerFlag[0];
       snrPeak = pol4_Ptr->snrPeakAfterFilter[0];
       //snrPeak = pol4_Ptr->snrPeakAfterFilter2[0];
       //snrPeak = pol4_Ptr->SNR_coherent[0];
       
       //snrPeak = pol4_Ptr->peakSNR_2[0];

       //thetaMap = thetaMap*180/3.14159;
       varnerFlag = pol4_Ptr->varnerFlag[0];
       didIFilter = pol4_Ptr->didIFilter[0];
       eventTracedFlag2 = pol4_Ptr->eventTracedFlag[0];
       varnerFlag2 = pol4_Ptr->varnerFlag2[0];
       double min_corr = min(pol4_Ptr->peakVal_box[2],pol4_Ptr->peakVal_box[3]);
       //peakVal = min_corr;
       peak_corr=0;
       for(int k=0;k<40;k++){
	  RMS_ants[0][k] = pol4_Ptr->RMS_ants[0][k];
	 for (int l=k+1;l<40;l++){
	  
	   max_correlation[0][k][k] = pol4_Ptr->max_correlation[0][k][l];
	   time_delay[0][k][l] = pol4_Ptr->time_delay[0][k][l];

	   if(time_delay[0][k][l] > peak_corr){
	     peak_corr = time_delay[0][k][l];
	     peak_corr_bin1 =k;
	     peak_corr_bin2 =l;
	     peak_corr_time = max_correlation[0][k][l];

	   }


	 }
       }
       //cout<<"peak_corr is "<<peak_corr<<"\n";

       if(eventNumber_cut != eventNumber){
	 cout<<"problem! trees are filled with different event! \n";
	 cout<<"eventNumber_cut, eventNumber are "<<eventNumber_cut<<" "<<eventNumber<<"\n";
       }
       for(int n=n_icemc;n<nevents_icemc;n++){
	 event_Tree->GetEvent(n);
	 //cout<<"n is "<<n<<"\n";
	 //cout<<"eventNumber_icemc is "<<eventNumber_icemc<<" weight is "<<weight<<"\n";
	 
	 if(eventNumber==eventNumber_icemc+1){
	   n_icemc = n;
	   break;
	 }
       }//n

        if(varnerFlag==1 || varnerFlag2==1){
	 continue;
       }


	

	/*	double theta1 =thetaMap*3.14159/180;
	double phi1 = phiMap*3.14159/180;

	double theta2 = secondThetaMap*3.14159/180;
	double phi2 = secondPhiMap*3.14159/180;
	*/
	//	cout<<"theta,phi, 2nd theta,phi, angle is "<<thetaMap<<" "<<phiMap<<" "<<secondThetaMap<<" "<<secondPhiMap<<" "<<cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2)<<"\n";

       if(peak2peak > peak2peak_max) peak2peak_max = peak2peak;
       // index = (int)peak2peak;
       // if(snrPeak>14) snrPeak=14;
       index = (int) snrPeak;

       // if(index ==6) cout<<"eventNumber "<<eventNumber_cut<<" has SNR of "<<snrPeak<<"\n";

       //	if(index==30) cout<<"eventNumber "<<eventNumber_cut<<" has SNR of "<<snrPeak<<"\n";
      
        if(index >= numbins){
	  index=numbins;
	 
	 index_ctr++;
	 
       }

      
       if(index < numbins){
	 volt_tracker[index] = volt_tracker[index]+weight; 
	 hnum_peak_vol->Fill(index);
       }
      

       float limit=0.9;
       // cout<<"distance from source is "<<distance_from_source<<"\n";
       hdistance->Fill(distance_from_source);
       hcross_corr->Fill(peakVal,weight);
       hcross_corr_SNR->Fill(peakVal,snrPeak,weight);
       hcross_corr_cut->Fill(peakVal,peakHilbertCoherent,weight);
       hcross_corr_cut_corrected->Fill(peakVal,peakHilbertCoherent*distance_from_source,weight);
       helevation_dist->Fill(thetaMap,weight);
       hcross_corr_cut_SNR->Fill(peakVal,snrPeak,weight);
       //double min_corr = min(pol4_Ptr->peakVal[2],pol4_Ptr->peakVal[3]);
       hcross_corr_LR->Fill(min_corr,peakHilbertCoherent,weight);
       hcross_corr_LR_corrected->Fill(min_corr,peakHilbertCoherent*distance_from_source,weight);
       hcross_corr_LR_SNR->Fill(min_corr,snrPeak,weight);
      

       hintegrandvsSNR->Fill(snrPeak,log10(peak_corr),weight);
       
       hdenomvsSNR->Fill(snrPeak,log10(RMS_ants[0][peak_corr_bin1]*RMS_ants[0][peak_corr_bin2]),weight);
       
       hCorrvstime->Fill(peak_corr_time,log10(peak_corr),weight);

       	double phiL = pol4_Ptr->phiMap[2];
	double phiR = pol4_Ptr->phiMap[3];

	if(abs(phiL-phiR) >180) {
	  if(phiL >phiR) phiL-=360;
	  else if(phiR > phiL) phiR-=360;
	  //cout<<"phiL and phiR are now "<<phiL<<" "<<phiR<<"\n";
	}


	//angle_cut = pow(phiL-phiR,2)+pow(cos((pol4_Ptr->thetaMap[2]-90)*3.14159/180)-cos((pol4_Ptr->thetaMap[3]-90)*3.14159/180),2);
		double theta1 =(pol4_Ptr->thetaMap[2]-90)*3.14159/180;
	double theta2 =(pol4_Ptr->thetaMap[3]-90)*3.14159/180;

	double phi1 = phiL*3.14159/180;
	double phi2 = phiR*3.14159/180;

	hdeltaTheta->Fill(snrPeak,abs(cos(theta1)-cos(theta2)));
	hdeltaPhi->Fill(snrPeak,abs(phiL-phiR));

	hpeakLvsR->Fill(pol4_Ptr->peakVal[3],pol4_Ptr->peakVal[2]);
       
       double ratioHV = pol4_Ptr->peakHilbertCoherent[1] / pol4_Ptr->peakHilbertCoherent[0];
       hpeaksratio->Fill(pol4_Ptr->snrPeakAfterFilter[0],ratioHV);

       ratioHV = pol4_Ptr->peakVal[1] / pol4_Ptr->peakVal[0];
       hcorrratio->Fill(pol4_Ptr->snrPeakAfterFilter[0],ratioHV);
       



       
       if (didIFilterHoriz>1 || didIFilter>1){
	 limit=.85;
       }
       else if (didIFilterAboveSatellite>0 || didIFilterAboveSatelliteHoriz>0){
	 limit=.85;
       }
       
       //do ctr for  last cuts////

       if(peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	 if(ratioFirstToSecondPeak<=(1/limit)){
	     ratio_last++;//ratio last
	   }
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	  if(peakVal<0.075){
	    peakVal_last++;//peakVal last
	   }
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075  && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	  if(peakHilbertCoherent <=15){
	    hilbert_last++;//hilbert last
	   }
       }
     if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && polFractionCoherent>0.3){
	 if(peakHilbertCoherent<-350*peakVal+57.14){
	   rotated_last++; //rotated last
	 }
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14)){
	 if(polFractionCoherent<=.3){
	   polfrac_last++;//rotated last	
       }
      
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	  if(thetaMap<=-35 || thetaMap>=0){
	    elevation_last++;//elevation last
       }
      
       }

       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	 if(eventTracedFlag2!=1){
	   traced_last++;//traced Last
	 }
       }

       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	 if(hwTriggerFlag==0){
	   //cout<<"eventnumber being cut by hwtrigger is "<<eventNumber<<" hwTriggerFlag is "<<hwTriggerFlag<<"\n";
	   triggered_last++;//triggered_last
	 }
       }
       //cout<<"eventNumber is "<<eventNumber<<"\n";
       //cout<<"lat,lon icemc are "<<lat_icemc<<" "<<lon_icemc<<" analysis source is "<<sourceLat<<" "<<sourceLon<<"\n";
       //cout<<"Arrival times are ";
       // for(int i =0;i<48;i++){
       // cout<<arrival_times[i]<<"\n";
       //}
      
       //cout<<"ratio,peakVal,hilbert are "<<ratioFirstToSecondPeak<<" "<<peakVal<<" "<<peakHilbertCoherent<<"\n";

       /////////////////////////////////////Stuff for plotting and in order cuts
       if(ratioFirstToSecondPeak <=(1/limit)){
	 //cout<<"eventNumber cut by ratio is "<<eventnumber2<<" and CW for this event was "<<CW<<" m is "<<m<<"\n";
	 //cout<<"event "<<eventNumber_cut<<" was cut by ratio!, SNR is "<<snrPeak<<" second peak is "<<peakSNR_2<<"\n";
	 if(snrPeak>10){
	   //cout<<"event "<<eventNumber_cut<<" was cut by ratio! \n";
	   //cout<<"ratioFirstoSecondPeak is "<<ratioFirstToSecondPeak<<"\n";
	 }
	 if(index < 1000) {
	   //cout<<"cut by ratio,index is "<<index<<" event is "<<eventNumber<<" weight is "<<weight<<" ratio is "<<ratioFirstToSecondPeak<<" limit is "<<1./limit<<" ("<<limit<<")\n";
	   //cout<<"lat,lon icemc are "<<lat_icemc<<" "<<lon_icemc<<" analysis source is "<<sourceLat<<" "<<sourceLon<<"\n";
	 }
	 // cout<<"event being cut is "<<eventNumber<<" SNRpeak, peak2peak, RMSnoise are "<<snrPeak<<" "<<peak2peak<<" "<<noiseBeforeFilter<<"\n";
	 ratiopeaksctr++;
	 continue;
       }
       
       ratiopeaks_vec[index]=ratiopeaks_vec[index]+=weight;
       
       // if(peakVal<.075){
	 if(peakVal<.03){
	 crosscorrctr++;
	 //cout<<"event being cut by crosscorr is "<<eventNumber<<" SNRpeak, peak2peak, RMSnoise are "<<snrPeak<<" "<<peak2peak<<" "<<noiseBeforeFilter<<" theta and phi are "<<thetaMap<<" "<<phiMap<<"\n";
	
	 continue;
       }
      	
       peakval_vec[index]=peakval_vec[index]+=weight;
       //cout<<"passed both cuts, event is "<<eventNumber<<"\n";
       //cout<<"icemc location is "<<lon_icemc<<" "<<lat_icemc<<" analysis loc is "<<sourceLon<<" "<<sourceLat<<"\n";
       if(peakHilbertCoherent <=15){
	 hilbertctr++;
	 continue;
       }
       	     
       hilbert_vec[index]=hilbert_vec[index]+=weight;
       /*
        if(polFractionCoherent<=.3){
	 polfractionctr++;
	 continue;
       }
      
       polfrac_vec[index]=polfrac_vec[index]+=weight;
       */
       
       // if(peakHilbertCoherent<-350*peakVal+57.14){
       if(peakHilbertCoherent<-100*peakVal+10){
       // if(peakHilbertCoherent<-117*peakVal+14){
	 rotatedctr++;
	 continue;
       }

      
       rotated_vec[index]=rotated_vec[index]+=weight;
       
       if(thetaMap<=-35 || thetaMap>=0){
	 elevationctr++;
	 continue;
       }
      
       elevation_vec[index]=elevation_vec[index]+=weight;
       
        if(eventTracedFlag2!=1){
	 tracedctr++;
	 //cout<<"elevatin of event cut is "<<thetaMap<<"\n";
	 helevation_traced->Fill(thetaMap,weight);
	 continue;
       }
      
       projected_vec[index]=projected_vec[index]+=weight;
      
       if(hwTriggerFlag==0){
	 triggeredctr++;
	 continue;
       }
      
       triggered_vec[index]=triggered_vec[index]+=weight;
       
        if(varnerFlag==1){// || varnerFlag2==1){
	 varnerevents++;
	 continue;
       }

       final_vec[index]=final_vec[index]+=weight;
      
       datactr++;
       // cout<<"eventnumber cut is "<<eventnumber2<<"\n";
      
     }//nevents0



    double binner;
    int j3=0;

    for(int j1 =0;j1<numbins;j1+=index_mult){
      for(int j2=0;j2<index_mult;j2++){
	//cout<<"j1 is "<<j1<<" adding bin "<<j1+j2<<" to "<<j3<<"\n";
	ratiopeaks_vec_sum[j3]+=ratiopeaks_vec[j1+j2];

	peakval_vec_sum[j3]+=peakval_vec[j1+j2];

	hilbert_vec_sum[j3]+=hilbert_vec[j1+j2];

	polfrac_vec_sum[j3]+=polfrac_vec[j1+j2];

	rotated_vec_sum[j3]+=rotated_vec[j1+j2];

	elevation_vec_sum[j3]+=elevation_vec[j1+j2];

	projected_vec_sum[j3]+=projected_vec[j1+j2];

	triggered_vec_sum[j3]+=triggered_vec[j1+j2];
	
	final_vec_sum[j3]+=triggered_vec[j1+j2];

	volt_tracker_sum[j3]+=volt_tracker[j1+j2];

	
      }//j2
      if(volt_tracker_sum[j3]==0.){
	  volt_tracker_sum[j3]=1;
      }
      j3++;
     
    }//j1
    // cout<<"j3 is "<<j3<<"\n";
    int j_mult;
   for(int j=0;j<numbins/index_mult;j++){
     j_mult = j*index_mult;
     // binner = (double) (j/index_mult);
     binner = (double) (j)+1.;
     /*   cout<<"ratiopeaks_vec_sum["<<j<<"] is "<<ratiopeaks_vec_sum[j]<<"\n";
     cout<<"ratiopeaks_vec_["<<j<<"] is "<<ratiopeaks_vec[j_mult]<<" "<<ratiopeaks_vec[j_mult+1]<<" "<<ratiopeaks_vec[j_mult+2]<<"\n";
     cout<<"volt_tracker_sum[j] is "<<volt_tracker_sum[j]<<"\n";
     cout<<"/////////////////////// \n";*/
     cout<<"j is "<<j<<" bin is "<<binner<<" number of events in bin is "<<volt_tracker_sum[j]<<"\n";
     cout<<"number passing cuts are "<<ratiopeaks_vec_sum[j]<<" "<<peakval_vec_sum[j]<<" "<<polfrac_vec_sum[j]<<" "<<rotated_vec_sum[j]<<" "<<elevation_vec_sum[j]<<" "<<projected_vec_sum[j]<<" "<<triggered_vec_sum[j]<<" "<<final_vec_sum[j]<<"\n";
     cout<<"ratiopeaks_vec_sum[j]/volt_tracker_sum[j] is "<<ratiopeaks_vec_sum[j]/volt_tracker_sum[j]<<"\n";
     hratiopeaks->SetBinContent(binner,ratiopeaks_vec_sum[j]/volt_tracker_sum[j]);
    
     hpeakVal->SetBinContent(binner,peakval_vec_sum[j]/volt_tracker_sum[j]);
     
     hhilbert->SetBinContent(binner,hilbert_vec_sum[j]/volt_tracker_sum[j]);
      
     hpolfrac->SetBinContent(binner,polfrac_vec_sum[j]/volt_tracker_sum[j]);
    
     hrotated->SetBinContent(binner,rotated_vec_sum[j]/volt_tracker_sum[j]);
     
     helevation->SetBinContent(binner,elevation_vec_sum[j]/volt_tracker_sum[j]);
    
     hprojected->SetBinContent(binner,projected_vec_sum[j]/volt_tracker_sum[j]);
     
     htriggered->SetBinContent(binner,triggered_vec_sum[j]/volt_tracker_sum[j]);
    
     hfinal->SetBinContent(binner,final_vec_sum[j]/volt_tracker_sum[j]);
   }
      

   cout<<"for filter "<<filter_name<<"\n";
   cout<<"number of events that got pointed is: "<<eventctr_total2<<"\n";
   cout<<"number of events that got skipped is "<<index_ctr<<"\n";
   cout<<"number of events in data tree is: "<<dataeventsctr<<"\n";
   cout<<"payload blast cut is: "<<payloadBlastctr<<"\n";
 

   cout<<"ratio of peaks cut is: "<<ratiopeaksctr<<" :"<<ratio_last<<" \n";
   cout<<"cross corr cut is: "<< crosscorrctr<<" :"<<peakVal_last<<"\n";
   cout<<"number of events with hilbert <15 is: "<<hilbertctr<<" :"<<hilbert_last<<"\n";
   cout<<"polfraction cut is: "<< polfractionctr<<" :"<<polfrac_last<<"\n";
   cout<<"rotated cross cut is: "<<rotatedctr<<" :"<<rotated_last<<" \n";
   cout<<"elevation cut is: "<<elevationctr<<" :"<<elevation_last<<"\n";
   cout<<"traced cut is: "<<tracedctr<<" :"<<traced_last<<"\n";
   cout<<"triggered cut is: "<<triggeredctr<<" :"<<triggered_last<<"\n";
   cout<<"varner events is: "<<varnerevents<<"\n";
   cout<<"number that pass all cuts is: "<<datactr<<"\n";

   cout<<"peak2peak_max is "<<peak2peak_max<<"\n";
     
char printer[256];
 TH2D  *haxes = new TH2D("axes","SNR peak After Filter; SNR; #frac{# passing cuts}{# passing trigger}",numbins/index_mult,0,16,15,0,1.1);
 TH2D  *haxes_S = new TH2D("axes","SNR peak After Filter; Peak-to-Peak; #frac{# passing cuts}{# passing trigger}",numbins/index_mult,0,numbins+2,15,0,1.1);
TH2D  *haxes_zoomed = new TH2D("axes z",";Cross Correlation Value;Peak Hilbert Envelope",numbins/index_mult,0,.2,100,0,100);
TH2D  *haxes_zoomed_SNR = new TH2D("axes z snr",";Cross Correlation Value;SNR",numbins/index_mult,0,.2,100,0,14);
TH2D  *haxes_zoomed_corrected = new TH2D("axes zc",";Cross Correlation Value;Peak Hilbert Envelope",numbins/index_mult,0,.2,100,0,40E6);

 TCanvas *c1 = new TCanvas("c1","c1",1000,800);
 haxes->Draw();
 //haxes_S->Draw();
 hratiopeaks->Draw("same");
 hpeakVal->Draw("same");
 hhilbert->Draw("same");
 hpolfrac->Draw("same");
 hrotated->Draw("same");
 helevation->Draw("same");
 hprojected->Draw("same");
 htriggered->Draw("same");
 hfinal->Draw("same");
 

 
 TLegend leg(.8,.6,1,1,"");
 leg.SetFillColor(0);
 
 leg.AddEntry(hratiopeaks,"ratio of peaks");
 leg.AddEntry(hpeakVal,"peak val of cross corr");
 leg.AddEntry(hhilbert,"hilbert");
 leg.AddEntry(hpolfrac,"pol frac");
 leg.AddEntry(hrotated,"rotated cross corr");
 leg.AddEntry(helevation,"elevation >-35 && <0");
 leg.AddEntry(hprojected,"traced back to ice");
 leg.AddEntry(htriggered,"triggered");
 leg.AddEntry(hfinal,"final efficiency");
 
 leg.DrawClone("same");

   
 sprintf(printer,"/home/dailey.110/analysis_oindree/efficiency_sim.png");
 c1->Print(printer);
 
 TCanvas *c2 = new TCanvas("c2","c2",800,800);
 helevation_dist->Draw();
 helevation_traced->Draw("same");


TLegend leg1(.7,.9,1,1,"");
 leg1.SetFillColor(0);
 
 leg1.AddEntry(helevation_dist,"theta for reconstruction");
 leg1.AddEntry(helevation_traced,"theta for event missing Antarctica");
 leg1.DrawClone("same");
 sprintf(printer,"/home/dailey.110/analysis_oindree/elevation_dist_sim.png");
 c2->Print(printer);

 TCanvas *c3 = new TCanvas("c3","c3",800,800);
 hcross_corr->Draw();

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_dist_sim.png");
 c3->Print(printer);
 /*
 TCanvas *c4 = new TCanvas("c4","c4",800,800);
 hcross_corr_SNR->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_SNR_sim.png");
 c4->Print(printer);
 */
  TCanvas *c5 = new TCanvas("c5","c5",800,800);
 hcross_corr_cut->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_cut_sim.png");
 c5->Print(printer);

 TCanvas *c6 = new TCanvas("c6","c6",800,800);
 hcross_corr_cut_corrected->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_cut_corrected_sim.png");
 c6->Print(printer);

TCanvas *c6a = new TCanvas("c6a","c6a",800,800);
 haxes_zoomed_corrected->Draw();
 hcross_corr_cut_corrected->Draw("same zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_cut_corrected_zoomed_sim.png");
 c6a->Print(printer);

 TCanvas *c5a = new TCanvas("c5a","c5a",800,800);
 haxes_zoomed->Draw();
 hcross_corr_cut->Draw("same zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_cut_zoomed_sim.png");
 c5a->Print(printer);



TCanvas *c7 = new TCanvas("c7","c7",800,800);
 hcross_corr_LR->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_LR_sim.png");
 c7->Print(printer);

TCanvas *c7a = new TCanvas("c7a","c7a",800,800);
 haxes_zoomed->Draw();
 hcross_corr_LR->Draw("same zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_LR_zoomed_sim.png");
 c7a->Print(printer);

TCanvas *c7b = new TCanvas("c7b","c7b",800,800);
 hcross_corr_LR_corrected->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_LR_corrected_sim.png");
 c7b->Print(printer);

TCanvas *c7c = new TCanvas("c7c","c7c",800,800);
 haxes_zoomed->Draw();
 hcross_corr_LR_corrected->Draw("same zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_LR_zoomed_corrected_sim.png");
 c7c->Print(printer);


TCanvas *c8 = new TCanvas("c8","c8",800,800);
 hdeltaTheta->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/deltaTheta_LR_sim.png");
 c8->Print(printer);

TCanvas *c9 = new TCanvas("c9","c9",800,800);
 hdeltaPhi->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/deltaPhi_LR_sim.png");
 c9->Print(printer);

TCanvas *c11 = new TCanvas("c11","c11",800,800);
 hcross_corr_cut_SNR->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_cut_SNR_sim.png");
 c11->Print(printer);

TCanvas *c11a = new TCanvas("c11a","c11a",800,800);
 haxes_zoomed_SNR->Draw();
 hcross_corr_cut_SNR->Draw("same zcol");

sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_cut_SNR_zoomed_sim.png");
 c11a->Print(printer);

TCanvas *c12 = new TCanvas("c12","c12",800,800);
 hcross_corr_LR_SNR->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_LR_SNR_sim.png");
 c12->Print(printer);

TCanvas *c12a = new TCanvas("c12a","c12a",800,800);
 haxes_zoomed_SNR->Draw();
 hcross_corr_LR_SNR->Draw("same zcol");

sprintf(printer,"/home/dailey.110/analysis_oindree/cross_corr_LR_SNR_zoomed_sim.png");
 c12a->Print(printer);

 TLine *line = new TLine(0,0,1,1);
 line->SetLineWidth(2);
TCanvas *c13 = new TCanvas("c13","c13",800,800);
 hpeakLvsR->Draw("zcol");
 line->Draw("same");
 sprintf(printer,"/home/dailey.110/analysis_oindree/peakLvsR_sim.png");
 c13->Print(printer);

TCanvas *c14 = new TCanvas("c14","c14",800,800);
 hpeaksratio->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/peaksratio_sim.png");
 c14->Print(printer);

TCanvas *c15 = new TCanvas("c15","c15",800,800);
 hcorrratio->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/corrratio_sim.png");
 c15->Print(printer);

TCanvas *c16 = new TCanvas("c16","c16",800,800);
//c16->SetLogy();
 hintegrandvsSNR->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/integrandSNR_sim.png");
 c16->Print(printer);

TCanvas *c17 = new TCanvas("c17","c17",800,800);
//c17->SetLogy();
 hdenomvsSNR->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/denomSNR_sim.png");
 c17->Print(printer);

TCanvas *c18 = new TCanvas("c18","c18",800,800);
// c18->SetLogy();
 hCorrvstime->Draw("zcol");

 sprintf(printer,"/home/dailey.110/analysis_oindree/Corrtime_sim.png");
 c18->Print(printer);



 delete c1;
 delete hpeakVal;
 delete hratiopeaks;
 delete hhilbert;
 delete hpolfrac;
 delete hprojected;
 delete htriggered;
 delete helevation;
 delete hrotated;
 delete hfinal;



  
}//main
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
  RootStyle->SetPadLeftMargin  (0.16);
  RootStyle->SetPadRightMargin (.2);
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
  RootStyle->SetStatX         (0);
  RootStyle->SetStatY         (0);
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
  RootStyle->SetLabelSize  ( 0.030,"X");
  RootStyle->SetLabelFont  ( 42   ,"X");

  RootStyle->SetTickLength ( 0.015,"Y");
  RootStyle->SetTitleSize  ( 0.050,"Y");
  RootStyle->SetTitleOffset( 1.4,"Y");
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

