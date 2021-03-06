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
#include "MyCorrelator.h"

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
analysis_info_4pol *simpol4_Ptr = 0;
//#endif
class MyCorrelator;
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
  MyCorrelator *magicPtr = new MyCorrelator("/data/anita/anitaplus/AnitaFlight0809RootData",12, WaveCalType::kDefault);
   gStyle=color;
   gStyle->SetPalette(1);
   gStyle -> SetOptStat(0000000000000);

   int index;
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
  
   float nadirFlag;
   
   double AnitaLat;
   double AnitaLon;
   double AnitaAlt;
   double heading;

   double simAnitaLat;
   double simAnitaLon;
   double simAnitaAlt;
   double simheading;

   int badNoiseFlag;

   
   float didIFilter;
   float didIFilterHoriz;
   float didIFilterAboveSatelliteHoriz;
   float didIFilterAboveSatellite;
    
   int payloadBlastctr=0;
   int hilbertctr=0;
   int ratiopeaksctr=0;
   int crosscorrctr=0;
   int polfractionctr=0;
   int rotatedctr=0;
   int elevationctr=0;
   int triggeredctr=0;
   int tracedctr=0;
   int badNoisectr=0;
   
   int varnerevents=0;
   int datactr=0;

     int ratio_last=0;
   int peakVal_last=0;
   int hilbert_last=0;
   int polfrac_last=0;
   int rotated_last=0;
   int elevation_last=0;
   int traced_last=0;
   int triggered_last=0;

   double snrPeak;

   double sourceLon;
   double sourceLat;
   double sourceAlt;
   int realtime;
   int false_neg_ctr=0;
   TH1D *hpeakVal;
   TH1D *hratiopeaks;
   TH1D *hhilbert;
   TH1D *hpolfrac;
   TH1D *helevation;
   TH1D *hrealTime;
   TH1D *hrealTime_sec;
  int numbins = 50;
   int index_mult = 1;
   // numbins=50;
   // index_mult=1;
   //create histograms that are filled at the end of the code
   hratiopeaks = new TH1D("ratiopeaks","; Ratio 2nd/1st Peak;Number of Events",100,0,1);
   hpeakVal = new TH1D("peakVal",";Peak Val of Cross Correlation;Number of Events",100,0,1);
   hhilbert = new TH1D("hilbert",";Peak of Coherently Summed Waveform;Number of Events",250,0,500);
   hpolfrac = new TH1D("polfrac","; Polarization Fraction;Number of Events",100,0,1);
   helevation = new TH1D("elevation",";Elevation Angle;Number of Events",180,-90,90);
   hrealTime = new TH1D("realTime",";RealTime of trigger;Number of Events",2497980/5,1229804070,1232350250);
   hrealTime->GetXaxis()->SetTimeDisplay(1);

   hrealTime_sec = new TH1D("realTime_sec",";RealTime of trigger;Number of Events",375,1231453575,1231453950);
   hrealTime_sec->GetXaxis()->SetTimeDisplay(1);
  

   double weight;
   UInt_t eventNumber_icemc;
   int n_icemc=0;
   double lat_icemc;
   double lon_icemc;
  
   double arrival_times[48];
   double peak2peak_max=0.;
  
  
  
  //read in root file
   
    int polarization =0;
   
    string temp = "~/analysis_oindree/Passed_events_90sample.root";
   
  
   char *rootfile = Form(temp.c_str());
  
 
  TChain *pol4_Tree = new TChain("passing_events");
  
  pol4_Tree->Add(rootfile);
  
     int mainrfcmflag;
     int bigenoughpeakflag;
     int dcoffsetflag;
     int shorttraceflag;
     int nadirrfcmflag;
     
     
     int n_old=0;
     int nevents0 = pol4_Tree->GetEntries();
     pol4_Tree->SetBranchAddress("pol4_Ptr",&pol4_Ptr);
     cout<<"nevents0 is "<<nevents0<<"\n";
    
    
     for(int m=0;m<nevents0;m++){
     //for(int m=1369290;m<nevents0;m++){
     //for(int m=0;m<1;m++){
       	
	//get events
      
       // cout<<"here \n";
       pol4_Tree->GetEvent(m);
       //cout<<"pol4_ptr->EventNumber is "<<pol4_Ptr->eventNumber<<"\n";
       
       //cout<<"Anita is at "<<AnitaLat<<" "<<AnitaLon<<" "<<AnitaAlt<<"\n";
       eventNumber_cut = pol4_Ptr->eventNumber;
       for(int pol=0;pol<2;pol++){
	 // cout<<"eventNumber_cut is "<<eventNumber_cut<<"\n";
	 peakHilbertCoherent=pol4_Ptr->peakHilbertCoherent[pol];
	 peakVal=pol4_Ptr->peakVal[pol];
	 ratioFirstToSecondPeak=pol4_Ptr->ratioFirstToSecondPeak[pol];
	 thetaMap=pol4_Ptr->thetaMap[pol];
	 snrcoherent = pol4_Ptr->SNR_coherent[pol];
	 polFractionCoherent=pol4_Ptr->polFractionCoherent[pol];
	 
	 //if(peakHilbertCoherent != peakHilbertCoherent) cout<<"eventnumber is "<<eventNumber_cut<<" snr is "<<snrcoherent<<"\n";
	 hwTriggerFlag=pol4_Ptr->hwTriggerFlag[pol];
	 snrPeak = pol4_Ptr->snrPeakAfterFilter[pol];
	 
	 varnerFlag = pol4_Ptr->varnerFlag[pol];
	 didIFilter = pol4_Ptr->didIFilter[pol];
	 eventTracedFlag2 = pol4_Ptr->eventTracedFlag[pol];
	 varnerFlag2 = pol4_Ptr->varnerFlag2[pol];
	 realtime = pol4_Ptr->realTime;
	 
	 mainrfcmflag=  pol4_Ptr->mainrfcmflag;
	 bigenoughpeakflag=pol4_Ptr->bigenoughpeakflag;
	 dcoffsetflag=pol4_Ptr->dcoffsetflag;
	 shorttraceflag=pol4_Ptr->shorttraceflag;
	 nadirrfcmflag=pol4_Ptr->nadirrfcmflag;
	 badNoiseFlag = pol4_Ptr->noiseFlag[pol];
	 if ( !mainrfcmflag || shorttraceflag !=0 || dcoffsetflag != 0 || bigenoughpeakflag !=1){
	   cout<<"quality cut event \n";
	   continue;
	 }
	 
	 
	 //PAYLOAD BLASTS
	 if(pol4_Ptr->payloadblastflag ==1){
	   payloadBlastctr++;
	   continue;
	 }
	 
	 index = (int) snrPeak;
	 
	 
	 if(index >= numbins){
	   index=numbins;
	   
	 }
	 
	 hpeakVal->Fill(peakVal);
	 hratiopeaks->Fill(1./ratioFirstToSecondPeak);
	 hhilbert->Fill(peakHilbertCoherent);
	 hpolfrac->Fill(polFractionCoherent);
	 helevation->Fill(thetaMap);
	 hrealTime->Fill(realtime);
	 hrealTime_sec->Fill(realtime);
	 cout<<"event is "<<eventNumber_cut<<" pol is "<<pol<<" peakHilbert is "<<peakHilbertCoherent<<" realtime is "<<realtime<<" lat,lon are "<<pol4_Ptr->sourceLat[pol]<<" "<<pol4_Ptr->sourceLon[pol]<<" thetaMap is "<<thetaMap<<"\n";
       }//pol
     }//nevents0

  

     
char printer[256];
 TH2D  *haxes = new TH2D("axes","; SNR; #frac{# passing cuts}{# passing trigger}",numbins/index_mult,0,16,15,0,1.1);
 TH2D  *haxes_S = new TH2D("axes","; Peak-to-Peak; #frac{# passing cuts}{# passing trigger}",numbins/index_mult,0,numbins+2,15,0,1.1);
TH2D  *haxes_zoomed = new TH2D("axes z",";Cross Correlation Value;Peak Hilbert Envelope",numbins/index_mult,0,.2,100,0,100);
TH2D  *haxes_zoomed_SNR = new TH2D("axes z snr",";Cross Correlation Value;SNR",numbins/index_mult,0,.2,100,0,14);
TH2D  *haxes_zoomed_corrected = new TH2D("axes zc",";Cross Correlation Value;Peak Hilbert Envelope",numbins/index_mult,0,.2,100,0,40E6);

 TCanvas *c1 = new TCanvas("c1","c1",880,800);
 c1->SetLogy();
 hratiopeaks->Draw();
 

 sprintf(printer,"RatioOfPeaks_passing.png");
 c1->Print(printer);

 TCanvas *c2 = new TCanvas("c2","c2",880,800);
 c2->SetLogy();
 hpeakVal->Draw();
 
 sprintf(printer,"PeakVal_passing.png");
 c2->Print(printer);

TCanvas *c3 = new TCanvas("c3","c3",880,800);
 c3->SetLogy();
 hhilbert->Draw();


 sprintf(printer,"PeakHilbert_passing.png");
 c3->Print(printer);

 TCanvas *c4 = new TCanvas("c4","c4",880,800);
 c4->SetLogy();
 hpolfrac->Draw();
 
 sprintf(printer,"PolFrac_passing.png");
 c4->Print(printer);

 TCanvas *c5 = new TCanvas("c5","c5",880,800);
 c5->SetLogy();
 helevation->Draw();


 sprintf(printer,"Elevation_passing.png");
 c5->Print(printer);

 TCanvas *c6 = new TCanvas("c6","c6",880,800);
 hrealTime->Draw();
 sprintf(printer,"RealTime_passing.png");
 c6->Print(printer);

  TCanvas *c6a = new TCanvas("c6a","c6a",880,800);
 hrealTime_sec->Draw();
 sprintf(printer,"RealTime_passing_sec.png");
 c6a->Print(printer);
 
 delete hpeakVal;
 delete hratiopeaks;
 delete hhilbert;
 delete hpolfrac;

 delete helevation;
 
 delete hrealTime;
  
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

