//tausgSystem->Reset();
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
#include "TNtuple.h"
#include "MyCorrelator.h"
#include "analysis_info_4pol.h"
#include "TMarker.h"

#include "healpix_base.h"
#include "chealpix.h"

using namespace std;

class MyCorrelator;

//void Pointsimdata() {
int k_value = 4;
long n_side =pow(2,k_value);
int n_pix_int=12*n_side*n_side; //Total number of pixels


void LatLon2phitheta(double lat,double lon, double &phi, double &theta);

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

int main() {
  ofstream myfile;
  myfile.open("Abby_output.txt");
  myfile<<"ratio, peakVal, hilbert, thetaMap, PolFrac, Trace, Trigger, Varner,SNR,pixel_num \n";
  TStyle* RootStyle();
  // TStyle *color=RootStyle();

  

  /*  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  gSystem->AddIncludePath("-I${EVENT_SIMULATION_DIR}");
  gSystem->AddIncludePath("-I${ANITA_ANALYSIS}/classes");
  gSystem->Load("libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");  
  gSystem->Load("libRootFftwWrapper.so");     	  
  gSystem->Load("libAnitaEvent.so");   	  
  gSystem->Load("libMyCorrelator.so");
  gSystem->Load("libAnitaCorrelator.so");
  gSystem->Load("libStephenAnalysisTools.so");
  */
  //int ctr_brian=0;
  int xCorPassFlag=0;
  int ratioOfPeaksPassFlag=0;
  int elevationAnglePassFlag=0;
  int peakCrossCorrFlag=0;
  int polFractionFlag=0;
  int peakHilbertFlag=0;
  int triggerflag=0;
  
  //int firstEvent;
  //int lastEvent;
  vector<int> vtimes_regular;
  vector<int> vtimes_regular1;
  vector<double> time_full;
  vector<double> time_spectrum;

  vector<int> eventNumber_vector;
  vector<int> eventNumber_tracker;
    long pixel_num;
  TMarker *eventPos[18];
  
 
    int run;
    int eventNumber;
    int whichPol=0;
     analysis_info_4pol* pol4_Ptr = 0;
     TH2D *hevents = new TH2D("events",";Peak Value;Peak Hilbder",1000,0,1,10000,0,500);
     
    int eventNumber_array[18] =  {4586631,8381355,9362397,10345208,11207106,12943715,13662401,15406954,19480638,20564174,24699044,25887362,15636066,16014510,14496361,14577785,21684774,27146983};
    int run_array[18] = {43,75,86,99,110,127,134,152,177,190,239,250,153,155,142,142,208,261};
    for (int i =0;i<18;i++){
      //  int i =0;
      eventNumber=eventNumber_array[i];
      run = run_array[i];
      whichPol=0;
      if (i > 13){
	whichPol=1;
      }
	 
      pol4_Ptr = new analysis_info_4pol();
      
      //MyCorrelator *magicPtr = new MyCorrelator("/u/home/agoodhue/icemcEventMaker/runs",run, WaveCalType::kDefault);
      MyCorrelator *magicPtr = new MyCorrelator("/data/anita/anitaplus/AnitaFlight0809RootData",run, WaveCalType::kDefault);
      //MyCorrelator *magicPtr = new MyCorrelator("/data/anita/btdailey/anita_data/sample10",run, WaveCalType::kDefault);
      // MyCorrelator *magicPtr = new MyCorrelator("/home/btdailey/icemc_svn/outputs",run, WaveCalType::kDefault);
      //MyCorrelator *magicPtr = new MyCorrelator("/u/osgstorage/anita/data/mcm08",run, WaveCalType::kDefault);
      cout<<"RUN NUMBER: "<<run<<endl;
      
      
      // pol4_Ptr = new analysis_info_4pol();
      //double deltaTheta, deltaPhi;
      double thetaMap,phiMap;
      double theta,phi;
      // double sourceLon, sourceLat, sourceAlt;
      double finaltheta;
      
      
      cout<<"eventNumber, run and pol are "<<eventNumber<<","<<run<<","<<whichPol<<"\n";
      
      magicPtr->pointThisEvent(eventNumber, 0, pol4_Ptr);
      
      //inputs are: event number, a flag that is 1 if you want to draw maps and 0 if you don't, first ntuple, 
      //                                       second ntuple, 
      //                                       third ntuple, fourth ntuple, where it pointed in theta, where it pointed in phi, 
      //                                       and polarization (0 is vertical pol, 1 h)
      cout<<"passed point this event \n\n\n";
      
      
      cout<<"ratio, peakVal, hilbert, thetaMap, PolFrac, Trace, Trigger, Varner are "<<pol4_Ptr->ratioFirstToSecondPeak[whichPol]<<" "<<pol4_Ptr->peakVal[whichPol]<<" "<<pol4_Ptr->peakHilbertCoherent[whichPol]<<" "<<pol4_Ptr->thetaMap[whichPol]<<" "<<pol4_Ptr->polFractionCoherent[whichPol]<<" "<<pol4_Ptr->eventTracedFlag[whichPol]<<" "<<pol4_Ptr->hwTriggerFlag[whichPol]<<" "<<pol4_Ptr->varnerFlag[whichPol]<<"\n";
      cout<<"peakVal, peakHilbert is "<<pol4_Ptr->peakVal[whichPol]<<" "<<pol4_Ptr->peakHilbertCoherent[whichPol]<<"\n";

      cout<<"lat, lon is "<<pol4_Ptr->sourceLat[whichPol]<<" "<<pol4_Ptr->sourceLon[whichPol]<<"\n";
      LatLon2phitheta(90-pol4_Ptr->sourceLat[whichPol],180+pol4_Ptr->sourceLon[whichPol],phi,theta);//90-Lats,180+Lons
	      
      pixel_num=ang2pix_ring(n_side,theta,phi);

      cout<<"pixel_num is "<<pixel_num<<"\n";
     
	myfile<<eventNumber<<"\t"<<pol4_Ptr->ratioFirstToSecondPeak[whichPol]<<"\t"<<pol4_Ptr->peakVal[whichPol]<<"\t"<<pol4_Ptr->peakHilbertCoherent[whichPol]<<"\t"<<pol4_Ptr->thetaMap[whichPol]<<"\t"<<pol4_Ptr->polFractionCoherent[whichPol]<<"\t"<<pol4_Ptr->eventTracedFlag[whichPol]<<"\t"<<pol4_Ptr->hwTriggerFlag[whichPol]<<"\t"<<pol4_Ptr->varnerFlag[whichPol]<<"\t"<<pol4_Ptr->SNR_coherent[whichPol]<<"\t"<<pixel_num<<"\n";

      hevents->Fill(pol4_Ptr->peakVal[whichPol],pol4_Ptr->peakHilbertCoherent[whichPol]);
      eventPos[i] = new TMarker(pol4_Ptr->peakVal[whichPol],pol4_Ptr->peakHilbertCoherent[whichPol],5);
      eventPos[i]->SetMarkerStyle(20);
      if(i ==12 || i==13){
	eventPos[i]->SetMarkerColor(kViolet-6);
      }
      if(whichPol==1){
	eventPos[i]->SetMarkerColor(kBlue);
      }
      delete pol4_Ptr;
      delete magicPtr;
				}//for loop
     // hevents->SetMarkerSize(5);
       TLine *line1 = new TLine(0.,57,.163,0.);
       
       line1->SetLineColor(kRed);
       line1->SetLineWidth(2);


       TCanvas *c0 = new TCanvas("c0","c0",880,800);
       hevents->Draw();
       for(int i=0;i<18;i++){
	 eventPos[i]->Draw("same");
       }
       line1->Draw("same");
       
       
       c0->Print("Abby_comparison_events.png");
       myfile.close();
}//main
//}
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
  RootStyle->SetPadLeftMargin  (0.13);
  RootStyle->SetPadRightMargin (.13);
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
  RootStyle->SetLabelSize  ( 0.040,"X");
  RootStyle->SetLabelFont  ( 42   ,"X");

  RootStyle->SetTickLength ( 0.015,"Y");
  RootStyle->SetTitleSize  ( 0.055,"Y");
  RootStyle->SetTitleOffset( 1.00,"Y");
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


