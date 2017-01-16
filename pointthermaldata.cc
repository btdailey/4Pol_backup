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
using namespace std;

class MyCorrelator;

//void Pointsimdata() {

int main() {
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

  /*double freq;
  double spec;
  double freq_cos;
  double freq_sin;
  double time;
  double time_spec;

  int times;
  int timens;
  //int eventNumber;
  int eventNumber_1;

  int startrun=12;//start of run 
  int endrun=224;//end of run 
  int basenumber=274;//base number 
  */
  // for (int run=startrun;run<=endrun;run++){
    int run=1;
    int eventNumber;

    //eventNumber = 604434;//CW1
    

    // eventNumber = 544001;
    //run =12;

    //eventNumber = 597021;
    //run =13;
    //eventNumber = 1518201;
    //run = 18;

    //eventNumber = 1301596;
    //run = 17;

    //eventNumber = 9885731;
    eventNumber = 9885751;
    run=92;
    int whichPol=0;


  //MyCorrelator *magicPtr = new MyCorrelator("/u/home/agoodhue/icemcEventMaker/runs",run, WaveCalType::kDefault);
  //  MyCorrelator *magicPtr = new MyCorrelator("/data/anita/btdailey/",run, WaveCalType::kDefault);
    MyCorrelator *magicPtr = new MyCorrelator("/data/anita/anitaplus/AnitaFlight0809RootData",run, WaveCalType::kDefault);
    // MyCorrelator *magicPtr = new MyCorrelator("/home/btdailey/icemc_svn/outputs",run, WaveCalType::kDefault);
  //MyCorrelator *magicPtr = new MyCorrelator("/u/osgstorage/anita/data/mcm08",run, WaveCalType::kDefault);
  cout<<"RUN NUMBER: "<<run<<endl;
   analysis_info_4pol* pol4_Ptr = 0;

  pol4_Ptr = new analysis_info_4pol();
  //double deltaTheta, deltaPhi;
  double thetaMap,phiMap;
  
  // double sourceLon, sourceLat, sourceAlt;
  double finaltheta;
   TNtuple *ndata = new TNtuple("ndata","stuff to plot",
         "eventNumber:deltaTheta:deltaPhi:deltamcmTheta:deltamcmPhi:thetaMap:phiMap:mapSNR:peakVal:ratioFirstToSecondPeak:snrCoherent:snrPeakAnt:maxSignalPeak:distanceTD:peakHilbertCoherent");
   TNtuple *ndata2 = new TNtuple("ndata2","stuff to plot 2","eventNumber:deltaTTD");
   TNtuple *ndata3 = new TNtuple("ndata2","stuff to plot 2","eventNumber");
   TNtuple *ndata4 = new TNtuple("ndata2","stuff to plot 2","eventNumber");

   double sourceLat =-80.4181;
   double sourceLon =149.857; 
   double sourceAlt =0.;
   std::string baseNames[1000];
   string baseName;


   //magicPtr->loopOverEvents(9648,9652,0,0,0,0);//neutrino search ish
   // magicPtr->processEventsFromAList(0, 0, 3, 0);
   // magicPtr->loopOverEvents(67455,67467,0,0,0,0);
   
   // magicPtr->loopOverEvents(0,1000000,0,0,0,0,ctr_brian);//do taylor dome events 1000000
   
   //magicPtr->loopOverEvents(57753,57757,0,0,0,0);//do thermal noise upward sample events
   //                                          inputs are:
   //                                          entry you want to start on in the tree, entry you want to end on in the tree, do
   //                                          you want to draw maps (0 is no), taylor dome flag (1 is do taylor dome events, 0 is don't)
   //                                          thermal noise flag, (3 is upward pointing thermal no25307044ise and 0 is don't do thermal noise),
   //                                          and polarization flag (0 is v, 1 is h)
   
   // magicPtr->pointThisEvent(eventNumber, 1, ndata, ndata2, ndata3, ndata4, thetaMap, phiMap,0,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerflag,finaltheta);

   // magicPtr->pointThisEvent(eventNumber, 1, ndata, ndata2, ndata3, ndata4, thetaMap, phiMap,whichPol,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerflag, finaltheta);
   magicPtr->pointThisEvent(eventNumber, 1, pol4_Ptr);
    //inputs are: event number, a flag that is 1 if you want to draw maps and 0 if you don't, first ntuple, 
   //                                       second ntuple, 
   //                                       third ntuple, fourth ntuple, where it pointed in theta, where it pointed in phi, 
   //                                       and polarization (0 is vertical pol, 1 h)

   /*
   int eventTracedFlag;
   cout<<"before thetaMap,phiMap,sourceLon,sourceLat,sourceAlt are "<<thetaMap<<" "<<phiMap<<" "<<sourceLon<<" "<<sourceLat<<" "<<sourceAlt<<"\n";
   eventTracedFlag=magicPtr->traceBackToContinent(eventNumber, thetaMap,phiMap, sourceLon, sourceLat, sourceAlt);
   cout<<"eventTracedFlag is "<<eventTracedFlag<<"\n";
   cout<<"after thetaMap,phiMap,sourceLon,sourceLat,sourceAlt are "<<thetaMap<<" "<<phiMap<<" "<<sourceLon<<" "<<sourceLat<<" "<<sourceAlt<<"\n";
   if (eventTracedFlag==1 || eventTracedFlag==2) magicPtr->checkIfNearAnyBase(eventNumber,thetaMap,phiMap, sourceLon, sourceLat, sourceAlt, baseName);
   cout<<"basename is "<<baseName<<"\n";
   */


   
   delete magicPtr;
   // }//for loop
   
  
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


