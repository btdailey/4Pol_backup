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

    int payloadBlastctr=0;
   int hilbertctr[2]=0;
   int ratiopeaksctr[2]=0;
   int crosscorrctr[2]=0;
   int polfractionctr[2]=0;
   int rotatedctr[2]=0;
   int elevationctr[2]=0;
   int triggeredctr[2]=0;
   int tracedctr[2]=0;
   int badNoisectr[2]=0;
   
   int varnerevents[2]=0;
   int datactr[2]=0;

   int ratio_last[2]=0;
   int peakVal_last[2]=0;
   int hilbert_last[2]=0;
   int polfrac_last[2]=0;
   int rotated_last[2]=0;
   int elevation_last[2]=0;
   int traced_last[2]=0;
   int triggered_last[2]=0;

   int ctrboth=0;
   
  //read in root file
   
  //read in root file
   
     run=12;
   temp = "/data/anita/btdailey/final_filter/10sample/geom_4pol_partial_0301/output%d_0.root";
  
  
   rootfile = Form(temp.c_str(),run);
  
 
  TChain *pol4_Tree = new TChain("analysis_info_4pol");
  TChain *Pointed_Tree =  new TChain("tdataTracing");
  pol4_Tree->Add(rootfile);
  Pointed_Tree->Add(rootfile);
  //extras
    for(int extra=13;extra<263;extra++){
	 
	 for(int delta=0;delta<200000;delta+=3000){
	   
	   sprintf(rootfile,"/data/anita/btdailey/final_filter/10sample/geom_4pol_partial_0301/output%d_%d.root",extra,delta);
	   
	 TFile *fpTest = TFile::Open(rootfile);
	 if(!fpTest){ 
	   cout<<"broke at extra "<<extra<<"_"<<delta<<"\n";
	   break;
	   //break;
	 }
	 else {
	   delete fpTest;
	   
	   
	   
	   pol4_Tree->Add(rootfile);
	   Pointed_Tree->Add(rootfile);
	 }
	  }//delta
       }//extra
  
   
      //////THESE TREES ARE FILLED FOR EVENTS PASSING ALL CUTS!////////////

    int eventNumber;
    double sourceLon, sourceLat, sourceHeight, thetaMap, phiMap, snrPeakAfterFilter, heading,
      altitude, latitude, longitude;
    float varnerFlag, hwTriggerFlag, triggerOrPhiMaskFlag;
    float varnerFlag2;
    int realTime, gpsBadFlag;
    float ratioFirstToSecondPeak;
    float polAngleCoherent;
    float peakVal;
    float peakHilbertCoherent;
    float polFractionCoherent;
    float didIFilter;
    float didIFilterHoriz;



    char filenameBig[150];
    sprintf(filenameBig,"/data/anita/btdailey/rootOutputs/allPointedEventsHandV_1.root");
    cout<<"outputting to file: "<<filenameBig<<endl;
    TFile *rootfileWrite = new TFile(filenameBig,"RECREATE");
    TTree *tdataPointedBig=new TTree("tdataPointedBig","all pointed events");
    tdataPointedBig->Branch("eventNumber",&eventNumber,"eventNumber/I");
    tdataPointedBig->Branch("sourceLon",&sourceLon,"sourceLon/D");
    tdataPointedBig->Branch("sourceLat",&sourceLat,"sourceLat/D");
    tdataPointedBig->Branch("sourceHeight",&sourceHeight,"sourceHeight/D");
    tdataPointedBig->Branch("heading",&heading,"heading/D");
    tdataPointedBig->Branch("latitude",&latitude,"latitude/D");
    tdataPointedBig->Branch("longitude",&longitude,"longitude/D");
    tdataPointedBig->Branch("altitude",&altitude,"altitude/D");
    tdataPointedBig->Branch("thetaMap",&thetaMap,"thetaMap/D");
    tdataPointedBig->Branch("phiMap",&phiMap,"phiMap/D");
    tdataPointedBig->Branch("snrPeakAfterFilter",&snrPeakAfterFilter,"snrPeakAfterFilter/D");
    tdataPointedBig->Branch("varnerFlag",&varnerFlag,"varnerFlag/F");
    tdataPointedBig->Branch("hwTriggerFlag",&hwTriggerFlag,"hwTriggerFlag/F");
    tdataPointedBig->Branch("triggerOrPhiMaskFlag",&triggerOrPhiMaskFlag,"triggerOrPhiMaskFlag/F");
    tdataPointedBig->Branch("varnerFlag2",&varnerFlag2,"varnerFlag2/F");
    tdataPointedBig->Branch("realTime",&realTime,"realTime/I");
    tdataPointedBig->Branch("gpsBadFlag",&gpsBadFlag,"gpsBadFlag/I");
    tdataPointedBig->Branch("ratioFirstToSecondPeak",&ratioFirstToSecondPeak,"ratioFirstToSecondPeak/F");
    tdataPointedBig->Branch("polAngleCoherent",&polAngleCoherent,"polAngleCoherent/F");
    tdataPointedBig->Branch("whichPolDidIUse",&whichPolDidIUse,"whichPolDidIUse/I");
    tdataPointedBig->Branch("peakVal",&peakVal,"peakVal/F");
    tdataPointedBig->Branch("peakHilbertCoherent",&peakHilbertCoherent,"peakHilbertCoherent/F");
    tdataPointedBig->Branch("thetaMap",&thetaMap,"thetaMap/F");
    tdataPointedBig->Branch("polFractionCoherent",&polFractionCoherent,"polFractionCoherent/F");
    tdataPointedBig->Branch("didIFilter",&didIFilter,"didIFilter/F");
    tdataPointedBig->Branch("didIFilterHoriz",&didIFilterHoriz,"didIFilterHoriz/F");
    

     int nevents0 = pol4_Tree->GetEntries();//dataTree->GetEntries();
     
     int nevents_pointed = Pointed_Tree->GetEntries();
    
     pol4_Tree->SetBranchAddress("pol4_Ptr",&pol4_Ptr);
    
     
     int whichPolDidIUse=0;
     
     // int polarization =0;
     for(int m=0;m<nevents0;m++){
     
       	if (m % (nevents0 / 100) == 0){
	  cout << m << " events. " <<(double(m)/double(nevents0)) * 100 << "% complete.\n";
	}
	//get events
      
	
      
       pol4_Tree->GetEvent(m);
       Pointed_Tree->GetEvent(m);
       


       //cout<<"Anita is at "<<AnitaLat<<" "<<AnitaLon<<" "<<AnitaAlt<<"\n";
       eventNumber_cut = pol4_Ptr->eventNumber;


        mainrfcmflag=  pol4_Ptr->mainrfcmflag;
	 bigenoughpeakflag=pol4_Ptr->bigenoughpeakflag;
	 dcoffsetflag=pol4_Ptr->dcoffsetflag;
	 shorttraceflag=pol4_Ptr->shorttraceflag;
	 nadirrfcmflag=pol4_Ptr->nadirrfcmflag;
	 badNoiseFlag = pol4_Ptr->noiseFlag[polarization];
	 if ( !mainrfcmflag || shorttraceflag !=0 || dcoffsetflag != 0 || bigenoughpeakflag !=1){
	   cout<<"quality cut event \n";
	   continue;
	 }
	 
	 
	 //PAYLOAD BLASTS
	 if(pol4_Ptr->payloadblastflag ==1){
	   payloadBlastctr++;
	   continue;
	 }

       for(int polarization=0;polarization<2;polarization++){
	 peakHilbertCoherent=pol4_Ptr->peakHilbertCoherent[polarization];
	 peakVal=pol4_Ptr->peakVal[polarization];
	 ratioFirstToSecondPeak=pol4_Ptr->ratioFirstToSecondPeak[polarization];
	 thetaMap=pol4_Ptr->thetaMap[polarization];
	 snrcoherent = pol4_Ptr->SNR_coherent[polarization];
	 polFractionCoherent=pol4_Ptr->polFractionCoherent[polarization];
	 
	
	 hwTriggerFlag=pol4_Ptr->hwTriggerFlag[polarization];
	 snrPeak = pol4_Ptr->snrPeakAfterFilter[polarization];
	 
	 varnerFlag = pol4_Ptr->varnerFlag[polarization];
	 didIFilter = pol4_Ptr->didIFilter[polarization];
	 eventTracedFlag2 = pol4_Ptr->eventTracedFlag[polarization];
	 varnerFlag2 = pol4_Ptr->varnerFlag2[polarization];
	 
	 
	
	 
	 
	 
	 float limit=0.9;
	 
	 
	 if (didIFilterHoriz>1 || didIFilter>1){
	   limit=.85;
	 }
	 
	 
	 //do ctr for  last cuts////
	 
	 if(peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	   if(ratioFirstToSecondPeak<=(1/limit)){
	     ratio_last[polarization]++;//ratio last
	   }
	 }
	 if(ratioFirstToSecondPeak > (1/limit) && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	   if(peakVal<0.075){
	     peakVal_last[polarization]++;//peakVal last
	   }
	 }
	 if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075  && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	   if(peakHilbertCoherent <=15){
	     hilbert_last[polarization]++;//hilbert last
	   }
	 }
	 if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && polFractionCoherent>0.3){
	   if(peakHilbertCoherent<-350*peakVal+57.14){
	     rotated_last[polarization]++; //rotated last
	   }
	 }
	 if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14)){
	   if(polFractionCoherent<=.3){
	     polfrac_last[polarization]++;//rotated last	
	   }
	   
	 }
	 if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	   if(thetaMap<=-35 || thetaMap>=0){
	     elevation_last[polarization]++;//elevation last
	   }
	   
	 }
	 
	 if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	   if(eventTracedFlag2!=1){
	     traced_last[polarization]++;//traced Last
	   }
	 }
	 
	 if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	   if(hwTriggerFlag==0){
	     //cout<<"eventnumber being cut by hwtrigger is "<<eventNumber<<" hwTriggerFlag is "<<hwTriggerFlag<<"\n";
	     triggered_last[polarization]++;//triggered_last
	   }
	 }
	 
	 
	 
	 
	 
	 /////////////////////////////////////Stuff for plotting and in order cuts
	 
	 
	 if(ratioFirstToSecondPeak <=(1/limit)){
	   ratiopeaksctr[polarization]++;
	   continue;
	 }
	 
	 
	 if(peakVal<.075){
	   //	 if(peakVal<.03){
	   crosscorrctr[polarization]++;
	   
	   continue;
	 }
	 
	 if(peakHilbertCoherent <=15){
	   hilbertctr[polarization]++;
	   continue;
	 }
	 
	 
	 if(polFractionCoherent<=.3){
	   polfractionctr[polarization]++;
	   continue;
	 }
	 
	 
	 
	 if(peakHilbertCoherent<-350*peakVal+57.14){
	   rotatedctr[polarization]++;
	   continue;
	 }
	 
	 
	 
	 if(thetaMap<=-35 || thetaMap>=0){
	   elevationctr[polarization]++;
	   continue;
	 }
	 
	 
	 if(eventTracedFlag2!=1){
	   tracedctr[polarization]++;
	   continue;
	 }
	 
	 
	 if(hwTriggerFlag==0){
	   triggeredctr[polarization]++;
	   continue;
	 }
	 
	 
	 if(varnerFlag==1){// || varnerFlag2==1){
	   varnerevents[polarization]++;
	   continue;
	 }
	 
	 
	 cout<<"eventnumber is "<<eventNumber_cut<<"\n";
	 
	 
	 datactr[polarization]++;
	 eventTree_out->Fill();
	 extraTree->Fill();
	 // cout<<"eventnumber cut is "<<eventnumber2<<"\n";
       }//polarization


    whichPolDidIUse=-1;
    
    if(peakValV >0 && peakValH > 0){
      ctrBoth++;
      if (peakValV>peakValH) whichPolDidIUse=0; //pick one based on polarization
      else whichPolDidIUse=1; //pick one based on polarization
    }
    else if( peakValV >0 && peakValH <0) whichPolDidIUse=0;
    else if( peakValV <0 && peakValH >0) whichPolDidIUse=1;

    eventNumber=pol4_Ptr->eventNumber;
    sourceLon=pol4_Ptr->sourceLon[whichPolDidIUse];
    sourceLat=pol4_Ptr->sourceLat[whichPolDidIUse];
    sourceHeight=pol4_Ptr->sourceAlt[whichPolDidIUse];
    thetaMap=pol4_Ptr->thetaMap[whichPolDidIUse];
    phiMap=pol4_Ptr->phiMap[whichPolDidIUse];
    heading=pol4_Ptr->heading;
    snrPeakAfterFilter=pol4_Ptr->snrPeakAfterFilter[whichPolDidIUse];
    altitude=pol4_Ptr->AnitaAlt;
    latitude=pol4_Ptr->AnitaLat;
    longitude=pol4_Ptr->AnitaLon;
    varnerFlag=pol4_Ptr->varnerFlag[whichPolDidIUse];
    varnerFlag2=pol4_Ptr->varnerFlag2[whichPolDidIUse];
    gpsBadFlag=0;
    hwTriggerFlag=pol4_Ptr->hwTriggerFlag[whichPolDidIUse];
    triggerOrPhiMaskFlag=pol4_Ptr->triggerOrPhiMaskFlag[whichPolDidIUse];
    realTime=pol4_Ptr->realTime;
    ratioFirstToSecondPeak=pol4_Ptr->ratioFirstToSecondPeak[whichPolDidIUse];
    polAngleCoherent=pol4_Ptr->polAngleCoherent[whichPolDidIUse];
    peakVal=pol4_Ptr->peakVal[whichPolDidIUse];
    peakHilbertCoherent=pol4_Ptr->peakHilbertCoherent[whichPolDidIUse];
    polFractionCoherent=pol4_Ptr->polFractionCoherent[whichPolDidIUse];
    didIFilter=pol4_Ptr->didIFilter[0];
    didIFilterHoriz=pol4_Ptr->didIFilterHoriz[1];


     }//nevents0

     cout<<"payload blast cut is: "<<payloadBlastctr<<"\n";
     for(int polarization=0;polarization<2;polarization++){
       cout<<"Abby Cuts pol: "<<polarization<<" /////////////// \n";
       cout<<"ratio of peaks cut is: "<<ratiopeaksctr[polarization]<<" :"<<ratio_last[polarization]<<" \n";
       cout<<"cross corr cut is: "<< crosscorrctr[polarization]<<" :"<<peakVal_last[polarization]<<"\n";
       cout<<"number of events with hilbert <15 is: "<<hilbertctr[polarization]<<" :"<<hilbert_last[polarization]<<"\n";
       cout<<"polfraction cut is: "<< polfractionctr[polarization]<<" :"<<polfrac_last[polarization]<<"\n";
       cout<<"rotated cross cut is: "<<rotatedctr[polarization]<<" :"<<rotated_last[polarization]<<" \n";
       cout<<"elevation cut is: "<<elevationctr[polarization]<<" :"<<elevation_last[polarization]<<"\n";
       cout<<"traced cut is: "<<tracedctr[polarization]<<" :"<<traced_last[polarization]<<"\n";
       cout<<"triggered cut is: "<<triggeredctr[polarization]<<" :"<<triggered_last[polarization]<<"\n";
       cout<<"varner events is: "<<varnerevents[polarization]<<"\n";
       cout<<"number that pass all cuts is: "<<datactr[polarization]<<"\n";
     }

   ////////////////////////////
   
   
   ////////////////////////////


   rootfile_out = eventTree_out->GetCurrentFile();
   rootfile_out->Write();

   rootfile_sim = eventTree_sim->GetCurrentFile();
   rootfile_sim->Write();




  
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

