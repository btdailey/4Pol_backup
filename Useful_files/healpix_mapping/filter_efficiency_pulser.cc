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
RawAnitaHeader *fHeadPtr = 0;
UsefulAnitaEvent *fUsefulEventPtr;
UsefulAnitaEvent *fUsefulEventPtr_CW;

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

   
   
   float peakHilbertCoherent=0;
   float peakVal=0;
   float ratioFirstToSecondPeak=0;
   float thetaMap=0;
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
   int eventTracedFlag2=0;
   int dataeventsctr=0;
   int varnerevents=0;

   int polarization=0;

   int eventnumber2=0;
   int realtime;

   int eventnumber2_old=0;
   int realtime_old=0;
   double noiseBeforeFilter=0.;
   double snrPeak;

   fHeadPtr = new RawAnitaHeader();
   
   
   int numbins = 20;
   int index_mult = 1;

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
  
  
    
   //create histograms that are filled at the end of the code
   hratiopeaks = new TH1D("ratiopeaks","SNR peak After Filter; SNR; #frac{# passing cuts}{# passing trigger}",numbins/index_mult,0,numbins);
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
   
   string temp = "/data/anita/btdailey/filter_study/pulser/%s/output13_0_0.root";
   int index=0;
   for(int runno=1;runno<=1;runno++){
         
     cout<<"runno is "<<runno<<"\n";
     char *rootfile = Form(temp.c_str(),filter_name.c_str());
       cout<<"rootfile is "<<rootfile<<"\n";
       TChain *dataTree = new TChain("ndata");
       TChain *dataTree2 = new TChain("ndata2");
       TChain *dataTree3 = new TChain("ndata3");
       TChain *dataTree4 = new TChain("ndata4");
       TChain *FlagTree = new TChain("ndataflags");
       TChain *TracingTree = new TChain("tdataTracing");
       TChain *PointedTree = new TChain("tdataPointed");
       
       dataTree->Add(rootfile);
       dataTree2->Add(rootfile);
       dataTree3->Add(rootfile);
       dataTree4->Add(rootfile);
       FlagTree->Add(rootfile);
       TracingTree->Add(rootfile);
       PointedTree->Add(rootfile);
       for(int extra=14;extra<262;extra++){
    
	 sprintf(rootfile,"/data/anita/btdailey/filter_study/pulser/%s/output%d_0_0.root",filter_name.c_str(),extra);
	 TFile *fpTest = TFile::Open(rootfile);
	 if(!fpTest){ 
	   cout<<"broke at extra "<<extra<<"\n";
	   continue;
	   //break;
	 }
	 else {
	   delete fpTest;
	   
	   dataTree->Add(rootfile);
	   dataTree2->Add(rootfile);
	   dataTree3->Add(rootfile);
	   dataTree4->Add(rootfile);
	   FlagTree->Add(rootfile);
	   TracingTree->Add(rootfile);
	   PointedTree->Add(rootfile);
	 }
       }//extra
    
     int nevents0 = dataTree->GetEntries();

     int nevents2 = PointedTree->GetEntries();

   
     dataeventsctr+=nevents0;
     eventctr_total2+=nevents2;
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
     
     TracingTree->SetBranchAddress("eventTracedFlag",&eventTracedFlag2);
     TracingTree->SetBranchAddress("eventNumber",&eventnumber2);
     TracingTree->SetBranchAddress("realTime",&realtime);
     TracingTree->SetBranchAddress("snrPeakAfterFilter",&snrPeak);
     TracingTree->SetBranchAddress("noiseBeforeFilter",&noiseBeforeFilter);

     double timer;
    
     cout<<"nevents0 is : "<<nevents0<<"\n";
    
      double max_volts=-1000;
      double min_volts = 1000;
       
      double peak2peak=0.;
      double SNRforindex=0.;

     for(int m=0;m<nevents0;m++){
      
	//get events
       dataTree->GetEvent(m);
       dataTree3->GetEvent(m);
       dataTree4->GetEvent(m);
       TracingTree->GetEvent(m);
       dataTree2->GetEvent(m);
      
       index = (int)snrPeak;
       
        if(index >= numbins){
	  index=numbins;
	 
	 index_ctr++;
	 //	 cout<<"indices are "<<index<<"\n";
	 // continue;
       }

      
       if(index < numbins){
	 volt_tracker[index] = volt_tracker[index]+1; 
	 hnum_peak_vol->Fill(index);
       }
      

       float limit=0.9;

       
      
       
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
	   triggered_last++;//triggered_last
	 }
       }

      
       /////////////////////////////////////Stuff for plotting and in order cuts
       if(ratioFirstToSecondPeak <=(1/limit)){
	 //cout<<"eventNumber cut by ratio is "<<eventnumber2<<" and CW for this event was "<<CW<<" m is "<<m<<"\n";
	 ratiopeaksctr++;
	 continue;
       }
       
       ratiopeaks_vec[index]=ratiopeaks_vec[index]++;
       
       if(peakVal<.075){
	 crosscorrctr++;
	 continue;
       }
      	
       peakval_vec[index]=peakval_vec[index]++;
      
       if(peakHilbertCoherent <=15){
	 hilbertctr++;
	 continue;
       }
       	     
       hilbert_vec[index]=hilbert_vec[index]++;
      
        if(polFractionCoherent<=.3){
	 polfractionctr++;
	 continue;
       }
      
       polfrac_vec[index]=polfrac_vec[index]++;
       
       
       if(peakHilbertCoherent<-350*peakVal+57.14){
	 rotatedctr++;
	 continue;
       }

      
       rotated_vec[index]=rotated_vec[index]++;
       
       if(thetaMap<=-35 || thetaMap>=0){
	 elevationctr++;
	 continue;
       }
      
       elevation_vec[index]=elevation_vec[index]++;
       
        if(eventTracedFlag2!=1){
	 tracedctr++;
	 continue;
       }
      
       projected_vec[index]=projected_vec[index]++;
      
       if(hwTriggerFlag==0){
	 triggeredctr++;
	 continue;
       }
      
       triggered_vec[index]=triggered_vec[index]++;
       
        if(varnerFlag==1){// || varnerFlag2==1){
	 varnerevents++;
	 continue;
       }

       final_vec[index]=final_vec[index]++;

       datactr++;
       // cout<<"eventnumber cut is "<<eventnumber2<<"\n";
      
     }//nevents0
   }//runno


    double binner;
    int j3=0;

    for(int j1 =0;j1<numbins;j1+=index_mult){
      for(int j2=0;j2<index_mult;j2++){
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
     binner = (double) (j_mult);
     /*   cout<<"ratiopeaks_vec_sum["<<j<<"] is "<<ratiopeaks_vec_sum[j]<<"\n";
     cout<<"ratiopeaks_vec_["<<j<<"] is "<<ratiopeaks_vec[j_mult]<<" "<<ratiopeaks_vec[j_mult+1]<<" "<<ratiopeaks_vec[j_mult+2]<<"\n";
     cout<<"volt_tracker_sum[j] is "<<volt_tracker_sum[j]<<"\n";
     cout<<"/////////////////////// \n";*/
     hratiopeaks->Fill(binner,ratiopeaks_vec_sum[j]/volt_tracker_sum[j]);
    
     hpeakVal->Fill(binner,peakval_vec_sum[j]/volt_tracker_sum[j]);
     
     hhilbert->Fill(binner,hilbert_vec_sum[j]/volt_tracker_sum[j]);
      
     hpolfrac->Fill(binner,polfrac_vec_sum[j]/volt_tracker_sum[j]);
    
     hrotated->Fill(binner,rotated_vec_sum[j]/volt_tracker_sum[j]);
     
     helevation->Fill(binner,elevation_vec_sum[j]/volt_tracker_sum[j]);
    
     hprojected->Fill(binner,projected_vec_sum[j]/volt_tracker_sum[j]);
     
     htriggered->Fill(binner,triggered_vec_sum[j]/volt_tracker_sum[j]);
    
     hfinal->Fill(binner,final_vec_sum[j]/volt_tracker_sum[j]);
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
     
char printer[256];
 TH2D  *haxes = new TH2D("axes","SNR peak After Filter; SNR; #frac{# passing cuts}{# passing trigger}",numbins/index_mult,0,numbins,15,0.85,1);
 TCanvas *c1 = new TCanvas("c1","c1",1000,800);
 // haxes->Draw();
 hratiopeaks->Draw();//"same");
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

   
 sprintf(printer,"/home/btdailey/analysis/filter_study/pulser/filter_efficiency/efficiency_%s_pulser.png",filter_name.c_str());
 c1->Print(printer);
 

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
  RootStyle->SetLabelSize  ( 0.040,"X");
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

