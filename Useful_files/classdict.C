//
// File generated by rootcint at Fri Nov 25 16:02:45 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME classdict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "classdict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
      #if !(defined(R__ACCESS_IN_SYMBOL) || defined(R__USE_SHADOW_CLASS))
      typedef ::analysis_info_4pol analysis_info_4pol;
      #else
      class analysis_info_4pol  {
         public:
         //friend XX;
         int eventNumber; //
         double anitaLat; //
         double anitaLon; //
         double anitaAlt; //
         double heading; //
         double deltaTheta[4]; //changes with map, needs 4 pols 
         double deltaPhi[4]; //4
         double deltamcmTheta[4]; //4
         double deltamcmPhi[4]; //4
         double mapSNR[4]; //4
         double peakVal[4]; //4
         double ratioFirstToSecondPeak[4]; //4
         double snrCoherent[4]; //4, changes with antennas you use which is based on map
         double snrPeakAnt[4]; //set really early, only based on vpol right now, but could go in getgraphsthisevent and change it
         double maxSignalPeak[4]; //only based on vpol right now, probably need for all 4 
         double peakHilbertCoherent[4]; //4
         double snrPeakAfterFilter[4]; //4
         double snrPeakAfterFilter2[4]; //4
         int didIFilter[4]; //4
         int triggerOrPhiMaskFlag[4]; //4
         double thetaMap[4]; //
         double phiMap[4]; //
         double secondThetaMap[4]; //4 
         double secondPhiMap[4]; //4 
         double headingOfThisEvent[4]; //4, outputs heading - peakphi of map, so changes with map 
         double tertiaryThetaMap[4]; //4,  CHANGED NAME SO ALSO CHANGE NAME IN POINTTHISEVENT
         double tertiaryPhiMap[4]; //4,    CHANGED NAME SO ALSO CHANGE NAME IN POINTTHISEVENT 
         int varnerFlag[4]; //4, SEE IF CAN THROW OUT    
         int varnerFlag2[4]; //4
         int phiMaskFlag[4]; //4, dependent on peak phi 
         int hwTriggerFlag[4]; //4
         double polAngleCoherent[4]; //4 
         double polFractionCoherent[4]; //4
         int didIFilterAboveSatellite[4]; //4
         double meanFreqVert[4]; //4, do we need all 4? 
         double peakThetaFinal[4]; //not needed 
         double peakPhiFinal[4]; //not needed 
         int eventPointedFlag[4]; //4 
         int eventTracedFlag[4]; //4
         double sourceLat[4]; //4 
         double sourceLon[4]; //4
         double sourceAlt[4]; //4
         double rmsNoiseCoherent[4]; //4 
         double noiseBeforeFilter[4]; //4  //we need to code this more
         double CWheight[4]; //4
         double SNR_ant[4]; //4, add in index 
         double SNR_ant_triggered[4]; //4, add in index 
         double SNR_ant_closest[4]; //4, add in index 
         double SNR_ant_coherent[4]; //4, add in index 
         double PowerCut[4]; //4, add in index, might not keep as array of 40 antennas
         int CoherentAnts[4]; //4, list of antennas 
         double distance_from_source[4]; //4 
         double peak2peak_signal[4]; //4, depend on after filtering which will change for each pol 
         double peakVoltage_2[4]; //
         double peakSNR_2[4]; //
         double peakVal_box[4]; //peakVal in box surrounding peak in Vpol
         double max_correlation[4][40][40]; //
         double RMS_ants[4][40]; //
         double time_delay[4][40][40]; //
         int num_bins_filtered; //
         int mainrfcmflag; //
         int bigenoughpeakflag; //
         int dcoffsetflag; //
         int shorttraceflag; //
         int nadirrfcmflag; //
         int payloadblastflag; //
         double SNR_coherent[4]; //
         double SNR_coherent2[4]; //for powerSNR min value==1
         int SNR_noise_bin[4]; //
      };
      #endif

   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void analysis_info_4pol_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void analysis_info_4pol_Dictionary();
   static void *new_analysis_info_4pol(void *p = 0);
   static void *newArray_analysis_info_4pol(Long_t size, void *p);
   static void delete_analysis_info_4pol(void *p);
   static void deleteArray_analysis_info_4pol(void *p);
   static void destruct_analysis_info_4pol(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::analysis_info_4pol*)
   {
      // Make sure the shadow class has the right sizeof
      R__ASSERT(sizeof(::analysis_info_4pol) == sizeof(::ROOT::Shadow::analysis_info_4pol));
      ::analysis_info_4pol *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::analysis_info_4pol),0);
      static ::ROOT::TGenericClassInfo 
         instance("analysis_info_4pol", "./analysis_info_4pol.h", 5,
                  typeid(::analysis_info_4pol), DefineBehavior(ptr, ptr),
                  &analysis_info_4pol_ShowMembers, &analysis_info_4pol_Dictionary, isa_proxy, 4,
                  sizeof(::analysis_info_4pol) );
      instance.SetNew(&new_analysis_info_4pol);
      instance.SetNewArray(&newArray_analysis_info_4pol);
      instance.SetDelete(&delete_analysis_info_4pol);
      instance.SetDeleteArray(&deleteArray_analysis_info_4pol);
      instance.SetDestructor(&destruct_analysis_info_4pol);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::analysis_info_4pol*)
   {
      return GenerateInitInstanceLocal((::analysis_info_4pol*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::analysis_info_4pol*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void analysis_info_4pol_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const ::analysis_info_4pol*)0x0)->GetClass();
   }

} // end of namespace ROOT

//______________________________________________________________________________
namespace ROOT {
   void analysis_info_4pol_ShowMembers(void *obj, TMemberInspector &R__insp)
   {
      // Inspect the data members of an object of class analysis_info_4pol.
      typedef ::ROOT::Shadow::analysis_info_4pol ShadowClass;
      ShadowClass *sobj = (ShadowClass*)obj;
      if (sobj) { } // Dummy usage just in case there is no datamember.

      TClass *R__cl  = ::ROOT::GenerateInitInstanceLocal((const ::analysis_info_4pol*)0x0)->GetClass();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "eventNumber", &sobj->eventNumber);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "anitaLat", &sobj->anitaLat);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "anitaLon", &sobj->anitaLon);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "anitaAlt", &sobj->anitaAlt);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "heading", &sobj->heading);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "deltaTheta[4]", sobj->deltaTheta);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "deltaPhi[4]", sobj->deltaPhi);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "deltamcmTheta[4]", sobj->deltamcmTheta);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "deltamcmPhi[4]", sobj->deltamcmPhi);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "mapSNR[4]", sobj->mapSNR);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "peakVal[4]", sobj->peakVal);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "ratioFirstToSecondPeak[4]", sobj->ratioFirstToSecondPeak);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "snrCoherent[4]", sobj->snrCoherent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "snrPeakAnt[4]", sobj->snrPeakAnt);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "maxSignalPeak[4]", sobj->maxSignalPeak);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "peakHilbertCoherent[4]", sobj->peakHilbertCoherent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "snrPeakAfterFilter[4]", sobj->snrPeakAfterFilter);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "snrPeakAfterFilter2[4]", sobj->snrPeakAfterFilter2);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "didIFilter[4]", sobj->didIFilter);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "triggerOrPhiMaskFlag[4]", sobj->triggerOrPhiMaskFlag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "thetaMap[4]", sobj->thetaMap);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "phiMap[4]", sobj->phiMap);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "secondThetaMap[4]", sobj->secondThetaMap);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "secondPhiMap[4]", sobj->secondPhiMap);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "headingOfThisEvent[4]", sobj->headingOfThisEvent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "tertiaryThetaMap[4]", sobj->tertiaryThetaMap);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "tertiaryPhiMap[4]", sobj->tertiaryPhiMap);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "varnerFlag[4]", sobj->varnerFlag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "varnerFlag2[4]", sobj->varnerFlag2);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "phiMaskFlag[4]", sobj->phiMaskFlag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "hwTriggerFlag[4]", sobj->hwTriggerFlag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "polAngleCoherent[4]", sobj->polAngleCoherent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "polFractionCoherent[4]", sobj->polFractionCoherent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "didIFilterAboveSatellite[4]", sobj->didIFilterAboveSatellite);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "meanFreqVert[4]", sobj->meanFreqVert);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "peakThetaFinal[4]", sobj->peakThetaFinal);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "peakPhiFinal[4]", sobj->peakPhiFinal);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "eventPointedFlag[4]", sobj->eventPointedFlag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "eventTracedFlag[4]", sobj->eventTracedFlag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "sourceLat[4]", sobj->sourceLat);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "sourceLon[4]", sobj->sourceLon);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "sourceAlt[4]", sobj->sourceAlt);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "rmsNoiseCoherent[4]", sobj->rmsNoiseCoherent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "noiseBeforeFilter[4]", sobj->noiseBeforeFilter);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "CWheight[4]", sobj->CWheight);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "SNR_ant[4]", sobj->SNR_ant);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "SNR_ant_triggered[4]", sobj->SNR_ant_triggered);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "SNR_ant_closest[4]", sobj->SNR_ant_closest);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "SNR_ant_coherent[4]", sobj->SNR_ant_coherent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "PowerCut[4]", sobj->PowerCut);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "CoherentAnts[4]", sobj->CoherentAnts);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "distance_from_source[4]", sobj->distance_from_source);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "peak2peak_signal[4]", sobj->peak2peak_signal);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "peakVoltage_2[4]", sobj->peakVoltage_2);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "peakSNR_2[4]", sobj->peakSNR_2);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "peakVal_box[4]", sobj->peakVal_box);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "max_correlation[4][40][40]", sobj->max_correlation);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "RMS_ants[4][40]", sobj->RMS_ants);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "time_delay[4][40][40]", sobj->time_delay);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "num_bins_filtered", &sobj->num_bins_filtered);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "mainrfcmflag", &sobj->mainrfcmflag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "bigenoughpeakflag", &sobj->bigenoughpeakflag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "dcoffsetflag", &sobj->dcoffsetflag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "shorttraceflag", &sobj->shorttraceflag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "nadirrfcmflag", &sobj->nadirrfcmflag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "payloadblastflag", &sobj->payloadblastflag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "SNR_coherent[4]", sobj->SNR_coherent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "SNR_coherent2[4]", sobj->SNR_coherent2);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "SNR_noise_bin[4]", sobj->SNR_noise_bin);
   }

}

namespace ROOT {
   // Wrappers around operator new
   static void *new_analysis_info_4pol(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::analysis_info_4pol : new ::analysis_info_4pol;
   }
   static void *newArray_analysis_info_4pol(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::analysis_info_4pol[nElements] : new ::analysis_info_4pol[nElements];
   }
   // Wrapper around operator delete
   static void delete_analysis_info_4pol(void *p) {
      delete ((::analysis_info_4pol*)p);
   }
   static void deleteArray_analysis_info_4pol(void *p) {
      delete [] ((::analysis_info_4pol*)p);
   }
   static void destruct_analysis_info_4pol(void *p) {
      typedef ::analysis_info_4pol current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::analysis_info_4pol

/********************************************************
* classdict.C
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableclassdict();

extern "C" void G__set_cpp_environmentclassdict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("analysis_info_4pol.h");
  G__cpp_reset_tagtableclassdict();
}
#include <new>
extern "C" int G__cpp_dllrevclassdict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* analysis_info_4pol */
// automatic default constructor
static int G__classdict_168_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   analysis_info_4pol *p;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new analysis_info_4pol[n];
     } else {
       p = new((void*) gvp) analysis_info_4pol[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new analysis_info_4pol;
     } else {
       p = new((void*) gvp) analysis_info_4pol;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__classdictLN_analysis_info_4pol));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__classdict_168_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   analysis_info_4pol* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new analysis_info_4pol(*(analysis_info_4pol*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__classdictLN_analysis_info_4pol));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef analysis_info_4pol G__Tanalysis_info_4pol;
static int G__classdict_168_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 0
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (analysis_info_4pol*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((analysis_info_4pol*) (soff+(sizeof(analysis_info_4pol)*i)))->~G__Tanalysis_info_4pol();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (analysis_info_4pol*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((analysis_info_4pol*) (soff))->~G__Tanalysis_info_4pol();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__classdict_168_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   analysis_info_4pol* dest = (analysis_info_4pol*) G__getstructoffset();
   *dest = *(analysis_info_4pol*) libp->para[0].ref;
   const analysis_info_4pol& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* analysis_info_4pol */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncclassdict {
 public:
  G__Sizep2memfuncclassdict(): p(&G__Sizep2memfuncclassdict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncclassdict::*p)();
};

size_t G__get_sizep2memfuncclassdict()
{
  G__Sizep2memfuncclassdict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceclassdict() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableclassdict() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__classdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__classdictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__classdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__classdictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__classdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__classdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__classdictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__classdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__classdictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__classdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* analysis_info_4pol */
static void G__setup_memvaranalysis_info_4pol(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__classdictLN_analysis_info_4pol));
   { analysis_info_4pol *p; p=(analysis_info_4pol*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->eventNumber)-(long)(p)),105,0,0,-1,-1,-1,1,"eventNumber=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->anitaLat)-(long)(p)),100,0,0,-1,-1,-1,1,"anitaLat=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->anitaLon)-(long)(p)),100,0,0,-1,-1,-1,1,"anitaLon=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->anitaAlt)-(long)(p)),100,0,0,-1,-1,-1,1,"anitaAlt=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->heading)-(long)(p)),100,0,0,-1,-1,-1,1,"heading=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->deltaTheta)-(long)(p)),100,0,0,-1,-1,-1,1,"deltaTheta[4]=",0,"changes with map, needs 4 pols ");
   G__memvar_setup((void*)((long)(&p->deltaPhi)-(long)(p)),100,0,0,-1,-1,-1,1,"deltaPhi[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->deltamcmTheta)-(long)(p)),100,0,0,-1,-1,-1,1,"deltamcmTheta[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->deltamcmPhi)-(long)(p)),100,0,0,-1,-1,-1,1,"deltamcmPhi[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->mapSNR)-(long)(p)),100,0,0,-1,-1,-1,1,"mapSNR[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->peakVal)-(long)(p)),100,0,0,-1,-1,-1,1,"peakVal[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->ratioFirstToSecondPeak)-(long)(p)),100,0,0,-1,-1,-1,1,"ratioFirstToSecondPeak[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->snrCoherent)-(long)(p)),100,0,0,-1,-1,-1,1,"snrCoherent[4]=",0,"4, changes with antennas you use which is based on map");
   G__memvar_setup((void*)((long)(&p->snrPeakAnt)-(long)(p)),100,0,0,-1,-1,-1,1,"snrPeakAnt[4]=",0,"set really early, only based on vpol right now, but could go in getgraphsthisevent and change it");
   G__memvar_setup((void*)((long)(&p->maxSignalPeak)-(long)(p)),100,0,0,-1,-1,-1,1,"maxSignalPeak[4]=",0,"only based on vpol right now, probably need for all 4 ");
   G__memvar_setup((void*)((long)(&p->peakHilbertCoherent)-(long)(p)),100,0,0,-1,-1,-1,1,"peakHilbertCoherent[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->snrPeakAfterFilter)-(long)(p)),100,0,0,-1,-1,-1,1,"snrPeakAfterFilter[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->snrPeakAfterFilter2)-(long)(p)),100,0,0,-1,-1,-1,1,"snrPeakAfterFilter2[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->didIFilter)-(long)(p)),105,0,0,-1,-1,-1,1,"didIFilter[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->triggerOrPhiMaskFlag)-(long)(p)),105,0,0,-1,-1,-1,1,"triggerOrPhiMaskFlag[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->thetaMap)-(long)(p)),100,0,0,-1,-1,-1,1,"thetaMap[4]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->phiMap)-(long)(p)),100,0,0,-1,-1,-1,1,"phiMap[4]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->secondThetaMap)-(long)(p)),100,0,0,-1,-1,-1,1,"secondThetaMap[4]=",0,"4 ");
   G__memvar_setup((void*)((long)(&p->secondPhiMap)-(long)(p)),100,0,0,-1,-1,-1,1,"secondPhiMap[4]=",0,"4 ");
   G__memvar_setup((void*)((long)(&p->headingOfThisEvent)-(long)(p)),100,0,0,-1,-1,-1,1,"headingOfThisEvent[4]=",0,"4, outputs heading - peakphi of map, so changes with map ");
   G__memvar_setup((void*)((long)(&p->tertiaryThetaMap)-(long)(p)),100,0,0,-1,-1,-1,1,"tertiaryThetaMap[4]=",0,"4,  CHANGED NAME SO ALSO CHANGE NAME IN POINTTHISEVENT");
   G__memvar_setup((void*)((long)(&p->tertiaryPhiMap)-(long)(p)),100,0,0,-1,-1,-1,1,"tertiaryPhiMap[4]=",0,"4,    CHANGED NAME SO ALSO CHANGE NAME IN POINTTHISEVENT ");
   G__memvar_setup((void*)((long)(&p->varnerFlag)-(long)(p)),105,0,0,-1,-1,-1,1,"varnerFlag[4]=",0,"4, SEE IF CAN THROW OUT    ");
   G__memvar_setup((void*)((long)(&p->varnerFlag2)-(long)(p)),105,0,0,-1,-1,-1,1,"varnerFlag2[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->phiMaskFlag)-(long)(p)),105,0,0,-1,-1,-1,1,"phiMaskFlag[4]=",0,"4, dependent on peak phi ");
   G__memvar_setup((void*)((long)(&p->hwTriggerFlag)-(long)(p)),105,0,0,-1,-1,-1,1,"hwTriggerFlag[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->polAngleCoherent)-(long)(p)),100,0,0,-1,-1,-1,1,"polAngleCoherent[4]=",0,"4 ");
   G__memvar_setup((void*)((long)(&p->polFractionCoherent)-(long)(p)),100,0,0,-1,-1,-1,1,"polFractionCoherent[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->didIFilterAboveSatellite)-(long)(p)),105,0,0,-1,-1,-1,1,"didIFilterAboveSatellite[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->meanFreqVert)-(long)(p)),100,0,0,-1,-1,-1,1,"meanFreqVert[4]=",0,"4, do we need all 4? ");
   G__memvar_setup((void*)((long)(&p->peakThetaFinal)-(long)(p)),100,0,0,-1,-1,-1,1,"peakThetaFinal[4]=",0,"not needed ");
   G__memvar_setup((void*)((long)(&p->peakPhiFinal)-(long)(p)),100,0,0,-1,-1,-1,1,"peakPhiFinal[4]=",0,"not needed ");
   G__memvar_setup((void*)((long)(&p->eventPointedFlag)-(long)(p)),105,0,0,-1,-1,-1,1,"eventPointedFlag[4]=",0,"4 ");
   G__memvar_setup((void*)((long)(&p->eventTracedFlag)-(long)(p)),105,0,0,-1,-1,-1,1,"eventTracedFlag[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->sourceLat)-(long)(p)),100,0,0,-1,-1,-1,1,"sourceLat[4]=",0,"4 ");
   G__memvar_setup((void*)((long)(&p->sourceLon)-(long)(p)),100,0,0,-1,-1,-1,1,"sourceLon[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->sourceAlt)-(long)(p)),100,0,0,-1,-1,-1,1,"sourceAlt[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->rmsNoiseCoherent)-(long)(p)),100,0,0,-1,-1,-1,1,"rmsNoiseCoherent[4]=",0,"4 ");
   G__memvar_setup((void*)((long)(&p->noiseBeforeFilter)-(long)(p)),100,0,0,-1,-1,-1,1,"noiseBeforeFilter[4]=",0,"4  //we need to code this more");
   G__memvar_setup((void*)((long)(&p->CWheight)-(long)(p)),100,0,0,-1,-1,-1,1,"CWheight[4]=",0,"4");
   G__memvar_setup((void*)((long)(&p->SNR_ant)-(long)(p)),100,0,0,-1,-1,-1,1,"SNR_ant[4]=",0,"4, add in index ");
   G__memvar_setup((void*)((long)(&p->SNR_ant_triggered)-(long)(p)),100,0,0,-1,-1,-1,1,"SNR_ant_triggered[4]=",0,"4, add in index ");
   G__memvar_setup((void*)((long)(&p->SNR_ant_closest)-(long)(p)),100,0,0,-1,-1,-1,1,"SNR_ant_closest[4]=",0,"4, add in index ");
   G__memvar_setup((void*)((long)(&p->SNR_ant_coherent)-(long)(p)),100,0,0,-1,-1,-1,1,"SNR_ant_coherent[4]=",0,"4, add in index ");
   G__memvar_setup((void*)((long)(&p->PowerCut)-(long)(p)),100,0,0,-1,-1,-1,1,"PowerCut[4]=",0,"4, add in index, might not keep as array of 40 antennas");
   G__memvar_setup((void*)((long)(&p->CoherentAnts)-(long)(p)),105,0,0,-1,-1,-1,1,"CoherentAnts[4]=",0,"4, list of antennas ");
   G__memvar_setup((void*)((long)(&p->distance_from_source)-(long)(p)),100,0,0,-1,-1,-1,1,"distance_from_source[4]=",0,"4 ");
   G__memvar_setup((void*)((long)(&p->peak2peak_signal)-(long)(p)),100,0,0,-1,-1,-1,1,"peak2peak_signal[4]=",0,"4, depend on after filtering which will change for each pol ");
   G__memvar_setup((void*)((long)(&p->peakVoltage_2)-(long)(p)),100,0,0,-1,-1,-1,1,"peakVoltage_2[4]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->peakSNR_2)-(long)(p)),100,0,0,-1,-1,-1,1,"peakSNR_2[4]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->peakVal_box)-(long)(p)),100,0,0,-1,-1,-1,1,"peakVal_box[4]=",0,"peakVal in box surrounding peak in Vpol");
   G__memvar_setup((void*)((long)(&p->max_correlation)-(long)(p)),100,0,0,-1,-1,-1,1,"max_correlation[4][40][40]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->RMS_ants)-(long)(p)),100,0,0,-1,-1,-1,1,"RMS_ants[4][40]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->time_delay)-(long)(p)),100,0,0,-1,-1,-1,1,"time_delay[4][40][40]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->num_bins_filtered)-(long)(p)),105,0,0,-1,-1,-1,1,"num_bins_filtered=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->mainrfcmflag)-(long)(p)),105,0,0,-1,-1,-1,1,"mainrfcmflag=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->bigenoughpeakflag)-(long)(p)),105,0,0,-1,-1,-1,1,"bigenoughpeakflag=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->dcoffsetflag)-(long)(p)),105,0,0,-1,-1,-1,1,"dcoffsetflag=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->shorttraceflag)-(long)(p)),105,0,0,-1,-1,-1,1,"shorttraceflag=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->nadirrfcmflag)-(long)(p)),105,0,0,-1,-1,-1,1,"nadirrfcmflag=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->payloadblastflag)-(long)(p)),105,0,0,-1,-1,-1,1,"payloadblastflag=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->SNR_coherent)-(long)(p)),100,0,0,-1,-1,-1,1,"SNR_coherent[4]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->SNR_coherent2)-(long)(p)),100,0,0,-1,-1,-1,1,"SNR_coherent2[4]=",0,"for powerSNR min value==1");
   G__memvar_setup((void*)((long)(&p->SNR_noise_bin)-(long)(p)),105,0,0,-1,-1,-1,1,"SNR_noise_bin[4]=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarclassdict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncanalysis_info_4pol(void) {
   /* analysis_info_4pol */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__classdictLN_analysis_info_4pol));
   // automatic default constructor
   G__memfunc_setup("analysis_info_4pol", 1869, G__classdict_168_0_1, (int) ('i'), G__get_linked_tagnum(&G__classdictLN_analysis_info_4pol), -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic copy constructor
   G__memfunc_setup("analysis_info_4pol", 1869, G__classdict_168_0_2, (int) ('i'), G__get_linked_tagnum(&G__classdictLN_analysis_info_4pol), -1, 0, 1, 1, 1, 0, "u 'analysis_info_4pol' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~analysis_info_4pol", 1995, G__classdict_168_0_3, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__classdict_168_0_4, (int) ('u'), G__get_linked_tagnum(&G__classdictLN_analysis_info_4pol), -1, 1, 1, 1, 1, 0, "u 'analysis_info_4pol' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncclassdict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalclassdict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcclassdict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__classdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__classdictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__classdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__classdictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__classdictLN_analysis_info_4pol = { "analysis_info_4pol" , 115 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableclassdict() {
  G__classdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__classdictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__classdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__classdictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__classdictLN_analysis_info_4pol.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableclassdict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__classdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__classdictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__classdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__classdictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__classdictLN_analysis_info_4pol),sizeof(analysis_info_4pol),-1,262144,(char*)NULL,G__setup_memvaranalysis_info_4pol,G__setup_memfuncanalysis_info_4pol);
}
extern "C" void G__cpp_setupclassdict(void) {
  G__check_setup_version(30051515,"G__cpp_setupclassdict()");
  G__set_cpp_environmentclassdict();
  G__cpp_setup_tagtableclassdict();

  G__cpp_setup_inheritanceclassdict();

  G__cpp_setup_typetableclassdict();

  G__cpp_setup_memvarclassdict();

  G__cpp_setup_memfuncclassdict();
  G__cpp_setup_globalclassdict();
  G__cpp_setup_funcclassdict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncclassdict();
  return;
}
class G__cpp_setup_initclassdict {
  public:
    G__cpp_setup_initclassdict() { G__add_setup_func("classdict",(G__incsetup)(&G__cpp_setupclassdict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initclassdict() { G__remove_setup_func("classdict"); }
};
G__cpp_setup_initclassdict G__cpp_setup_initializerclassdict;

