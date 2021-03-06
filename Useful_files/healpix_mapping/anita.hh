////////////////////////////////////////////////////////////////////////////////////////////////
//class Anita:
////////////////////////////////////////////////////////////////////////////////////////////////
class TGraph;
class TFile;
class TTree;
class TH1F;
class Vector;
class Position;
class Balloon;
class IceModel;
class TF1;
class Settings;
class TRandom3;

using std::ofstream;
using std::vector;


class Anita {
    
private:
    
    TRandom Rand3;
    
    double GaintoHeight(double gain,double freq,double nmedium_receiver);
    TGraph *gshort[4];
    void setTrigRequirement(int WHICH);
    
public:
    
    static const int NBANDS_MAX=100; // max number of bands 
    static const int NFREQ=128;  // number of frequency bins
    //const int NFREQ=4096; 
    
    
    static const int NTRIG=5;
    // these should all be in anita
    static const int NANTENNAS_MAX=2000;
    static const int NLAYERS_MAX=5; // max number of layers (in smex design, it's 4)
    static const int NTRIGGERLAYERS_MAX=3;
    static const int NPHI_MAX=400; // max number of antennas around in phi (in smex, 16)
    Vector ANTENNA_POSITION_START[NLAYERS_MAX][NPHI_MAX]; // antenna positions from Kurt's measurements
    double ANTENNA_DOWN[NLAYERS_MAX][NPHI_MAX]; // down angles of antennas from Kurt's measurements
    
    
    
    int NRX_PHI[NLAYERS_MAX]; // number of antennas around in each layer. (radians)
    double PHI_EACHLAYER[NLAYERS_MAX][NPHI_MAX];// phi of the center of each antenna on each layer 
    //before correcting for offset for the layer. 
    //only used if it is cylindrically symmetric (radians)
    double PHI_OFFSET[NLAYERS_MAX]; // antenna offset in phi for each layer (radians)
    double THETA_ZENITH[NLAYERS_MAX]; // how the antenna is tilted in theta (in radians with 0=up)
    // 0=horizontal,+90=down
    
    // what the payload looks like
    
    double LAYER_VPOSITION[Anita::NLAYERS_MAX];  // position of layers in z relative to vertical center of the payload
    // anita proposal "says that the separation between upper and lower
    // 2 layers of antennas is just under 4m.
    // for anita hill, consider the positions of the "layers" of the "payload" (the stations) to be sitting on the horizontal grid defined by polar coordinates
    double LAYER_HPOSITION[Anita::NLAYERS_MAX]; // distance in horizontal plane between center axis of the "payload" and each "layer".
    double LAYER_PHIPOSITION[Anita::NLAYERS_MAX];//phi corresponding to the position of each "layer" on the "payload"
    double RRX[Anita::NLAYERS_MAX]; // radius that the antenna sits from the axis of the payload (feedpoint)
    
    Anita(); // constructor
    ~Anita();
    void Initialize(Settings *settings1,ofstream &foutput,int inu); // initialize a bunch of stuff
    // takes arrays that span NFREQ and turn them into arrays that span HALFNFOUR
    void MakeArraysforFFT(double *vsignalarray_e,double *vsignalarray_h,double *vsignal_e_forfft,double *vsignal_h_forfft);
    int Match(int ilayer,int ifold,int rx_minarrivaltime);
    int getLabAttn(int NPOINTS_LAB,double *freqlab,double *labattn);
    
    void labAttn(double *vhz);
    void SetNoise(Settings *settings1,Balloon *bn1,IceModel *antarctica);  
    
    TFile *fnoise;
    TTree *tdiode;
    
    static const int NFOUR=1024;
    static const int HALFNFOUR=512;
    
    // these are used for the satellite thing
    int NBANDS; // number of frequency sub-bands (not counting full band)
    int PERCENTBW; // percent bandwidth
    
    // these variables are for filling the tsignals tree
    double signal_vpol_inanita[5][HALFNFOUR]; // this is the signal waveform in the vertical polarization, before converting to LCP, RCP where applicable
    //double noise_vpol_inanita[5][HALFNFOUR]; // this is the noise waveform in the vertical polarization, before converting to LCP, RCP where applicable
    double total_vpol_inanita[5][HALFNFOUR];// this is the sum of the signal and noise in the vertical polarization, before converting to LCP, RCP where applicable
    
    double timedomainsignal_rfcm[HALFNFOUR];
    double timedomainsignal_lab[HALFNFOUR];
    
    
    
    double total_diodeinput_1_inanita[5][HALFNFOUR]; // this is the waveform that is input to the tunnel diode in the first (LCP or vertical) polarization
    double total_diodeinput_2_inanita[5][HALFNFOUR]; // this is the waveform that is input to the tunnel diode in the second (RCP or horizontal) polarization
    
    double timedomain_output_1_inanita[5][HALFNFOUR]; // this is just for writing out to the following tree
    double timedomain_output_2_inanita[5][HALFNFOUR]; // this is just for writing out to the following tree
    
    int flag_e_inanita[5][HALFNFOUR];
    int flag_h_inanita[5][HALFNFOUR];
    double dangle_inanita,emfrac_inanita,hadfrac_inanita;
    double ston[5]; // signal to noise;
    
    double peak_v_banding_rfcm_e[5]; // peak V in e polarization after rfcm's and banding
    double peak_v_banding_rfcm_h[5]; // peak V in h polarization
    double peak_e_rx_signalonly; // peak voltage in e polarization received by the antenna
    double peak_h_rx_signalonly; // peak voltage in h polarization received by the antenna
    double peak_e_rx_rfcm; // peak voltage in e polarization received by the antenna
    double peak_h_rx_rfcm; // peak voltage in h polarization received by the antenna
    
    double peak_e_rx_rfcm_signalonly; // peak voltage in e polarization received by the antenna
    double peak_h_rx_rfcm_signalonly; // peak voltage in h polarization received by the antenna
    
    double peak_e_rx_rfcm_lab; // peaks of the previous arrays
    double peak_h_rx_rfcm_lab;
    //I think this is the numerator of the vertical axis on Matt's plot
    
    
    
    
    int channels_passing_e[5]; // channels passing.  This is reset for every antenna for every event
    int channels_passing_h[5]; // channels passing 
    int l1_passing; // l1 passing
    
    double volts_rx_rfcm_lab_e[HALFNFOUR]; // time domain voltage vs. time after rx, rfcm's and lab
    double volts_rx_rfcm_lab_h[HALFNFOUR]; // for one specific antenna
    double volts_rx_rfcm_lab_e_all[48][HALFNFOUR]; // time domain voltage vs. time after rx, rfcm's and lab
    double volts_rx_rfcm_lab_h_all[48][HALFNFOUR]; // for all antennas
    
    int irx;
    void BoxAverageComplex(double *array,const int n,int navg);
    void BoxAverage(double *array,const int n,int navg);
    int GetRx(int ilayer, int ifold); // get antenna number based on which layer and position it is
    
    double avgfreq_rfcm[NFREQ];
    double avgfreq_rfcm_lab[NFREQ];
    
    
    double volts_rx_rfcm_e[HALFNFOUR]; // time domain voltage vs. time after rx, rfcm's and lab
    double volts_rx_rfcm_h[HALFNFOUR];
    
    
    
    double vmmhz_banding[NFREQ]; // V/m/MHz after banding 
    double vmmhz_banding_rfcm[NFREQ]; // V/m/MHz after banding and rfcms
    
    double rms_rfcm_e_single_event; // This is in Volts, not mV!
    
    // Note: The following 4 RMS noise variables are for all antennas of all events.
    // In fact, they don't represent RMS until after all events are finished!
    double rms_rfcm_e; // rms noise just after rfcm's
    double rms_lab_e;// rms noise at lab chip
    double rms_rfcm_h; // rms noise just after rfcm's
    double rms_lab_h;// rms noise at lab chip
    
    
    TFile *fsignals;
    TTree *tsignals;
    
    TFile *fdata;
    TTree *tdata; // writing data out for the analysers
    
    TTree *tglob;
    
    TH1F *hsignals[5]; // s/n (max diode output/mean diode output) for vertical polarization in each band
    
    int PULSER;
    
    double f_pulser[NFOUR/4];
    double f_phases[NFOUR/4];
    double f_noise[NFOUR/4];
    double v_pulser[NFOUR/4];
    double v_phases[NFOUR/4];
    double v_noise[NFOUR/4];
    
    
    
    double cumulat_prob[9];
    double cumulat_prob_plus1[9];
    
    
    // for filling tsignals tree
    double timedomainnoise_rfcm_banding_e[5][HALFNFOUR];
    double timedomainnoise_rfcm_banding_h[5][HALFNFOUR];
    double timedomainnoise_rfcm_e[HALFNFOUR];
    double timedomainnoise_rfcm_h[HALFNFOUR];
    double timedomainnoise_lab_e[HALFNFOUR];
    double timedomainnoise_lab_h[HALFNFOUR];
    
    double phases[5][HALFNFOUR];
    
    // for filling tglob
    int passglobtrig;
    double integral_vmmhz_foranita;
    
    
    int nnoiseevents; // total number of noise events we're choosing from
    int noiseeventcounter; // counts which event we're on so we go in order
    
    double FREQ_LOW; // lowest frequency
    double FREQ_HIGH; // highest frequency
    
    double NOTCH_MIN; // low edge of notch filter.  This is set in the input file
    double NOTCH_MAX; // upper edge of notch filter
    
    int BANDING;// set in the input file
    // whether or not you set your own banding (1)
    // or use anita-1 banding
    
    double freq[NFREQ];  // frequency for each bin
    double freq_forfft[NFOUR]; // frequencies for taking fft of signal
    double freq_forplotting[NFOUR/4]; // just one entry for frequency, unlike the above.
    double time[NFOUR/2];
    double time_long[NFOUR];
    
    double time_centered[NFOUR/2];
    double freqdomain_rfcm_banding[5][HALFNFOUR/2]; // average noise in frequency domain
    double freqdomain_rfcm[HALFNFOUR/2]; // average noise in frequency domain
    double freqdomain_rfcm_theory[HALFNFOUR/2]; // average noise in frequency domain
    double avgfreqdomain_lab[HALFNFOUR/2]; // average noise in frequency domain
    double phases_rfcm_banding_e[5][HALFNFOUR/2];
    double phases_rfcm_banding_h[5][HALFNFOUR/2];
    
    double phases_rfcm_e[HALFNFOUR/2];
    double phases_rfcm_h[HALFNFOUR/2];
    double phases_lab_e[HALFNFOUR/2];
    double phases_lab_h[HALFNFOUR/2];
    
    
    // this goes from 0 to Fmax, and represents both real and imaginary components
    
    
    
    
    void getDiodeModel();
    TF1 fdiode;
    double maxt_diode;
    int idelaybeforepeak[5];
    int iwindow[5];
    double diode_real[5][NFOUR]; // This is the time domain of the diode response. (actually NFOUR/2 array is used.)
    double fdiode_real[5][NFOUR]; // This is the fft of the diode response. (use all NFOUR array. This is for doubling array size for zero padding)
    
    
    void myconvlv(double *timedomain_forconvl,const int NFOUR,double *fdiode,double &maxdiodeconvl,double &onediodeconvl,double *power_noise,double *diodeconv);
    
    
    static int SurfChanneltoBand(int isurf);
    int AntennaWaveformtoSurf(int ilayer,int ifold); // find surf that generates this antenna's waveform
    static int AntennaNumbertoSurfNumber(int ilayer,int ifold); // find surf where this antenna is triggered
    static int GetAntennaNumber(int ilayer,int ifold); // given icemc indices ilayer, ifold, find antenna number as defined officially on anita
    static int GetLayer(int rx);
    static int GetIfold(int rx);
    static int GetSurfChannel(int antenna, int ibw,int ipol); // which channel on the surf this channel on this antenna corresponds to.
    static int WhichBand(int ibw,int ipol); // which band, 1-8, in order as they are on the surf
    void Banding(int j,double *freq_noise,double *powerperfreq,int NPOINTS_NOISE);
    void Banding(int iband,double *vmmhz);
    void RFCMs(int ilayer,int ifold,double *freq_noise,double *vmmhz,int NPOINTS_NOISE);
    void normalize_for_nsamples(double *spectrum, double nsamples, double nsamp);
    void convert_power_spectrum_to_voltage_spectrum_for_fft(double *spectrum, double domain[], double phase[]);
    void GetNoiseWaveforms(); // make time domain noise waveform based on avgnoise being the v^2
    //void GetNoiseWaveform(int iband); // make time domain noise waveform based on avgnoise being the v^2
    void GetPhases();
    int count_getnoisewaveforms; //count how many times we're in GetNoiseWaveforms for calculating rms voltages
    
    // each of the above graphs has 601 bins in it
    static const int NPOINTS_BANDS=601;
    
    double freq_bands[5][NPOINTS_BANDS]; // a frequency array for each of the four bands
    double attn_bands[5][NPOINTS_BANDS]; // attn array for each of the four bands in dB
    double bandsattn[5][NPOINTS_BANDS]; // as a fraction
    //double correl[4][NPOINTS_BANDS]; // correlation between each band and the fullband
    double correl_banding[5][NPOINTS_BANDS]; // correlation between each band and the fullband
    double correl_lab[NPOINTS_BANDS]; // correlation between each band and the fullband
    //double correl[5][NPOINTS_BANDS]; // correlation between each band and the fullband
    
    
    static const int NPOINTS_AMPL=58;// bins in amplification
    double freq_ampl[NANTENNAS_MAX][NPOINTS_AMPL]; // frequencies for each amp bin
    double ampl[NANTENNAS_MAX][NPOINTS_AMPL]; // amplification
    double ampl_notdb[NANTENNAS_MAX][NPOINTS_AMPL];// amplification again, but as a fraction
    double noisetemp[NANTENNAS_MAX][NPOINTS_AMPL]; // noise temp each amp bin
    
    
    static const int NPOINTS_NOISE=2000;
    
    
    //double bwslice_thresholds[4]={2.319,2.308,2.300,2.290}; // this allows you to set different thresholds for each band  
    double bwslice_vnoise[NLAYERS_MAX][5]; // expected noise voltage for antenna layer and 
    //each slice in bandwidth
    
    double probability[5];
    double bwslice_enoise[5]; // average integrated power in each band
    double bwslice_fwhmnoise[5]; // 1/2 of fwhm of hnoise
    double bwslice_rmsdiode[5]; // average rms diode output across noise waveforms in each band 
    double bwslice_meandiode[5]; // mean diode output across all samples in a sample of noise waveforms generated for each band
    double bwslice_vrms[5]; // rms noise voltage for this bandwidth slice 
    double freq_noise[5][NPOINTS_NOISE]; // frequency array that goes with vnoise array
    
    
    double impedence;
    double phase;
    double powerthreshold[5];
    double powerthreshold_nadir[5];
    
    double l1window; // time window where we require coincidences at L1
    
    double minsignalstrength; // minimum signal strength (measured as output of the diode) that a signal has to be for it to be worth adding to noise and performing the diode integration (each time we do this is uses up a noise waveform)
    
    double INTEGRATIONTIME; // integration time of the trigger diode
    static const int nsamp=100; // number of samples that were used to measure the noise data
    double TIMESTEP; // time step between samples
    double DEADTIME;
    
    double energythreshold; // relative to expected energy from noise
    
    
    
    int NTRIGGERLAYERS; // number of layers considered by the trigger.  may be different from nlayers, the number of physical layers on the payload.
    // In Anita 1 and Anita 2, the number of physical layers were 3 while the number of trigger layers were 2.
    int PHITRIG[NLAYERS_MAX]; // number of positions in phi for each trigger layer
    int REQUIRE_CENTRE; // require centre antenna in clump to be one of those hit
    static const int NTRIGPHISECTORS=16; // number of phi sectors in the trigger
    
    int GAINS;// whether to use constant gains as entered in GetBeamWidths (0) or to use Ped's measurements as entered in ReadGains (1)
    static const int NPOINTS_GAIN =131; // number of pointqs within bandwidth that gain is measured. from 200 MHz to 1.5 GHz with step size of 10 MHz
    double gainv_measured[NPOINTS_GAIN]; // may way of making the program use 3994760 fewer bytes than if these five arrays had still had 100000 elements
    double gainh_measured[NPOINTS_GAIN];
    double gainhv_measured[NPOINTS_GAIN];
    double gainvh_measured[NPOINTS_GAIN];
    double frequency_forgain_measured[NPOINTS_GAIN];
    
    double gain_angle[4][NPOINTS_GAIN][7]; /* first term: 0 = v polarization channel, a angle
					    1 = h polarization channel, a angle
					    2 = h polarization channel, e angle
					    3 = v polarization channel, e angle
					    second term: frequency
					    third term: angle */
    
    // frequency binning
    // anita proposal says frequency range is 0.2-1.2 GHz.
    // specs for the quad ridge antenna Model 0312-810 say 0.3-1.5 GHz
    double flare[4][NFREQ];  // for coarse antenna models:  beam width: e-plane: vp/hp, h-plane: vp/hp
    double gain[2][NFREQ];   // for coarse antenna models:  gain vert pol,h pol
    
    int GetBeamWidths(Settings *settings1); // for getting beam widths using coarse models (horn specs or simple model for EeVA)
    void Set_gain_angle(Settings *settings1,double nmedium_receiver);
    double Get_gain_angle(int gain_type, int k, double hitangle);
    void ReadGains();
    
    void AntennaGain(Settings *settings1,double hitangle_e,double hitangle_h,double e_component,double h_component,int k,double &vsignalarray_e,double &vsignalarray_h);
    
    
    
    double reference_angle[7]; // reference angles for finding gains of antenna
    
    double inv_angle_bin_size[6];
    int whichbin[NFREQ]; // these are for finding gains as a function of frequency
    double scalef2[NFREQ], scalef1[NFREQ]; // they are set in Set_gain_angle
    double vvGaintoHeight[NFREQ], hhGaintoHeight[NFREQ], hvGaintoHeight[NFREQ], vhGaintoHeight[NFREQ]; // holds results of the function double GaintoHeight
    
    double diffraction[2][89][NFREQ];
    void SetDiffraction();
    double GetDiffraction(int ilayer, double zenith_angle, int ifreq);
    
    
    
    static const int NPOINTS_LAB=272; // from note 137
    
    double freqlab[NPOINTS_LAB]; // frequency for each lab attn. bin
    
    double labattn[NPOINTS_LAB]; // lab attenuation
    
    
    double VNOISE[NLAYERS_MAX]; // noise calculated for each antenna layer depending on cant angle- this is only used right now for the chance in hell cuts
    
    
    int trigRequirements[NLAYERS_MAX];//  0th element - L1 - how many channels per antenna should pass
    // 1st element- L2 - how many antennas on a layer
    // 2nd element - L3 - how many L2 triggers should be coincident
    
    int antennatosurf[32];
    
    double maxthreshold;
    double bwslice_thresholds[5]; // thresholds for each band -- this is just an initialization- this is set in the input file
    int bwslice_allowed[5]; // these bands are allowed to contribute to the trigger sum -- this is set in the input file
    int bwslice_required[5]; // these bands are required to be among the channels that pass -- this is set in the input file
    int pol_allowed[2];// which polarisations are allowed to have channels that fire (V,H)
    
    
    
    double bwslice_center[5]; // center frequencies
    double bwslice_width[5]; // 3 dB bandwidths, without overlap
    
    
    double bwslice_min[5]; //minimum of each bandwidth slice
    
    double bwslice_max[5]; //minimum of each bandwidth slice
    double bwmin; // minimum width of any allowed bandwidth slice
    
    TFile* coherent_datafile;
    TTree* coherent_waveform_sum_tree;
    static const unsigned int NUM_COHERENT_ANTENNAS = 9;
    unsigned cwst_event_number;
    double cwst_rms_noise;
    double cwst_actual_rms;
    double cwst_threshold;
    unsigned cwst_window_start;
    unsigned cwst_window_end;
    double cwst_deg_theta;
    double cwst_deg_phi;
    double cwst_timesteps[HALFNFOUR];
    vector <double>* cwst_whole_wfms[NUM_COHERENT_ANTENNAS];
    vector <double>* cwst_wfms[NUM_COHERENT_ANTENNAS];
    vector <double>* cwst_aligned_wfms[NUM_COHERENT_ANTENNAS];
    vector <double> cwst_summed_wfm;
    vector <double> cwst_power_of_summed_wfm;
    double cwst_power;
    void fill_coherent_waveform_sum_tree(unsigned inu, Settings* settings1, double rms_noise, double actual_rms, unsigned window_start, unsigned window_end, double deg_theta, double deg_phi, double timesteps[], vector <vector <vector <double> > >&, vector <vector <vector <double> > >&, vector <vector <double > >&, vector <double>&, vector <double>&, double power);
    void GetPayload(Settings*, Balloon*);
    double VNOISE_ANITALITE[NPHI_MAX]; // noise for each antenna, for the anita-lite trigger configuration.
    double INCLINE_TOPTHREE; // cant angle of top three layers of antennas
    double INCLINE_NADIR; // cant angle of nadir (bottom) layer of antennas
    double INCLUDE_NADIRONLY; // cant angle of nadir (bottom) layer of antennas
    double LIVETIME;
    
    double SIGMA_THETA; // resolution on the polar angle of the signal
}; //class Anita

/*// stole this little number from Ryan's AnitaConventions.h
namespace AnitaPol {
    typedef enum EAnitaPol {
	kLeft=0,
	kRight=1
    } AnitaPol_t;
    }*/

