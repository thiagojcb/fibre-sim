// Input constants

Int_t     trials = 1000;
Int_t     mySeed = 0;
TRandom3 *rand1;

Double_t posZ      = 0.;                          // translation from centre (z=0)
Double_t fibre_dt  = 7.0;                         // decay time WLS fibre
Double_t TTS       = 0.4;                         // fibre transit time spread / m
Double_t timeFibre = 6.26;                        // average time, in ns, for photon to travel 1m of fibre
Double_t Att_leng  = 5.;                          // fibre attenuation length, in metres
Double_t factor    = 1.00;                        // factor to boost or decrease LY of Susie sim
Double_t myEff     = 100*(factor*0.108*0.4*0.9);  // trapping, QE, coupling
                                                  // with all these numbers + Susie 2 MeV electron file: 400 PE/MeV
Int_t sampling_bins = 160;                        // number of bins to sample the waveform (SMAPIC style)

Double_t jitter    = 0.0001;                         // electronics jitter (gauss sigma), ns

//input file
TFile *fMyFile;
TTree *myTree;
Double_t Hit_Z, Time_ns; // branches from Susie G4 output
Int_t    fibreNumber;

map<int, vector<double>> channelHits_front;
map<int, vector<double>> channelHits_back;

TH1D *hTime_front;
TH1D *hTime_front_0;
TH1D *hTime_front_1;
TH1D *hTime_front_2;
TH1D *hTime_back;
TH1D *hTime_back_0;
TH1D *hTime_back_1;
TH1D *hTime_back_2;
TH1D *hTime_total;

TH1D *hTime_Max_bin_cont;

TH1D *hRecoZ;
Double_t first_hit_f, first_hit_b;
TGraph *gBF;

TSpline3 *spline;

//TH1F channelWF[];
//TH1D hTime_front_map[];

map<int, TH1F> channelWF;
map<int, TH1D> hTime_front_map;

TH1F* channelWF_ref;
TAxis *xaxis;

TCanvas* cWF;

Double_t maxAmp;
Int_t maxChC = 0;
Int_t maxCh  = 0;

TH1D* hMaxAmp;
