// Input constants

Int_t     trials = 1000;
Int_t     mySeed = 0;
TRandom3 *rand1;

Double_t posZ      = 1.89;                        // translation from centre (z=0)
Double_t fibre_dt  = 2.;                          // decay time WLS fibre
Double_t TTS       = 0.4;                         // fibre transit time spread / m
Double_t timeFibre = 6.26;                        // average time, in ns, for photon to travel 1m of fibre
Double_t Att_leng  = 5.;                          // fibre attenuation length, in metres
Double_t factor    = 3.41;                        // factor to boost or decrease LY of Susie sim
Double_t myEff     = 100*(factor*0.108*0.4*0.9);  // trapping, QE, coupling
                                                    // with all these numbers + Susie 2 MeV electron file: 400 PE/MeV
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

TH1F** channelWF;
TH1D** hTime_front_map;

TH1F* channelWF_ref;
TAxis *xaxis;
