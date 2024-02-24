
Double_t posZ      = 1.89;                        // translation from centre (z=0)
Double_t fibre_dt  = 2.;                          // decay time WLS fibre
Double_t TTS       = 0.4;                         // fibre transit time spread / m
Double_t timeFibre = 6.26;                        // average time, in ns, for photon to travel 1m of fibre
Double_t Att_leng  = 5.;                          // fibre attenuation length, in metres
Double_t factor    = 3.41;                        // factor to boost or decrease LY of Susie sim
Double_t myEff     = 100*(factor*0.108*0.4*0.9);  // trapping, QE, coupling
                                                    // with all these numbers + Susie 2 MeV electron file: 400 PE/MeV
