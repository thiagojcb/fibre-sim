{
  TFile *fMyFile = new TFile("electron_2MeV.root");
  TTree *myTree = (TTree*)fMyFile->Get("Hits");
Double_t Hit_Z, Time_ns;
  myTree->SetBranchAddress("Hit_Z"  , &Hit_Z);
  myTree->SetBranchAddress("Time_ns", &Time_ns);
  Int_t nEvt = myTree->GetEntries("Event_Number==0");

  Double_t posZ = 0.;                 // translation to add "by hand". Susie generated events at Z=0cm, detector centre
  Double_t fibre_dt = 2.8;            // decay time WLS fibre
  Double_t myEff=100*(0.108*0.4*0.9); // trapping, QE, optical coupling, in %
  Double_t TTS = 0.4;                 // fibre transit time spread / metre
  Double_t timeFibre = 6.26;          // average time, in ns, for photon to travel 1m of fibre
  Double_t Att_leng = 5.;             // fibre attenuation length metres
  
  Int_t mySeed = 2;
  TRandom3 *rand1 = new TRandom3(mySeed);

  for (int i = 0; i<nEvt; ++i) {
      myTree->GetEntry(i);
      Float_t detected = rand1->Rndm(); // to choose if this hit makes a signal
      Float_t side = rand1->Rndm();     // to choose if this hit goes to detector's back or front side
      if (side <= 0.5) { //front electronics
	
    	Double_t distTravel = 2. - posZ  - Hit_Z / 1e3; //distance travelled from hit position to SiPM.
	    Double_t totalEff = myEff * TMath::Exp(- distTravel / Att_leng);

    	if(detected < totalEff/100.) { // we have a hit!
          Double_t WLStime = rand1->Exp(fibre_dt);                                 //decay time fibre
	      Double_t spreadT = (rand1->Rndm() - 0.5 ) * distTravel * TTS;            //time spread, flat around average.
          Double_t hitTime = Time_ns + WLStime + spreadT + distTravel * timeFibre; //hit total time since G4 t0, smeared by fibre properties
	  cout<<hitTime<<endl;
        }//detected-iff
      }//side-if
   }//for-loop
}
