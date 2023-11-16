//macro to get a time smeared distributions for CLOUD DAQ

{
  TFile *fMyFile = new TFile("positron_2MeV.root");
  TTree *myTree = (TTree*)fMyFile->Get("Hits");

  Double_t posZ = 1.5;
  //Double_t posZ = 0.;
  Double_t fibre_dt = 2.8; // decay time WLS fibre
  Double_t myEff=100*(0.108*0.4*0.9);
  //Double_t myEff=100.;
  Double_t TTS = 0.4;
  Double_t timeFibre = 6.24; // average time, in ns, for photon to travel 1m of fibre
  Double_t Att_leng = 5.; // fibre attenuation length metres
    
  Double_t             fibreX;
  Double_t             fibreY;
  Double_t             Hit_X;
  Double_t             Hit_Y;
  Double_t             Hit_Z;
  Double_t             Time_ns;
  
  myTree->SetBranchAddress("fibreX" , &fibreX);
  myTree->SetBranchAddress("fibreY" , &fibreY);
  myTree->SetBranchAddress("Hit_X"  , &Hit_X);
  myTree->SetBranchAddress("Hit_Y"  , &Hit_Y);
  myTree->SetBranchAddress("Hit_Z"  , &Hit_Z);
  myTree->SetBranchAddress("Time_ns", &Time_ns);

  TH1D *hTime_front = new TH1D("hTft","Front electronics",501,-0.5,501.5);  
  TH1D *hTime_front_0 = new TH1D("hTf0","Scint. Time + Rndm. Walk",501,-0.5,501.5);  
  TH1D *hTime_front_1 = new TH1D("hTf1","+ WLS",501,-0.5,501.5);  
  TH1D *hTime_front_2 = new TH1D("hTf2","+ Ttransit Time",501,-0.5,501.5);  
  TH1D *hTime_back  = new TH1D("hTb","Back electronics",501,-0.5,501.5);
  TH1D *hTime_back_0  = new TH1D("hTb0","Scint. Time + Rndm. Walk",501,-0.5,501.5);
  TH1D *hTime_back_1  = new TH1D("hTb1","+ WLS",501,-0.5,501.5);
  TH1D *hTime_back_2  = new TH1D("hTb2","+ Transity Time",501,-0.5,501.5);
  TH1D *hTime_total = new TH1D("hTtotal","Total",501,-0.5,501.5);  

  hTime_front->GetYaxis()->SetTitle("Entries / 1 ns");
  hTime_front->GetXaxis()->SetTitle("Hit time (ns)");
  hTime_front->SetLineWidth(4); hTime_front->SetFillColor(kBlue+2); hTime_front->SetFillStyle(3005);
  hTime_back->SetLineWidth(4);
  hTime_back->SetLineColor(kRed+2);  hTime_back->SetMarkerColor(kRed+2);
  hTime_back->SetFillColor(kRed+2); hTime_back->SetFillStyle(3004);

  hTime_front_0->GetYaxis()->SetTitle("Entries / 1 ns");
  hTime_front_0->GetXaxis()->SetTitle("Hit time (ns)");
  hTime_front_0->SetLineWidth(4);
  hTime_front_1->SetLineWidth(4);
  hTime_front_2->SetLineWidth(4);
  hTime_back_0->SetLineWidth(4);
  hTime_front_1->SetLineColor(kRed+2);   hTime_front_1->SetMarkerColor(kRed+2); 
  hTime_front_2->SetLineColor(kViolet-3);   hTime_front_2->SetMarkerColor(kViolet-3);
  hTime_back_0->SetLineColor(kRed+2);   hTime_back_0->SetMarkerColor(kRed+2);

  hTime_total->SetLineWidth(2);
  hTime_total->SetLineColor(kBlack);

  Int_t mySeed = 0;
  TRandom3 *rand1 = new TRandom3(mySeed);
  

  Int_t nEvt=0, iEvt=0;
  Int_t pastEvt = nEvt;
  nEvt += myTree->GetEntries(Form("Event_Number==%i",iEvt));

  for (int i = pastEvt; i<nEvt; ++i) {
    myTree->GetEntry(i);
    Float_t detected = rand1->Rndm(); // choose if this hit makes a signal
      Float_t side = rand1->Rndm(); // choose if this hit goes to detector's back or front
      if (side <= 0.5) {

	Double_t distTravel = 2. - posZ  - Hit_Z / 1e3;
	Double_t totalEff = myEff * TMath::Exp(- distTravel / Att_leng);

	if(detected < totalEff/100.) {

	  if(distTravel < -4 || distTravel > 4) {
	    cout<<distTravel<<" "<< Hit_Z<<endl;
	  }
	  Double_t WLStime = rand1->Exp(fibre_dt); //decay time fibre
	  Double_t spreadT = (rand1->Rndm() - 0.5 ) * distTravel * TTS; //time spread, flat around average.
	  Double_t hitTime = Time_ns + WLStime + spreadT + distTravel * timeFibre;
	  hTime_front_0->Fill(Time_ns);
	  hTime_front_1->Fill(Time_ns + WLStime);
	  hTime_front_2->Fill(Time_ns + WLStime + distTravel * timeFibre + spreadT);
	  hTime_front->Fill(hitTime);
	  hTime_total->Fill(hitTime);
	  //}
	  //hTime->Fill(hitTime, myEff/100.);
	}
      } else {
	Double_t distTravel = 2. + posZ  + Hit_Z / 1e3;
	if(distTravel < -4 || distTravel > 4) {
	  cout<<distTravel<<" "<< Hit_Z<<endl;
	}
	Double_t totalEff = myEff * TMath::Exp(- distTravel / Att_leng);
	if(detected < totalEff/100.) {
	  Double_t WLStime = rand1->Exp(fibre_dt); //decay time fibre
	  Double_t spreadT = (rand1->Rndm() - 0.5 ) * distTravel * TTS;
	  Double_t hitTime = Time_ns + WLStime + spreadT + distTravel * timeFibre;
	  hTime_back_0->Fill(Time_ns);
	  hTime_back_1->Fill(Time_ns + WLStime);
	  hTime_back_2->Fill(Time_ns + WLStime + distTravel * timeFibre + spreadT);
	  hTime_back->Fill(hitTime);
	  hTime_total->Fill(hitTime);
	  //}
	  //hTime->Fill(hitTime, myEff/100.);
	}
      }
  }

  hTime_front->SetStats(kFALSE);
  hTime_front_0->SetStats(kFALSE);
  hTime_back->SetStats(kFALSE);

  Double_t maxCount = hTime_front->GetMaximum();
  if(hTime_back->GetMaximum() > maxCount)
    maxCount = hTime_back->GetMaximum();
  
  hTime_front->Draw();
  hTime_front->GetYaxis()->SetRangeUser(0.1,maxCount*1.1);
  //hTime_front_0->Draw();
  //hTime_front_1->Draw("same");
  //hTime_front_2->Draw("same");
  hTime_back->Draw("same");
  //hTime_total->Draw("same");
  //hTime_back_0->Draw("same");
  //hTime_front_0->GetXaxis()->SetRangeUser(0,80);
  hTime_front->GetXaxis()->SetRangeUser(0,100);

  //gPad->SetLogy();
  //gPad->BuildLegend(39.2097,1.78702,80.0083,3.10247);


  hTime_front->SetTitle(Form("Front Channels (%2.1fm away)", 2.0-posZ));
  hTime_back->SetTitle(Form("Back Channels (%2.1fm away)", 2.0+posZ));
  gPad->BuildLegend(0.491404,0.555789,0.893983,0.890526);
  //hTime_front_0->SetTitle(Form("%2.0fm fibre",2. - posZ));
  hTime_front->SetTitle(Form("Hit time distribution for 2 MeV positron at Z = %2.1fm",posZ));

  
    /// debug
  /*
  TH1D *hRatio = (TH1D*)hTime_front->Clone();
  hRatio->Sumw2();
  hTime_back->Sumw2();
  hRatio->Divide(hTime_back);
  hRatio->GetXaxis()->SetRangeUser(0,31);
  new TCanvas;
  hRatio->SetLineWidth(3);
  hRatio->Draw();
  */

  gPad->SetGridx();
  gPad->SetGridy();

  
  
}
