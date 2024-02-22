//macro to get a time smeared distributions for CLOUD DAQ

{
  TFile *fMyFile = new TFile("electron_2MeV.root");
  TTree *myTree = (TTree*)fMyFile->Get("Hits");

  Double_t posZ = 1.89;
  Double_t fibre_dt = 2.; // decay time WLS fibre
  Double_t myEff=2*100*(0.108*0.4*0.9); // trapping, QE, coupling
  //Double_t myEff=100;
  Double_t TTS = 0.4; // fibre transit time spread / m
  Double_t timeFibre = 6.26; // average time, in ns, for photon to travel 1m of fibre
  Double_t Att_leng = 5.; // fibre attenuation length metres
  
  Double_t Hit_Z, Time_ns; // branches from Susie G4 output
  Int_t    fibreNumber;
  myTree->SetBranchAddress("Hit_Z"      , &Hit_Z);
  myTree->SetBranchAddress("Time_ns"    , &Time_ns);
  myTree->SetBranchAddress("fibreNumber", &fibreNumber);

  map<int, vector<double>> channelHits_front;
  map<int, vector<double>> channelHits_back;
  
  TH1D *hTime_front = new TH1D("hTft","Front electronics",501,-0.5,501.5);  
  TH1D *hTime_front_0 = new TH1D("hTf0","Scint. Time + Rndm. Walk",501,-0.5,501.5);  
  TH1D *hTime_front_1 = new TH1D("hTf1","+ WLS",501,-0.5,501.5);  
  TH1D *hTime_front_2 = new TH1D("hTf2","+ Ttransit Time",501,-0.5,501.5);  
  TH1D *hTime_back  = new TH1D("hTb","Back electronics",501,-0.5,501.5);
  TH1D *hTime_back_0  = new TH1D("hTb0","Scint. Time + Rndm. Walk",501,-0.5,501.5);
  TH1D *hTime_back_1  = new TH1D("hTb1","+ WLS",501,-0.5,501.5);
  TH1D *hTime_back_2  = new TH1D("hTb2","+ Transity Time",501,-0.5,501.5);
  TH1D *hTime_total = new TH1D("hTtotal","Total",501,-0.5,501.5);

  TH1D *hTime_Max_bin_cont = new TH1D("hTime_Max_bin_cont","",101,-0.5,100.5);
  
  hTime_front->GetYaxis()->SetTitle("Entries / 1 ns");
  hTime_front->GetXaxis()->SetTitle("Hit time (ns)");
  hTime_front->SetLineWidth(4); hTime_front->SetFillColor(kBlue+2); hTime_front->SetFillStyle(3005);
  hTime_back->SetLineWidth(4);
  hTime_back->SetLineColor(kRed+2);  hTime_back->SetMarkerColor(kRed+2);
  hTime_back->SetFillColor(kRed+2); hTime_back->SetFillStyle(3004);

  hTime_front_0->GetYaxis()->SetTitle("Entries / 1 ns");
  hTime_front_0->GetXaxis()->SetTitle("Hit time (ns)");
  hTime_front_0->SetLineWidth(4);
  hTime_front_0->SetLineWidth(4); hTime_front_0->SetFillColor(kBlue+2); hTime_front_0->SetFillStyle(3005);
  hTime_front_1->SetLineWidth(4);
  hTime_front_2->SetLineWidth(4);
  hTime_back_0->SetLineWidth(4);
  hTime_front_1->SetLineColor(kRed+2);   hTime_front_1->SetMarkerColor(kRed+2);
  hTime_front_1->SetFillStyle(3004);   hTime_front_1->SetFillColor(kRed+2);
  hTime_front_2->SetLineColor(kViolet+2);   hTime_front_2->SetMarkerColor(kViolet+2);
  hTime_front_2->SetFillStyle(3003);   hTime_front_2->SetFillColor(kViolet+2);
  hTime_back_0->SetLineColor(kRed+2);   hTime_back_0->SetMarkerColor(kRed+2);

  hTime_total->SetLineWidth(2);
  hTime_total->SetLineColor(kBlack);

  TH1D *hRecoZ = new TH1D("hRecoZ","Reco Z resolution",100,-20,20);
  hRecoZ->GetXaxis()->SetTitle("Reco Z - True Z (cm)");
  hRecoZ->GetYaxis()->SetTitle("Entries");
  
  Int_t mySeed = 0;
  TRandom3 *rand1 = new TRandom3(mySeed);
  
  Int_t nEvt=0, iEvt=1, pastEvt=0;
  for(int i=0;i<=iEvt;++i){
    pastEvt=nEvt;
    nEvt += myTree->GetEntries(Form("Event_Number==%i",i));
    //cout<<i<<" "<<pastEvt<<" "<<nEvt<<endl;
  }

  Int_t trials = 1;
  Double_t first_hit_f, first_hit_b;

  TGraph *gBF = new TGraph();

  TNtuple t1("t1","","time:volt");
  //t1.ReadFile("SPE_WF.txt"); //positive pulse 
  t1.ReadFile("ampSim.txt"); //negative pulse
  t1.Draw("volt:time","","*");
  auto gWF = (TGraph*)gPad->GetPrimitive("Graph");
  TSpline3 *spline = new TSpline3("spline", gWF);
  spline->SetLineColor(kRed);
  spline->Draw("csame");
  
  for (int j=0; j<trials; ++j){ //randomization loop, for same event
    hTime_front_0->Reset();
    hTime_front_1->Reset();
    hTime_front_2->Reset();
    hTime_front->Reset();
    hTime_total->Reset();
    hTime_back_0->Reset();
    hTime_back_1->Reset();
    hTime_back_2->Reset();
    hTime_back->Reset();
    first_hit_f = 100;
    first_hit_b = 100;
    Double_t hitF, hitB, zF, zB;
  
    for (int i = pastEvt; i<nEvt; ++i) {
      myTree->GetEntry(i);

      Float_t detected = rand1->Rndm(); // to choose if this hit makes a signal
      Float_t side = rand1->Rndm();     // to choose if this hit goes to detector's back or front
      if (side <= 0.5) {
	
    	Double_t distTravel = 2. - posZ  - Hit_Z / 1e3;
	    Double_t totalEff = myEff * TMath::Exp(- distTravel / Att_leng);

    	if(detected < totalEff/100.) {
	  
	      if(distTravel < -4 || distTravel > 4) {//checking if translated hit is outside detector edge
	        cout<<distTravel<<" "<< Hit_Z<<endl;
          }
          Double_t WLStime = rand1->Exp(fibre_dt); //decay time fibre
	  Double_t spreadT = (rand1->Rndm() - 0.5 ) * distTravel * TTS; //time spread, flat around average.
          Double_t hitTime = Time_ns + WLStime + spreadT + distTravel * timeFibre; //total time, since G4 t0
    	  hTime_front_0->Fill(Time_ns);
          hTime_front_1->Fill(Time_ns + WLStime);
    	  hTime_front_2->Fill(Time_ns + WLStime + distTravel * timeFibre + spreadT);
          hTime_front->Fill(hitTime);
    	  hTime_total->Fill(hitTime);
	  //cout<<fibreNumber<<endl;
	  channelHits_front[fibreNumber].push_back(hitTime);

	  if(hitTime<first_hit_f){
	    first_hit_f = hitTime;
	    hitF = Time_ns;
	    zF = Hit_Z;
	    //cout<<"Front: "<<Hit_Z<<" "<<Time_ns<<" "<<WLStime<<" "<<spreadT<<" "<<distTravel<<" "<<hitTime<<endl;
          }
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

	  channelHits_back[fibreNumber].push_back(hitTime);
	  
	  if(hitTime<first_hit_b){
	    first_hit_b = hitTime;
	    hitB = Time_ns;
	    zB = Hit_Z;
	    //cout<<"Back: "<<Hit_Z<<" "<<Time_ns<<" "<<WLStime<<" "<<spreadT<<" "<<distTravel<<" "<<hitTime<<endl;
	  }
        }
      }
    }//hits on a event loop
    Double_t recoPosZ = (first_hit_b-first_hit_f)/(2*timeFibre);
    /*
      cout<<"Trial "<<j<<" , Front: "<<first_hit_f<<" , Back: "<<first_hit_b<<" , recoZ = "<<recoPosZ*100<<endl;
      cout<<"HitF = "<<hitF<<" . HitB = "<<hitB;
      cout<<" zF = "<<zF<<" . zB = "<<zB<<endl<<endl;
    */
    //gBF->AddPoint(zF,zB);
    gBF->AddPoint(hTime_total->Integral(),100*(recoPosZ - posZ));
    hRecoZ->Fill(100*(recoPosZ - posZ)); //in cm

    hTime_Max_bin_cont->Fill(hTime_front->GetMaximum());
    
  }//rndm loop

  /// plotting

  hTime_front->SetStats(kFALSE);
  hTime_front_0->SetStats(kFALSE);
  hTime_back->SetStats(kFALSE);

  Double_t maxCount = hTime_front->GetMaximum();
  if(hTime_back->GetMaximum() > maxCount)
    maxCount = hTime_back->GetMaximum();

  new TCanvas();
  hTime_front->Draw();
  hTime_front->GetYaxis()->SetRangeUser(0.1,maxCount*1.1);

  hTime_back->Draw("same");
  //hTime_total->Draw("same");
  //hTime_front->GetXaxis()->SetRangeUser(0,100);

  hTime_front->SetTitle(Form("Front Channels (%2.1fm away)", 2.0-posZ));
  hTime_back->SetTitle(Form("Back Channels (%2.1fm away)", 2.0+posZ));
  gPad->BuildLegend(0.491404,0.555789,0.893983,0.890526);
  hTime_front->SetTitle(Form("Hit time distribution for 2 MeV electron at Z = %2.1fm",posZ));
  gPad->SetGridx();
  gPad->SetGridy();
  
  new TCanvas();
  hRecoZ->Draw();
  gPad->SetGridx();
  gPad->SetGridy();
  
  new TCanvas();
  hTime_front_0->Draw();
  hTime_front_1->Draw("same");
  hTime_front_2->Draw("same");
  //hTime_back_0->Draw("same");
  gPad->SetLogy();
  hTime_front_0->GetXaxis()->SetRangeUser(0,80);
  gPad->BuildLegend(0.6689,0.7094,0.9971,0.995);
  hTime_front_0->SetTitle(Form("%2.1fm fibre",2. - posZ));
  
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

  new TCanvas;
  hTime_Max_bin_cont->Draw();
  
  gPad->SetGridx();
  gPad->SetGridy();
  
  TH1F** channelWF = new TH1F*[channelHits_front.size()];
  TH1D** hTime_front_map = new TH1D*[channelHits_front.size()];
  Int_t i=1;

  TH1F* channelWF_ref = new TH1F("channelWF_ref","channelWF_ref",160,0,100);// SAMPIC bin size: 6.2500000e-10 ns (64 samples giving 40ns)
  TAxis *xaxis = channelWF_ref->GetXaxis();
  
  for (const auto& [key, value] : channelHits_front){
    cout<<i<<" : "<<key<<" : ";
    channelWF[i] = new TH1F(Form("h%d",key),Form("Fibre %d, front ch",key),160,0,100);// SAMPIC bin size: 6.2500000e-10 ns (64 samples giving 40ns)
    channelWF[i]->GetXaxis()->SetTitle("tick time (ns)");
    channelWF[i]->GetYaxis()->SetTitle("Voltage (V)");
    hTime_front_map[i] = new TH1D(Form("hTfm%d",key),Form("Fibre %d, front ch",key),160,0,100);
    hTime_front_map[i]->GetXaxis()->SetTitle("Hit time (ns)");
    hTime_front_map[i]->GetYaxis()->SetTitle("Entries / 0.625 ns");

    for (const auto& n : value){
      hTime_front_map[i]->Fill(n);
      Int_t    iBin    = xaxis->FindBin(n);
      Double_t binW    = hTime_front_map[i]->GetBinWidth(iBin);
      Double_t iCenter = hTime_front_map[i]->GetBinCenter(iBin);
      Double_t tDif    = n - iCenter;
      Double_t t0      = tDif;
      if(tDif>0){
	iBin++; //start from the next bin, since this one will have 0 volts
	t0 = binW - tDif; //conpensate for the shift
      }else{
	t0 = t0*(-1.0); //otherwise will get region where no pulse info
      }
      
      for(int j=iBin; j<160;++j){
	//Double_t pulse_i = spline->Eval(1e-9 +(t0)*1e-9 + (j-iBin)*binW*1e-9); //positive pulse starts at 1ns
	Double_t pulse_i = spline->Eval(10e-9 +(t0)*1e-9 + (j-iBin)*binW*1e-9); //negative pulse starts at 10ns
	
	Double_t voltage = channelWF[i]->GetBinContent(j) + pulse_i;
	channelWF[i]->SetBinContent(j,voltage);
      }
      cout<<n<<" ";
    }
    cout<<endl;
    ++i;
  }
  TCanvas* c2 = new TCanvas;
  c2->Divide(2,1);
  c2->cd(1);
  hTime_front_map[17]->Draw();
  c2->cd(2);
  channelWF[17]->Draw();

  //hTime_front_map->Draw();
  //hTime_front->Draw("same");  

}
