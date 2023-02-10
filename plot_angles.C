//macro to plot the angles using root

{
  TNtuple t1("t1","","angle:angle2:angle3:isCore:isClad1:isClad2:isLoss");
  t1.ReadFile("Fibre_toy_MC_1.0.txt");
  //t1.ReadFile("Fibre_toy_MC_1.33.txt");

  TH1D *hCore    = new TH1D("hCore"   ,"Core Light",900,0,90);
  TH1D *hCladIn  = new TH1D("hCladIn" ,"Inner Clad Light",900,0,90);
  TH1D *hCladOut = new TH1D("hCladOut","Outter Clad Light",900,0,90);
  TH1D *hLoss    = new TH1D("hLoss"   ,"Loss Light",900,0,90);

  t1.Draw("angle>>hCore"   ,"isCore");
  t1.Draw("angle>>hCladIn" ,"isClad1","same");
  t1.Draw("angle>>hCladOut","isClad2","same");
  t1.Draw("angle>>hLoss"   ,"isLoss","same");

  hCore->SetLineColor(kBlue);
  hCladIn->SetLineColor(kOrange);
  hCladOut->SetLineColor(kOrange+2);
  hLoss->SetLineColor(kRed);

  hCore->SetMarkerColor(kBlue);
  hCladIn->SetMarkerColor(kOrange);
  hCladOut->SetMarkerColor(kOrange+2);
  hLoss->SetMarkerColor(kRed);

  hCore->SetFillColor(kBlue);
  hCladIn->SetFillColor(kOrange);
  hCladOut->SetFillColor(kOrange+2);
  hLoss->SetFillColor(kRed);

  gPad->BuildLegend();
  hCore->SetStats(0);
  hCore->SetTitle("Air as outter media");
  //hCore->SetTitle("Water as outter media");
  hCore->GetXaxis()->SetTitle("Light angle [deg]");
  hCore->GetYaxis()->SetTitle("Entries / 0.1 deg");
}
