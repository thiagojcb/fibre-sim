//macro to draw the attenuation length of the fibre (2 expo)

{
  TF1 f1("f1","[0]*exp(-x/[1])+[2]*exp(-x/[3])",0,30);
  f1.SetParameters(20,30,10,500);
  f1.SetLineColor(kBlack);
  f1.SetLineWidth(3);
  f1.SetTitle("Total");

  TF1 f2("f2","[0]*exp(-x/[1])",0,30);
  f2.SetParameters(20,30);
  f2.SetLineColor(kBlue);
  f2.SetLineWidth(3);
  f2.SetTitle("Short Att. Length, 30 cm");

  TF1 f3("f3","[0]*exp(-x/[1])",0,30);
  f3.SetParameters(10,500);
  f3.SetLineColor(kRed);
  f3.SetLineWidth(3);
  f3.SetTitle("Long Att. Length, 500 cm");

  f1.Draw();
  f1.GetXaxis()->SetTitle("Distance (cm)");
  f1.GetYaxis()->SetTitle("Light Intensity (A.U.)");
  f1.GetYaxis()->SetRangeUser(0,35);

  f2.Draw("same");
  f3.Draw("same");

  gPad->BuildLegend();

  gPad->SetGridy();
  gPad->SetGridx();
  gPad->SetTicky();
  gPad->SetTickx();
  

}
