//macro to draw the attenuation length of the fibre (2 expo)

{
  TF1 f1("f1","[0]*exp(-x/[1])+[2]*exp(-x/[3])",0,100);
  f1.SetParameters(20,30,10,500);
  f1.SetLineColor(kBlack);
  f1.SetLineWidth(3);
  f1.SetTitle("Total");

  TF1 f2("f2","[0]*exp(-x/[1])",0,100);
  f2.SetParameters(20,30);
  f2.SetLineColor(kBlue);
  f2.SetLineWidth(3);
  f2.SetTitle("Short Att. Length, 30 cm");

  TF1 f3("f3","[0]*exp(-x/[1])",0,100);
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

  cout<<"Total light at 65cm: "<< f1.Eval(65) <<endl;
  cout<<"Short light component at 65cm: "<< f2.Eval(65) <<endl;
  cout<<"Long light component at 65cm: "<< f3.Eval(65) <<endl;

  TLine l1(65,0,65,35);
  l1.SetLineWidth(3);
  l1.SetLineStyle(7);
  l1.Draw();
    
  
}
