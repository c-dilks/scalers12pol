// draws beam intensity vs. day plots to compare the two

void DrawPlots()
{
  TTree * blue_tr = new TTree();
  TTree * yell_tr = new TTree();
  blue_tr->ReadFile("blue-wcm-wcmbeamm.dat","t/D:d/D:i/D");
  yell_tr->ReadFile("yell-wcm-wcmbeamm.dat","t/D:d/D:i/D");
  Double_t t_blue,d_blue,i_blue;
  Double_t t_yell,d_yell,i_yell;
  blue_tr->SetBranchAddress("t",&t_blue);
  blue_tr->SetBranchAddress("d",&d_blue);
  blue_tr->SetBranchAddress("i",&i_blue);
  yell_tr->SetBranchAddress("t",&t_yell);
  yell_tr->SetBranchAddress("d",&d_yell);
  yell_tr->SetBranchAddress("i",&i_yell);

  TGraph * blue_gr = new TGraph();
  TGraph * yell_gr = new TGraph();
  Int_t blue_pt=0;
  Int_t yell_pt=0;


  for(Int_t x=0; x<blue_tr->GetEntries(); x++)
  {
    blue_tr->GetEntry(x);
    blue_gr->SetPoint(blue_pt,d_blue,i_blue);
    blue_pt++;
  };
  for(Int_t x=0; x<yell_tr->GetEntries(); x++)
  {
    yell_tr->GetEntry(x);
    yell_gr->SetPoint(yell_pt,d_yell,i_yell);
    yell_pt++;
  };

  blue_gr->SetMarkerColor(kBlue);
  yell_gr->SetMarkerColor(kOrange);
  blue_gr->SetMarkerStyle(kFullCircle);
  yell_gr->SetMarkerStyle(kFullCircle);

  TCanvas * cc = new TCanvas("cc","cc",1500,600);
  TMultiGraph * mg = new TMultiGraph();
  mg->Add(blue_gr);
  mg->Add(yell_gr);
  mg->Draw("ape");
}
