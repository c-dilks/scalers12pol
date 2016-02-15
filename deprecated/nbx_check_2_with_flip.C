// plots nbx vs bXing ... what is this structure?
// -- RUN IN BACKGROUND!
// -- outputs nbx_vs_bxing.pdf
// -- SYMMETRY TEST (AKIO'S IDEA)

void nbx_check_2_with_flip(const char * filename="counts.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("sca");

  Int_t NRUNS_tmp = tr->GetMaximum("i");
  const Int_t NRUNS = NRUNS_tmp;
  TH1D * h[NRUNS];
  TH1D * s[NRUNS]; // symmetry test
  char h_n[NRUNS][256];
  char s_n[NRUNS][256];
  char cut[NRUNS][256];

  TCanvas * c = new TCanvas("c","c",1400,1000);

  Double_t max;

  gStyle->SetOptStat(0);


  for(Int_t i=0; i<NRUNS; i++)
  {
    printf("i=%d\n",i);

    sprintf(h_n[i],"N_{bx} vs. bXing for i=%d",i+1);
    sprintf(s_n[i],"N_{bx} vs. bXing for i=%d sym",i+1);
    h[i] = new TH1D(h_n[i],h_n[i],120,0,120);
    s[i] = new TH1D(s_n[i],s_n[i],120,0,120);
    sprintf(cut[i],"tot_bx*(i==%d)",i+1);
    tr->Project(h_n[i],"bx",cut[i]);
    tr->Project(s_n[i],"119-bx",cut[i]);


    s[i]->SetLineColor(kRed);

    max = h[i]->GetMaximum();
    h[i]->Scale(1/max);
    s[i]->Scale(1/max);

    c->Clear();
    c->SetGrid(1,1);
    h[i]->Draw();
    s[i]->Draw("same");

    if(i==0) c->Print("nbx_vs_bxing.pdf(","pdf");
    else if(i+1==NRUNS) c->Print("nbx_vs_bxing.pdf)","pdf");
    else c->Print("nbx_vs_bxing.pdf");
  };
};
