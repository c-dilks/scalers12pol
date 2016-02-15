// plots nbx vs bXing for each run ... what is this structure?
// -- RUN IN BACKGROUND!
// -- outputs nbx_vs_bxing.pdf

void nbx_check_2(const char * filename="counts.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("sca");

  Int_t NRUNS_tmp = tr->GetMaximum("i");
  const Int_t NRUNS = NRUNS_tmp;
  TH1D * h[NRUNS];
  char h_n[NRUNS][256];
  char h_t[NRUNS][256];
  char cut[NRUNS][256];

  TCanvas * c = new TCanvas("c","c",1400,1000);

  Double_t max;

  for(Int_t i=0; i<NRUNS; i++)
  {
    printf("i=%d\n",i);

    sprintf(h_n[i],"nbx_bx_%d",i+1);
    sprintf(h_t[i],"N_{bx} vs. bXing for i=%d",i+1);
    h[i] = new TH1D(h_n[i],h_n[i],120,0,120);
    sprintf(cut[i],"tot_bx*(i==%d)",i+1);
    tr->Project(h_n[i],"bx",cut[i]);

    max = h[i]->GetMaximum();
    h[i]->Scale(1/max);

    c->Clear();
    c->SetGrid(1,1);
    h[i]->Draw();

    if(i==0) c->Print("nbx_check/nbx_vs_bxing.pdf(","pdf");
    else if(i+1==NRUNS) c->Print("nbx_check/nbx_vs_bxing.pdf)","pdf");
    else c->Print("nbx_check/nbx_vs_bxing.pdf");
  };
};
