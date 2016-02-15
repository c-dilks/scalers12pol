// plots nbx vs bXing ... what is this structure?
// -- RUN IN BACKGROUND!
//
// computes difference between bin n+1 and n

void nbx_step_check(const char * filename="counts.root")
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

  Double_t step[NRUNS][120];
  Double_t step_e[NRUNS][120];
  Double_t bx_arr[120];
  Double_t bx_arr_e[120];
  for(Int_t b=0; b<120; b++)
  {
    bx_arr[b]=b;
    bx_arr_e[b]=0;
  };
  TGraphErrors * tg[NRUNS];

  for(Int_t i=0; i<NRUNS; i++)
  {
    printf("i=%d\n",i);

    sprintf(h_n[i],"nbx_bx_%d",i+1);
    sprintf(h_t[i],"N_{bx} vs. bXing for i=%d",i+1);
    h[i] = new TH1D(h_n[i],h_n[i],120,0,120);
    sprintf(cut[i],"tot_bx*(i==%d)",i+1);
    tr->Project(h_n[i],"bx",cut[i]);

    max = h[i]->GetMaximum();
    //h[i]->Scale(1/max);




    for(Int_t b=1; b<=119; b++)
    {
      step[i][b] = h[i]->GetBinContent(b) - h[i]->GetBinContent(b+1);
      step_e[i][b] = sqrt(step[i][b]);
    };
    tg[i] = new TGraphErrors(120,bx_arr,step[i],bx_arr_e,step_e[i]);

    c->Clear();
    c->SetGrid(1,1);
    tg[i]->SetMarkerStyle(kFullCircle); 
    tg[i]->Draw("ape");

    if(i==0) c->Print("nbx_step_check.pdf(","pdf");
    else if(i+1==NRUNS) c->Print("nbx_step_check.pdf)","pdf");
    else c->Print("nbx_step_check.pdf");

  };

};
