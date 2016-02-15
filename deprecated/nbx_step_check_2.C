// plots nbx vs bXing ... what is this structure?
// -- RUN IN BACKGROUND!
//
// -- for bin n, this fits the data from bins n-1 to n+1 to
// a straight line and then computes the difference between 
// bin n and where the fit says bin n should be

void nbx_step_check_2(const char * filename="counts.root")
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

  Double_t step[NRUNS][118];
  Double_t bx_arr[118];
  for(Int_t b=1; b<119; b++) bx_arr[b]=b;
  TGraph * tg[NRUNS];
  TF1 * fitf[118];
  char fitf_n[118][32];
  for(Int_t b=2; b<=119; b++) 
  {
    sprintf(fitf_n[b-2],"fit_%d",b-2);
    fitf[b-2] = new TF1(fitf_n[b-2],"pol1",b-2,b);
  };
  Double_t y;


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


    for(Int_t b=2; b<=119; b++)
    {
      h[i]->Fit(fitf[b-2],"Q","",b-2,b);
      y = fitf[b-2]->GetParameter(1) * (b-1) + fitf[b-2]->GetParameter(0);
      step[i][b-2] = h[i]->GetBinContent(b-1) - y;
    };
    tg[i] = new TGraph(120,bx_arr,step[i]);

    c->Clear();
    c->SetGrid(1,1);
    tg[i]->SetMarkerStyle(kFullCircle); 
    tg[i]->Draw("ape");

    if(i==0) c->Print("nbx_step_check_2.pdf(","pdf");
    else if(i+1==NRUNS) c->Print("nbx_step_check_2.pdf)","pdf");
    else c->Print("nbx_step_check_2.pdf");

  };

};
