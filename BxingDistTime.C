// reads spin patterns and fills total bXing distribution (and four spin bits), 
// where each bin is weighted by run time (to account for shorter vs. longer runs);
//
// this was written to see if there is any structure in the bXing distributions

void BxingDistTime()
{
  TFile * tf = new TFile("counts.root","READ");
  TTree * tr = (TTree*) tf->Get("sca");
  TH1D * bxing_dist[5];
  char bxing_dist_n[5][64];
  char spin_cut[5][128];
  strcpy(spin_cut[0],"t*(!kicked && blue==-1 && yell==-1)");
  strcpy(spin_cut[1],"t*(!kicked && blue==-1 && yell==1)");
  strcpy(spin_cut[2],"t*(!kicked && blue==1 && yell==-1)");
  strcpy(spin_cut[3],"t*(!kicked && blue==1 && yell==1)");
  strcpy(spin_cut[4],"t*(!kicked && blue!=0 && yell!=0)");
  for(Int_t s=0; s<5; s++)
  {
    sprintf(bxing_dist_n[s],"bxing_dist_%d",s);
    bxing_dist[s] = new TH1D(bxing_dist_n[s],
                             "raw bXing dist * time (Grn:-- Orn:-+ Red:+- Blue:++ Blk:all)",
                             120,0,120);
    tr->Project(bxing_dist_n[s],"bx",spin_cut[s]);
  };
  bxing_dist[0]->SetLineColor(kGreen+2);
  bxing_dist[1]->SetLineColor(kOrange+7);
  bxing_dist[2]->SetLineColor(kRed);
  bxing_dist[3]->SetLineColor(kBlue);
  bxing_dist[4]->SetLineColor(kBlack);
  TCanvas * canv = new TCanvas("canv","canv",1100,500);
  bxing_dist[4]->Draw();
  for(Int_t s=0; s<4; s++) bxing_dist[s]->Draw("same");
};
