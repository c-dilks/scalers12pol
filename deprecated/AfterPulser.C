// afterpulsing algorithm

void AfterPulser(Int_t sigma_offset=1)
{
  // --- control variables ---
  Int_t scalers=1e6; // total number of bXings (possible scales)
  // random number threshholds: prob must be greater than threshhold for each case
  Double_t nonabort_th = 5e-2; // non-abort scaler
  //Double_t abort_th = 1-0.14e-3; // abort scaler
  Double_t abort_th = 0.99; // abort scaler
  Double_t nonabort_afterpulse_th = 0.1; // afterpulsing from nonabort bXing
  Double_t abort_afterpulse_th = 0.9; // afterpulsing from abort bXing

  // define distributions
  TH1D * base = new TH1D("base","base",120,0,120); // base bXing dist
  TH1D * pulsed = new TH1D("pulsed","pulsed",120,0,120); // after-pulsed bXing dist
  TH1D * offset_dist = new TH1D("offset_dist","offset_dist",100,0,sigma_offset*3);

  // bXing loop
  Int_t bx=0;
  TRandom * rng = new TRandom();
  Double_t prob; // random number - uniform
  Double_t prob_gaus; // random number -- gaussian around zero
  Int_t offset;
  for(Int_t s=0; s<scalers; s++)
  {
    prob = rng->Uniform(1);
    prob_gaus = rng->Gaus(0,sigma_offset);
    //prob_gaus = rng->Uniform(sigma_offset);
    if(prob_gaus<0) prob_gaus*=-1;
    offset = (Int_t) prob_gaus; // truncate
    
    // add count for non-aborts
    if(!((bx>=31 && bx<=39) || (bx>=111 && bx<=120)) && prob > nonabort_th)
    {
      base->Fill(bx);
      pulsed->Fill(bx);
      if(prob>nonabort_afterpulse_th) 
      {
        pulsed->Fill((bx+offset)%120);
        offset_dist->Fill(offset);
      };
    }
    // add count for aborts
    else if(((bx>=31 && bx<=39) || (bx>=111 && bx<=120)) && prob > abort_th)
    {
      base->Fill(bx);
      pulsed->Fill(bx);
      if(prob>abort_afterpulse_th) 
      {
        pulsed->Fill((bx+offset)%120);
        offset_dist->Fill(offset);
      };
    };

    bx=(bx+1)%120;
  };

  // normalise pulsed distribution
  pulsed->Scale(1.0/pulsed->GetMaximum());

  // draw
  gStyle->SetOptStat(0);
  TCanvas * base_canv = new TCanvas("base_canv","base_canv",700,500);
  base->Draw();
  TCanvas * offset_dist_canv = new TCanvas("offset_dist_canv","offset_dist_canv",700,500);
  offset_dist->Draw();



  // read bXing distribution from matx
  TFile * matx_file = new TFile("rootfiles/all.root","READ");
  TTree * matx = (TTree*) matx_file->Get("matx");
  TH1D * matx_dist = new TH1D("matx_dist","matx_dist",120,0,120);
  matx->Project("matx_dist","bx-1","mul*(i>700 && tbit==1 && cbit==1)");
  matx_dist->Scale(1.0/matx_dist->GetMaximum());

  // compare dists
  TH1D * div = new TH1D("div","sim/data",120,0,120);
  div->Divide(pulsed,matx_dist,1.0,1.0);

  // set colors
  matx_dist->SetLineColor(kBlue);
  pulsed->SetLineColor(kRed);


  // draw comparison
  TCanvas * pulsed_canv = new TCanvas("pulsed_canv","pulsed_canv",700,500);
  pulsed->Draw();
  matx_dist->Draw("same");

  TCanvas * div_canv = new TCanvas("div_canv","div_canv",700,500);
  div->Draw();
}
