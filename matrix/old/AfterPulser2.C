void AfterPulser2(Int_t ncycles=1e4)
{
  const Double_t PHD_max = 100;
  const Double_t TH_base = 15; 
  Double_t tau = PHD_max/5;
  Double_t tau_noise = PHD_max/150;


  TRandom * RNG = new TRandom();


  TH1D * PHD = new TH1D("PHD","pulse height",PHD_max,0,PHD_max);
  TH1D * NHD = new TH1D("NHD","noise",PHD_max,0,PHD_max);
  TH1D * BXD = new TH1D("BXD","baseline bXing",120,0,120);
  Double_t pulse_height;
  Double_t noise_height;
  Int_t n_above_thresh=0;

  Double_t threshold[120];

  for(Int_t i=0; i<ncycles; i++)
  {
    for(Int_t bx=0; bx<120; bx++)
    {
      pulse_height = RNG->Exp(tau);
      noise_height = RNG->Exp(tau_noise);

      if((bx>=31 && bx<=39) || (bx>=111 && bx<=120)) 
      {
        pulse_height = noise_height; // set pulse equal to background noise
        NHD->Fill(noise_height);
      }
      else 
      {
        pulse_height += noise_height; // add background noise to pulse
        PHD->Fill(pulse_height);
        if(pulse_height>TH_base) n_above_thresh++;
      };
      if(pulse_height>TH_base) BXD->Fill(bx);
    };
  };

  TCanvas * PHD_canv = new TCanvas("PHD_canv","PHD_canv",700,500);
  PHD->Draw();
  TCanvas * NHD_canv = new TCanvas("NHD_canv","NHD_canv",700,500);
  NHD_canv->SetLogy();
  NHD->Draw();
  printf("%.1f%% events above baseline threshold\n",((Double_t)n_above_thresh)/((Double_t)(120*ncycles))*100);
  TCanvas * BXD_canv = new TCanvas("BXD_canv","BXD_canv",700,500);
  BXD->Draw();
}
