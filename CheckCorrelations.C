// checks the value of Pearson's correlation coefficients on a run-by-run basis
// for correlations of raw counts

void CheckCorrelations(const char * filename="counts.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("sca");

  Int_t IMAX_tmp = tr->GetMaximum("i");
  const Int_t IMAX = IMAX_tmp;

  // tbit
  // - 0 = bbc
  // - 1 = zdc
  // - 2 = vpd
  // cbit
  // - 0 = E
  // - 1 = W
  // - 2 = X
  // corrbit
  // - 0 = X.E
  // - 1 = X.W
  // - 2 = W.E
  char tbit[3][8];
  strcpy(tbit[0],"bbc");
  strcpy(tbit[1],"zdc");
  strcpy(tbit[2],"vpd");
  char cbit[3][4];
  strcpy(cbit[0],"e");
  strcpy(cbit[1],"w");
  strcpy(cbit[2],"x");


  // get minima & maxima of counts
  Double_t maxim[3][3]; // [tbit] [cbit]
  Double_t minim[3][3];
  char tcbit[3][3][16];
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      sprintf(tcbit[t][c],"%s%s",tbit[t],cbit[c]);
      maxim[t][c]=tr->GetMaximum(tcbit[t][c]);
      minim[t][c]=tr->GetMinimum(tcbit[t][c]);
    };
  };


  // define distribution of correlation coefficients
  TH1D * corr_dist[3][3]; // [tbit] [corrbit]
  char corr_dist_n[3][3][128];
  char corr_dist_t[3][3][128];

  // define correlation vs. run index plots
  TH1D * corr_vs_run[3][3];
  char corr_vs_run_n[3][3][128];
  char corr_vs_run_t[3][3][128];


  // 2d-histograms for computing correlation coefficients for each run
  // - since plot minima are sensitive to abort gap counts, I used a very fine binning 
  const Int_t NBINS = 1000;
  TH2D * corr[3][3][IMAX]; // [tbit] [corrbit] [run index]
  char corr_n[3][3][IMAX][64];
  char corr_t[3][3][IMAX][64];
  Int_t ccx,ccy;
  char project[128];
  char cut[256];
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t x=0; x<3; x++)
    {
      if(x==0)
      {
        ccy=2; //X
        ccx=0; //E
      }
      else if(x==1)
      {
        ccy=2; //X
        ccx=1; //W
      }
      else if(x==2)
      {
        ccy=1; //W
        ccy=0; //E
      };

      sprintf(corr_dist_n[t][x],"corr_dist_%s_%s%s",tbit[t],cbit[ccy],cbit[ccx]);
      sprintf(corr_dist_t[t][x],"#rho(%s%s,%s%s) distribution",tbit[t],cbit[ccy],tbit[t],cbit[ccx]);
      corr_dist[t][x] = new TH1D(corr_dist_n[t][x],corr_dist_t[t][x],1000,0.9,1.1);

      sprintf(corr_vs_run_n[t][x],"corr_vs_run_%s_%s%s",tbit[t],cbit[ccy],cbit[ccx]);
      sprintf(corr_vs_run_t[t][x],"#rho(%s%s,%s%s) vs. run index",tbit[t],cbit[ccy],tbit[t],cbit[ccx]);
      corr_vs_run[t][x] = new TH1D(corr_vs_run_n[t][x],corr_vs_run_t[t][x],IMAX,1,IMAX+1);

      for(Int_t i=0; i<IMAX; i++)
      {
        sprintf(corr_n[t][x][i],"corr_t%d_x%d_i%d",t,x,i);
        sprintf(corr_t[t][x][i],"corr_t%d_x%d_i%d",t,x,i);
        corr[t][x][i] = new TH2D(corr_n[t][x][i],corr_t[t][x][i],NBINS,minim[t][ccx],maxim[t][ccx],NBINS,minim[t][ccy],maxim[t][ccy]);
        sprintf(project,"%s%s:%s%s",tbit[t],cbit[ccy],tbit[t],cbit[ccx]);
        sprintf(cut,"!(bx>=31 && bx<40) && !(bx>=111) && !kicked && blue!=0 && yell!=0 && i==%d",i+1);
        tr->Project(corr_n[t][x][i],project,cut);
        corr_dist[t][x]->Fill(corr[t][x][i]->GetCorrelationFactor());
        corr_vs_run[t][x]->SetBinContent(i+1,corr[t][x][i]->GetCorrelationFactor());
        delete corr[t][x][i];
        printf("t=%d x=%d i=%d\n",t,x,i);
      };
    };
  };


  // write output
  TFile * outfile = new TFile("corr.root","RECREATE");
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t x=0; x<3; x++)
    {
      corr_dist[t][x]->Write();
    };
  };
  for(Int_t t=0; t<3; t++)
  {
    for(Int_t x=0; x<3; x++)
    {
      corr_vs_run[t][x]->Write();
    };
  };
};
