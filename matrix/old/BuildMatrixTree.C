void BuildMatrixTree(const char * filename="matrix_out")
{
  // sort the matrix
  char sorter[256];
  char filename_sorted[64];
  sprintf(sorter,".! cat %s | sort -n > %s_sorted",filename,filename);
  sprintf(filename_sorted,"%s_sorted",filename);
  gROOT->ProcessLine(sorter);


  // read tree
  TTree * tr = new TTree();
  char branches[1024];
  tr->ReadFile(filename_sorted,"i/I:fi/I:t/D:tbit/I:cbit/I:bx/I:raw/D:acc/D:mul/D:fac/D:R1/D:R2/D:R3/D:R4/D:R5/D:R6/D:R7/D:R8/D:R9/D");


  // draw histogram
  gStyle->SetOptStat(0);
  Int_t max_i = tr->GetMaximum("i");
  TH2F * mat = new TH2F("mat","mat",max_i,1,max_i+1,120,1,121); // FIX THIS!!!
  tr->Project("mat","bx:i","mul/t*(tbit==1 && cbit==2)");
  mat->SetMinimum(500);
  TCanvas * cc = new TCanvas("cc","cc",1000,800);
  mat->Draw("colz");

  // draw fill lines
  TLine * fill_line[300];
  Int_t f_index,index;
  Int_t f_index_tmp;
  tr->SetBranchAddress("fi",&f_index);
  tr->SetBranchAddress("i",&index);
  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);
    if(f_index!=f_index_tmp)
    {
      fill_line[f_index-1] = new TLine(index,1,index,120);
      fill_line[f_index-1]->Draw();
      f_index_tmp = f_index;
    };
  };


  // write output
  TFile * outfile = new TFile("matrix.root","RECREATE");
  tr->Write("tr");
  mat->Write("mat");
  printf("matrix.root written\n");


  // print matrix for mathematica
  Int_t bn;
  Float_t bc;
  gSystem->RedirectOutput("for_mathematica","w");
  for(Int_t y=1; y<=mat->GetNbinsY(); y++)
  {
    for(Int_t x=1; x<=mat->GetNbinsX(); x++)
    {
      bn = mat->GetBin(x,y);
      bc = mat->GetBinContent(bn);
      if(x==mat->GetNbinsX()) printf("%f\n",bc);
      else printf("%f ",bc);
    };
  };
  gSystem->RedirectOutput(0);
  // -- old
  /*
  Int_t bn;
  Float_t bc;
  gSystem->RedirectOutput("for_mathematica");
  printf("{\n");
  for(Int_t y=1; y<=mat->GetNbinsY(); y++)
  {
    printf("{\n");
    for(Int_t x=1; x<=mat->GetNbinsX(); x++)
    {
      bn = mat->GetBin(x,y);
      bc = mat->GetBinContent(bn);
      if(x==mat->GetNbinsX()) printf("%f\n",bc);
      else printf("%f,",bc);
    };
    if(y==mat->GetNbinsY()) printf("}\n");
    else printf("},\n");
  };
  printf("}\n");
  gSystem->RedirectOutput(0);
  */




};
