// before running this, hadd all rootfiles together to rootfiles/all.root
//
// this script takes the total matx tree and draws the requested bx:i matrix, weighted
// by parameter "var" for the specified detector & bit combo (tbit & cbit)

void DrawMatrix(const char * var="mul", Int_t tbit0=1, Int_t cbit0=2, Bool_t printPNG=0)
{
  TFile * infile = new TFile("rootfiles/all.root","READ");
  TFile * outfile = new TFile("matrix.root","RECREATE");
  TTree * matx = (TTree*) infile->Get("matx");

  // strings
  char tbit_str[3][16];
  char cbit_str[3][16];
  strcpy(tbit_str[0],"bbc");
  strcpy(tbit_str[1],"zdc");
  strcpy(tbit_str[2],"vpd");
  strcpy(cbit_str[0],"e");
  strcpy(cbit_str[1],"w");
  strcpy(cbit_str[2],"x");
  char mat_title[128];
  sprintf(mat_title,"%s%s %s rate over bx:i",tbit_str[tbit0],cbit_str[cbit0],var);

  
  // check to see if matx tree is sorted by run index from hadd
  // -- if it's not, you'll need to write a sorting routine
  printf("checking matx...\n");
  Int_t index,findex;
  Int_t index_tmp=0;
  matx->SetBranchAddress("i",&index);
  matx->SetBranchAddress("fi",&findex);
  for(Int_t i=0; i<matx->GetEntries(); i++)
  {
    matx->GetEntry(i);
    if(index != index_tmp)
    {
      if(index == index_tmp+1) index_tmp = index;
      else 
      {
        fprintf(stderr,"ERROR: hadded matx tree not correctly sorted\n");
        return;
      };
    };
  };


  // draw 2d histogram for matrix
  printf("projecting...\n");
  gStyle->SetOptStat(0);
  Int_t max_i = matx->GetMaximum("i");
  TH2F * mat = new TH2F("mat",mat_title,max_i,1,max_i+1,120,1,121);
  mat->GetXaxis()->SetTitle("run index");
  mat->GetYaxis()->SetTitle("bXing");
  char project_cut[128];
  if(!strcmp(var,"fac")) sprintf(project_cut,"fac*(tbit==%d && cbit==%d)",tbit0,cbit0);
  else sprintf(project_cut,"%s/t*(tbit==%d && cbit==%d)",var,tbit0,cbit0);
  matx->Project("mat","bx:i",project_cut);
  if(strcmp(var,"fac")) mat->SetMinimum(500); // sets minimum if var != "fac"
  else mat->SetMinimum(0); // MINUS TWO ISSUE ?????????????????????????????????????????????????????????
  Float_t sf;
  if(printPNG) sf=3; else sf=1;
  TCanvas * cc = new TCanvas("cc","cc",sf*1000,sf*800);
  //cc->SetLogz();
  mat->Draw("colz");


  // set fill lines
  TLine * fill_line[500]; // assumes max 500 fills
  Int_t findex_tmp=0;
  for(Int_t i=0; i<matx->GetEntries(); i++)
  {
    matx->GetEntry(i);
    if(findex != findex_tmp)
    {
      fill_line[findex-1] = new TLine(index,1,index,120);
      fill_line[findex-1]->Draw();
      findex_tmp = findex;
    };
  };


  // write output
  //matx->Write("tr");
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
  printf("text form matrix written.\n");
  if(printPNG) cc->Print("cc.png","png");
}
