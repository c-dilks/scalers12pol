// builds scaler tree from acc file
// -- this also allows for empty bunches documented in "pathologies.dat" to be manually omitted
//    from relative luminosity computation ("CLEAN UP PROCEDURE")

void mk_tree(const char * acc_file="datfiles/acc.dat")
{
  // read acc file into tree
  TTree * acc = new TTree("acc","counts tree from acc.dat");
  char cols[2048];
  char bbc_cols[256];
  char zdc_cols[256];
  char vpd_cols[256];
  for(Int_t i=0; i<=7; i++)
  {
    if(i==0) 
    {
      sprintf(bbc_cols,"bbc_%d/D",i);
      sprintf(zdc_cols,"zdc_%d/D",i);
      sprintf(vpd_cols,"vpd_%d/D",i);
    }
    else
    {
      sprintf(bbc_cols,"%s:bbc_%d/D",bbc_cols,i);
      sprintf(zdc_cols,"%s:zdc_%d/D",zdc_cols,i);
      sprintf(vpd_cols,"%s:vpd_%d/D",vpd_cols,i);
    };
  };
  sprintf(cols,"i/I:runnum/I:fi/I:fill/I:t/D:bx/I:%s:%s:%s:tot_bx/D:blue/I:yell/I",bbc_cols,zdc_cols,vpd_cols);
  printf("%s\n",cols);
  acc->ReadFile(acc_file,cols);

  acc->Print();

  // set branch addresses to read through acc tree
  Int_t index,runnum,fill_index,fill,bx;
  Double_t bbc[8];
  Double_t zdc[8];
  Double_t vpd[8];
  Double_t time;
  Double_t tot_bx;
  Int_t blue,yell;
  acc->SetBranchAddress("i",&index);
  acc->SetBranchAddress("runnum",&runnum);
  acc->SetBranchAddress("fi",&fill_index);
  acc->SetBranchAddress("fill",&fill);
  acc->SetBranchAddress("t",&time);
  acc->SetBranchAddress("bx",&bx);
  char str[16];
  for(Int_t i=0; i<8; i++) { sprintf(str,"bbc_%d",i); acc->SetBranchAddress(str,&bbc[i]); };
  for(Int_t i=0; i<8; i++) { sprintf(str,"zdc_%d",i); acc->SetBranchAddress(str,&zdc[i]); };
  for(Int_t i=0; i<8; i++) { sprintf(str,"vpd_%d",i); acc->SetBranchAddress(str,&vpd[i]); };
  acc->SetBranchAddress("tot_bx",&tot_bx);
  acc->SetBranchAddress("blue",&blue);
  acc->SetBranchAddress("yell",&yell);
  

  // restructure tree into one suitable for analysis
  TTree * sca = new TTree("sca","restructured tree");
  Double_t bbce,bbcw,bbcx; // e=east, w=west, x=coincidence
  Double_t zdce,zdcw,zdcx;
  Double_t vpde,vpdw,vpdx;
  Bool_t okEntry,kicked;
  sca->Branch("i",&index,"i/I");
  sca->Branch("runnum",&runnum,"runnum/I");
  sca->Branch("fi",&fill_index,"fi/I");
  sca->Branch("fill",&fill,"fill/I");
  sca->Branch("t",&time,"t/D");
  sca->Branch("bx",&bx,"bx/I");
  sca->Branch("bbce",&bbce,"bbce/D");
  sca->Branch("bbcw",&bbcw,"bbcw/D");
  sca->Branch("bbcx",&bbcx,"bbcx/D");
  sca->Branch("zdce",&zdce,"zdce/D");
  sca->Branch("zdcw",&zdcw,"zdcw/D");
  sca->Branch("zdcx",&zdcx,"zdcx/D");
  sca->Branch("vpde",&vpde,"vpde/D");
  sca->Branch("vpdw",&vpdw,"vpdw/D");
  sca->Branch("vpdx",&vpdx,"vpdx/D");
  sca->Branch("tot_bx",&tot_bx,"tot_bx/D");
  sca->Branch("blue",&blue,"blue/I");
  sca->Branch("yell",&yell,"yell/I");
  sca->Branch("kicked",&kicked,"kicked/O");


  // read kicked bunches tree from "kicked" file
  TTree * kicked_tr = new TTree();
  kicked_tr->ReadFile("kicked","fill/I:bx/I:spinbit/I");
  Int_t kicked_fill,kicked_bx,kicked_spinbit;
  kicked_tr->SetBranchAddress("fill",&kicked_fill);
  kicked_tr->SetBranchAddress("bx",&kicked_bx);
  kicked_tr->SetBranchAddress("spinbit",&kicked_spinbit);
  
  for(Int_t q=0; q<acc->GetEntries(); q++)
  {
    acc->GetEntry(q);

    // -- see doc for bit details
    // BBC, ZDC, VPD bits: [ x w e ]
    bbce = bbc[1] + bbc[3] + bbc[5] + bbc[7]; // e + we + xe + xwe
    bbcw = bbc[2] + bbc[3] + bbc[6] + bbc[7]; // w + we + xw + xwe
    bbcx = bbc[3] + bbc[7]; // we + xwe

    zdce = zdc[1] + zdc[3] + zdc[5] + zdc[7]; // e + we + xe + xwe
    zdcw = zdc[2] + zdc[3] + zdc[6] + zdc[7]; // w + we + xw + xwe
    zdcx = zdc[3] + zdc[7]; // we + xwe
    
    vpde = vpd[1] + vpd[3] + vpd[5] + vpd[7]; // e + we + xe + xwe
    vpdw = vpd[2] + vpd[3] + vpd[6] + vpd[7]; // w + we + xw + xwe
    vpdx = vpd[3] + vpd[7]; // we + xwe


    // KICKED BUNCHES
    // manually omit empty bunches documented in pathologies.dat -- CLEAN UP PROCEDURE
    // (see 09.01.14 log entry)
    okEntry=true;
    // kicked bunches (presumably empty) 
    /*
    if(fill==17384 && (bx==29 || bx==30 || bx==117)) okEntry=false;
    if(fill==17416 && bx==79) okEntry=false;
    if(fill==17491 && bx==105) okEntry=false;
    if(fill==17519 && (bx==94 || bx==109)) okEntry=false;
    if(fill==17520 && bx==0) okEntry=false;
    if(fill==17529 && bx==97) okEntry=false;
    if(fill==17534 && bx==112) okEntry=false;
    if(fill==17553 && bx==73) okEntry=false;
    if(fill==17554 && (bx==7 || bx==14)) okEntry=false;
    if(fill==17555 && bx==61) okEntry=false;
    if(fill==17576 && bx==94) okEntry=false;
    // afterpulse-like bunches -- remove 1st 2 bunches after abort gaps 
    //if(fill==17512 && (bx>=40 && bx<=59)) okEntry=false;
    //if((fill>=17513 && fill<=17520) && ((bx>=0 && bx<=19) || (bx>=40 && bx<=59))) okEntry=false;
    */
    for(Int_t kk=0; kk<kicked_tr->GetEntries(); kk++)
    {
      kicked_tr->GetEntry(kk);
      if(fill==kicked_fill && bx==kicked_bx) okEntry=false;
    };
    
    kicked=!okEntry; // cleaned up analysis
    //kicked=0; // take all bXings

    sca->Fill(); // new kick method
  };

  TFile * outfile = new TFile("counts.root","RECREATE");
  acc->Write("acc");
  sca->Write("sca");
  printf("counts.root written\n");
};
      

