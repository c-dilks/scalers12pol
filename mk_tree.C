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

  Int_t IMAX_tmp = acc->GetMaximum("i");
  const Int_t IMAX = IMAX_tmp;

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


  // build arrays for restructuring; arrays are needed so that
  // we can implement bXing shift corrections
  Double_t bbce_arr[IMAX][120];
  Double_t bbcw_arr[IMAX][120];
  Double_t bbcx_arr[IMAX][120];
  Double_t zdce_arr[IMAX][120];
  Double_t zdcw_arr[IMAX][120];
  Double_t zdcx_arr[IMAX][120];
  Double_t vpde_arr[IMAX][120];
  Double_t vpdw_arr[IMAX][120];
  Double_t vpdx_arr[IMAX][120];
  Int_t runnum_arr[IMAX];
  Int_t fi_arr[IMAX];
  Int_t fill_arr[IMAX];
  Double_t time_arr[IMAX];
  Double_t tot_bx_arr[IMAX][120];
  Int_t blue_arr[IMAX][120];
  Int_t yell_arr[IMAX][120];
  Bool_t kicked_arr[IMAX][120];


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
    

    // store data into arrays, implementing bXing shift corrections on scalers
    if(fill==16570)
    {
      bbce_arr[index-1][(bx+113)%120] = bbce; // shift down 7 bXings
      bbcw_arr[index-1][(bx+113)%120] = bbcw;
      bbcx_arr[index-1][(bx+113)%120] = bbcx;
      zdce_arr[index-1][(bx+113)%120] = zdce; // shift down 7 bXings
      zdcw_arr[index-1][(bx+113)%120] = zdcw;
      zdcx_arr[index-1][(bx+113)%120] = zdcx;
      vpde_arr[index-1][(bx+113)%120] = vpde; // shift down 7 bXings
      vpdw_arr[index-1][(bx+113)%120] = vpdw;
      vpdx_arr[index-1][(bx+113)%120] = vpdx;
    }
    else if(fill == 16582 ||
            fill == 16586 ||
            fill == 16587 || 
            fill == 16592 || 
            fill == 16593 || 
            fill == 16594 || 
            fill == 16597 || 
            fill == 16602)
    {
      bbce_arr[index-1][bx] = bbce; // no shift
      bbcw_arr[index-1][bx] = bbcw;
      bbcx_arr[index-1][bx] = bbcx;
      zdce_arr[index-1][bx] = zdce; // no shift
      zdcw_arr[index-1][bx] = zdcw;
      zdcx_arr[index-1][bx] = zdcx;
      vpde_arr[index-1][(bx+1)%120] = vpde; // shift up 1 bXings
      vpdw_arr[index-1][(bx+1)%120] = vpdw;
      vpdx_arr[index-1][(bx+1)%120] = vpdx;
    }
    else
    {
      bbce_arr[index-1][bx] = bbce;
      bbcw_arr[index-1][bx] = bbcw;
      bbcx_arr[index-1][bx] = bbcx;
      zdce_arr[index-1][bx] = zdce;
      zdcw_arr[index-1][bx] = zdcw;
      zdcx_arr[index-1][bx] = zdcx;
      vpde_arr[index-1][bx] = vpde;
      vpdw_arr[index-1][bx] = vpdw;
      vpdx_arr[index-1][bx] = vpdx;
    };

    runnum_arr[index-1] = runnum;
    fi_arr[index-1] = fill_index;
    fill_arr[index-1] = fill;
    time_arr[index-1] = time;
    tot_bx_arr[index-1][bx] = tot_bx;
    blue_arr[index-1][bx] = blue;
    yell_arr[index-1][bx] = yell;
    kicked_arr[index-1][bx] = kicked;
  };


  // fill restructured tree
  for(Int_t i=0; i<IMAX; i++)
  {
    index = i+1;
    runnum = runnum_arr[i];
    fill_index = fi_arr[i];
    fill = fill_arr[i];
    time = time_arr[i];
    for(Int_t b=0; b<120; b++)
    {
      bx = b;
      bbce = bbce_arr[i][b];
      bbcw = bbcw_arr[i][b];
      bbcx = bbcx_arr[i][b];
      zdce = zdce_arr[i][b];
      zdcw = zdcw_arr[i][b];
      zdcx = zdcx_arr[i][b];
      vpde = vpde_arr[i][b];
      vpdw = vpdw_arr[i][b];
      vpdx = vpdx_arr[i][b];
      tot_bx = tot_bx_arr[i][b];
      blue = blue_arr[i][b];
      yell = yell_arr[i][b];
      kicked = kicked_arr[i][b];
      sca->Fill();
    };
  };

  TFile * outfile = new TFile("counts.root","RECREATE");
  acc->Write("acc");
  sca->Write("sca");
  printf("counts.root written\n");
};
      

