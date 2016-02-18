// builds scaler tree from acc file
// -- this also allows for empty bunches documented in "pathologies.dat" to be manually omitted
//    from relative luminosity computation ("CLEAN UP PROCEDURE")

void mk_tree(const char * acc_file="datfiles/acc.dat", Int_t board)
{
  // read acc file into tree
  TString outfile_n = Form("counts_bd%d.root",board);
  TFile * outfile = new TFile(outfile_n.Data(),"RECREATE");
  TTree * acc = new TTree("acc","counts tree from acc.dat");
  char cols[2048];
  char hor_cols[256];
  char ver_cols[256];
  char trun_cols[256];
  for(Int_t i=0; i<=7; i++)
  {
    if(i==0) 
    {
      sprintf(hor_cols,"hor%d/D",i);
      sprintf(ver_cols,"ver%d/D",i);
      sprintf(trun_cols,"trun%d/D",i);
    }
    else
    {
      sprintf(hor_cols,"%s:hor%d/D",hor_cols,i);
      sprintf(ver_cols,"%s:ver%d/D",ver_cols,i);
      sprintf(trun_cols,"%s:trun%d/D",trun_cols,i);
    };
  };
  sprintf(cols,"i/I:runnum/I:fi/I:fill/I:t/D:bx/I:%s:%s:%s:front/D:back/D:goodTAC/D:tot_bx/D:blue/I:yell/I",hor_cols,ver_cols,trun_cols);
  printf("%s\n",cols);
  acc->ReadFile(acc_file,cols);

  acc->Print();

  Int_t IMAX_tmp = acc->GetMaximum("i");
  const Int_t IMAX = IMAX_tmp;

  // set branch addresses to read through acc tree
  Int_t index,runnum,fill_index,fill,bx;
  Double_t hor[8];
  Double_t ver[8];
  Double_t trun[8];
  Double_t front,back,goodTAC;
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
  for(Int_t i=0; i<8; i++) { sprintf(str,"hor%d",i); acc->SetBranchAddress(str,&hor[i]); };
  for(Int_t i=0; i<8; i++) { sprintf(str,"ver%d",i); acc->SetBranchAddress(str,&ver[i]); };
  for(Int_t i=0; i<8; i++) { sprintf(str,"trun%d",i); acc->SetBranchAddress(str,&trun[i]); };
  acc->SetBranchAddress("front",&front);
  acc->SetBranchAddress("back",&back);
  acc->SetBranchAddress("goodTAC",&goodTAC);
  acc->SetBranchAddress("tot_bx",&tot_bx);
  acc->SetBranchAddress("blue",&blue);
  acc->SetBranchAddress("yell",&yell);


  // build arrays for restructuring; arrays are needed so that
  // we can implement bXing shift corrections

  Double_t hor_arr[8][IMAX][128];
  Double_t ver_arr[8][IMAX][128];
  Double_t trun_arr[8][IMAX][128];
  Double_t front_arr[IMAX][128];
  Double_t back_arr[IMAX][128];
  Double_t goodTAC_arr[IMAX][128];

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
  Bool_t okEntry,kicked;
  Double_t hor_shift[8];
  Double_t ver_shift[8];
  Double_t trun_shift[8];
  Double_t front_shift,back_shift,goodTAC_shift;
  Int_t bxq,fq;
  sca->Branch("i",&index,"i/I");
  sca->Branch("runnum",&runnum,"runnum/I");
  sca->Branch("fi",&fill_index,"fi/I");
  sca->Branch("fill",&fill,"fill/I");
  sca->Branch("t",&time,"t/D");
  sca->Branch("bx",&bx,"bx/I");
  char str2[16];
  for(Int_t i=0; i<8; i++) { sprintf(str,"hor%d",i); sprintf(str2,"%s/D",str); sca->Branch(str,&hor_shift[i],str2); };
  for(Int_t i=0; i<8; i++) { sprintf(str,"ver%d",i); sprintf(str2,"%s/D",str); sca->Branch(str,&ver_shift[i],str2); };
  for(Int_t i=0; i<8; i++) { sprintf(str,"trun%d",i); sprintf(str2,"%s/D",str); sca->Branch(str,&trun_shift[i],str2); };
  sca->Branch("front",&front_shift,"front/D");
  sca->Branch("back",&back_shift,"back/D");
  sca->Branch("goodTAC",&goodTAC_shift,"goodTAC/D");
  sca->Branch("tot_bx",&tot_bx,"tot_bx/D");
  sca->Branch("blue",&blue,"blue/I");
  sca->Branch("yell",&yell,"yell/I");
  //sca->Branch("kicked",&kicked,"kicked/O"); // deprecated
  sca->Branch("bxq",&bxq,"bxq/I"); // bx quality [0=good, 1=first bad, 2=afterpulse after bad (for npulse bx's)]
  sca->Branch("fq",&fq,"fq/I"); // fill quality [0=good, 1=modulation, 2=some slats out, 3=modulation+slats out, 4=very bad]


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
    /*
    bbce = bbc[1] + bbc[3] + bbc[5] + bbc[7]; // e + we + xe + xwe
    bbcw = bbc[2] + bbc[3] + bbc[6] + bbc[7]; // w + we + xw + xwe
    bbcx = bbc[3] + bbc[7]; // we + xwe

    zdce = zdc[1] + zdc[3] + zdc[5] + zdc[7]; // e + we + xe + xwe
    zdcw = zdc[2] + zdc[3] + zdc[6] + zdc[7]; // w + we + xw + xwe
    zdcx = zdc[3] + zdc[7]; // we + xwe
    
    vpde = vpd[1] + vpd[3] + vpd[5] + vpd[7]; // e + we + xe + xwe
    vpdw = vpd[2] + vpd[3] + vpd[6] + vpd[7]; // w + we + xw + xwe
    vpdx = vpd[3] + vpd[7]; // we + xwe
    */


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
      for(int x=0; x<8; x++)
      {
        /*
        hor_arr[x][index-1][(bx+113)%128] = hor[x];     // shift down 7 bXings (not needed in ZDC-SMD??)
        ver_arr[x][index-1][(bx+113)%128] = ver[x];
        trun_arr[x][index-1][(bx+113)%128] = trun[x];
        front_arr[index-1][(bx+113)%128] = front;
        back_arr[index-1][(bx+113)%128] = back;
        goodTAC_arr[index-1][(bx+113)%128] = goodTAC;
        */
        hor_arr[x][index-1][bx] = hor[x];     // no shift
        ver_arr[x][index-1][bx] = ver[x];
        trun_arr[x][index-1][bx] = trun[x];
        front_arr[index-1][bx] = front;
        back_arr[index-1][bx] = back;
        goodTAC_arr[index-1][bx] = goodTAC;
      };
    }
    /*
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
    */
    else
    {
      for(int x=0; x<8; x++)
      {
        hor_arr[x][index-1][bx] = hor[x];     // no shift
        ver_arr[x][index-1][bx] = ver[x];
        trun_arr[x][index-1][bx] = trun[x];
        front_arr[index-1][bx] = front;
        back_arr[index-1][bx] = back;
        goodTAC_arr[index-1][bx] = goodTAC;
      };
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


  // BXING QA
  // ---------------------------------------
  Int_t NF_tmp = acc->GetMaximum("fi") + 1;
  const Int_t NF = NF_tmp;
  Int_t bxq_arr[NF][120]; 
  Int_t fq_arr[NF];
  for(int f=0; f<NF; f++) {
    fq_arr[f]=0;
    for(int g=0; g<120; g++) {
      bxq_arr[f][g]=0;
    };
  };

  // bad bXings
  bxq_arr[1][22] = 1; // [fi] [bx]
  bxq_arr[2][61] = 1;
  bxq_arr[2][78] = 1;
  bxq_arr[6][9] = 1;
  bxq_arr[8][27] = 1;
  bxq_arr[9][56] = 1;
  bxq_arr[14][23] = 1;
  bxq_arr[15][12] = 1;
  bxq_arr[15][19] = 1;
  bxq_arr[18][104] = 1;
  bxq_arr[21][75] = 1;
  bxq_arr[22][15] = 1;
  bxq_arr[26][26] = 1;
  bxq_arr[27][26] = 1;
  bxq_arr[28][96] = 1;
  bxq_arr[29][70] = 1;
  bxq_arr[30][23] = 1;
  bxq_arr[31][17] = 1;
  bxq_arr[32][108] = 1;
  bxq_arr[35][21] = 1;
  bxq_arr[35][28] = 1;
  bxq_arr[36][90] = 1;
  bxq_arr[36][101] = 1;
  bxq_arr[37][51] = 1;
  bxq_arr[38][78] = 1;
  bxq_arr[40][7] = 1;
  bxq_arr[44][78] = 1;
  bxq_arr[46][60] = 1;

  // bad fills
  fq_arr[1] = 2; // [fi]
  fq_arr[9] = 1;
  fq_arr[10] = 1;
  fq_arr[13] = 1;
  fq_arr[14] = 1;
  fq_arr[15] = 1;
  fq_arr[18] = 1;
  fq_arr[20] = 4;
  fq_arr[23] = 1;
  fq_arr[25] = 1;
  fq_arr[26] = 1;
  fq_arr[27] = 1;
  fq_arr[28] = 1;
  fq_arr[29] = 3;
  fq_arr[31] = 1;
  fq_arr[32] = 1;
  fq_arr[33] = 1;
  fq_arr[35] = 1;
  fq_arr[37] = 1;
  fq_arr[38] = 1;
  fq_arr[39] = 1;
  fq_arr[40] = 1;
  fq_arr[41] = 1;
  fq_arr[42] = 1;
  fq_arr[44] = 1;
  fq_arr[46] = 1;
  fq_arr[47] = 1;
  fq_arr[48] = 4;

  // mark few bXings after a bad bXing
  Int_t npulse=3; // number of additional bXings to mark
  for(int f=0; f<NF; f++) {
    for(int g=0; g<120; g++) {
      if(bxq_arr[f][g]==1) {
        for(int nn=1; nn<=npulse; nn++) {
          bxq_arr[f][(g+nn)%120] = 2;
        };
      };
    };
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
      for(int x=0; x<8; x++)
      {
        hor_shift[x] = hor_arr[x][i][b];
        ver_shift[x] = ver_arr[x][i][b];
        trun_shift[x] = trun_arr[x][i][b];
        front_shift = front_arr[i][b];
        back_shift = front_arr[i][b];
        goodTAC_shift = goodTAC_arr[i][b];
      };
      tot_bx = tot_bx_arr[i][b];
      blue = blue_arr[i][b];
      yell = yell_arr[i][b];
      kicked = kicked_arr[i][b];
      bxq = bxq_arr[fill_index][bx];
      fq = fq_arr[fill_index];
      sca->Fill();
    };
  };

  acc->Write("acc");
  sca->Write("sca");
  printf("%s written\n",outfile_n.Data());
};
      

