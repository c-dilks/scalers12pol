// builds sums.root, which is scaler counts summed
// for each run number
//
// --determines spin pattern and adds 3 branches to the tree which
//   state the spin patterns collided for each run
//
// --zeroAborts will filter out counts in abort gaps and
//   bunches listed as empty by cdev; this causes slope for
//   total bXings / bXing rate vs. run time to be less than unity
//   and therefore should be left as zeroAborts=false

void sumTree(const char * filename="counts.root")
{
  Bool_t zeroAborts = false;

  TFile * infile = new TFile(filename,"READ");
  TTree * str = (TTree*) infile->Get("sca");

  // read counts.root tree
  Int_t i,runnum,fi,fill,bx,blue,yell;
  Double_t t,bbce,bbcw,bbcx,zdce,zdcw,zdcx,vpde,vpdw,vpdx,tot_bx;
  str->SetBranchAddress("i",&i);
  str->SetBranchAddress("runnum",&runnum);
  str->SetBranchAddress("fi",&fi);
  str->SetBranchAddress("fill",&fill);
  str->SetBranchAddress("t",&t);
  str->SetBranchAddress("bx",&bx);
  str->SetBranchAddress("bbce",&bbce);
  str->SetBranchAddress("bbcw",&bbcw);
  str->SetBranchAddress("bbcx",&bbcx);
  str->SetBranchAddress("zdce",&zdce);
  str->SetBranchAddress("zdcw",&zdcw);
  str->SetBranchAddress("zdcx",&zdcx);
  str->SetBranchAddress("vpde",&vpde);
  str->SetBranchAddress("vpdw",&vpdw);
  str->SetBranchAddress("vpdx",&vpdx);
  str->SetBranchAddress("tot_bx",&tot_bx);
  str->SetBranchAddress("blue",&blue);
  str->SetBranchAddress("yell",&yell);


  // get number of runs in each fill (used for averaging things over runs in a fill)
  Int_t fi_max_tmp = str->GetMaximum("fi");
  Int_t i_max_tmp = str->GetMaximum("i");
  const Int_t fi_max = fi_max_tmp;
  const Int_t i_max = i_max_tmp;
  Int_t num_runs[fi_max];
  for(Int_t n=0; n<fi_max; n++) num_runs[n]=0;
  for(Int_t n=0; n<str->GetEntries(); n++)
  {
    str->GetEntry(n);
    if(bx==0) num_runs[fi-1]+=1;
  };
  //for(Int_t n=0; n<fi_max; n++) printf("%d - %d\n",n,num_runs[n]);



  // spin pattern definitions
  // (overall pattern no.) = 10 * (blue pattern no.) + (yell pattern no.)
  Int_t pattern[4][8]; // [pattern no.] [bXing]

  pattern[0][0] =  1; // pattern 1 + - + - - + - + 
  pattern[0][1] = -1;
  pattern[0][2] =  1;
  pattern[0][3] = -1;
  pattern[0][4] = -1;
  pattern[0][5] =  1;
  pattern[0][6] = -1;
  pattern[0][7] =  1;

  pattern[1][0] = -1; // pattern 2 - + - + + - + - 
  pattern[1][1] =  1;
  pattern[1][2] = -1;
  pattern[1][3] =  1;
  pattern[1][4] =  1;
  pattern[1][5] = -1;
  pattern[1][6] =  1;
  pattern[1][7] = -1;

  pattern[2][0] =  1; // pattern 3 + + - - + + - - 
  pattern[2][1] =  1;
  pattern[2][2] = -1;
  pattern[2][3] = -1;
  pattern[2][4] =  1;
  pattern[2][5] =  1;
  pattern[2][6] = -1;
  pattern[2][7] = -1;

  pattern[3][0] = -1; // pattern 4 - - + + - - + + 
  pattern[3][1] = -1;
  pattern[3][2] =  1;
  pattern[3][3] =  1;
  pattern[3][4] = -1;
  pattern[3][5] = -1;
  pattern[3][6] =  1;
  pattern[3][7] =  1;



  // fill sum tree
  Int_t num_runs_cur;
  Double_t bbce_tot,bbcw_tot,bbcx_tot,zdce_tot,zdcw_tot;
  Double_t zdcx_tot,vpde_tot,vpdw_tot,vpdx_tot,tot_bx_tot;
  Double_t tau; // total bXings / bXing rate
  Int_t pattern_no; 
  Int_t cnt_blue_pos,cnt_blue_neg;
  Int_t cnt_yell_pos,cnt_yell_neg;
  Int_t cnt_0,cnt_1,cnt_2,cnt_3;
  cnt_blue_pos=cnt_blue_neg=cnt_yell_pos=cnt_yell_neg=0;
  cnt_0=cnt_1=cnt_2=cnt_3=0;
  TFile * outfile = new TFile("sums.root","RECREATE");
  TTree * sum = new TTree("sum","sum");
  sum->Branch("i",&i,"i/I");
  sum->Branch("runnum",&runnum,"runnum/I");
  sum->Branch("fi",&fi,"fi/I");
  sum->Branch("fill",&fill,"fill/I");
  sum->Branch("t",&t,"t/D");
  sum->Branch("tau",&tau,"tau/D");
  sum->Branch("num_runs",&num_runs_cur,"num_runs/I"); // total no. runs in a fill
  sum->Branch("bbce",&bbce_tot,"bbce/D");
  sum->Branch("bbcw",&bbcw_tot,"bbcw/D");
  sum->Branch("bbcx",&bbcx_tot,"bbcx/D");
  sum->Branch("zdce",&zdce_tot,"zdce/D");
  sum->Branch("zdcw",&zdcw_tot,"zdcw/D");
  sum->Branch("zdcx",&zdcx_tot,"zdcx/D");
  sum->Branch("vpde",&vpde_tot,"vpde/D");
  sum->Branch("vpdw",&vpdw_tot,"vpdw/D");
  sum->Branch("vpdx",&vpdx_tot,"vpdx/D");
  sum->Branch("tot_bx",&tot_bx_tot,"tot_bx/D");
  sum->Branch("pattern",&pattern_no,"pattern/I"); // overall pattern no.
  sum->Branch("cnt_blue_pos",&cnt_blue_pos,"cnt_blue_pos/I"); // no. blue=1
  sum->Branch("cnt_blue_neg",&cnt_blue_neg,"cnt_blue_neg/I"); // no. blue=-1
  sum->Branch("cnt_yell_pos",&cnt_yell_pos,"cnt_yell_pos/I"); // no. yell=1
  sum->Branch("cnt_yell_neg",&cnt_yell_neg,"cnt_yell_neg/I"); // no. yell=-1
  sum->Branch("cnt_0",&cnt_0,"cnt_0/I"); // no. - - collisions
  sum->Branch("cnt_1",&cnt_1,"cnt_1/I"); // no. - + collisions
  sum->Branch("cnt_2",&cnt_2,"cnt_2/I"); // no. + - collisions
  sum->Branch("cnt_3",&cnt_3,"cnt_3/I"); // no. + + collisions

  // tree loop; this assumes tree from counts.root is ordered by
  // run number and further by bXing number
  Int_t blue_pattern_cnt[4]; // patern check counters
  Int_t yell_pattern_cnt[4];
  Int_t blue_mul,yell_mul; // multipliers to determine pattern no.
  Int_t blue_err,yell_err; // error catching
  for(Int_t n=0; n<str->GetEntries(); n++)
  {
    str->GetEntry(n);
    num_runs_cur = num_runs[fi-1];

    if(bx==0)
    {
      // reset pattern check counters
      pattern_no=-1;
      for(Int_t p=0; p<4; p++)
      {
        blue_pattern_cnt[p] = 0;
        yell_pattern_cnt[p] = 0;
      };
    };

    if(bx<8)
    {
      // iterate pattern check counters
      for(Int_t p=0; p<4; p++)
      {
        if(blue == pattern[p][bx]) blue_pattern_cnt[p]++;
        if(yell == pattern[p][bx]) yell_pattern_cnt[p]++;
      };
    };

    // determine pattern numbers
    if(bx==8)
    {
      blue_mul=0;
      yell_mul=0;
      blue_err=0;
      yell_err=0;
      for(Int_t p=0; p<4; p++)
      {
        if(blue_pattern_cnt[p] == 8) 
        {
          blue_mul = p+1;
          blue_err++;
        };
        if(yell_pattern_cnt[p] == 8) 
        {
          yell_mul = p+1;
          yell_err++;
        };
      };
      if(blue_err>1 || yell_err>1) pattern_no=99; // more than one pattern matched
      else if(blue_err==0 || yell_err==0) pattern_no=0; // no pattern matched
      else pattern_no = 10*blue_mul+yell_mul;
    };

    // increment counters
    if(blue==-1 && yell==-1)
    {
      cnt_blue_neg++;
      cnt_yell_neg++;
      cnt_0++;
    }
    else if(blue==-1 && yell==1)
    {
      cnt_blue_neg++;
      cnt_yell_pos++;
      cnt_1++;
    }
    else if(blue==1 && yell==-1)
    {
      cnt_blue_pos++;
      cnt_yell_neg++;
      cnt_2++;
    }
    else if(blue==1 && yell==1)
    {
      cnt_blue_pos++;
      cnt_yell_pos++;
      cnt_3++;
    };


    // increment scalers
    if(!zeroAborts || (blue!=0 && yell!=0))
    {
      bbce_tot += bbce;
      bbcw_tot += bbcw;
      bbcx_tot += bbcx;
      zdce_tot += zdce;
      zdcw_tot += zdcw;
      zdcx_tot += zdcx;
      vpde_tot += vpde;
      vpdw_tot += vpdw;
      vpdx_tot += vpdx;
      tot_bx_tot += tot_bx;
    };

    // fill sum tree
    if(bx==119)
    {
      tau = tot_bx_tot/(9.382e6);
      sum->Fill();
      bbce_tot = 0;
      bbcw_tot = 0;
      bbcx_tot = 0;
      zdce_tot = 0;
      zdcw_tot = 0;
      zdcx_tot = 0;
      vpde_tot = 0;
      vpdw_tot = 0;
      vpdx_tot = 0;
      tot_bx_tot = 0;
      cnt_blue_pos=cnt_blue_neg=cnt_yell_pos=cnt_yell_neg=0;
      cnt_0=cnt_1=cnt_2=cnt_3=0;
    };
  };

  sum->Write();
  //sum->Scan("fi:i:num_runs");
  sum->Print();
};


    


