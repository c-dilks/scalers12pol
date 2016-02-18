// SCALER READER FOR ZDC POLARIMETRY

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sca_read.h"

int main(int argc, char * argv[])
{
  int debug = 0;

  // open scaler file
  FILE * fp;
  char filename[64];
  if(argc!=2)
  {
    fprintf(stderr,"Usage: sca_read_bin.o [*.sca filename]\n");
    return 0;
  }
  else
  {
    strcpy(filename,argv[1]);
    if((fp=fopen(filename,"r"))==NULL)
    {
      fprintf(stderr,"ERROR: not able to open %s\n",filename);
      return 0;
    };
  };
  int runnum, boardnum, unixtime;
  sscanf(filename,"sca2012/run%d_%d_%d.sca",&runnum,&boardnum,&unixtime);
  printf("reading %s\n",filename);
  printf("runnum=%d\nboard=%d\n",runnum,boardnum);


  // initialise scaler counters
  unsigned long long HorZDC[8][128]; // [3 bits --> integer between 0-7]
  unsigned long long VerZDC[8][128]; // [3 bits --> integer between 0-7]
  unsigned long long TrunZDC[8][128]; // [3 bits --> integer between 0-7]
  unsigned long long FrontZDC[128];
  unsigned long long BackZDC[128];
  unsigned long long GoodTAC[128];
  unsigned long long bx_cnt[120];
  int b, d;
  for(b=0; b<120; b++)
  {
    for(d=0; d<8; d++)
    {
      HorZDC[d][b]=0;
      VerZDC[d][b]=0;
      TrunZDC[d][b]=0;
    };
    bx_cnt[b]=0;
    FrontZDC[b]=0;
    BackZDC[b]=0;
    GoodTAC[b]=0;
  };

  
  // read scaler file and increment counters
  // -- see 24bit.pdf
  int num, chn, bx;
  int cHor, cVer, cTrun;
  int j;
  unsigned long long val;
  while(sca_read_bin(0, fp, &num, &chn, &val)==1)
  {
    if(!(num==1 && chn==0))
    {
      bx = (chn>>17) & 0x7F; // mask bits 17-23
      chn = chn & 0xFFF;   // mask bits 0-11


      if(bx>=0 && bx<=119)
      {
        HorZDC[ chn & 0x7 ][bx] += val; // bits 0-2
        VerZDC[ (chn>>3) & 0x7 ][bx] += val; // bits 3-5
        TrunZDC[ (chn>>6) & 0x7 ][bx] += val; // bits 6-8
        FrontZDC[bx] += ((chn>>9) & 0x1)==1 ? val:0; // bit 9
        BackZDC[bx] += ((chn>>10) & 0x1)==1 ? val:0; // ...
        GoodTAC[bx] = ((chn>>11) & 0x1)==1 ? 1:0; // ...

        if(debug)
          printf("[DEBUG] bx=%d val=%lld tac=%lld\n",bx,val,GoodTAC[bx]);

        bx_cnt[bx]+=val;
      };
    };
  };

  fclose(fp);


  // print scaler counts to datfile
  char datfile[128];
  sprintf(datfile,"datfiles/run%d_%d.dat",runnum,boardnum);
  FILE * fo;
  fo=fopen(datfile,"w");
  for(b=0; b<120; b++)
  {
    fprintf(fo,"%d",b);
    for(d=0; d<8; d++) fprintf(fo," %lld",HorZDC[d][b]);
    for(d=0; d<8; d++) fprintf(fo," %lld",VerZDC[d][b]);
    for(d=0; d<8; d++) fprintf(fo," %lld",TrunZDC[d][b]);
    fprintf(fo," %lld",FrontZDC[b]);
    fprintf(fo," %lld",BackZDC[b]);
    fprintf(fo," %lld",GoodTAC[b]);
    fprintf(fo," %lld\n",bx_cnt[b]);
  };
  fclose(fo);


  return 1;
}
