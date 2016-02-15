#include<stdio.h>
#include<stdlib.h>

#include"./myscaler.h"

#define NVAULT 120

int sca_read_bin(int debug, FILE *infile, int *num, int *chn, unsigned long long *val);

int try_scaler(char *filename, struct sca_run *scap)
{
  int run,board,unixTime;
  int list;

  int num, chn;
  int chnBBC,chnZDC,chnFZDC,chnVPD; 
  int clk;

  int bbc,zdc,fzdc,vpd;

  unsigned long long val;
  unsigned long long strobeCrt;

  int i,j;
  int cntr=0;
  FILE *fp;


  if((fp=fopen(filename,"r"))==NULL)
  {
    printf("can't open %s\n",filename);
    return 0;
  }
  printf("%s\n",filename);

  //    sscanf(filename,"../../pp500_production_2012/pp500_production_2012_list_%d/run%d_%d_%d.sca",&list,&run,&board,&unixTime);
  sscanf(filename,"../pp500_production_2012/run%d_%d_%d.sca",&run,&board,&unixTime);

  scap->run=run;
  scap->board=board;
  scap->unixTime=unixTime;

  if(board==12)
  {
    for(i=0;i<NVAULT;i++)
    {
      for(j=0;j<8;j++)
      {
        scap->scabx[i].valbbc[j]=0;
        scap->scabx[i].valzdc[j]=0;
        scap->scabx[i].valfzdc[j]=0;
        scap->scabx[i].valvpd[j]=0;
      }
      scap->scabx[i].val=0;
    }

    strobeCrt=0;
    cntr=0;

    while(sca_read_bin(0,fp,&num,&chn,&val)==1)
    {
      if(num==1&&chn==0)
      {
        strobeCrt=val+strobeCrt;
        scap->strobe=strobeCrt;
        cntr++;
      }
      else
      {
        clk=(chn>>17)&0x7F;
        chn=(chn&0x1FFFF);

        chnBBC=(chn&0x1C000);

        chnZDC=(chn&0xE00);

        chnFZDC=(chn&0x3800);

        chnVPD=(chn&0x7);

        switch(chnBBC)
        {
          case 0x0:
            bbc=0;
            break;
          case 0x04000:
            bbc=1;
            break;
          case 0x08000:
            bbc=2;
            break;
          case 0x0C000:
            bbc=3;
            break;
          case 0x10000:
            bbc=4;
            break;
          case 0x14000:
            bbc=5;
            break;
          case 0x18000:
            bbc=6;
            break;
          case 0x1C000:
            bbc=7;
            break;
        }

        switch(chnZDC)
        {
          case 0x0 :
            zdc=0;
            break;
          case 0x200:
            zdc=1;
            break;
          case 0x400:
            zdc=2;
            break;
          case 0x600:
            zdc=3;
            break;
          case 0x800:
            zdc=4;
            break;
          case 0xA00:
            zdc=5;
            break;
          case 0xC00:
            zdc=6;
            break;
          case 0xE00:
            zdc=7;
            break;
        }

        switch(chnFZDC)
        {
          case 0x0 :
            fzdc=0;
            break;
          case 0x1000:
            fzdc=1;
            break;
          case 0x2000:
            fzdc=2;
            break;
          case 0x3000:
            fzdc=3;
            break;
          case 0x800:
            fzdc=4;
            break;
          case 0x1800:
            fzdc=5;
            break;
          case 0x2800:
            fzdc=6;
            break;
          case 0x3800:
            fzdc=7;
            break;
        }

        switch(chnVPD)
        {
          case 0x0:
            vpd=0;
            break;
          case 0x1:
            vpd=1;
            break;
          case 0x2:
            vpd=2;
            break;
          case 0x3:
            vpd=3;
            break;
          case 0x4:
            vpd=4;
            break;
          case 0x5:
            vpd=5;
            break;
          case 0x6:
            vpd=6;
            break;
          case 0x7:
            vpd=7;
            break;
        }
        if(clk>=0&&clk<=119)
        {
          scap->scabx[clk].valbbc[bbc]+=val;
          scap->scabx[clk].valzdc[zdc]+=val;
          scap->scabx[clk].valfzdc[fzdc]+=val;
          scap->scabx[clk].valvpd[vpd]+=val;
          scap->scabx[clk].val+=val;
        }
      }
    }

    fclose(fp);

    //      printf("%llu\n",scap->scabx[117].valvpd[1]);

    printf("%s is generated, strobecounter=%llu, counter=%d\n",filename,strobeCrt,cntr);
  }

  return 1;
}
