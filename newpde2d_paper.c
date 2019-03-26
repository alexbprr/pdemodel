#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "sds.h"
#include "sdsalloc.h"
#include "hashmap.h"

#define KEY_MAX_LENGTH (100)
#define KEY_COUNT (100)

typedef struct data_struct_s
{
    char key_string[KEY_MAX_LENGTH];
    float value;
} data_struct_t;

map_t pmap; //Map parameter name to its value

float deltaT  = 0.;
float days = 0;
int T = 0;
float deltaY = 0.;
float deltaX = 0.;
//#define size_ (int)(10/deltaX) + 1
int size_  = 0;
float maxDensity = 0;

//Bacteria parameters
float r_b = 0; // 1.3// 0.8
float m_b = 0.;// 0.
float b_kdif = 0.;// 0.1 //0.5, 0.1
float b_nexp = 0;// 5
float D_b = 0.; //G diffusion
float D2_b = 0.; //Nonlinear diffusion
float D3_b = 0.; //Classic diffusion
float D_f = 0.; //G diffusion
float D2_f = 0.; //Nonlinear diffusion
float D3_f = 0.; //Classic diffusion

float B0 = 0;
float B02 = 0.;
float F0 = 0;
float N0 = 0;

float k = 0;
//Coa, vWbp and fibrin parameters
float p = 0;//1.5; //p=1, wbfib=1.5 -> eliminacao da bacteria
float m_f = 0.;
float X_f = 0.;
float f_kdif = 0.; //0.5, 0.1
float f_nexp = 0.;
//float lambda_ = 4.;
//float phi = 1.;

float X_fib = 0.;

//Neutrophil parameters
float s = 0;
float l = 0;
float m_n = 0.;
float D_n = 0.;
float X_n = 0.;
float s2 = 0.;
float s3 = 0.;
float alpha = 0.;

//Toxin parameters
float betato = 0.;
float m_to = 0.;
float D_to = 0.;

//Dead neutrophil parameters
float alphato = 0.;
float gammand =0.;

//Weights
//Influence factors for bacteria
float wbb = 0.;
float wfb = 0.;//4.
float wnb = 0;
float wfibb = 0.;

//Influence factors for coa/vWbp
float wbf =0;//1.5
float wff =0.;
float wnf =0.;
float wfibf = 0;
//Influence factors for fibrin
float wbfib = 0.;
float wffib =0.;
float wnfib =0.;
float wfibfib =0.;
//Influence factors for neutrophil
float wbn =0;
float wfn =0.;
float wnn =0.;
float wfibn =0.;

//Influence factors for toxin
float wbto =0.;
float wfto =0.;
float wnto =0.;
float wfibto =0.;
float wtoto =0.;

//Influence factors diffusion
float wbb_dif = 0;
float wfb_dif = 0;
float wnb_dif = 0;
float wfibb_dif = 0;

float wff_dif = 0;
float wbf_dif = 0;
float wnf_dif = 0;
float wfibf_dif = 0;

float wfn_dif = 0;
float wbn_dif = 0;
float wnn_dif = 0;
float wfibn_dif = 0;

//mr
float wbmr = 0;
float wfmr = 0;
float wnmr = 0;
float wfibmr = 0;
float wmrmr = 0;
float wbmr_dif = 0;
float wfibmr_dif = 0;
float wndmr_dif = 0;
float smr = 0;
float D_mr = 0;
float lmr = 0;
float X_mr = 0;
//ma
float wbma = 0;
float wfma = 0;
float wnma = 0;
float wfibma = 0;
float wmama = 0;
float wbma_dif = 0;
float wfibma_dif = 0;
float wndma_dif = 0;
float D_ma = 0;
float lma = 0;
float X_ma = 0;
//damage tissue
float b_dmg = 0;
float nd_dmg = 0;
float gamma_dmg = 0;
float wdn = 0;
float wdn_dif = 0;

//Migration due to inflammation
float s_inf = 0;


float verifyMaxDensity(float value)
{
     if (value > maxDensity)
         value = maxDensity;
     return value;
}

float verifyDensity(float value)
{
    if (value < 0.) return 0.;
    else return value;
}

float g(float v)
{
    return 1.0 - (v/(float)maxDensity);
}

float gDiffusion(float v)
{
    return (1.0 - (v/(float)maxDensity)) < 0.? 0.: (1.0 - (v/(float)maxDensity));
}

float diffusionX(float*** u, int x, int y, int tAntigo)
{
    float resX;
    if (x == 0)
    {
        resX = (u[tAntigo][x+1][y] - u[tAntigo][x][y]) /(deltaX*deltaX);
    }
    else if (x == (size_ -1))
    {
        resX = (u[tAntigo][x-1][y] - u[tAntigo][x][y])/(deltaX*deltaX);
    }
    else
    {
        resX = (u[tAntigo][x+1][y] - 2*u[tAntigo][x][y] + u[tAntigo][x-1][y])/(deltaX*deltaX); 	//dentro do dominio mas fora da extremidade
    }
    return resX;
}

float diffusionY(float*** u, int x, int y, int tAntigo)
{
    float resY;
    if (y == 0)
    {
        resY = (u[tAntigo][x][y+1] - u[tAntigo][x][y]) /(deltaY*deltaY);
    }
    else if (y == (size_ -1))
    {
        resY = (u[tAntigo][x][y-1] - u[tAntigo][x][y])/(deltaY*deltaY);
    }
    else
    {
        resY = (u[tAntigo][x][y+1] - 2*u[tAntigo][x][y] + u[tAntigo][x][y-1])/(deltaY*deltaY); 	//dentro do dominio mas fora da extremidade
    }
    return resY;
}

float diffusion(float*** u, int x, int y, int tAntigo)
{
  return diffusionX(u,x,y,tAntigo) + diffusionY(u,x,y,tAntigo);
}

void imprimeDados(FILE* arq, float*** u, int tAtual)
{
    for (int x = 0; x < size_; x++)
    {
      for (int y = 0; y < size_; y++)
        fprintf(arq, "%f ", u[tAtual][x][y]);
      fprintf(arq,"\n");
    }
    fprintf(arq,"--");//A cada vez que imprimeDados é chamada, o separador -- é gravado após a matriz do tempo atual
    fprintf(arq,"\n");
}

float chemotaxis(float*** u, float ***ch, int x, int y, float*** w, int tAntigo)
{
    float resX = 0, resY = 0, flux_left = 0, flux_right = 0;
    if (x > 0)
    {
        if((ch[tAntigo][x][y] -  ch[tAntigo][x-1][y])>0)
        {
            flux_left = -((ch[tAntigo][x][y] - ch[tAntigo][x-1][y]) * u[tAntigo][x-1][y] *gDiffusion(w[tAntigo][x-1][y]))/ deltaX;
        }
        else
        {
            flux_left = -((ch[tAntigo][x][y] - ch[tAntigo][x-1][y]) * u[tAntigo][x][y] *gDiffusion(w[tAntigo][x][y]))/ deltaX;
        }
    }
    if(x < (size_-1))
    {
        if((ch[tAntigo][x+1][y] - ch[tAntigo][x][y]) > 0)
        {
            flux_right = ((ch[tAntigo][x+1][y] - ch[tAntigo][x][y]) * u[tAntigo][x][y]*gDiffusion(w[tAntigo][x][y])) / deltaX;
        }
        else
        {
            flux_right = ((ch[tAntigo][x+1][y] - ch[tAntigo][x][y]) * u[tAntigo][x+1][y]*gDiffusion(w[tAntigo][x+1][y])) / deltaX;
        }
    }
    resX = (flux_left + flux_right)/deltaX;

    if (y > 0)
    {
        if((ch[tAntigo][x][y] -  ch[tAntigo][x][y-1])>0)
        {
            flux_left = -((ch[tAntigo][x][y] - ch[tAntigo][x][y-1]) * u[tAntigo][x][y-1] *gDiffusion(w[tAntigo][x][y-1]))/ deltaY;
        }
        else
        {
            flux_left = -((ch[tAntigo][x][y] - ch[tAntigo][x][y-1]) * u[tAntigo][x][y] *gDiffusion(w[tAntigo][x][y]))/ deltaY;
        }
    }
    if(y < (size_-1))
    {
        if((ch[tAntigo][x][y+1] - ch[tAntigo][x][y]) > 0)
        {
            flux_right = ((ch[tAntigo][x][y+1] - ch[tAntigo][x][y]) * u[tAntigo][x][y]*gDiffusion(w[tAntigo][x][y])) / deltaY;
        }
        else
        {
            flux_right = ((ch[tAntigo][x][y+1] - ch[tAntigo][x][y]) * u[tAntigo][x][y+1]*gDiffusion(w[tAntigo][x][y+1])) / deltaY;
        }
    }
    resY = (flux_left + flux_right)/deltaY;
    return resX + resY;
}

float gRightBorderMediumX(int x, int y, float*** w, int t)
{
    return gDiffusion((w[t][x+1][y] + w[t][x][y])/2);
}

float gRightBorderMediumY(int x, int y, float*** w, int t)
{
    return gDiffusion((w[t][x][y+1] + w[t][x][y])/2);
}

float gLeftBorderMediumX(int x, int y, float*** w, int t)
{
    return gDiffusion((w[t][x][y] + w[t][x-1][y])/2);
}

float gLeftBorderMediumY(int x, int y, float*** w, int t)
{
    return gDiffusion((w[t][x][y] + w[t][x][y-1])/2);
}

float uRightBorderX(float*** u, int x, int y, int t)
{
    return u[t][x+1][y] - u[t][x][y];
}

float uRightBorderY(float*** u, int x, int y, int t)
{
    return u[t][x][y+1] - u[t][x][y];
}

float uLeftBorderX(float*** u, int x, int y, int t)
{
    return u[t][x][y] - u[t][x-1][y];
}

float uLeftBorderY(float*** u, int x, int y, int t)
{
    return u[t][x][y] - u[t][x][y-1];
}

float uRightBorderMediumX(float*** u, int x, int y, int t)
{
    return (u[t][x+1][y] + u[t][x][y])/2;
}

float uRightBorderMediumY(float*** u, int x, int y, int t)
{
    return (u[t][x][y+1] + u[t][x][y])/2;
}

float uLeftBorderMediumX(float*** u, int x, int y, int t)
{
    return (u[t][x][y] + u[t][x-1][y])/2;
}

float uLeftBorderMediumY(float*** u, int x, int y, int t)
{
    return (u[t][x][y] + u[t][x][y-1])/2;
}

float localAverageX(float*** u, int x, int y, float*** w, int t)
{
    if(x == 0)
        return (gRightBorderMediumX(x, y, w, t)*(uRightBorderX(u, x, y, t))) /(deltaX*deltaX);
    else if(x == size_-1)
        return - (gLeftBorderMediumX(x, y, w, t)*(uLeftBorderX(u, x, y, t))) /(deltaX*deltaX);
    else
        return (gRightBorderMediumX(x, y, w, t)*uRightBorderX(u, x, y, t) -  gLeftBorderMediumX(x, y, w, t)*
        uLeftBorderX(u, x, y, t) ) /(deltaX*deltaX);
}

float localAverageY(float*** u, int x, int y, float*** w, int t)
{
    if(y == 0)
        return (gRightBorderMediumY(x, y, w, t)*(uRightBorderY(u, x, y, t))) /(deltaY*deltaY);
    else if(y == size_-1)
        return - (gLeftBorderMediumY(x, y, w, t)*(uLeftBorderY(u, x, y, t))) /(deltaY*deltaY);
    else
        return (gRightBorderMediumY(x, y, w, t)*uRightBorderY(u, x, y, t) -  gLeftBorderMediumY(x, y, w, t)*
        uLeftBorderY(u, x, y, t) ) /(deltaY*deltaY);
}

float localAverage(float*** u, int x, int y, float*** w, int t)
{
  return localAverageX(u,x,y,w,t) + localAverageY(u,x,y,w,t);
}

float expURightBorderMediumX(float*** u, int x, int y, int t, int nexp)
{
    return pow(uRightBorderMediumX(u, x, y, t), nexp);
}

float expURightBorderMediumY(float*** u, int x, int y, int t, int nexp)
{
    return pow(uRightBorderMediumY(u, x, y, t), nexp);
}

float expULeftBorderMediumX(float*** u, int x, int y, int t, int nexp)
{
    return pow(uLeftBorderMediumX(u, x, y, t), nexp);
}

float expULeftBorderMediumY(float*** u, int x, int y, int t, int nexp)
{
    return pow(uLeftBorderMediumY(u, x, y, t), nexp);
}

float nonLinearFuncRightX(float*** v, int x, int y, int t, float kdif, int nexp)
{
    return (kdif + 1)*expURightBorderMediumX(v, x, y, t, nexp)*(1/(kdif + expURightBorderMediumX(v, x, y, t, nexp)));
}

float nonLinearFuncRightY(float*** v, int x, int y, int t, float kdif, int nexp)
{
    return (kdif + 1)*expURightBorderMediumY(v, x, y, t, nexp)*(1/(kdif + expURightBorderMediumY(v, x, y, t, nexp)));
}

float nonLinearFuncLeftX(float*** v, int x, int y, int t, float kdif, int nexp)
{
    return (kdif + 1)*expULeftBorderMediumX(v, x, y, t, nexp)*(1/(kdif + expULeftBorderMediumX(v, x, y, t, nexp)));
}

float nonLinearFuncLeftY(float*** v, int x, int y, int t, float kdif, int nexp)
{
    return (kdif + 1)*expULeftBorderMediumY(v, x, y, t, nexp)*(1/(kdif + expULeftBorderMediumY(v, x, y, t, nexp)));
}

float flocalAverageX(float*** u, int x, int y, float*** w, int t, float*** b, float kdif, int nexp)
{
    if(x == 0)
        return (gRightBorderMediumX(x, y, w, t)*nonLinearFuncRightX(b, x, y, t, kdif, nexp)*(uRightBorderX(u, x, y, t))) /(deltaX*deltaX);
    else if(x == size_-1)
        return - (gLeftBorderMediumX(x, y, w, t)*nonLinearFuncLeftX(b, x, y, t, kdif, nexp)*(uLeftBorderX(u, x, y, t))) /(deltaX*deltaX);
    else
        return ( (gRightBorderMediumX(x, y, w, t)*nonLinearFuncRightX(b, x, y, t, kdif, nexp)*(uRightBorderX(u, x, y, t))) -
                 (gLeftBorderMediumX(x, y, w, t)*nonLinearFuncLeftX(b, x, y, t, kdif, nexp)*(uLeftBorderX(u, x, y, t))) ) /(deltaX*deltaX);
}

float flocalAverageY(float*** u, int x, int y, float*** w, int t, float*** b, float kdif, int nexp)
{
    if(y == 0)
        return (gRightBorderMediumY(x, y, w, t)*nonLinearFuncRightY(b, x, y, t, kdif, nexp)*(uRightBorderY(u, x, y, t))) /(deltaY*deltaY);
    else if(y == size_-1)
        return - (gLeftBorderMediumY(x, y, w, t)*nonLinearFuncLeftY(b, x, y, t, kdif, nexp)*(uLeftBorderY(u, x, y, t))) /(deltaY*deltaY);
    else
        return ( (gRightBorderMediumY(x, y, w, t)*nonLinearFuncRightY(b, x, y, t, kdif, nexp)*(uRightBorderY(u, x, y, t))) -
                 (gLeftBorderMediumY(x, y, w, t)*nonLinearFuncLeftY(b, x, y, t, kdif, nexp)*(uLeftBorderY(u, x, y, t))) ) /(deltaY*deltaY);
}

float flocalAverage(float*** u, int x, int y, float*** w, int t, float*** b, float kdif, int nexp)
{
  return flocalAverageX(u,x,y,w,t,b,kdif,nexp) + flocalAverageY(u,x,y,w,t,b,kdif,nexp);
}

void updateGWReaction(float*** b, float*** f, float*** n, float*** fib, float*** wb, float*** wf, float*** wn, float*** wfib,
  float*** gwb, float*** gwf, float*** gwn, float*** gwfib, float*** to, float*** wto, float*** gwto, float*** wmr, float*** gwmr,
  float*** wma, float*** gwma , float*** ma, float*** mr,float*** dmt,int x, int y, int t)
{
    wb[t][x][y] = wbb*b[t][x][y] + wfb*f[t][x][y] + wnb*n[t][x][y] + wfibb*fib[t][x][y]; //*b[t][x][y];
    gwb[t][x][y] = g(wb[t][x][y]);

    wf[t][x][y] = wbf*b[t][x][y] + wff*f[t][x][y] + wnf*n[t][x][y] + wfibf*fib[t][x][y];
    gwf[t][x][y] = g(wf[t][x][y]);

    wfib[t][x][y] = wbfib*b[t][x][y] + wffib*f[t][x][y] + wnfib*n[t][x][y] + wfibfib*fib[t][x][y];
    gwfib[t][x][y] = g(wfib[t][x][y]);

    wn[t][x][y] = wbn*b[t][x][y] + wfn*f[t][x][y] + wnn*n[t][x][y] + wfibn*fib[t][x][y] + wdn*dmt[t][x][y];
    gwn[t][x][y] = g(wn[t][x][y]);

    wto[t][x][y] = wbto*b[t][x][y] + wfto*f[t][x][y] + wnto*n[t][x][y] + wfibto*fib[t][x][y] + wtoto*to[t][x][y];
    gwto[t][x][y] = g(wto[t][x][y]);
//---
    wmr[t][x][y] = wbmr*b[t][x][y] + wfmr*f[t][x][y] + wnmr*n[t][x][y] + wfibmr*fib[t][x][y] + wmrmr*mr[t][x][y];
    gwmr[t][x][y] = g(wmr[t][x][y]);

    wma[t][x][y] = wbma*b[t][x][y] + wfma*f[t][x][y] + wnma*n[t][x][y] + wfibma*fib[t][x][y] + wmama*ma[t][x][y];
    gwma[t][x][y] = g(wma[t][x][y]);

}


void updateWDiffusion(float*** b, float*** f, float*** n, float*** fib, float*** wbdif, float*** wfdif, float*** wndif, float*** wfibdif,
  float*** to, float*** wtodif, float*** wmrdif, float*** wmadif, float*** dmt, int x, int y, int t)
{
    wbdif[t][x][y] = wbb_dif*b[t][x][y] + wfb_dif*f[t][x][y] + wnb_dif*n[t][x][y] + wfibb_dif*fib[t][x][y];
    wbdif[t][x][y] = verifyMaxDensity(wbdif[t][x][y]);

    wfdif[t][x][y] = wbf_dif*b[t][x][y] + wff_dif*f[t][x][y] + wnf_dif*n[t][x][y] + wfibf_dif*fib[t][x][y];
    wfdif[t][x][y] = verifyMaxDensity(wfdif[t][x][y]);

    // wfibdif[t][x][y] = wbfib*b[t][x][y] + wnfib*n[t][x][y] + wfibfib*fib[t][x][y];
    //
    wndif[t][x][y] = wbn_dif*b[t][x][y] + wfn_dif*f[t][x][y] + wnn_dif*n[t][x][y] + wfibn_dif*fib[t][x][y] + wdn_dif*dmt[t][x][y];
    wndif[t][x][y] = verifyMaxDensity(wndif[t][x][y]);
    //
    // wtodif[t][x][y] = wbto*b[t][x][y] + wfto*f[t][x][y] + wnto*n[t][x][y] + wfibto*fib[t][x][y] + wtoto*to[t][x][y];
//---
    wmrdif[t][x][y] = wbmr_dif*b[t][x][y] + wfibmr_dif*fib[t][x][y]; // + wndmr_dif*nd[t][x][y];
    wmrdif[t][x][y] = verifyMaxDensity(wmrdif[t][x][y]);

    wmadif[t][x][y] = wbma_dif*b[t][x][y] + wfibma_dif*fib[t][x][y]; // + wndma_dif*nd[t][x][y];
    wmadif[t][x][y] = verifyMaxDensity(wmadif[t][x][y]);
}

void readAndSetParameters(FILE* fp)
{
  sds* tokens;
  int count;
  data_struct_t* hash_entry;
  char c[100];
  float fvalue;
  pmap = hashmap_new();
  //while (fscanf(fp, "%m[^=]=%ms", &key, &value) == 2) {
  while ((fgets (c, sizeof(c), fp)) != NULL)
  {
    printf("%s\n", c);
    sds line = sdsnew(c);
    sdstrim(line, " ");
    printf("%s\n", line);
    if (sdscmp(line,sdsnew("--\n")) == 0)
        break;
    tokens = sdssplitlen(line,sdslen(line),"=",1,&count);
    //printf("key: %s\n", tokens[0]);
    //printf("value: %s\n", tokens[1]);
    fflush(stdout);
    hash_entry = malloc(sizeof(data_struct_t));
    sprintf(hash_entry->key_string,"%s",tokens[0]);
    fvalue = atof(tokens[1]);
    hash_entry->value = fvalue;
    hashmap_put(pmap, hash_entry->key_string, hash_entry);
  }
}

void readAndSetInitialCondition(FILE* fp, float ***b, float ***f, float ***fib, float ***n, float ***mr, float ***ma, int t)
{
  int x, y, option, xinf, xsup, yinf, ysup;
  sds *tokens;
  int count, j, num;
  char c[100];
  long pos = 0;
  float value;
  while ((fgets (c, sizeof(c), fp)) != NULL)
  {
    printf("%s\n", c);
    sds line = sdsnew(c);
    sdstrim(line, " ");
    //tokens = sdssplitlen(line,sdslen(line),"\n",1,&count);
    tokens = sdssplitlen(line,sdslen(line),":",1,&count);
    //printf("tokens: %s\n", tokens[0]);
    //printf("count: %d\n", count);
    sds celltype = tokens[0];
    if (count >= 2)
    {
        tokens = sdssplitlen(tokens[1],sdslen(tokens[1]),",",1,&count);
        sscanf(tokens[0],"%d",&option);
        value = atof(tokens[1]);
        //printf("option: %d\n", option);
        if (option == 0)
        {
          //value = atof(tokens[1]);
          for (j = 2; j < count; j++)
          {
            sdstrim(tokens[j], "(");
            sdstrim(tokens[j], ")");
            if (j % 2 == 0)
            {
              sscanf(tokens[j], "%d", &x);
            }
            else
            {
              sscanf(tokens[j], "%d", &y);
              if (sdscmp(celltype,sdsnew("b0")) == 0)
                b[t][y][x] = value;
              else if (sdscmp(celltype,sdsnew("n0")) == 0)
                n[t][y][x] = value;
              else if (sdscmp(celltype,sdsnew("coa0")) == 0)
                f[t][y][x] = value;
              else if (sdscmp(celltype,sdsnew("fib0")) == 0)
                fib[t][y][x] = value;
              else if (sdscmp(celltype,sdsnew("mr0")) == 0)
                mr[t][y][x] = value;
              else if (sdscmp(celltype,sdsnew("ma0")) == 0)
                ma[t][y][x] = value;

              //printf("saída: (%d,%d)\n",x,y);
            }
          }
        }
        else if (option == 1)
        {

          for (j = 2; j < count; j++)
          {
            sdstrim(tokens[j], "(");
            sdstrim(tokens[j], ")");
            if (j % 4 == 2)
            {
              sscanf(tokens[j], "%d", &xinf);
            }
            else if (j % 4 == 3)
            {
              sscanf(tokens[j], "%d", &xsup);
            }
            else if (j % 4 == 0)
            {
              sscanf(tokens[j], "%d", &yinf);
            }
            else if (j % 4 == 1)
            {
              sscanf(tokens[j], "%d", &ysup);
              printf("x limits: %d - %d\n",xinf,xsup);
              printf("y limits: %d - %d\n",yinf,ysup);
              if (sdscmp(celltype,sdsnew("b0")) == 0)
                for (x = xinf; x <= xsup; x++)
                  for (y = yinf; y <= ysup; y++)
                    b[t][y][x] = value;
              else if (sdscmp(celltype,sdsnew("n0")) == 0)
                for (x = xinf; x <= xsup; x++)
                  for (y = yinf; y <= ysup; y++)
                    n[t][y][x] = value;
              else if (sdscmp(celltype,sdsnew("coa0")) == 0)
                for (x = xinf; x <= xsup; x++)
                  for (y = yinf; y <= ysup; y++)
                    f[t][y][x] = value;
              else if (sdscmp(celltype,sdsnew("fib0")) == 0)
                for (x = xinf; x <= xsup; x++)
                  for (y = yinf; y <= ysup; y++)
                    fib[t][y][x] = value;
              else if (sdscmp(celltype,sdsnew("mr0")) == 0)
                for (x = xinf; x <= xsup; x++)
                  for (y = yinf; y <= ysup; y++)
                    mr[t][y][x] = value;
              else if (sdscmp(celltype,sdsnew("ma0")) == 0)
                for (x = xinf; x <= xsup; x++)
                  for (y = yinf; y <= ysup; y++)
                    ma[t][y][x] = value;
            }
          }
        }
    }
    sdsfree(line);
  }
}

int main(int argc, char* argv[])
{
    int x=0, y=0, t=0, i=0, j=0;
    int numIteracoesGravadas = 10;
    FILE* entryfile;
    if ((entryfile = fopen(argv[1],"r")) == NULL)
    {
      printf("Error on file opening! \n");
      exit(1);
    }
    readAndSetParameters(entryfile);
    data_struct_t* hash_entry_temp = malloc(sizeof(data_struct_t));
    hashmap_get(pmap, "deltaT", (void**)(&hash_entry_temp));
    deltaT = hash_entry_temp->value;
    hashmap_get(pmap, "days", (void**)(&hash_entry_temp));
    days = hash_entry_temp->value;
    T = (float)days*((float)(1./deltaT));
    hashmap_get(pmap, "days", (void**)(&hash_entry_temp));
    days = hash_entry_temp->value;
    hashmap_get(pmap, "deltaX", (void**)(&hash_entry_temp));
    deltaX = hash_entry_temp->value;
    hashmap_get(pmap, "deltaY", (void**)(&hash_entry_temp));
    deltaY = hash_entry_temp->value;
    hashmap_get(pmap, "size_", (void**)(&hash_entry_temp));
    size_ = hash_entry_temp->value;
    hashmap_get(pmap, "maxDensity", (void**)(&hash_entry_temp));
    maxDensity = hash_entry_temp->value;
    hashmap_get(pmap, "r_b", (void**)(&hash_entry_temp));
    r_b = hash_entry_temp->value;
    hashmap_get(pmap, "b_kdif", (void**)(&hash_entry_temp));
    b_kdif = hash_entry_temp->value;
    hashmap_get(pmap, "b_nexp", (void**)(&hash_entry_temp));
    b_nexp = hash_entry_temp->value;
    hashmap_get(pmap, "D_b", (void**)(&hash_entry_temp));
    D_b = hash_entry_temp->value;
    hashmap_get(pmap, "D2_b", (void**)(&hash_entry_temp));
    D2_b = hash_entry_temp->value;
    hashmap_get(pmap, "k", (void**)(&hash_entry_temp));
    k = hash_entry_temp->value;
    hashmap_get(pmap, "p", (void**)(&hash_entry_temp));
    p = hash_entry_temp->value;
    hashmap_get(pmap, "m_f", (void**)(&hash_entry_temp));
    m_f = hash_entry_temp->value;
    hashmap_get(pmap, "X_f", (void**)(&hash_entry_temp));
    X_f = hash_entry_temp->value;
    hashmap_get(pmap, "f_kdif", (void**)(&hash_entry_temp));
    f_kdif = hash_entry_temp->value;
    hashmap_get(pmap, "f_nexp", (void**)(&hash_entry_temp));
    f_nexp = hash_entry_temp->value;
    hashmap_get(pmap, "D_f", (void**)(&hash_entry_temp));
    D_f = hash_entry_temp->value;
    hashmap_get(pmap, "D2_f", (void**)(&hash_entry_temp));
    D2_f = hash_entry_temp->value;
    hashmap_get(pmap, "s", (void**)(&hash_entry_temp));
    s = hash_entry_temp->value;
    hashmap_get(pmap, "s2", (void**)(&hash_entry_temp));
    s2 = hash_entry_temp->value;
    hashmap_get(pmap, "s3", (void**)(&hash_entry_temp));
    s3 = hash_entry_temp->value;
    hashmap_get(pmap, "l", (void**)(&hash_entry_temp));
    l = hash_entry_temp->value;
    hashmap_get(pmap, "m_n", (void**)(&hash_entry_temp));
    m_n = hash_entry_temp->value;
    hashmap_get(pmap, "D_n", (void**)(&hash_entry_temp));
    D_n = hash_entry_temp->value;
    hashmap_get(pmap, "X_n", (void**)(&hash_entry_temp));
    X_n = hash_entry_temp->value;
    hashmap_get(pmap, "alpha", (void**)(&hash_entry_temp));
    alpha = hash_entry_temp->value;
    hashmap_get(pmap, "betato", (void**)(&hash_entry_temp));
    betato = hash_entry_temp->value;
    hashmap_get(pmap, "alphato", (void**)(&hash_entry_temp));
    alphato = hash_entry_temp->value;
    hashmap_get(pmap, "m_to", (void**)(&hash_entry_temp));
    m_to = hash_entry_temp->value;
    hashmap_get(pmap, "D_to", (void**)(&hash_entry_temp));
    D_to = hash_entry_temp->value;
    hashmap_get(pmap, "gammand", (void**)(&hash_entry_temp));
    gammand = hash_entry_temp->value;
    hashmap_get(pmap, "wbb", (void**)(&hash_entry_temp));
    wbb = hash_entry_temp->value;
    hashmap_get(pmap, "wfb", (void**)(&hash_entry_temp));
    wfb = hash_entry_temp->value;
    hashmap_get(pmap, "wnb", (void**)(&hash_entry_temp));
    wnb = hash_entry_temp->value;
    hashmap_get(pmap, "wfibb", (void**)(&hash_entry_temp));
    wfibb = hash_entry_temp->value;
    hashmap_get(pmap, "wbf", (void**)(&hash_entry_temp));
    wbf = hash_entry_temp->value;
    hashmap_get(pmap, "wff", (void**)(&hash_entry_temp));
    wff = hash_entry_temp->value;
    hashmap_get(pmap, "wnf", (void**)(&hash_entry_temp));
    wnf = hash_entry_temp->value;
    hashmap_get(pmap, "wfibf", (void**)(&hash_entry_temp));
    wfibf = hash_entry_temp->value;
    hashmap_get(pmap, "wbfib", (void**)(&hash_entry_temp));
    wbfib = hash_entry_temp->value;
    hashmap_get(pmap, "wffib", (void**)(&hash_entry_temp));
    wffib = hash_entry_temp->value;
    hashmap_get(pmap, "wnfib", (void**)(&hash_entry_temp));
    wnfib = hash_entry_temp->value;
    hashmap_get(pmap, "wfibfib", (void**)(&hash_entry_temp));
    wfibfib = hash_entry_temp->value;
    hashmap_get(pmap, "wbn", (void**)(&hash_entry_temp));
    wbn = hash_entry_temp->value;
    hashmap_get(pmap, "wfn", (void**)(&hash_entry_temp));
    wfn = hash_entry_temp->value;
    hashmap_get(pmap, "wnn", (void**)(&hash_entry_temp));
    wnn = hash_entry_temp->value;
    hashmap_get(pmap, "wfibn", (void**)(&hash_entry_temp));
    wfibn = hash_entry_temp->value;
    hashmap_get(pmap, "wbto", (void**)(&hash_entry_temp));
    wbto = hash_entry_temp->value;
    hashmap_get(pmap, "wfto", (void**)(&hash_entry_temp));
    wfto = hash_entry_temp->value;
    hashmap_get(pmap, "wnto", (void**)(&hash_entry_temp));
    wnto = hash_entry_temp->value;
    hashmap_get(pmap, "wfibto", (void**)(&hash_entry_temp));
    wfibto = hash_entry_temp->value;
    hashmap_get(pmap, "wtoto", (void**)(&hash_entry_temp));
    wtoto = hash_entry_temp->value;

    hashmap_get(pmap, "wbb_dif", (void**)(&hash_entry_temp));
    wbb_dif = hash_entry_temp->value;
    hashmap_get(pmap, "wfb_dif", (void**)(&hash_entry_temp));
    wfb_dif = hash_entry_temp->value;
    hashmap_get(pmap, "wff_dif", (void**)(&hash_entry_temp));
    wff_dif = hash_entry_temp->value;
    hashmap_get(pmap, "wbf_dif", (void**)(&hash_entry_temp));
    wbf_dif = hash_entry_temp->value;
    hashmap_get(pmap, "wfibb_dif", (void**)(&hash_entry_temp));
    wfibb_dif = hash_entry_temp->value;
    hashmap_get(pmap, "wfibf_dif", (void**)(&hash_entry_temp));
    wfibf_dif = hash_entry_temp->value;
    hashmap_get(pmap, "wfibn_dif", (void**)(&hash_entry_temp));
    wfibn_dif = hash_entry_temp->value;

//--
    hashmap_get(pmap, "wbma", (void**)(&hash_entry_temp));
    wbma = hash_entry_temp->value;
    hashmap_get(pmap, "wfma", (void**)(&hash_entry_temp));
    wfma = hash_entry_temp->value;
    hashmap_get(pmap, "wnma", (void**)(&hash_entry_temp));
    wnma = hash_entry_temp->value;
    hashmap_get(pmap, "wfibma", (void**)(&hash_entry_temp));
    wfibma = hash_entry_temp->value;
    hashmap_get(pmap, "wmama", (void**)(&hash_entry_temp));
    wmama = hash_entry_temp->value;
    hashmap_get(pmap, "wbma_dif", (void**)(&hash_entry_temp));
    wbma_dif = hash_entry_temp->value;
    hashmap_get(pmap, "wfibma_dif", (void**)(&hash_entry_temp));
    wfibma_dif = hash_entry_temp->value;
    hashmap_get(pmap, "wndma_dif", (void**)(&hash_entry_temp));
    wndma_dif = hash_entry_temp->value;

    hashmap_get(pmap, "wbmr", (void**)(&hash_entry_temp));
    wbmr = hash_entry_temp->value;
    hashmap_get(pmap, "wfmr", (void**)(&hash_entry_temp));
    wfmr = hash_entry_temp->value;
    hashmap_get(pmap, "wnmr", (void**)(&hash_entry_temp));
    wnmr = hash_entry_temp->value;
    hashmap_get(pmap, "wfibmr", (void**)(&hash_entry_temp));
    wfibmr = hash_entry_temp->value;
    hashmap_get(pmap, "wmrmr", (void**)(&hash_entry_temp));
    wmrmr = hash_entry_temp->value;
    hashmap_get(pmap, "wbmr_dif", (void**)(&hash_entry_temp));
    wbmr_dif = hash_entry_temp->value;
    hashmap_get(pmap, "wfibmr_dif", (void**)(&hash_entry_temp));
    wfibmr_dif = hash_entry_temp->value;
    hashmap_get(pmap, "wndmr_dif", (void**)(&hash_entry_temp));
    wndmr_dif = hash_entry_temp->value;
    hashmap_get(pmap, "smr", (void**)(&hash_entry_temp));
    smr = hash_entry_temp->value;
    hashmap_get(pmap, "D_mr", (void**)(&hash_entry_temp));
    D_mr = hash_entry_temp->value;
    hashmap_get(pmap, "D_ma", (void**)(&hash_entry_temp));
    D_ma = hash_entry_temp->value;
    hashmap_get(pmap, "lmr", (void**)(&hash_entry_temp));
    lmr = hash_entry_temp->value;
    hashmap_get(pmap, "lma", (void**)(&hash_entry_temp));
    lma = hash_entry_temp->value;
    hashmap_get(pmap, "X_mr", (void**)(&hash_entry_temp));
    X_mr = hash_entry_temp->value;
    hashmap_get(pmap, "X_ma", (void**)(&hash_entry_temp));
    X_ma = hash_entry_temp->value;
//-- at
    hashmap_get(pmap, "wdn", (void**)(&hash_entry_temp));
    wdn = hash_entry_temp->value;
    hashmap_get(pmap, "wdn_dif", (void**)(&hash_entry_temp));
    wdn_dif = hash_entry_temp->value;
    hashmap_get(pmap, "b_dmg", (void**)(&hash_entry_temp));
    b_dmg = hash_entry_temp->value;
    hashmap_get(pmap, "nd_dmg", (void**)(&hash_entry_temp));
    nd_dmg = hash_entry_temp->value;
    hashmap_get(pmap, "gamma_dmg", (void**)(&hash_entry_temp));
    gamma_dmg = hash_entry_temp->value;
    hashmap_get(pmap, "s_inf", (void**)(&hash_entry_temp));
    s_inf = hash_entry_temp->value;

    int interval = (int)T/numIteracoesGravadas;
    printf("Size: %d\n", size_);
    printf("Interval: %d\n", interval);

    float ***b = (float***)malloc(2*sizeof(float**)); // bacteria
    float ***f = (float***)malloc(2*sizeof(float**));
    float ***n = (float***)malloc(2*sizeof(float**));
    float ***fib = (float***)malloc(2*sizeof(float**));
    float ***to = (float***)malloc(2*sizeof(float**));
    float ***nd = (float***)malloc(2*sizeof(float**));
    float ***mr = (float***)malloc(2*sizeof(float**));
    float ***ma = (float***)malloc(2*sizeof(float**));
    float ***dmt = (float***)malloc(2*sizeof(float**));

    float ***wb = (float***)malloc(2*sizeof(float**));
    float ***wf = (float***)malloc(2*sizeof(float**));
    float ***wn = (float***)malloc(2*sizeof(float**));
    float ***wfib = (float***)malloc(2*sizeof(float**));
    float ***wto = (float***)malloc(2*sizeof(float**));

    float ***wbdif = (float***)malloc(2*sizeof(float**));
    float ***wfdif = (float***)malloc(2*sizeof(float**));
    float ***wndif = (float***)malloc(2*sizeof(float**));
    float ***wfibdif = (float***)malloc(2*sizeof(float**));
    float ***wtodif = (float***)malloc(2*sizeof(float**));

    float ***gwb = (float***)malloc(2*sizeof(float**));
    float ***gwf = (float***)malloc(2*sizeof(float**));
    float ***gwn = (float***)malloc(2*sizeof(float**));
    float ***gwfib = (float***)malloc(2*sizeof(float**));
    float ***gwto = (float***)malloc(2*sizeof(float**));

    float ***gwone = (float***)malloc(2*sizeof(float**));

    float ***wmrdif = (float***)malloc(2*sizeof(float**));
    float ***wmadif = (float***)malloc(2*sizeof(float**));
    float ***wmr = (float***)malloc(2*sizeof(float**));
    float ***gwmr = (float***)malloc(2*sizeof(float**));
    float ***wma = (float***)malloc(2*sizeof(float**));
    float ***gwma = (float***)malloc(2*sizeof(float**));

    for(t=0; t < 2; t++)
    {
        b[t] = (float**)malloc(size_*sizeof(float*));
        f[t] = (float**)malloc(size_*sizeof(float*));
        n[t] = (float**)malloc(size_*sizeof(float*));
        fib[t] = (float**)malloc(size_*sizeof(float*));
        to[t] = (float**)malloc(size_*sizeof(float*));
        nd[t] = (float**)malloc(size_*sizeof(float*));
        mr[t] = (float**)malloc(size_*sizeof(float*));
        ma[t] = (float**)malloc(size_*sizeof(float*));
        dmt[t] = (float**)malloc(size_*sizeof(float*));

        wb[t] = (float**)malloc(size_*sizeof(float*));
        wf[t] = (float**)malloc(size_*sizeof(float*));
        wn[t] = (float**)malloc(size_*sizeof(float*));
        wfib[t] = (float**)malloc(size_*sizeof(float*));
        wto[t] = (float**)malloc(size_*sizeof(float*));

        wbdif[t] = (float**)malloc(size_*sizeof(float*));
        wfdif[t] = (float**)malloc(size_*sizeof(float*));
        wndif[t] = (float**)malloc(size_*sizeof(float*));
        wfibdif[t] = (float**)malloc(size_*sizeof(float*));
        wtodif[t] = (float**)malloc(size_*sizeof(float*));

        gwb[t] = (float**)malloc(size_*sizeof(float*));
        gwf[t] = (float**)malloc(size_*sizeof(float*));
        gwn[t] = (float**)malloc(size_*sizeof(float*));
        gwfib[t] = (float**)malloc(size_*sizeof(float*));
        gwto[t] = (float**)malloc(size_*sizeof(float*));

        gwone[t] = (float**)malloc(size_*sizeof(float*));
//----
        wmrdif[t] = (float**)malloc(size_*sizeof(float*));
        wmadif[t] = (float**)malloc(size_*sizeof(float*));
        wmr[t] = (float**)malloc(size_*sizeof(float*));
        gwmr[t] = (float**)malloc(size_*sizeof(float*));
        wma[t] = (float**)malloc(size_*sizeof(float*));
        gwma[t] = (float**)malloc(size_*sizeof(float*));

        for(x=0; x < size_; x++)
        {
          b[t][x] = (float*)malloc(size_*sizeof(float));
          f[t][x] = (float*)malloc(size_*sizeof(float));
          n[t][x] = (float*)malloc(size_*sizeof(float));
          fib[t][x] = (float*)malloc(size_*sizeof(float));
          to[t][x] = (float*)malloc(size_*sizeof(float));
          nd[t][x] = (float*)malloc(size_*sizeof(float));
          mr[t][x] = (float*)malloc(size_*sizeof(float));
          ma[t][x] = (float*)malloc(size_*sizeof(float));
          dmt[t][x] = (float*)malloc(size_*sizeof(float));


          wb[t][x] = (float*)malloc(size_*sizeof(float));
          wf[t][x] = (float*)malloc(size_*sizeof(float));
          wn[t][x] = (float*)malloc(size_*sizeof(float));
          wfib[t][x] = (float*)malloc(size_*sizeof(float));
          wto[t][x] = (float*)malloc(size_*sizeof(float));

          wbdif[t][x] = (float*)malloc(size_*sizeof(float));
          wfdif[t][x] = (float*)malloc(size_*sizeof(float));
          wndif[t][x] = (float*)malloc(size_*sizeof(float));
          wfibdif[t][x] = (float*)malloc(size_*sizeof(float));
          wtodif[t][x] = (float*)malloc(size_*sizeof(float));

          gwb[t][x] = (float*)malloc(size_*sizeof(float));
          gwf[t][x] = (float*)malloc(size_*sizeof(float));
          gwn[t][x] = (float*)malloc(size_*sizeof(float));
          gwfib[t][x] = (float*)malloc(size_*sizeof(float));
          gwto[t][x] = (float*)malloc(size_*sizeof(float));

          gwone[t][x] = (float*)malloc(size_*sizeof(float));

          wmrdif[t][x] = (float*)malloc(size_*sizeof(float));
          wmadif[t][x] = (float*)malloc(size_*sizeof(float));
          wmr[t][x] = (float*)malloc(size_*sizeof(float));
          gwmr[t][x] = (float*)malloc(size_*sizeof(float));
          wma[t][x] = (float*)malloc(size_*sizeof(float));
          gwma[t][x] = (float*)malloc(size_*sizeof(float));

          for(y=0; y < size_; y++)
            gwone[t][x][y] = 1;
            b[t][x][y] =  0;
            f[t][x][y] =  0;
            n[t][x][y] =  0;
            mr[t][x][y] = 0;
            ma[t][x][y] = 0;
            fib[t][x][y] = 0.;
            to[t][x][y] = 0.;
            nd[t][x][y] = 0.;
            dmt[t][x][y] = 0.;
            wmrdif[t][x][y] = 0;
            wmadif[t][x][y] = 0;
            wmr[t][x][y] = 0;
            gwmr[t][x][y] = 0;
            wma[t][x][y] = 0;
            gwma[t][x][y] = 0;
        }
    }
//----
    //write space file
    FILE* espaco = fopen("x.dat", "w");
    float inc = 0.0;
    for (x=0; x < size_; x++)
    {
        fprintf(espaco, "%.1f\n", inc);
        inc = inc + (float)deltaX;
    }
    fclose(espaco);

    espaco = fopen("y.dat", "w");
    inc = 0.0;
    for (y=0; y < size_; y++)
    {
        fprintf(espaco, "%.1f\n", inc);
        inc = inc + (float)deltaY;
    }
    fclose(espaco);
    //write time file
    float time = 0.0;
    FILE* tempo = fopen("t.dat", "w");
    for (t=0; t <= T; t++)
    {
        if ((t % interval) == 0)
        {
            fprintf(tempo, "%.0f \n", time);
        }
        time = time + (float)deltaT;//(float)deltaT;
    }
    fclose(tempo);


    FILE* bfile = fopen("bacteria.dat", "w");
    FILE* ffile = fopen("coa_vwbp.dat", "w");
    FILE* nfile = fopen("neutrophil.dat", "w");
    FILE* fibfile = fopen("fibrin.dat", "w");
    FILE* tofile = fopen("toxin.dat", "w");
    FILE* ndfile = fopen("nd.dat", "w");
    FILE* mrfile = fopen("mr.dat", "w");
    FILE* mafile = fopen("ma.dat", "w");
    FILE* dmt_file = fopen("dmt.dat", "w");

    int tAntigo = 0;
    int tAtual = 1;
    int totalBacterias = 0, totalNeutrofilos = 0;

    // for(x=0; x < size_; x++)
    // {
    //   for(y=0; y < size_; y++)
    //   {
    //
    //   }
    // }

    fseek(entryfile,0,SEEK_SET);
    readAndSetInitialCondition(entryfile,(float***)b,(float***)f,(float***)fib,(float***)n,(float***)mr,(float***)ma,tAntigo);
    for(x=0; x < size_; x++)
    {
      for(y=0; y < size_; y++)
      {
        totalBacterias += b[tAntigo][x][y];
        totalNeutrofilos += n[tAntigo][x][y];

        updateGWReaction((float***)b,(float***)f,(float***)n,(float***)fib,(float***)wb,(float***)wf,(float***)wn,(float***)wfib,
            (float***)gwb,(float***)gwf,(float***)gwn,(float***)gwfib,(float***)to, (float***)wto,(float***)gwto,(float***)wmr,
            (float***)gwmr,(float***)wma,(float***)gwma, (float***)ma, (float***)mr ,(float***)dmt,x, y, tAntigo);
        updateWDiffusion((float***)b,(float***)f,(float***)n,(float***)fib,(float***)wbdif,(float***)wfdif,(float***)wndif,(float***)wfibdif,
                    (float***)to, (float***)wtodif, (float***)wmrdif, (float***)wmadif,(float***)dmt, x, y, tAntigo);
      }
    }
    imprimeDados(bfile,(float***)b,tAntigo);
    imprimeDados(ffile,(float***)f,tAntigo);
    imprimeDados(nfile,(float***)n,tAntigo);
    imprimeDados(fibfile,(float***)fib,tAntigo);
    imprimeDados(tofile,(float***)to,tAntigo);
    imprimeDados(ndfile,(float***)nd,tAntigo);
    imprimeDados(mrfile,(float***)mr,tAntigo);
    imprimeDados(mafile,(float***)ma,tAntigo);
    imprimeDados(dmt_file,(float***)dmt,tAntigo);

    float b_, f_, n_, fib_, to_, nd_, mr_, ma_, d_;
    float gwb_, gwf_, gwfib_, gwn_, gwmr_, gwma_;
    float newb, newf, newn, newfib, newto, newnd, newmr, newma;

    for(int step = 1; step < T+1; step++)
    {
        if(step % 2 == 0)
        {
            tAtual = 0;
            tAntigo = 1;
        }
        else
        {
            tAtual = 1;
            tAntigo = 0;
        }
        //Criar para MR, MA, DMT, TO, ND - Gravar em um arquivo
        totalBacterias = 0;
        totalNeutrofilos = 0;
        #pragma omp parallel for num_threads(4) private(x,y,b_,f_,n_,fib_,nd_,to_,mr_,ma_,newb)
        for (x=0; x < size_; x++)
        {
          for (y=0; y < size_; y++)
          {

              #pragma omp critital
              {
                  totalBacterias += b[tAtual][x][y];
                  totalNeutrofilos += n[tAtual][x][y];
              }
            //Bacteria pde
            //b[tAtual][x][y] = ( (r_b - l*n[tAntigo][x][y])*b[tAntigo][x][y]*gwb[tAntigo][x][y] //- l*n[tAntigo][x][y]*b[tAntigo][x][y]*gwn[tAntigo][x][y]
            //+ D_b*localAverage((float***)b,x,y,(float***)wbdif,tAntigo)
            //+ D2_b*flocalAverage((float***)b,x,y,(float***)wbdif,tAntigo,(float***)b,b_kdif,b_nexp)
            //)*deltaT + b[tAntigo][x][y];

            b_ = b[tAntigo][x][y];
            f_ = f[tAntigo][x][y];
            n_ = n[tAntigo][x][y];
            fib_ = fib[tAntigo][x][y];
            to_ = to[tAntigo][x][y];
            nd_ = nd[tAntigo][x][y];
            //dmt_ = dmt[tAntigo][x][y];
            mr_ = mr[tAntigo][x][y];
            ma_ = ma[tAntigo][x][y];

            b[tAtual][x][y] = (((r_b - l*n[tAntigo][x][y]
                - lmr*mr[tAntigo][x][y] - lma*ma[tAntigo][x][y])*b[tAntigo][x][y])* gwb[tAntigo][x][y] + D_b*localAverage(b,x,y,wbdif,tAntigo)
            )*deltaT + b[tAntigo][x][y];
            b[tAtual][x][y] = verifyDensity(b[tAtual][x][y]);

            //Coa/vWbp pde
            f[tAtual][x][y] = ( k*b[tAntigo][x][y]*gwf[tAntigo][x][y]    //k*f[tAntigo][x][y]*b[tAntigo][x][y]*(1 - f[tAntigo][x][y]) - p*f[tAntigo][x][y]*b[tAntigo][x][y]*gwfib[tAntigo][x][y] //*gwf[tAntigo][x][y]   //((1 + b_kdif)*pow(b[tAntigo][x][y],b_nexp)/(b_kdif + pow(b[tAntigo][x][y],b_nexp)))*gwf[tAntigo][x][y]    //*b[tAntigo][x][y]*gwf[tAntigo][x][y]
            - m_f*f[tAntigo][x][y]
            + D_f*localAverage((float***)f,x,y,(float***)wfdif,tAntigo)
            + D2_f*flocalAverage((float***)f,x,y,(float***)wfdif,tAntigo,(float***)b,f_kdif,f_nexp)
            - X_f*chemotaxis((float***)f,(float***)b,x,y,(float***)wfdif,tAntigo)
            )*deltaT + f[tAntigo][x][y];
            f[tAtual][x][y] = verifyDensity(f[tAtual][x][y]);

            //Fibrin pde
            fib[tAtual][x][y] = (p*b[tAntigo][x][y]*f[tAntigo][x][y]*gwfib[tAntigo][x][y])*deltaT + fib[tAntigo][x][y];;
            //fib[tAtual][x][y] = (p*b[tAntigo][x][y]*f[tAntigo][x][y]*gwfib[tAntigo][x][y])*deltaT + fib[tAntigo][x][y]; //*gwfib[tAntigo][x][y])*deltaT + fib[tAntigo][x][y];
            fib[tAtual][x][y] = verifyDensity(fib[tAtual][x][y]);

            //Neutrophil PDE
           // n[tAtual][x][y] =  ( (s*b[tAntigo][x][y]*n[tAntigo][x][y] + s2*b[tAntigo][x][y]+ s3)*gwn[tAntigo][x][y]
           // - m_n*n[tAntigo][x][y] - alphato*to[tAntigo][x][y]*n[tAntigo][x][y] -alpha*b[tAntigo][x][y]*n[tAntigo][x][y]
           // + D_n*localAverage((float***)n,x,y,(float***)wndif,tAntigo)
            //- X_n*chemotaxis((float***)n, (float***)b, x, y, (float***)wndif, tAntigo)
           // )*deltaT + n[tAntigo][x][y];

            n[tAtual][x][y] = (s_inf*n[tAntigo][x][y]*dmt[tAntigo][x][y]*gwn[tAntigo][x][y] + s*ma[tAntigo][x][y]*gwn[tAntigo][x][y] - alpha*to[tAntigo][x][y]*n[tAntigo][x][y]
                - m_n*n[tAntigo][x][y]
            + D_n*localAverage(n,x,y,wndif,tAntigo)- X_n*chemotaxis((float***)n, (float***)b, x, y, (float***)wndif, tAntigo)
            )*deltaT + n[tAntigo][x][y];
            n[tAtual][x][y] = verifyDensity(n[tAtual][x][y]);

            to[tAtual][x][y] = (betato*b[tAntigo][x][y]*gwto[tAntigo][x][y] - m_to*to[tAntigo][x][y]
              + D_to*diffusion((float***)to,x,y,tAntigo)
            )*deltaT + to[tAntigo][x][y];
            //*(1 - to[tAntigo][x])
//---
            mr[tAtual][x][y] = (s_inf*mr[tAntigo][x][y]*dmt[tAntigo][x][y]*gwmr[tAntigo][x][y] + (smr - lmr)*b[tAntigo][x][y]*mr[tAntigo][x][y]*gwmr[tAntigo][x][y]
            - alpha*to[tAntigo][x][y]*mr[tAntigo][x][y]
             + D_mr*localAverage(mr,x,y,wmrdif,tAntigo)
             - X_mr*chemotaxis((float***)mr, (float***)b, x, y, (float***)wmrdif, tAntigo))*deltaT + mr[tAntigo][x][y];
            mr[tAtual][x][y] = verifyDensity(mr[tAtual][x][y]);

            //Ajustar taxas de fagocitose. Estudar efeito toxina

            ma[tAtual][x][y] = (lmr*mr[tAntigo][x][y]*b[tAntigo][x][y]*gwmr[tAntigo][x][y] - alpha*to[tAntigo][x][y]*ma[tAntigo][x][y]
             + D_ma*localAverage(ma,x,y,wmadif,tAntigo)
             - X_ma*chemotaxis((float***)ma, (float***)b, x, y, (float***)wmadif, tAntigo))*deltaT + ma[tAntigo][x][y];
            //alterar B[t][x][y] = (r*B - ln*N*B - lmr*MR*B - lma*MA*B) * GB(WB) + D_B*localAverage(B,x,y,wbdif,tAntigo);
            //alterar N[t][x][y] = B2*MA*GN(WN) - alpha*to*N + Dn*localAverage(N,x,y,wndif,tAntigo);
            ma[tAtual][x][y] = verifyDensity(ma[tAtual][x][y]);

/* //--            nd[tAtual][x][y] = (
                alphato*to[tAntigo][x][y]*ma[tAntigo][x][y] +
                alphato*to[tAntigo][x][y]*mr[tAntigo][x][y] +
                alphato*to[tAntigo][x][y]*n[tAntigo][x][y] - gammand*nd[tAntigo][x][y]*ma[tAntigo][x][y]
                + m_n*n[tAntigo][x][y]
            )*deltaT
            + nd[tAntigo][x][y];
*/
            nd[tAtual][x][y] = (
              alphato*to[tAntigo][x][y]*n[tAntigo][x][y] +
              m_n*n[tAntigo][x][y] - gammand*nd[tAntigo][x][y] * (ma[tAntigo][x][y] + mr[tAntigo][x][y])
            )*deltaT + nd[tAntigo][x][y];

            nd[tAtual][x][y] = verifyDensity(nd[tAtual][x][y]);

            dmt[tAtual][x][y] = ( (b_dmg * b[tAntigo][x][y] + nd_dmg * nd[tAntigo][x][y])*(1 - dmt[tAntigo][x][y])
            - gamma_dmg * (mr[tAntigo][x][y] + ma[tAntigo][x][y])
            )*deltaT + dmt[tAntigo][x][y];

            updateGWReaction((float***)b,(float***)f,(float***)n,(float***)fib,(float***)wb,(float***)wf,(float***)wn,(float***)wfib,
                      (float***)gwb,(float***)gwf,(float***)gwn,(float***)gwfib,(float***)to, (float***)wto, (float***)gwto,
                      (float***)wmr,(float***)gwmr,(float***)wma, (float***)gwma, (float***)ma, (float***)mr,(float***)dmt, x,y, tAtual);
            updateWDiffusion((float***)b,(float***)f,(float***)n,(float***)fib,(float***)wbdif,(float***)wfdif,(float***)wndif,
            (float***)wfibdif, (float***)to, (float***)wtodif, (float***)wmrdif, (float***)wmadif ,(float***)dmt,x, y, tAntigo);
          }
        }
        //print('Total de bacterias para o passo '+ str(step) + ' = ' + str(totalBacterias))
        if((step % interval) == 0){
            imprimeDados(bfile,(float***)b,tAtual);
            imprimeDados(ffile,(float***)f,tAtual);
            imprimeDados(nfile,(float***)n,tAtual);
            imprimeDados(fibfile,(float***)fib,tAtual);
            imprimeDados(tofile,(float***)to,tAtual);
            imprimeDados(ndfile,(float***)nd,tAtual);
            imprimeDados(mrfile,(float***)mr,tAtual);
            imprimeDados(mafile,(float***)ma,tAtual);
            imprimeDados(dmt_file,(float***)dmt,tAtual);
            // fflush(bfile);
            // fflush(ffile);
            // fflush(nfile);
            // fflush(fibfile);
            // fflush(tofile);
            // fflush(ndfile);
        }
    }
    for(t=0; t < 2; t++)
    {
      for(x=0; x < size_; x++)
      {
        free(b[t][x]);
        free(f[t][x]);
        free(n[t][x]);
        free(fib[t][x]);
        free(to[t][x]);
        free(nd[t][x]);
        free(dmt[t][x]);
        free(mr[t][x]);
        free(ma[t][x]);
        free(wb[t][x]);
        free(wf[t][x]);
        free(wn[t][x]);
        free(wto[t][x]);
        free(wfib[t][x]);
        free(wbdif[t][x]);
        free(wfdif[t][x]);
        free(wndif[t][x]);
        free(wfibdif[t][x]);
        free(wtodif[t][x]);
        free(gwb[t][x]);
        free(gwf[t][x]);
        free(gwn[t][x]);
        free(gwfib[t][x]);
        free(gwto[t][x]);
        free(wmrdif[t][x]);
        free(wmadif[t][x]);
        free(wmr[t][x]);
        free(gwmr[t][x]);
        free(wma[t][x]);
        free(gwma[t][x]);
      }
      free(b[t]);
      free(f[t]);
      free(n[t]);
      free(fib[t]);
      free(to[t]);
      free(nd[t]);
      free(dmt[t]);
      free(mr[t]);
      free(ma[t]);
      free(wb[t]);
      free(wf[t]);
      free(wn[t]);
      free(wfib[t]);
      free(wto[t]);
      free(wbdif[t]);
      free(wfdif[t]);
      free(wndif[t]);
      free(wfibdif[t]);
      free(wtodif[t]);
      free(gwb[t]);
      free(gwf[t]);
      free(gwn[t]);
      free(gwfib[t]);
      free(gwto[t]);
      free(wmrdif[t]);
      free(wmadif[t]);
      free(wmr[t]);
      free(gwmr[t]);
      free(wma[t]);
      free(gwma[t]);
    }
    free(b);
    free(f);
    free(n);
    free(fib);
    free(to);
    free(nd);
    free(dmt);
    free(mr);
    free(ma);
    free(wb);
    free(wf);
    free(wn);
    free(wfib);
    free(wto);
    free(wbdif);
    free(wfdif);
    free(wndif);
    free(wfibdif);
    free(wtodif);
    free(gwb);
    free(gwf);
    free(gwn);
    free(gwfib);
    free(gwto);
    free(wmrdif);
    free(wmadif);
    free(wmr);
    free(gwmr);
    free(wma);
    free(gwma);

    fclose(entryfile);
    fclose(bfile);
    fclose(ffile);
    fclose(nfile);
    fclose(fibfile);
    fclose(tofile);
    fclose(ndfile);
    fclose(mrfile);
    fclose(mafile);
    fclose(dmt_file);

    hashmap_free(pmap);

    return 0;
}
