//consec_corr not matching due to possible cummulative effect of rounding off

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct cmplx
{
    int real;
    int imag;
};

struct maxPar
{
    int maxValue;
    int maxIndex;
};

typedef struct cmplx complex;
typedef struct maxPar maxParam;

#define DEFAULT_DFT_SIZE 512
#define DFT_SIZE 128
#define SUBCARRIER_NUMBER 11
#define SYMBOL_NUMBER 11
#define pi 3.1416
#define TD_RX_LENGTH 1920
#define SF 3


complex X[DEFAULT_DFT_SIZE];
complex Xx[DEFAULT_DFT_SIZE];
complex corrOut[DEFAULT_DFT_SIZE];
complex referenceGrid[SUBCARRIER_NUMBER][SYMBOL_NUMBER];
complex z[SUBCARRIER_NUMBER][SYMBOL_NUMBER];
complex zIDFT[DFT_SIZE][SYMBOL_NUMBER];
complex rxGrid[SUBCARRIER_NUMBER][SYMBOL_NUMBER];
complex rxGrid_corrected[SUBCARRIER_NUMBER][SYMBOL_NUMBER];
int absSquare_idft[DFT_SIZE]={0};
complex one;
complex negativeOne;
complex refNPSSgrid_freqDomain[SUBCARRIER_NUMBER][SYMBOL_NUMBER];
complex refNPSSgrid_timeDomain[DFT_SIZE][SYMBOL_NUMBER];
complex rxNPSS_timeDomain[TD_RX_LENGTH];
complex rxNPSS_timeDomain_corrected[TD_RX_LENGTH];
complex rxNPSS_timeDomain1[DFT_SIZE][SYMBOL_NUMBER+3];
complex rxNPSS_timeDomain1_corrected[DFT_SIZE][SYMBOL_NUMBER+3];
complex rxNPSSgrid_timeDomain1[DFT_SIZE][SYMBOL_NUMBER];
complex rxNPSSgrid_freqDomain[SUBCARRIER_NUMBER][SYMBOL_NUMBER];
complex colRxGridtd[128]={0};

complex fnCmplxAdd(complex a, complex b){
    complex c;
    c.real = a.real + b.real;
    c.imag = a.imag + b.imag;
    return c;
};

complex fnCmplxSub(complex a, complex b){
    complex c;
    c.real = a.real - b.real;
    c.imag = a.imag - b.imag;
    return c;
};

complex fnCmplxMult(complex a, complex b){
    complex c;
    long int temp_real,temp_imag;
    temp_real = a.real * b.real - a.imag * b.imag;
    temp_imag = a.real * b.imag + a.imag * b.real;

    temp_real >>= SF;
    temp_imag >>= SF;

    c.real = temp_real;
    c.imag = temp_imag;
    return c;
};


complex fnConj(complex a){
    complex c;
    c.real = a.real;
    c.imag = 0 - a.imag;
    return c;
};

void fnIDFT(complex in[], float N ){

    int n,k;
    float Xx1[512];
    float Xx2[512];


     for(n=0;n<N;n++)
     {
        Xx1[n]=0;
        Xx2[n]=0;
        for(k=0;k<N;k++)
        {
            Xx1[n] = Xx1[n]+in[k].real*cos((2*pi*k*n)/N) - in[k].imag*sin((2*pi*k*n)/N);
            Xx2[n] = Xx2[n]+in[k].imag*cos((2*pi*k*n)/N) + in[k].real*sin((2*pi*k*n)/N);

        }
        X[n].real=round(Xx1[n]/N);
        X[n].imag=round(Xx2[n]/N);
     }

};


void correlator(complex in1[], complex in2[], int M, int N){
    int m;
    complex  in2_conj[11];

    for(m=0;m<M;m++)
    {
        in2_conj[m].real = in2[m].real;
        in2_conj[m].imag = -in2[m].imag;

        corrOut[m].real = fnCmplxMult(in1[m],in2_conj[m]).real;
        corrOut[m].imag = fnCmplxMult(in1[m],in2_conj[m]).imag;
    }

};

int absoluteSquare(complex in1){

    int abs;
    long int temp;

    temp = in1.real*in1.real + in1.imag*in1.imag;
    temp >>= SF;
    abs = temp;
    return abs;

};


maxParam maxParameters(int a[], int length){
    maxParam m ;
    m.maxValue = a[0];
    m.maxIndex = 0;

    int i;
    for(i=0;i<length;i++){
        if(a[i]>= m.maxValue){
            m.maxValue = a[i] ;
            m.maxIndex = i;
        }
    }
    return m;
    }


void fnDFT(complex in_new[], int N ){

    int n,k;
    float Xx_real[512];
    float Xx_imag[512];


     for(n=0;n<N;n++)
     {
        Xx_real[n]=0;
        Xx_imag[n]=0;

      for(k=0;k<N;k++)
      {
        Xx_real[n] = Xx_real[n]+in_new[k].real*cos((2*pi*k*n)/N) + in_new[k].imag*sin((2*pi*k*n)/N);
        Xx_imag[n] = Xx_imag[n]+in_new[k].imag*cos((2*pi*k*n)/N) - in_new[k].real*sin((2*pi*k*n)/N);

      }
     // printf("%d + i %d\n", Xx_real[n],Xx_imag[n] );

     Xx[n].real = round(Xx_real[n]);
     Xx[n].imag = round(Xx_imag[n]);
     //printf("%d + i %d\n", Xx[n].real,Xx[n].imag);
     }
};

int main()
{

    int td_refNPSSgrid[2816];
    int td_rxNPSSgrid[3840];
    complex td_refNPSSgrid_vector[1408];
    int i,j,m,n,l,s,q,r;


////////////////////////////////////////////////////////////////////
    FILE *rm2;      //importing td_refNPSSgrid1
        int ss2;
        float buftd[1000];
        rm2 = fopen("td_refNPSSgrid_fixed.dat", "r");

        while (!feof(rm2)) {
            for(ss2=0;ss2<2816;ss2++)
            {
                fscanf(rm2,"%f",&buftd[ss2]);
            }
        }
    fclose(rm2);

    for(m=0;m<2816;m++){    //reassigning

        td_refNPSSgrid[m] = buftd[m];

    }

    for(i=0;i<sizeof(td_refNPSSgrid)/sizeof(td_refNPSSgrid[0]);i++){     //orgainzing into real and complex

        if(i%2 == 0){
            td_refNPSSgrid_vector[i/2].real = td_refNPSSgrid[i];
        }

        if(i%2 == 1){
            td_refNPSSgrid_vector[(i-1)/2].imag = td_refNPSSgrid[i];
        }
    }

    for(n=0;n<128 ;n++){        //creating refNPSSgrid_timeDomain grid
        for(l=0;l<11 ;l++){
            refNPSSgrid_timeDomain[n][l].real = td_refNPSSgrid_vector[n + 128*l].real;
            refNPSSgrid_timeDomain[n][l].imag = td_refNPSSgrid_vector[n + 128*l].imag;
        }
    }

//// //   display TD_reference grid
//    printf("time domain REFERENCE GRID \n\n");
//    for(l=0;l<11 ;l++){
//        for(n=0;n<128 ;n++){
//
//            printf("%d + %di \n ",refNPSSgrid_timeDomain[n][l].real,refNPSSgrid_timeDomain[n][l].imag);
//        }
//    printf("\n");
//    }
//

    complex colRefGridtd[128]={0};
    for(q=0;q<11;q++){      //dft of the imported time-domain grid refNPSSgrid_timeDomain to get freq-domain refNPSSgrid

        for(r=0;r<128 ;r++){


            colRefGridtd[r].real = refNPSSgrid_timeDomain[r][q].real;
            colRefGridtd[r].imag = refNPSSgrid_timeDomain[r][q].imag;
            //printf("%d + i%d\n", colRefGridtd[r].real, colRefGridtd[r].imag);

        }
            //printf("\n");

        fnDFT(colRefGridtd, 128);
        for(s=0;s<11;s++){
            refNPSSgrid_freqDomain[s][q].real = Xx[s].real;
            refNPSSgrid_freqDomain[s][q].imag = Xx[s].imag;
           // printf("%d + i%d\n", refNPSSgrid_freqDomain[s][q].real, refNPSSgrid_freqDomain[s][q].imag);
        }
    }


////display refNPSSgrid_freqDomain grid
//    for(q=0;q<11;q++){
//        for(r=0;r<11 ;r++){
//            printf("%d + i %d \n",refNPSSgrid_freqDomain[r][q].real, refNPSSgrid_freqDomain[r][q].imag);
//        }
//            printf("\n");
//    }


/////////////////////////////////////////////////////////////////////////////////////////////
    FILE *rm3;//importing rx_td_vector
        int ss3;
        float buftd1[10000];
        rm3 = fopen("rx_td_vector_fixed.dat", "r");

        while (!feof(rm3)) {
            for(ss3=0;ss3<3840;ss3++)
            {
                fscanf(rm3,"%f",&buftd1[ss3]);
            }
        }
    fclose(rm3);

    for(m=0;m<3840;m++){    //reassigning

        td_rxNPSSgrid[m] = buftd1[m];
       // printf("%d\n", td_rxNPSSgrid[m]);
    }

    for(i=0;i<sizeof(td_rxNPSSgrid)/sizeof(td_rxNPSSgrid[0]);i++){  //orgainzing into real and complex

        if(i%2 == 0){
            rxNPSS_timeDomain[i/2].real = td_rxNPSSgrid[i];
        }

        if(i%2 == 1){
            rxNPSS_timeDomain[(i-1)/2].imag = td_rxNPSSgrid[i];
        }
    }

////    display received grid
//    printf("\n\n RECEIVED time domain waveform \n\n");
//    for(l=0;l<1920 ;l++){
//
//
//            printf("%d + %di ",rxNPSS_timeDomain[l].real,rxNPSS_timeDomain[l].imag);
//
//    printf("\n");
//    }
//
//

///////////////////////////

    //remove cyclic prefix from 1920 length txWaveform1

    int var = 1;

    i = 0;
    for(var=0;var<14;var++){ //remove cyclic prefix from 1920 length rxNPSS_timeDomain
        if(var==0||var==7){

            i=i+10;
            for(q=0;q<128;q++){  //10 samples for 1st and 8th symbol and 9 for the rest
                rxNPSS_timeDomain1[q][var].real = rxNPSS_timeDomain[i+q].real;
                rxNPSS_timeDomain1[q][var].imag = rxNPSS_timeDomain[i+q].imag;


            }

        }
        else{

            i=i+9;
            for(q=0;q<128;q++){
                rxNPSS_timeDomain1[q][var].real = rxNPSS_timeDomain[i+q].real;
                rxNPSS_timeDomain1[q][var].imag = rxNPSS_timeDomain[i+q].imag;
            }

        }
        i = i+128;
    }

////    //display time domain received grid without cp
//    printf("\n\n RECEIVED time domain waveform without cp \n\n");
//    for(n=3;n<4;n++){
//        for(l=0;l<128 ;l++){
//
//            printf("%d + %di \n",rxNPSS_timeDomain1[l][n].real,rxNPSS_timeDomain1[l][n].imag);
//}
//
//    printf("\n");
//    }



    for(q=3;q<14;q++){      //dft of the imported time-domain grid rxNPSSgrid_timeDomain to get freq domain received signal

        for(r=0;r<128 ;r++){


            colRxGridtd[r].real = rxNPSS_timeDomain1[r][q].real;
            colRxGridtd[r].imag = rxNPSS_timeDomain1[r][q].imag;

        }

        fnDFT(colRxGridtd, 128);
        for(s=0;s<11;s++){
            rxGrid[s][q-3].real = Xx[s].real;    /////////////////SAME AS rxGrid but shifted
            rxGrid[s][q-3].imag = Xx[s].imag;
          // printf(" %d %d + %di\n ",s,Xx[s].real,Xx[s].imag);
        }
        //printf("\n");
    }



////display refNPSSgrid_freqDomain grid----------------fd_rx_our
//    for(q=0;q<11;q++){
//        for(r=0;r<11 ;r++){
//            printf("%d + i %d \n",refNPSSgrid_freqDomain[r][q].real, refNPSSgrid_freqDomain[r][q].imag);
//        }
//            printf("\n");
//    }
////////////////////////////////////DOPPLER ESTIMATION/////////////////////////////////////////

    for(j=0;j<11;j++){  //creating z
        for(i=0;i<11;i++){


            if(j==4||j==5||j==9){

                z[i][j].real= -1*rxGrid[i][j].real;
                z[i][j].imag= -1*rxGrid[i][j].imag;
            }
            else{
                z[i][j].real= rxGrid[i][j].real;
                z[i][j].imag= rxGrid[i][j].imag;
            }
        }
    }

    complex colRxGrid1[128]={0};
    for(q=0;q<11;q++){  //take idft of z

        for(r=0;r<11 ;r++){


            colRxGrid1[r].real = z[r][q].real;
            colRxGrid1[r].imag = z[r][q].imag;

        }

        fnIDFT(colRxGrid1, 128);

        for(s=0;s<128;s++){


            zIDFT[s][q].real = X[s].real;
            zIDFT[s][q].imag = X[s].imag;
        }
    }


////     display zIDFT grid
//    printf("\n\n zIDFT GRID \n\n");
//    for(l=5;l<6 ;l++){
//        for(n=0;n<128 ;n++){  //n is supposed to be till 128
//
//            printf("%d + %di ",zIDFT[n][l].real,zIDFT[n][l].imag);
//            printf("\n");
//        }
//    printf("\n");
//    }


    printf("\n\n\n");
    complex consec_corr[128];
    complex consec_corr_sum[11] = {0};
    complex zConj;
    float d_est[11];
    float d_est_angle[11];
    float d_est_sum;
    float dopplerEstimate={0};

    for(i=0;i<1;i++){      //estimate doppler
        consec_corr_sum[i].real = 0;
        consec_corr_sum[i].imag = 0;
        for(j=0;j<12;j++){
            zConj = fnConj(zIDFT[j][i]);

            consec_corr[j] = fnCmplxMult(zIDFT[j][i+1] , zConj);
            printf("%d: consec_corr_sum = %d + i%d\n",i,consec_corr[j].real,consec_corr[j].imag );
            consec_corr_sum[i] = fnCmplxAdd(consec_corr_sum[i] ,consec_corr[j]);

        }

           // printf("%d: consec_corr_sum = %d + i%d\n",i,consec_corr_sum[i].real,consec_corr_sum[i].imag );
        d_est_angle[i] = atan((float)consec_corr_sum[i].imag/(float)consec_corr_sum[i].real);
        //printf("d_est angle = %f \n",d_est_angle[i] );
        d_est[i] = d_est_angle[i]* 1920000/(2*pi*128);
        // printf("d_est = %f \n",d_est[i] );
    }

    for(i=0;i<11;i++){
        d_est_sum += d_est[i];
    }
    dopplerEstimate = d_est_sum/10;
    printf(" dopplerEstimate = %f \n",dopplerEstimate);

    printf("\n\n\n");


    /////////////////////////////DOPPLER CORRECTION///////////////////////////



   // %Correct the received signal for Doppler shift
     //       rxWaveform_corrected = rxWaveform.*exp(2*pi*(d_est1(iter)/fs).*[0:length(rxWaveform)-1]');cos((2*pi*k*n)/N)
    for(i=0;i<1920;i++){

        rxNPSS_timeDomain_corrected[i].real = rxNPSS_timeDomain[i].real*cos(2*pi*dopplerEstimate/1920000) + rxNPSS_timeDomain[i].imag*sin(2*pi*dopplerEstimate/1920000);
        rxNPSS_timeDomain_corrected[i].imag = rxNPSS_timeDomain[i].imag*cos(2*pi*dopplerEstimate/1920000) - rxNPSS_timeDomain[i].real*sin(2*pi*dopplerEstimate/1920000);
    //printf("%d +%di \n", rxNPSS_timeDomain_corrected[i].real, rxNPSS_timeDomain_corrected[i].imag);
    }


 //remove cyclic prefix from 1920 length txWaveform1

//    int var = 1;

    i = 0;
    for(var=0;var<14;var++){ //remove cyclic prefix from 1920 length txWaveform1

        if(var==0||var==7){

            i=i+10;
            for(q=0;q<128;q++){  //10 samples for 1st and 8th symbol and 9 for the rest
                rxNPSS_timeDomain1_corrected[q][var].real = rxNPSS_timeDomain_corrected[i+q].real;
                rxNPSS_timeDomain1_corrected[q][var].imag = rxNPSS_timeDomain_corrected[i+q].imag;


            }

        }
        else{

            i=i+9;
            for(q=0;q<128;q++){
                rxNPSS_timeDomain1_corrected[q][var].real = rxNPSS_timeDomain_corrected[i+q].real;
                rxNPSS_timeDomain1_corrected[q][var].imag = rxNPSS_timeDomain_corrected[i+q].imag;
            }

        }
        i = i+128;
    }

////    //display time domain received grid
//    printf("\n\n RECEIVED time domain waveform without cp \n\n");
//     for(n=0;n<14;n++){
//        for(l=0;l<128 ;l++){
//
//
//
//            printf("%d + %di \n",rxNPSS_timeDomain1[l][n].real,rxNPSS_timeDomain1[l][n].imag);
//}
//
//    printf("\n");
//    }


    complex colRxGridtd_corrected[128]={0};
    for(q=3;q<14;q++){      //dft of the corrected time-domain grid rxNPSSgrid_timeDomain_coorected to get freq domain corrected received signal

        for(r=0;r<128 ;r++){


            colRxGridtd_corrected[r].real = rxNPSS_timeDomain1_corrected[r][q].real;
            colRxGridtd_corrected[r].imag = rxNPSS_timeDomain1_corrected[r][q].imag;

        }

        fnDFT(colRxGridtd_corrected, 128);
        for(s=0;s<11;s++){
            rxGrid_corrected[s][q-3].real = Xx[s].real;    /////////////////SAME AS corrected rxGrid but shifted
            rxGrid_corrected[s][q-3].imag = Xx[s].imag;
          //  printf(" %d %d + %di\n ",s,Xx[s].real,Xx[s].imag);
        }
    }
//////////////////////////////////TIME OFFSET ESTIMATION////////////////////////////////////////
//////////////////////////////////////CORRELATOR//////////////////////////////////////////
    complex colRefGrid[128]={0};
    complex colRxGrid[128]={0};

    for(q=0;q<11 ;q++){
        for(r=0;r<11 ;r++){

            colRxGrid[r].real = rxGrid_corrected[r][q].real ;
            colRxGrid[r].imag = rxGrid_corrected[r][q].imag ;

            colRefGrid[r].real = refNPSSgrid_freqDomain[r][q].real;
            colRefGrid[r].imag = refNPSSgrid_freqDomain[r][q].imag;

          //  printf("%d:  %d + i%d , %d + i%d\n", r, colRxGrid[r].real, colRxGrid[r].imag , colRefGrid[r].real, colRefGrid[r].imag  );
        }
        correlator(colRxGrid,colRefGrid,11,DFT_SIZE);
        fnIDFT(corrOut, 128);

        for(i=0;i<128 ;i++){
        //printf("%d:  %f + i%f \n", i, corrOut[i].real, corrOut[i].imag   );

            absSquare_idft[i] = absSquare_idft[i] + absoluteSquare(X[i]);

            }
    }


    maxParam maxABS = maxParameters(absSquare_idft, DFT_SIZE);
    printf("maxValue = %d , maxIndex = %d \n",maxABS.maxValue, maxABS.maxIndex);




    return 0;
}


