// radix4_256.c
// 256 point FFT radix 4 Algorithm
// 12.10.2016 by Dmitry Chistyakov 
//
// Updates:
//      12.10.2016 - release.
//      27.04.2017 - optimized using of dragonfly's, deleted math functions.
// Time of calculation (FFT stages + sort stage + calculation amplitudes)
// is 2.4 ms (for Clock = 160MHz (RISC)).
//*****************************************************************************
//*****************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <stdint.h>

// Arrays for Real and Imaginary parts of spectrum
float fft_r[256];    
float fft_j[256];    

// Arrays for Real and Imaginary parts of complex exponents
float exp_r[320];
float exp_j[320];

// Size of input data 
int framesize = 64; 

// Intermediate results
float   rimd0, jimd0,
        rimd1, jimd1,
        rimd2, jimd2,
        rimd3, jimd3;

// Variables for cycles
int i, kk, n, step, count, dis;
int d = 0;
int shift;
int mult;

// Indexes
int index0, index1, index2, index3;

// Variables for jumping to another frame size.
int jumpto16 = 0;
int jumpto64 = 0;

// Function for getting table of complex exponents
void fft_expcalculation(){
    float arg;
    float pi = 3.1415926535897932384626433832795;
    int n;
    int m = 1;
    // Calculate exponents for 16 and 64 point stage
    for (n=0; n<16; n++){
        arg = 0;
        exp_r[0+n] = 1;
        exp_j[0+n] = 0;

        arg = (float)(2*pi*n)/64.0;
        exp_r[16+n] = cos(arg); 
        exp_j[16+n] = sin(-arg);

        arg = arg * 2;
        exp_r[32+n] = cos(arg); 
        exp_j[32+n] = sin(-arg);

        arg = arg * 3 / 2;
        exp_r[48+n] = cos(arg); 
        exp_j[48+n] = sin(-arg);
    }
    // Calculate exponents for 256 point stage
    for (n=0; n<64; n++){
        arg = 0;
        exp_r[64+n] = 1;
        exp_j[64+n] = 0;

        arg = (float)(2*pi*n)/256.0;
        exp_r[128+n] = cos(arg); 
        exp_j[128+n] = sin(-arg);

        arg = arg * 2;
        exp_r[192+n] = cos(arg); 
        exp_j[192+n] = sin(-arg);

        arg = arg * 3 / 2;
        exp_r[256+n] = cos(arg); 
        exp_j[256+n] = sin(-arg);
    }
}

int main()
{   

    int quarter = framesize / 4;

    if (framesize == 64) {
        jumpto64 = 1000;
    }
    else if (framesize == 16){
        jumpto16 = 1000;
        jumpto64 = 1000;
    }

    // Get table of exponents
    fft_expcalculation();

    // Frame forming stage (getting data for analysis)
    for (i=0; i<framesize; i++){
        fft_r[i]=i; // Put sample into array
        fft_j[i]=0;
    }

uint32_t high, low;
__asm__ __volatile__ (
"rdtsc\n"
"movl %%edx, %0\n"
"movl %%eax, %1\n"
: "=r" (high), "=r" (low)
:: "%rax", "%rbx", "%rcx", "%rdx"
);
uint64_t ticks = ((uint64_t)high << 32) | low;

    for (n=0; n<64; n++){

        if (framesize!=256) continue;
        index0 = n;
        index1 = index0+64;
        index2 = index1+64;
        index3 = index2+64;

		rimd0 = fft_r[index0];
		jimd0 = fft_j[index0];
		rimd1 = fft_r[index1];
		jimd1 = fft_j[index1];
		rimd2 = fft_r[index2];
		jimd2 = fft_j[index2];
		rimd3 = fft_r[index3];
		jimd3 = fft_j[index3];

        // "Dragonfly" operations
        fft_r[index0] = rimd0 + rimd1 + rimd2 + rimd3;
        fft_j[index0] = jimd0 + jimd1 + jimd2 + jimd3;

        fft_r[index1] = rimd0 + jimd1 - rimd2 - jimd3;
        fft_j[index1] = jimd0 - rimd1 - jimd2 + rimd3;

        fft_r[index2] = rimd0 - rimd1 + rimd2 - rimd3;
        fft_j[index2] = jimd0 - jimd1 + jimd2 - jimd3;

        fft_r[index3] = rimd0 - jimd1 - rimd2 + jimd3;
        fft_j[index3] = jimd0 + rimd1 - jimd2 - rimd3;

                
        rimd0 = fft_r[index0];
        jimd0 = fft_j[index0];
        rimd1 = fft_r[index1];
        jimd1 = fft_j[index1];
        rimd2 = fft_r[index2];
        jimd2 = fft_j[index2];
        rimd3 = fft_r[index3];
        jimd3 = fft_j[index3];

        // Multiply on complex expornents
        fft_r[index0] = rimd0 * exp_r[index0+64] - jimd0 * exp_j[index0+64];
        fft_j[index0] = rimd0 * exp_j[index0+64] + jimd0 * exp_r[index0+64];
        fft_r[index1] = rimd1 * exp_r[index1+64] - jimd1 * exp_j[index1+64];
        fft_j[index1] = rimd1 * exp_j[index1+64] + jimd1 * exp_r[index1+64];
        fft_r[index2] = rimd2 * exp_r[index2+64] - jimd2 * exp_j[index2+64];
        fft_j[index2] = rimd2 * exp_j[index2+64] + jimd2 * exp_r[index2+64];
        fft_r[index3] = rimd3 * exp_r[index3+64] - jimd3 * exp_j[index3+64];
        fft_j[index3] = rimd3 * exp_j[index3+64] + jimd3 * exp_r[index3+64];
    }

    for (shift=0; shift<256; shift+=64+jumpto64){
        // Sequently perform 64 point fft for 8 times
        for (count=4; count>=1; count-=3){
            // Perform second stage when count = 4.
            for (kk=0; kk<64/count; kk+=16+jumpto16){
                // Perform third stage when count = 1, step = 1.
                for (step=1; step<=4/count; step+=3){
                    // Perform fourth stage when count = 1, step = 4.
                    mult = count*4/step;
                    for (n=0; n<4*count; n++){

                            dis=n * step + kk + shift;

                            index0 = dis;
                            index1 = index0+mult;
                            index2 = index1+mult;
                            index3 = index2+mult;

                            rimd0 = fft_r[index0];
                            jimd0 = fft_j[index0];
                            rimd1 = fft_r[index1];
                            jimd1 = fft_j[index1];
                            rimd2 = fft_r[index2];
                            jimd2 = fft_j[index2];
                            rimd3 = fft_r[index3];
                            jimd3 = fft_j[index3];

                            // "Dragonfly" operations
                            fft_r[index0] = rimd0 + rimd1 + rimd2 + rimd3;
                            fft_j[index0] = jimd0 + jimd1 + jimd2 + jimd3;

                            fft_r[index1] = rimd0 + jimd1 - rimd2 - jimd3;
                            fft_j[index1] = jimd0 - rimd1 - jimd2 + rimd3;

                            fft_r[index2] = rimd0 - rimd1 + rimd2 - rimd3;
                            fft_j[index2] = jimd0 - jimd1 + jimd2 - jimd3;

                            fft_r[index3] = rimd0 - jimd1 - rimd2 + jimd3;
                            fft_j[index3] = jimd0 + rimd1 - jimd2 - rimd3;

                            // Multiply on complex expornents
                            if (step!=4){
                                    
                                rimd0 = fft_r[index0];
                                jimd0 = fft_j[index0];
                                rimd1 = fft_r[index1];
                                jimd1 = fft_j[index1];
                                rimd2 = fft_r[index2];
                                jimd2 = fft_j[index2];
                                rimd3 = fft_r[index3];
                                jimd3 = fft_j[index3];

                                if (count==4){ 
                                    dis = 0; 
                                }
                                else{ 
                                    dis=d-kk; 
                                }
                                dis-=shift;

                                fft_r[index0] = rimd0 * exp_r[index0+dis] - jimd0 * exp_j[index0+dis];
                                fft_j[index0] = rimd0 * exp_j[index0+dis] + jimd0 * exp_r[index0+dis];
                                fft_r[index1] = rimd1 * exp_r[index1+dis] - jimd1 * exp_j[index1+dis];
                                fft_j[index1] = rimd1 * exp_j[index1+dis] + jimd1 * exp_r[index1+dis];
                                fft_r[index2] = rimd2 * exp_r[index2+dis] - jimd2 * exp_j[index2+dis];
                                fft_j[index2] = rimd2 * exp_j[index2+dis] + jimd2 * exp_r[index2+dis];
                                fft_r[index3] = rimd3 * exp_r[index3+dis] - jimd3 * exp_j[index3+dis];
                                fft_j[index3] = rimd3 * exp_j[index3+dis] + jimd3 * exp_r[index3+dis];
                                d+=15;
                            }
                    }
                    d=0;
                }
            }
        }
    }

    // fft end.
    // Calculate amplitude spectrum.
   
    for  (i=0; i<framesize; i++){
    	fft_j[i] = sqrt(fft_r[i] * fft_r[i] + fft_j[i] * fft_j[i]);
    }

    // Sort stage.
    for (kk=0; kk<64; kk+=16){
        for (shift=0; shift<256; shift+=64+jumpto64){
        	for (i=0; i<4; i++){
                index0 = 4*i+shift+kk;
                index1 = i*quarter/4+d;
        		fft_r[index0]=fft_j[index1];
        		fft_r[index0+1]=fft_j[index1+quarter];
        		fft_r[index0+2]=fft_j[index1+2*quarter];
        		fft_r[index0+3]=fft_j[index1+3*quarter];
            }
            d++;
        }
    }

__asm__ __volatile__ (
"rdtsc\n"
"movl %%edx, %0\n"
"movl %%eax, %1\n"
: "=r" (high), "=r" (low)
:: "%rax", "%rbx", "%rcx", "%rdx"
);
ticks = (((uint64_t)high << 32) | low) - ticks;
float timeest = ticks * 1.3 / 30501;
printf("Elapsed ticks: %i, Time estimate (ms): %f\n", ticks, timeest);

 

    // See results
    for (i=0; i<framesize; i++){
        printf("  %f\n", fft_r[i]);
    }

    
    
    return 0;
}
