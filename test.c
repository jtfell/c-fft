#include "fft.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    complex * input1 = (complex*) malloc(sizeof(struct complex_t) * 30);
    complex * input2 = (complex*) malloc(sizeof(struct complex_t) * 30);
    complex * result1, * result2;
    
    /* Init inputs */
    for (int i=0; i < 30; i++) {
        input1[i].re = (double) i;
        input1[i].im = 0.0;
        input2[i].re = (double) i;
        input2[i].im = 0.0;
    }
    
    /* Do FFT */
    result1 = FFT_CooleyTukey(input1, 30, 6, 5);
    result2 = FFT_GoodThomas(input2, 30, 6, 5);
    
    /* Compare results */
    printf("Index \t Cooley-Tukey Output \t \t Good-Thomas Output \n");
    for (int i=0; i < 30; i++) {
        printf("%d: \t %f + %fi \t %f + %fi \n", i, result1[i].re, result1[i].im, 
                result2[i].re, result2[i].im);
    }
  
    return 0;
}

