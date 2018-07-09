#include"dgemm.h"
#include <immintrin.h>

void dgemm(float *a, float *b, float *c, int n)
{
	// NOTE: n%8 == 0
	__m256 vec1, vec2, curAcc, m;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
        	curAcc = _mm256_setzero_ps();
            for(int k = 0; k < n; k+=8){
                // c[i * n + j] += a[i * n  + k] * b[j * n  + k];
                vec1 = _mm256_loadu_ps(a+(i*n+k));
                vec2 = _mm256_loadu_ps(b+(j*n+k));
                m = _mm256_mul_ps(vec1, vec2);
                curAcc = _mm256_add_ps(curAcc, m);
            }
            c[i*n + j] = curAcc[0] + curAcc[1] + curAcc[2] + curAcc[3] + curAcc[4] + curAcc[5] + curAcc[6] + curAcc[7];
        }
    }
}
