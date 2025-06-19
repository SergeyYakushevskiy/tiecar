package dstu.csae.utils;

import dstu.csae.galois.extended.GaloisFieldExtension;
import dstu.csae.polynomial.FieldPolynomial;

import java.util.Arrays;

public class BerlekampMassey {

    public static FieldPolynomial findErrorLocatorPolynomial(GaloisFieldExtension field, int[] syndrome) {
        int N = syndrome.length;
        int[] C = new int[N];
        int[] B = new int[N];
        C[0] = 1;
        B[0] = 1;
        int L = 0;
        int m = 1;
        int b = 1;
        for(int n = 0; n < N; n ++){
            int d = syndrome[n];
            for(int j = 1; j <= L; j ++){
                if(n - j < 0){
                    continue;
                }
                d = field.add(d, field.multiply(syndrome[n - j], C[j]));
            }
            if(d == 0){
                m++;
                continue;
            }

            int[] T = Arrays.copyOf(C, C.length);
            int coefficient = field.divide(d, b);
            int[] scaled = new int[N];
            for(int i = 0; i < B.length; i ++) {
                if(i + m >= N){
                    continue;
                }
                scaled[i + m] = field.multiply(coefficient, B[i]);
            }
            for(int i = 0; i < C.length; i ++) {
                C[i] = field.subtract(C[i], scaled[i]);
            }

            if (2 * L <= n) {
                L = n + 1 - L;
                B = T;
                b = d;
                m = 1;
            } else {
                m++;
            }
        }
        return new FieldPolynomial(field, Arrays.copyOf(C, L + 1));
    }
}
