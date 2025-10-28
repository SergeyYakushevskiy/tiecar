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
        for(int n = 0; n < N; n ++) {
            int d = syndrome[n];
            for(int j = 1; j <= L; j ++){
                int monomial = field.multiply(C[j], syndrome[n - j]);
                d = field.add(d, monomial);
            }
            if (d == 0) {
                m++;
                continue;
            }
            int[] T = Arrays.copyOf(C, C.length);
            int coefficient = field.divide(d, b);
            int[] scale = new int[B.length + m];
            for(int j = 0; j < B.length; j ++){
                scale[j + m] = field.inverseOfAddition(field.multiply(coefficient, B[j]));
            }
            for(int j = 0; j < C.length; j ++){
                scale[j] = field.add(C[j], scale[j]);
            }
            C = scale;
            if (2 * L <= n) {
                B = T;
                L = n + 1 - L;
                b = d;
                m = 1;
            }else {
                m ++;
            }
        }
        return new FieldPolynomial(field, Arrays.copyOf(C, L + 1));
    }
}
