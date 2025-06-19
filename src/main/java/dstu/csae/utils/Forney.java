package dstu.csae.utils;

import dstu.csae.code.Code;
import dstu.csae.galois.extended.GaloisFieldExtension;
import dstu.csae.polynomial.FieldPolynomial;
import lombok.NonNull;

public class Forney {

    public static int[] forneyAlgorithm(@NonNull Code code,
                                        @NonNull FieldPolynomial syndromes,
                                        @NonNull FieldPolynomial locator,
                                        @NonNull int[] errorPositions) {
        GaloisFieldExtension field = syndromes.getField();
        int invAlpha = field.inverseOfMultiplication(code.getAlpha());
        int[] errorValues = new int[errorPositions.length];
        FieldPolynomial omega = generateOmega(code, syndromes, locator);
        FieldPolynomial locatorDerivative = locator.derivative();
        for (int i = 0; i < errorPositions.length; i++) {
            int X = field.powMod(invAlpha, errorPositions[i]);
            int invX = field.inverseOfMultiplication(X);
            int numerator = omega.evaluate(X);
            int denominator = locatorDerivative.evaluate(X);
            errorValues[i] = field.inverseOfAddition(field.divide(numerator, denominator));
        }
        return errorValues;
    }

    public static FieldPolynomial generateOmega(Code code, FieldPolynomial syndromes, FieldPolynomial locator){
        GaloisFieldExtension field = code.getField();
        int r = code.getR();
        int[] modC = new int[r + 1];
        modC[r] = 1;
        FieldPolynomial mod = new FieldPolynomial(field, modC);
        return syndromes.multiply(locator).mod(mod);
    }
}
