package dstu.csae.utils;

import dstu.csae.code.Code;
import dstu.csae.galois.extended.GaloisFieldExtension;
import dstu.csae.polynomial.FieldPolynomial;
import lombok.NonNull;

import java.util.HashSet;

public class Chien {

    public static int[] findErrorPositions(@NonNull Code code,
                                           @NonNull FieldPolynomial locator) {
        GaloisFieldExtension field = locator.getField();
        int invAlpha = field.inverseOfMultiplication(code.getAlpha());
        HashSet<Integer> errorPositions = new HashSet<>();
        // Важно: getCharacteristic() возвращает размер поля (q), а не базовую характеристику (p)
        int l = field.getCharacteristic() - 1;
        for (int i = 0; i < l; i++) {
            int eval = locator.evaluate(field.powMod(invAlpha, i));
            if (eval == 0) {
                errorPositions.add(i);
            }
        }
        return errorPositions.stream()
                .mapToInt(Integer::valueOf)
                .toArray();
    }

}
