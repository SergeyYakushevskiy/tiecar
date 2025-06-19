import dstu.csae.galois.GaloisField;
import dstu.csae.galois.extended.GaloisFieldExtension;
import dstu.csae.polynomial.FieldPolynomial;
import dstu.csae.polynomial.Polynomial;
import dstu.csae.rs_code.RSCode;
import dstu.csae.rs_code.RSFCoder;
import dstu.csae.utils.FieldPolynomialGenerator;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class CodecTest {

    private static final RSFCoder coder = getGF9Coder();
    private static final FieldPolynomialGenerator generator = new FieldPolynomialGenerator(
            coder.getCode().getField(), coder.getCode().getK());

    public static Stream<FieldPolynomial> testGF9(){
        System.out.println(coder);
        List<FieldPolynomial> polynomialList = new ArrayList<>();
        while (generator.hasNext()){
            polynomialList.add(generator.next());
        }
        Collections.sort(polynomialList);
        return polynomialList.stream();
    }

    @Test
    public void testAlpha(){
        GaloisFieldExtension field = coder.getCode().getField();
        int n = coder.getCode().getN();
        int k = coder.getCode().getK();
        int alpha = coder.getCode().getAlpha();
        HashSet<Integer> values = new HashSet<>();
        for(int i = 0; i < n; i ++){
            int alphaInDeg = field.powMod(alpha, i);
            values.add(alphaInDeg);
        }
        assertEquals(field.getCharacteristic() - 1, values.size());
    }

    @ParameterizedTest
    @MethodSource("testGF9")
    public void testFieldPolynomial(FieldPolynomial p){
        FieldPolynomial encoded = coder.encode(p);
        FieldPolynomial decoded = coder.decode(encoded);
        assertEquals(p, decoded);
    }


    private static RSFCoder getGF9Coder(){
        GaloisField field = new GaloisField(3);
        GaloisFieldExtension fieldExtension = new GaloisFieldExtension(field, new Polynomial(new int[]{1, 0, 1}));
        int n = fieldExtension.getCharacteristic() - 1;
        int k = n / 2;
        RSCode code = new RSCode(fieldExtension, n, k);
        return new RSFCoder(code);
    }

}
