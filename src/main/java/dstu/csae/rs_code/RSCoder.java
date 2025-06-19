package dstu.csae.rs_code;

import dstu.csae.codec.Codec;
import dstu.csae.galois.extended.GaloisFieldExtension;
import dstu.csae.polynomial.FieldPolynomial;
import dstu.csae.utils.BerlekampMassey;
import dstu.csae.utils.Chien;
import dstu.csae.utils.Forney;
import lombok.NonNull;

import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

public class RSCoder extends Codec<RSCode, FieldPolynomial, FieldPolynomial> {

    public RSCoder(RSCode code){
        super(code);
    }

    @Override
    public FieldPolynomial encode(@NonNull FieldPolynomial input)
            throws IllegalArgumentException{
        int k = code.getK();
        if(input.getDegree() >= k){
            throw new IllegalArgumentException(String.format(
                    "Степень полинома %s должна быть меньше значения k = %d",
                    input,k
            ));
        }
        int r = code.getR();
        GaloisFieldExtension field = code.getField();
        FieldPolynomial generator = new FieldPolynomial(field, code.getG().getCoefficients());

        FieldPolynomial shifted = new FieldPolynomial(field, input.getDegree() + r);
        for (int i = 0; i <= input.getDegree(); i++) {
            shifted.set(i + r, input.get(i));
        }
        FieldPolynomial remainder = shifted.mod(generator);
        FieldPolynomial multiple = shifted.add(remainder);
        return multiple;
    }

    @Override
    public FieldPolynomial decode(@NonNull FieldPolynomial input)
            throws IllegalArgumentException{
        GaloisFieldExtension field = code.getField();
        int n = code.getN();
        if(input.getDegree() >= n){
            throw new IllegalArgumentException(String.format(
                    "Степень полинома %s долна быть меньше n = %d", input, n
            ));
        }
        int[] syndromes = computeSyndromes(input);
        boolean allZero = true;
        for (int s : syndromes) {
            if (s != 0) {
                allZero = false;
                break;
            }
        }

        if (allZero) {
            return getInformationalCoefficients(input);
        }
        FieldPolynomial errorLocatorPolynomial = BerlekampMassey.findErrorLocatorPolynomial(field, syndromes);
        int[] errorPositions = Chien.findErrorPositions(code, errorLocatorPolynomial);
        int[] errorValues = Forney.forneyAlgorithm(
                code,
                new FieldPolynomial(field, syndromes),
                errorLocatorPolynomial,
                errorPositions
        );
        input = correctErrors(input, errorPositions, errorValues);
        return getInformationalCoefficients(input);
    }

    private FieldPolynomial getInformationalCoefficients(FieldPolynomial received){
        GaloisFieldExtension field = received.getField();
        int k = code.getK();
        int n = code.getN();
        FieldPolynomial decoded = new FieldPolynomial(field, k);
        int[] informationCoefficients = Arrays.copyOfRange(received.getCoefficients(), n - k, n);
        IntStream.range(0, informationCoefficients.length)
                .forEach(index -> decoded.set(index, informationCoefficients[index]));
        return decoded;
    }

    private int[] computeSyndromes(FieldPolynomial received) {
        GaloisFieldExtension field = received.getField();
        int r = code.getR();
        int[] syndromes = new int[r];
        int alpha = code.getAlpha();

        for (int i = 1; i <= r; i++) {
            int power = field.powMod(alpha, i);
            syndromes[i - 1] = received.evaluate(power);
        }

        return syndromes;
    }

    public FieldPolynomial correctErrors(FieldPolynomial received,
                                         int[] errorPositions,
                                         int[] errorValues) {
        GaloisFieldExtension field = code.getField();
        int[] corrected = new int[code.getN()];
        IntStream.range(0, received.getDegree() + 1)
                .forEach(i -> corrected[i] = received.get(i));
        for (int i = 0; i < errorPositions.length; i++) {
            int position = errorPositions[i];
            if(position < 0 || position >= corrected.length){
                continue;
            }
            int errorValue = errorValues[i];
            corrected[position] = field.subtract(corrected[position], errorValue);
        }
        return new FieldPolynomial(field, corrected);
    }

    @Override
    public String toString() {
        return String.join("\n", List.of(
                "RSCoder",
                code.toString(),
                "α(x) = " + code.getField().get(code.getAlpha()),
                "g(x) = " + code.getG(),
                "h(x) = " + code.getH()
        ));
    }
}
