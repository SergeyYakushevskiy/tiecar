package dstu.csae.rs_code;

import dstu.csae.codec.Codec;
import dstu.csae.galois.GaloisField;
import dstu.csae.galois.extended.GaloisFieldExtension;
import dstu.csae.mathutils.MathUtils;
import dstu.csae.polynomial.FieldPolynomial;
import dstu.csae.utils.BerlekampMassey;
import dstu.csae.utils.Chien;
import dstu.csae.utils.Forney;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;


public class RSFCoder extends Codec<RSCode, FieldPolynomial, FieldPolynomial> {

    private static final Logger logger = LogManager.getLogger(RSFCoder.class);

    public RSFCoder(RSCode code) {
        super(code);
    }

    @Override
    public FieldPolynomial encode(FieldPolynomial input)
            throws IllegalArgumentException{
        int k = code.getK();
        if(input.getDegree() >= k){
            throw new IllegalArgumentException(String.format(
                "Степень полинома %s должна быть меньше k = %d", input, k
            ));
        }
        logger.info("КОДИРОВАНИЕ");
        int n = code.getN();
        GaloisFieldExtension field = code.getField();
        int alpha = code.getAlpha();
        int[] encoded = new int[n];
        logger.info("Принятый полином: {}", input);
        for(int i = 0; i < n; i ++){
            int power = field.powMod(alpha, i);
            encoded[i] = input.evaluate(power);
        }
        FieldPolynomial output = new FieldPolynomial(field, encoded);
        logger.info("Закодированный полином: {}", output);
        return output;
    }

    @Override
    public FieldPolynomial decode(FieldPolynomial input)
            throws IllegalArgumentException{
        int n = code.getN();
        if(input.getDegree() >= n){
            throw new IllegalArgumentException(String.format(
                    "Степень полинома %s долна быть меньше n = %d", input, n
            ));
        }
        logger.info("ДЕКОДИРОВАНИЕ");
        logger.info("Принятый полином: {}", input);
        GaloisFieldExtension field = input.getField();
        int[] syndromes = computeSyndromes(input);
        boolean allZero = true;
        for (int s : syndromes) {
            if (s != 0) {
                allZero = false;
                break;
            }
        }
        logger.info("Синдромы: {}", Arrays.toString(syndromes));
        if (allZero) {
            FieldPolynomial output = getInformationalCoefficients(input);
            logger.info("Восстановленный полином: {}", output);
            return output;
        }
        FieldPolynomial errorLocatorPolynomial = BerlekampMassey.findErrorLocatorPolynomial(field, syndromes);
        logger.info("Локатор ошибок: {}", errorLocatorPolynomial);
        int[] errorPositions = Chien.findErrorPositions(code, errorLocatorPolynomial);
        logger.info("Позиции ошибок: {}", Arrays.toString(errorPositions));
        int[] errorValues = Forney.forneyAlgorithm(
                code,
                new FieldPolynomial(field, syndromes),
                errorLocatorPolynomial,
                errorPositions
        );
        logger.info("Значения ошибок: {}", Arrays.toString(errorValues));
        input = correctErrors(input, errorPositions, errorValues);
        logger.info("Полином после исправления ошибок: {}", input);
        FieldPolynomial output = getInformationalCoefficients(input);
        logger.info("Восстановленный полином: {}", output);
        return output;
    }

    private FieldPolynomial getInformationalCoefficients(FieldPolynomial received){
        GaloisFieldExtension fieldExtension = received.getField();
        GaloisField field = fieldExtension.getGaloisField();
        int n = code.getN();
        int k = code.getK();
        int alpha = code.getAlpha();
        int invAlpha = fieldExtension.inverseOfMultiplication(alpha);
        int invN = field.inverseOfMultiplication(field.bringToField(n));
        logger.info(invN);
        int[] decoded = new int[n];

        for(int i = 0; i < k; i ++){
            for(int j = 0; j < n; j ++){
                int alphaInDeg = fieldExtension.powMod(invAlpha, i * j);
                decoded[i] = fieldExtension.add(decoded[i], fieldExtension.multiply(received.get(j), alphaInDeg));
            }
            decoded[i] = fieldExtension.multiply(decoded[i], invN);
        }
        return new FieldPolynomial(fieldExtension, decoded);
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

    private FieldPolynomial correctErrors(FieldPolynomial received,
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
                "RSFCoder",
                code.toString(),
                "α(x) = " + code.getField().get(code.getAlpha()),
                "g(x) = " + code.getG(),
                "h(x) = " + code.getH()
        ));
    }
}
