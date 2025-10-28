import dstu.csae.galois.GaloisField;
import dstu.csae.galois.extended.GaloisFieldExtension;
import dstu.csae.polynomial.FieldPolynomial;
import dstu.csae.polynomial.Polynomial;
import dstu.csae.rscode.RSCode;
import dstu.csae.rscode.RSFCoder;
import dstu.csae.utils.FieldPolynomialGenerator;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

public class CodecTest {

    private static final RSFCoder coder = getGF9Coder();
    private static final FieldPolynomialGenerator generator = new FieldPolynomialGenerator(
            coder.getCode().getField(), coder.getCode().getK());

    // Воспроизводимый Random с фиксированным зерном
    private static final Random random = new Random(12345L); // Зерно можно изменить


    public static Stream<FieldPolynomial> testGF9(){
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

    // Новый тест: проверка корректности декодирования с одной случайной ошибкой
    @ParameterizedTest
    @MethodSource("testGF9")
    public void testFieldPolynomialWithSingleError(FieldPolynomial p) {
        System.out.println(coder);
        FieldPolynomial encoded = coder.encode(p);
        // Вносим одну случайную ошибку в закодированное слово
        FieldPolynomial encodedWithError = addSingleRandomError(encoded, random);

        // Проверяем, что ошибка действительно внесена
        assertNotEquals(encoded, encodedWithError, "Ошибка не была внесена в полином");

        // Декодируем слово с ошибкой
        FieldPolynomial decoded = coder.decode(encodedWithError);

        // Проверяем, что результат декодирования совпадает с исходным информационным полиномом
        assertEquals(p, decoded, "Декодер не смог исправить внесённую ошибку");
    }


    private static RSFCoder getGF9Coder(){
        GaloisField field = new GaloisField(3);
        GaloisFieldExtension fieldExtension = new GaloisFieldExtension(field, new Polynomial(new int[]{1, 0, 1}));
        int n = fieldExtension.getCharacteristic() - 1;
        int k = n / 2;
        RSCode code = new RSCode(fieldExtension, n, k);
        return new RSFCoder(code);
    }

    public static RSFCoder getGF8Coder(){
        GaloisField field = new GaloisField(2);
        GaloisFieldExtension fieldExtension = new GaloisFieldExtension(field, new Polynomial(new int[]{1, 0, 1, 1}));
        int n = fieldExtension.getCharacteristic() - 1;
        int k = n / 2;
        RSCode code = new RSCode(fieldExtension, n, k);
        return new RSFCoder(code);
    }

    public static FieldPolynomial addSingleRandomError(FieldPolynomial original, Random random) {
        GaloisFieldExtension field = original.getField();
        int characteristic = field.getCharacteristic();

        if (characteristic < 2) {
            throw new IllegalArgumentException("Характеристика поля должна быть >= 2 для внесения ошибки.");
        }

        int degree = original.getDegree();
        if (degree < 0) {
            throw new IllegalArgumentException("Нельзя внести ошибку в нулевой полином (степень < 0).");
        }

        int n = degree + 1; // Длина полинома (количество коэффициентов)

        // 1. Выбираем случайную позицию для ошибки (от 0 до n-1)
        int errorPosition = random.nextInt(n);

        // 2. Генерируем случайное ненулевое значение ошибки
        int errorValue;
        do {
            errorValue = random.nextInt(characteristic);
        } while (errorValue == 0);

        // 3. Получаем текущий коэффициент в позиции ошибки
        int originalCoefficient = original.get(errorPosition);

        // 4. Вычисляем новый коэффициент: original + error_value (в поле Галуа)
        int newCoefficient = field.add(originalCoefficient, errorValue);

        // 5. Создаем новый массив коэффициентов
        int[] newCoeffs = new int[n];
        for (int i = 0; i < n; i++) {
            newCoeffs[i] = original.get(i);
        }
        newCoeffs[errorPosition] = newCoefficient;

        // 6. Создаем и возвращаем новый полином с ошибкой
        return new FieldPolynomial(field, newCoeffs);
    }
}