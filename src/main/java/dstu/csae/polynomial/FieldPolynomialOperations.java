package dstu.csae.polynomial;

import dstu.csae.galois.extended.GaloisFieldExtension;
import lombok.NonNull;

import java.util.Arrays;

public class FieldPolynomialOperations {

    static FieldPolynomial addition(@NonNull FieldPolynomial first,
                                            @NonNull FieldPolynomial second){
        if(!first.getField().equals(second.getField())){
            return null;
        }
        GaloisFieldExtension field = first.getField();
        int[] firstC = first.getCoefficients();
        int[] secondC = second.getCoefficients();
        return new FieldPolynomial(field, addition(field, firstC, secondC));
    }

    static FieldPolynomial substraction(@NonNull FieldPolynomial reduced,
                                                  @NonNull FieldPolynomial subtracted){
        if(!reduced.getField().equals(subtracted.getField())){
            return null;
        }
        GaloisFieldExtension field = reduced.getField();
        int[] reducedC = reduced.getCoefficients();
        int[] subtractedC = subtracted.getCoefficients();
        return new FieldPolynomial(field, subtraction(field, reducedC, subtractedC));
    }

    static FieldPolynomial multiplication(@NonNull FieldPolynomial first,
                                                    @NonNull FieldPolynomial second){
        if(!first.getField().equals(second.getField())){
            return null;
        }
        GaloisFieldExtension field = first.getField();
        int[] firstC = first.getCoefficients();
        int[] secondC = second.getCoefficients();
        return new FieldPolynomial(field, multiplication(field, firstC, secondC));
    }

    static FieldPolynomial division(@NonNull FieldPolynomial divisible,
                                    @NonNull FieldPolynomial divisor){
        if(!divisible.getField().equals(divisor.getField())){
            return null;
        }
        GaloisFieldExtension field = divisible.getField();
        divisible.removeLastZero();
        divisor.removeLastZero();
        int[] divisibleC = divisible.getCoefficients();
        int[] divisorC = divisor.getCoefficients();
        return new FieldPolynomial(field, division(field, divisibleC, divisorC));
    }

    static FieldPolynomial mod(@NonNull FieldPolynomial divisible,
                               @NonNull FieldPolynomial divisor){
        if(!divisible.getField().equals(divisor.getField())){
            return null;
        }
        FieldPolynomial division = division(divisible, divisor);
        FieldPolynomial multiplication = multiplication(division, divisor);
        FieldPolynomial mod = substraction(divisible, multiplication);
        return mod;
    }

    static FieldPolynomial derivative(@NonNull FieldPolynomial polynomial)
            throws RuntimeException{
        GaloisFieldExtension field = polynomial.getField();
        int degree = polynomial.getDegree();
        if (degree == 0) {
            return new FieldPolynomial(field, new int[]{0});
        }
        int[] derivativeCoefficients = new int[degree];
        for (int i = 1; i <= degree; i++) {
            int coefficient = field.multiply(polynomial.get(i), field.getGaloisField().bringToField(i));
            derivativeCoefficients[i - 1] = coefficient;
        }
        return new FieldPolynomial(field, derivativeCoefficients);
    }

    static int evaluation(@NonNull FieldPolynomial polynomial, int x){
        GaloisFieldExtension field = polynomial.getField();
        return evaluation(field, polynomial.coefficients, x);
    }

    static FieldPolynomial interpolation(@NonNull int[] x, @NonNull FieldPolynomial y){
        GaloisFieldExtension field = y.getField();
        return new FieldPolynomial(field, interpolation(field, x, y.coefficients));
    }

    protected static  int[] addition(GaloisFieldExtension field, int[] first, int[] second){
        int maxLength = Math.max(first.length, second.length);
        int[] addition = new int[maxLength];
        System.arraycopy(first, 0, addition, 0, first.length);
        for(int i = 0; i < second.length; i ++){
            addition[i] = field.add(addition[i], second[i]);
        }
        return addition;
    }

    protected static int[] subtraction(GaloisFieldExtension field, int[] reduced, int[] subtracted){
        subtracted = Arrays.stream(subtracted).map(field::inverseOfAddition).toArray();
        return addition(field, reduced, subtracted);
    }

    protected static int[] multiplication(GaloisFieldExtension field, int[] first, int[] second){
        int[] multiplication = new int[first.length + second.length - 1];
        for(int i = 0; i < first.length; i ++){
            for(int j = 0; j < second.length; j ++){
                int multiplied = field.multiply(first[i], second[j]);
                multiplication[i + j] = field.add(multiplication[i + j], multiplied);
            }
        }
        return multiplication;
    }

    private static int[] division(GaloisFieldExtension field, int[] divisible, int[] divisor)
            throws ArithmeticException{
        if(divisible.length < divisor.length){
            return new int[]{0};
        }
        if(Arrays.stream(divisor).allMatch(x -> x == 0)){
            throw new ArithmeticException("Попытка деления на нулевой полином");
        }
        int[] division = new int[divisible.length - divisor.length + 1];
        int[] current;
        int divisibleDegree = divisible.length - 1;
        int divisorDegree = divisor.length - 1;
        while(divisibleDegree >= divisorDegree){
            current = new int[division.length];
            int divided = field.divide(divisible[divisibleDegree], divisor[divisorDegree]);
            division[divisibleDegree - divisorDegree] = divided;
            current[divisibleDegree - divisorDegree] = divided;
            divisible = subtraction(field, divisible, multiplication(field, current, divisor));
            divisibleDegree--;
        }
        return division;
    }

    private static int evaluation(GaloisFieldExtension field, int[] coefficients, int x){
        if(x == 0){
            return coefficients[0];
        }
        int result = 0;
        for(int i = 0; i < coefficients.length; i++){
            int monomial = field.powMod(x, i);
            monomial = field.multiply(coefficients[i], monomial);
            result = field.add(result, monomial);
        }
        return result;
    }

    private static int[] interpolation(GaloisFieldExtension field, int[] x, int[] y)
            throws IllegalArgumentException, ArithmeticException{
        if(x.length != y.length){
            throw new IllegalArgumentException(String.format(
                    "Количество значений x не совпадает с количеством значений y: %d != %d",
                    x.length, y.length
            ));
        }
        int[] interpolated = new int[]{0};
        int pointCount = x.length;
        for(int i = 0; i < pointCount; i ++){
            int[] multiplier = new int[]{1};
            for(int j = 0; j < pointCount; j ++){
                if(i == j) {
                    continue;
                }
                int[] numerator = new int[]{field.inverseOfAddition(x[j]), 1};
                int[] denominator = new int[]{field.inverseOfMultiplication(field.subtract(x[i], x[j]))};
                multiplier = multiplication(field, multiplier, numerator);
                multiplier = multiplication(field, multiplier, denominator);
            }
            multiplier = multiplication(field, multiplier, new int[]{y[i]});
            interpolated = addition(field, interpolated, multiplier);
        }
        return interpolated;
    }

}
