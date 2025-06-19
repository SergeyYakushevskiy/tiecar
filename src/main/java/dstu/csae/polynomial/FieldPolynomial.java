package dstu.csae.polynomial;

import dstu.csae.exceptions.EmptyCoefficientsException;
import dstu.csae.galois.extended.GaloisFieldExtension;
import lombok.Getter;
import lombok.NonNull;

import java.util.Arrays;

public class FieldPolynomial extends Polynomial {

    private static final String POLYNOMIAL_IS_NOT_NORMALIZED = "Коэффициенты %s не нормированы над полем %s";

    @Getter
    private GaloisFieldExtension field;


    public FieldPolynomial(@NonNull GaloisFieldExtension field)
            throws EmptyCoefficientsException {
        this(field, 0);
    }

    public FieldPolynomial(@NonNull GaloisFieldExtension field,
                           int degree)
            throws EmptyCoefficientsException {
        this(field, new int[degree + 1]);
    }

    public FieldPolynomial(@NonNull GaloisFieldExtension field,
                           @NonNull  int[] coefficients)
            throws EmptyCoefficientsException {
        super(coefficients);
        if(!field.isNormalized(this)){
            throw new IllegalArgumentException(String.format(
                    POLYNOMIAL_IS_NOT_NORMALIZED,
                    Arrays.toString(coefficients),
                    field));
        }
        this.field = field;
    }

    public int evaluate(int x){
        return FieldPolynomialOperations.evaluation(this, x);
    }

    public FieldPolynomial add(@NonNull FieldPolynomial polynomial){
        return FieldPolynomialOperations.addition(this, polynomial);
    }

    public FieldPolynomial subtract(@NonNull FieldPolynomial polynomial){
        return FieldPolynomialOperations.substraction(this, polynomial);
    }

    public FieldPolynomial multiply(int number){
        return multiply(new FieldPolynomial(this.field, new int[]{number}));
    }

    public FieldPolynomial multiply(@NonNull FieldPolynomial polynomial){
        return FieldPolynomialOperations.multiplication(this, polynomial);
    }

    public FieldPolynomial divide(@NonNull FieldPolynomial polynomial){
        return FieldPolynomialOperations.division(this, polynomial);
    }

    public FieldPolynomial mod(@NonNull FieldPolynomial polynomial){
        return FieldPolynomialOperations.mod(this, polynomial);
    }

    public FieldPolynomial derivative(){
        return FieldPolynomialOperations.derivative(this);
    }

    public FieldPolynomial interpolate(@NonNull int[] x){
        return FieldPolynomialOperations.interpolation(x, this);
    }
}
