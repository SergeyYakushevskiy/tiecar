package dstu.csae.utils;

import dstu.csae.galois.extended.GaloisFieldExtension;
import dstu.csae.polynomial.FieldPolynomial;
import lombok.NonNull;

public class FieldPolynomialGenerator {

    private final GaloisFieldExtension field;
    private final int fieldPower;
    private final int size;
    private long counter;
    private long count;

    public FieldPolynomialGenerator(@NonNull GaloisFieldExtension field, int size){
        assert size > 1;
        this.field = field;
        this.fieldPower = field.getCharacteristic();
        this.size = size;
        this.count = (long) Math.pow(field.getCharacteristic(), size);
        this.counter = 0;
    }

    public boolean hasNext(){
        return counter < count;
    }

    public FieldPolynomial next(){
        if(!hasNext()){
            return null;
        }
        return new FieldPolynomial(field, getState(counter ++));
    }

    private int[] getState(long current){
        int[] state = new int[size];
        for(int i = size - 1; i >= 0; i --){
            state[i] = (int)(current % fieldPower);
            current /= fieldPower;
        }
        return state;
    }

}
