package dstu.csae.rs_code;

import dstu.csae.code.Code;
import dstu.csae.galois.extended.GaloisFieldExtension;
import dstu.csae.polynomial.FieldPolynomial;
import dstu.csae.polynomial.Polynomial;
import lombok.Getter;
import lombok.NonNull;

@Getter
public class RSCode implements Code<FieldPolynomial> {
    private final GaloisFieldExtension field;
    private final int n;
    private final int k;
    private final int d;
    private final int r;
    private final int t;
    private final int alpha;
    private final Polynomial g;
    private final Polynomial h;

    public RSCode(@NonNull GaloisFieldExtension field, int n, int k)
            throws IllegalArgumentException{
        if(n >= field.getCharacteristic()){
            throw new IllegalArgumentException("Значение n должно быть не больше q - 1");
        }
        if(k >= n){
            throw new IllegalArgumentException("Значение k должно быть меньше n: k < n");
        }
        this.field = field;
        this.n = n;
        this.k = k;
        this.r = n - k;
        this.d = r + 1;
        this.t = (d - 1) / 2;
        alpha = field.findFirstIntPrimitive();
        g = setG(field, r, alpha);
        h = setH(field, r, alpha);
    }

    public RSCode(@NonNull GaloisFieldExtension field, int d){
        this.field = field;
        this.d = d;
        n = field.getCharacteristic() - 1;
        r = d - 1;
        k = n - r;
        t = (d - 1) / 2;
        alpha = field.findFirstIntPrimitive();
        g = setG(field, r, alpha);
        h = setH(field, r, alpha);
    }

    private Polynomial setG(GaloisFieldExtension field, int r, Integer alpha){
        FieldPolynomial g = new FieldPolynomial(field, r + 1);
        g.set(0, 1);
        for(int i = 1; i <= r ; i++){
            FieldPolynomial multiplier = new FieldPolynomial(field, 1);
            multiplier.set(1, 1);
            int currentAlpha = field.powMod(alpha, i);
            multiplier.set(0, field.inverseOfAddition(currentAlpha));
            g = g.multiply(multiplier);
        }
        return g;
    }

    private Polynomial setH(GaloisFieldExtension field, int r, Integer alpha){
        FieldPolynomial h = new FieldPolynomial(field, n - r);
        h.set(0, 1);
        for(int i = 0; i <= n - 1; i ++){
            if(i > 0 && i <= r){
                continue;
            }
            FieldPolynomial multiplier = new FieldPolynomial(field, 1);
            multiplier.set(1, 1);
            int currentAlpha = field.powMod(alpha, i);
            multiplier.set(0, field.inverseOfAddition(currentAlpha));
            h = h.multiply(multiplier);
        }
        return h;
    }

    @Override
    public String toString() {
        return "RSCode(" + n + "," + k + "), " + field;
    }
}