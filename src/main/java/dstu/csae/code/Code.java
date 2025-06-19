package dstu.csae.code;

import dstu.csae.galois.extended.GaloisFieldExtension;
import dstu.csae.polynomial.Polynomial;

public interface Code<T>{
    GaloisFieldExtension getField();
    int getN();
    int getK();
    int getD();
    int getT();
    int getR();
    int getAlpha();
    Polynomial getG();
    Polynomial getH();


}
