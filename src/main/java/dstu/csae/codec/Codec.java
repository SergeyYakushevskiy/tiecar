package dstu.csae.codec;

import dstu.csae.code.Code;
import lombok.Getter;
import lombok.Setter;

public abstract class Codec <C extends Code, T, V> implements Encoder<T, V>, Decoder<V, T>{

    @Getter
    @Setter
    protected C code;

    public Codec(C code){
        this.code = code;
    }

}
