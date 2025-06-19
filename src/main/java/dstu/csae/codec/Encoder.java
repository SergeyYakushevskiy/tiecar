package dstu.csae.codec;

public interface Encoder <T, V>{

    V encode(T input);

}
