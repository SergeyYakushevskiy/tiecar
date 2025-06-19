package dstu.csae.codec;

public interface Decoder <T, V>{

    V decode(T input);

}
