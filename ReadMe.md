A Rust CLI that prints out a newton fractal for a given polynomial.

Polynomial can be define through either coefficients: 


```bash
cargo run --release --coeffs -1,0,0,0,0,0,1,0
```
    creates the wada.png files, containing the newton fractal for

    P(X) = (-1 + 0i) + (0 + 0i)*X + (0 + 0i)*X² + (1 + 0i)X³

    eg.

    P(X) = X³ - 1


or through roots: 

```bash
cargo run --release --roots 1,0,-0.5,-0.87,-0.5,0.87
```
    creates roughly the same fractal as those are approximation of 1, -1/2 - i*sqrt(3)/2 and -1/2 + i*sqrt(3)/2, the 3rd roots of 1.
