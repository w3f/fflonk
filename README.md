# FFlonk

FFlonk: a Fast-Fourier inspired [verifier efficient version of PlonK](https://eprint.iacr.org/2021/1167). Ariel Gabizon and Zachary J. Williamson suggest reducing opening multiple polynomials at a single point x, to opening a single polynomial at many points via an ``FFT-like identity''.

## Building 

```
cargo build --release
```

## Tests

```
cargo test --all
```
