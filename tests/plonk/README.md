> cargo test test_vanilla_plonk_batch_kzg_opening --release --features "parallel print-trace" -- --nocapture

```
Start:   domain_size = 65536,  curve = ark_bls12_381
··Start:   Setup
····Start:   Computing 196606 scalars powers
····End:     Computing 196606 scalars powers .......................................14.267ms
····Start:   196606-scalar mul in G1
····End:     196606-scalar mul in G1 ...............................................1.278s
····Start:   2-scalar mul in G1
····End:     2-scalar mul in G1 ....................................................5.610ms
··End:     Setup ...................................................................1.308s
··Start:   Preprocessing
····Start:   Committing to batch of 8 polynomials
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................461.575ms
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................429.677ms
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................422.081ms
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................429.792ms
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................418.964ms
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................444.253ms
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................445.711ms
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................415.064ms
····End:     Committing to batch of 8 polynomials ..................................3.468s
··End:     Preprocessing ...........................................................3.469s
··Start:   Proving
····Start:   Committing to batch of 3 polynomials
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................421.111ms
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................415.927ms
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................419.734ms
····End:     Committing to batch of 3 polynomials ..................................1.257s
····Start:   Committing to degree 65535 polynomials
····End:     Committing to degree 65535 polynomials ................................388.513ms
····Start:   Committing to degree 196604 polynomials
····End:     Committing to degree 196604 polynomials ...............................1.169s
····Start:   Extra: commiting to the linearization polynomial
······Start:   Committing to degree 65535 polynomials
······End:     Committing to degree 65535 polynomials ..............................496.133ms
····End:     Extra: commiting to the linearization polynomial ......................497.160ms
··End:     Proving .................................................................4.221s
··Start:   Verifying
····Start:   Reconstructing the commitment to the linearization polynomial: 7-multiexp
····End:     Reconstructing the commitment to the linearization polynomial: 7-multiexp 1.277ms
····Start:   KZG batch verification
······Start:   aggregate evaluation claims at zeta
······End:     aggregate evaluation claims at zeta .................................509.800µs
······Start:   batched KZG openning
······End:     batched KZG openning ................................................3.967ms
····End:     KZG batch verification ................................................5.080ms
··End:     Verifying ...............................................................7.252ms
End:     domain_size = 65536,  curve = ark_bls12_381 ...............................9.010s
```

> cargo test test_vanilla_plonk_with_fflonk_opening --release --features "parallel print-trace" -- --nocapture

```
Start:   domain_size = 256,  curve = ark_bls12_381
··Start:   Setup
····Start:   Computing 3060 scalars powers
····End:     Computing 3060 scalars powers .........................................125.100µs
····Start:   3060-scalar mul in G1
····End:     3060-scalar mul in G1 .................................................20.446ms
····Start:   2-scalar mul in G1
····End:     2-scalar mul in G1 ....................................................3.521ms
··End:     Setup ...................................................................27.504ms
··Start:   Preprocessing
····Start:   Committing to combination #0
······Start:   combining 8 polynomials: t = 8, max_degree = 255
······End:     combining 8 polynomials: t = 8, max_degree = 255 ....................23.400µs
······Start:   committing to the combined polynomial: degree = 2047
······End:     committing to the combined polynomial: degree = 2047 ................14.400ms
····End:     Committing to combination #0 ..........................................15.314ms
··End:     Preprocessing ...........................................................15.729ms
··Start:   Proving
····Start:   Committing to 2 proof polynomials
······Start:   Committing to combination #1
········Start:   combining 4 polynomials: t = 4, max_degree = 509
········End:     combining 4 polynomials: t = 4, max_degree = 509 ..................8.600µs
········Start:   committing to the combined polynomial: degree = 2039
········End:     committing to the combined polynomial: degree = 2039 ..............10.475ms
······End:     Committing to combination #1 ........................................10.859ms
······Start:   Committing to combination #2
········Start:   combining 4 polynomials: t = 4, max_degree = 764
········End:     combining 4 polynomials: t = 4, max_degree = 764 ..................17.800µs
········Start:   committing to the combined polynomial: degree = 3058
········End:     committing to the combined polynomial: degree = 3058 ..............12.998ms
······End:     Committing to combination #2 ........................................13.867ms
····End:     Committing to 2 proof polynomials .....................................25.572ms
····Start:   Opening
······Start:   polynomial divisions
······End:     polynomial divisions ................................................1.747ms
······Start:   commitment to a degree-3050 polynomial
······End:     commitment to a degree-3050 polynomial ..............................17.542ms
······Start:   linear combination of polynomials
······End:     linear combination of polynomials ...................................437.600µs
····End:     Opening ...............................................................46.147ms
··End:     Proving .................................................................72.668ms
··Start:   Verifying
····Start:   barycentric evaluations
····End:     barycentric evaluations ...............................................95.800µs
····Start:   multiexp
····End:     multiexp ..............................................................435.900µs
··End:     Verifying ...............................................................2.399ms
End:     domain_size = 256,  curve = ark_bls12_381 .................................119.557ms
```