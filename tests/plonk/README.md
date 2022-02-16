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
Start:   domain_size = 65536,  curve = ark_bls12_381
··Start:   Setup
····Start:   Computing 786420 scalars powers
····End:     Computing 786420 scalars powers .......................................41.347ms
····Start:   786420-scalar mul in G1
····End:     786420-scalar mul in G1 ...............................................4.923s
····Start:   2-scalar mul in G1
····End:     2-scalar mul in G1 ....................................................5.045ms
··End:     Setup ...................................................................5.008s
··Start:   Preprocessing
····Start:   Committing to combination #0
······Start:   combining 8 polynomials: t = 8, max_degree = 65535
······End:     combining 8 polynomials: t = 8, max_degree = 65535 ..................10.839ms
······Start:   committing to the combined polynomial: degree = 524287
······End:     committing to the combined polynomial: degree = 524287 ..............2.984s
····End:     Committing to combination #0 ..........................................2.996s
··End:     Preprocessing ...........................................................2.998s
··Start:   Proving
····Start:   Committing to 2 proof polynomials
······Start:   Committing to combination #1
········Start:   combining 4 polynomials: t = 4, max_degree = 131069
········End:     combining 4 polynomials: t = 4, max_degree = 131069 ...............5.757ms
········Start:   committing to the combined polynomial: degree = 524279
········End:     committing to the combined polynomial: degree = 524279 ............1.914s
······End:     Committing to combination #1 ........................................1.921s
······Start:   Committing to combination #2
········Start:   combining 4 polynomials: t = 4, max_degree = 196604
········End:     combining 4 polynomials: t = 4, max_degree = 196604 ...............6.881ms
········Start:   committing to the combined polynomial: degree = 786418
········End:     committing to the combined polynomial: degree = 786418 ............1.866s
······End:     Committing to combination #2 ........................................1.873s
····End:     Committing to 2 proof polynomials .....................................3.799s
····Start:   Opening
······Start:   polynomial divisions
······End:     polynomial divisions ................................................445.258ms
······Start:   commitment to a degree-786410 polynomial
······End:     commitment to a degree-786410 polynomial ............................3.334s
······Start:   linear combination of polynomials
······End:     linear combination of polynomials ...................................93.415ms
····End:     Opening ...............................................................8.324s
··End:     Proving .................................................................12.148s
··Start:   Verifying
····Start:   barycentric evaluations
····End:     barycentric evaluations ...............................................135.500µs
····Start:   multiexp
····End:     multiexp ..............................................................505.100µs
··End:     Verifying ...............................................................3.324ms
End:     domain_size = 65536,  curve = ark_bls12_381 ...............................20.161s
```