pub mod ec;
pub(crate) mod poly;



use ark_ff::Field;



pub fn powers<F: Field>(base: F) -> impl Iterator<Item=F> {
    ark_std::iter::successors(Some(F::one()), move |power| Some(base * power))
}



#[cfg(test)]
pub fn curve_name<E: ark_ec::PairingEngine>() -> &'static str {
    // ark_ec::models::bw6::BW6<ark_bw6_761::curves::Parameters>
    let full_name = std::any::type_name::<E>();
    full_name.split_once("<").unwrap().1
        .split_once(":").unwrap().0
}