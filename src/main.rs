use std::time::Instant;
use perf_comp::*;


const PI: [u32;32] = [1,2,4,6,11,18,31,54,
                97,172,309,564,1028,1900,3512,6542,
        12251,23000,43390,82025,155611,295947,564163,1077871,
2063689,3957809,7603553,14630843,28192750,54400028,105097565,203280221];


// test aux sieve
pub fn main() {

    let start = Instant::now();

    let pattern = init_pattern();
    let mut counter1 = 0;
    let mut counter2 = 0;
    let mut x = 1;
    let mut c = 1;
    let mut l0 = 0;
    let mut l1 = 1;
    eprintln!("testing aux_sieve: ");
    let mut aux_base: u32 = 0;
    let (mut aux_sieve, aux_primes) = init_aux_primes();
    while l0 < 32 {
        counter1 += 1;
        update_aux_sieve(aux_base, &mut aux_sieve, &aux_primes, &pattern);
        let mut i = 0;
        while i < _AUX_SIEVE_WORDS_ {
            let mut j = (l1 - x) / (16 * _POINTER_SIZE_) ;
            if j > 0 {
                if j > _AUX_SIEVE_WORDS_ - i {
                    j = _AUX_SIEVE_WORDS_ - i;
                }
                counter2 += 1;
                c += count_zero_bits(&aux_sieve, i, j);
                i += j;
                x += j * (16 * _POINTER_SIZE_);
            } else {
                for j in 0..(8*_POINTER_SIZE_) {
                    if aux_sieve[i as usize] & MARK_MASK[j as usize] == 0 { c+=1; }
                    if x == l1 {
                        assert_eq!(c, PI[l0]);
                            //,"bad pi(2^{})={} != {}\n",l0 + 1,pi[l0],c);
                        l0 += 1;
                        println!("l0: {}",l0);
                        l1 = 2 * l1 + 1;
                    }
                    x += 2;
                    assert!(x <= l1);
                }
                i += 1;
            }

        }
        aux_base += _AUX_SIEVE_SPAN_;
    }

    let end = Instant::now();
    let duration = end.duration_since(start);
    eprintln!("good, {:?} {} {} \n", duration, counter1, counter2);

}
