//use crossbeam_utils::CachePadded;
use std::time::Instant;
use std::cmp;

//const _MARK_MASK_: [u64; 64] = array_init::array_init(|i| 1u64 << i); // come on Const Generics
pub const MARK_MASK: [u64;64] = {
    let mut res = [0; 64];
    let mut i = 0;
    while i < 64 {
        res[i] = 1 << i;
        i += 1;
    }
    res
};

// base-2 logarithm of the size (in bytes) of the bucket structure
const _BUCKET_SIZE_LOG2_ : usize = 12; 
const _SIEVE_BITS_LOG2_ : u32 = 23;
const _PRIME_THRESHOLD_: u64 = 1 << _SIEVE_BITS_LOG2_;
const _PRIMES_PER_BUCKET_ : usize = (1 << (_BUCKET_SIZE_LOG2_ - 3)) - 1;
const _BUCKET_ADDR_MASK_  : usize = (1 << _BUCKET_SIZE_LOG2_) - 1;

pub const _NUMBER_OF_AUX_PRIMES_: usize = 6536;
pub const _AUX_SIEVE_SPAN_: u32 = 1 << 19;
pub const _SIEVE_SPAN_: u64 = 2 << _SIEVE_BITS_LOG2_;
pub const _SIEVE_WORDS_: u32 = 1 << (_SIEVE_BITS_LOG2_ - 6);
pub const _POINTER_SIZE_: u32 = 8; // please just use 64Bit / 8 Byte. 
pub const _AUX_SIEVE_WORDS_: u32 = _AUX_SIEVE_SPAN_ / (2 * 8 * _POINTER_SIZE_);

#[inline(always)]
fn mark_2(arr: &mut [u64], index: usize) { arr[index >> 6] |= MARK_MASK[index & 63usize]; }
#[inline]
const fn test_2(arr: &[u64], index: usize) -> bool { arr[index >> 6] & MARK_MASK[index & 63usize] == 0 }


#[inline(always)]
pub fn count_zero_bits(arr: &[u64], start_at: u32, count_entries: u32) -> u32 {
    let a = start_at as usize;
    let b = a + count_entries as usize;
    let slice = &arr[a..b];
    slice.iter().copied().map(u64::count_zeros).sum()
}

//static pattern: mut[u64; 3 * 5 * 7 * 11 * 13] = [0; 3 * 5 * 7 * 11 * 13];

pub fn init_pattern() -> [u64; 3*5*7*11*13] {
    println!("initializing pattern");
    const LIMIT: usize = 3 * 5 * 7 * 11 * 13;
    let limit2 = LIMIT * 8 * _POINTER_SIZE_ as usize;
    let mut res: [u64; LIMIT] = [0; LIMIT];
    for i in ((3>>1)..limit2).step_by(3) { mark_2(&mut res, i); }
    for i in ((5>>1)..limit2).step_by(5) { mark_2(&mut res, i); }
    for i in ((7>>1)..limit2).step_by(7) { mark_2(&mut res, i); }
    for i in ((11>>1)..limit2).step_by(11) { mark_2(&mut res, i); }
    for i in ((13>>1)..limit2).step_by(13) { mark_2(&mut res, i); }
    res
}


#[inline]
pub fn update_aux_sieve(aux_base: &mut u32, aux_sieve: &mut[u64]
                        , aux_primes: &[u32], pattern: &[u64]) {
    
    assert_eq!((*aux_base & (_AUX_SIEVE_SPAN_ -1)), 0, "update_aux_sieve: illegal aux_base!");
    
    let first_primes = 3 * 5 * 7 * 11  * 13;

    // now initialize the aux_sieve array

    let o = *aux_base % first_primes;
    let offset = (o + ((o * 105) & 127) * first_primes) >> 7;// 105 = -1/15015 mod 128

    let mut i = cmp::min(_AUX_SIEVE_WORDS_, first_primes - offset);
    for j in 0..i {
        aux_sieve[j as usize] = pattern[(offset + j) as usize];
    }


    while i < _AUX_SIEVE_WORDS_ {
        let k = cmp::min(_AUX_SIEVE_WORDS_ - i, first_primes);
        for j in 0..k {
            aux_sieve[(i + j) as usize] = pattern[j as usize ];
        }
        i += k;
    }
    
    if *aux_base == 0
    {
        println!("aux_base is 0");
        // mark 1 as not prime, and mark 3, 5, 7, 11, and 13 as prime
        aux_sieve[0] |= MARK_MASK[0]; 
        aux_sieve[0] &= !(MARK_MASK[1] | MARK_MASK[2] 
                          | MARK_MASK[3] | MARK_MASK[5] | MARK_MASK[6]);
    }

    for i in 0.._NUMBER_OF_AUX_PRIMES_ {
        let mut j = aux_primes[i] * aux_primes[i];
        if j > *aux_base + (_AUX_SIEVE_SPAN_ - 1) { break; }
        if j > *aux_base { 
            j = ( j - *aux_base ) >> 1; 
        } else {
            j = aux_primes[i] - *aux_base % aux_primes[i];
            if (j & 1) == 0 {
                j += aux_primes[i];
            }
            j >>= 1;
        }
        while j < (_AUX_SIEVE_SPAN_ / 2) {
            mark_2(aux_sieve, j as usize);
            j += aux_primes[i];
        }
    }
}



pub fn init_aux_primes() -> ([u64;_AUX_SIEVE_WORDS_ as usize], [u32;_NUMBER_OF_AUX_PRIMES_]) {
    //aux_sieve: &mut [u64], aux_primes: &mut [u32]) {
    let mut aux_sieve:  [u64;_AUX_SIEVE_WORDS_ as usize] = [0;_AUX_SIEVE_WORDS_ as usize];
    let mut aux_primes: [u32;_NUMBER_OF_AUX_PRIMES_] = [0;_NUMBER_OF_AUX_PRIMES_];
    let start = Instant::now();
    println!("Initializing aux primes");
    for i in (3..256).step_by(2) {
        if test_2(&mut aux_sieve, i >> 1usize) {
            for j in (((i*i) >> 1)..32768).step_by(i) {
                mark_2(&mut aux_sieve, j as usize);
            }
        }
    }
    let mut j = 0;
    for i in (17>>1)..32768 {  // start at 17
        if test_2(&mut aux_sieve, i) {
            aux_primes[j] = 2 * (i as u32) + 1;
            j += 1;
        }
    }
    assert!( j == _NUMBER_OF_AUX_PRIMES_ && aux_primes[_NUMBER_OF_AUX_PRIMES_ - 1] == 65521
        ,"update_aux_sieve: initialization of the aux_primes array failed\n {} {} {} {}"
        , j, _NUMBER_OF_AUX_PRIMES_, aux_primes[_NUMBER_OF_AUX_PRIMES_ -1], 65521);



    let end = Instant::now();
    let duration = end.duration_since(start);
    eprintln!("init_aux_primes took, {:?}\n", duration);
    (aux_sieve, aux_primes)
}
