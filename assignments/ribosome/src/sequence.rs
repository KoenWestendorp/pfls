use std::{borrow::Borrow, marker::PhantomData};

pub(crate) type Codon = [RNA; 3];

const fn fits_in<T>(bits: usize) -> usize {
    const WORDSIZE: usize = 8;
    (std::mem::size_of::<T>() * WORDSIZE) / bits
}

#[derive(Debug, Clone)]
pub(crate) struct Sequence {
    len: usize,
    pub chunks: Vec<usize>,
}

impl Sequence {
    pub fn new() -> Self {
        Self {
            len: 0,
            chunks: Vec::new(),
        }
    }

    pub fn from_vec(vec: Vec<DNA>) -> Self {
        let mut chunks = Vec::with_capacity(vec.len() / fits_in::<usize>(2));
        for (index, nt) in vec.iter().enumerate() {
            let nt = *nt as u8 as usize;
            let shift = index % fits_in::<usize>(2) * 2;
            if shift == 0 {
                chunks.push(0usize)
            }
            // This is a safe unwrap, because in the first iteration of the loop, shift is
            // guaranteed to be 0, since 0 % n == 0 for all n.
            let last = chunks.last_mut().unwrap();
            let shifted_nt = nt << shift;
            *last |= shifted_nt;
        }

        Self {
            len: vec.len(),
            chunks,
        }
    }

    #[inline]
    fn chunk(&self, nt_index: usize) -> usize {
        let chunk_index = nt_index / fits_in::<usize>(2);
        self.chunks[chunk_index]
    }

    pub fn to_dna_vec(&self) -> Vec<DNA> {
        let mut dna_vec = Vec::with_capacity(self.chunks.len() * fits_in::<usize>(2));
        let mut nt_index = 0;
        while nt_index < self.len {
            // We can optimize this.
            let chunk = self.chunk(nt_index);
            let nt_bits = chunk >> nt_index % fits_in::<usize>(2) * 2 & 3;
            // Safety: we can safely cast this u8 to DNA because we know it has been created by
            // casting from DNA to u8 before.
            let nt = unsafe { std::mem::transmute(nt_bits as u8) };
            dna_vec.push(nt);
            nt_index += 1;
        }

        dna_vec
    }

    pub fn len(&self) -> usize {
        self.len
    }

    // NOTE: As far as I know, there is not a nice way to impl Index<usize> for Sequence. This
    // would require we can return a reference to the inner data, but this is not really possible,
    // since we get our return value by bit fiddling.
    pub fn get(&self, index: usize) -> DNA {
        // If index is out of bounds, panic here, rather than in our internal functions.
        {
            let len = self.chunks.len() * fits_in::<usize>(2);
            if index >= len {
                panic!("index out of bounds: the len is {len} but the index is {index}")
            }
        }
        // Get the chunk by nt index.
        let chunk = self.chunk(index);
        let nt_bits = chunk >> index % fits_in::<usize>(2) * 2 & 3;
        // Safety: we can safely cast this u8 to DNA because we know it has been created by
        // casting from DNA to u8 before.
        let nt = unsafe { std::mem::transmute(nt_bits as u8) };
        nt
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub(crate) enum DNA {
    A = 0b00,
    T = 0b01,
    C = 0b10,
    G = 0b11,
}

impl From<char> for DNA {
    fn from(value: char) -> Self {
        match value {
            'A' => Self::A,
            'T' => Self::T,
            'C' => Self::C,
            'G' => Self::G,
            unknown => panic!("unknown nucleotide '{unknown}'"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub(crate) enum RNA {
    A,
    U,
    C,
    G,
}

impl From<char> for RNA {
    fn from(value: char) -> Self {
        match value {
            'A' => Self::A,
            'U' => Self::U,
            'C' => Self::C,
            'G' => Self::G,
            unknown => panic!("unknown nucleotide '{unknown}'"),
        }
    }
}

impl From<DNA> for RNA {
    fn from(value: DNA) -> Self {
        match value {
            DNA::A => Self::A,
            DNA::T => Self::U,
            DNA::C => Self::C,
            DNA::G => Self::G,
        }
    }
}

pub(crate) type Seq<T> = Vec<T>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum AminoAcid {
    /// Alanine 	    'A'
    Ala,
    /// Arginine 	    'R'
    Arg,
    /// Asparagine 	    'N'
    Asn,
    /// Aspartate 	    'D'
    Asp,
    /// Cysteine 	    'C'
    Cys,
    /// Glutamine 	    'Q'
    Gln,
    /// Glutamate 	    'E'
    Glu,
    /// Glycine 	    'G'
    Gly,
    /// Histidine 	    'H'
    His,
    /// Isoleucine 	    'I'
    Ile,
    /// Leucine 	    'L'
    Leu,
    /// Lysine 	        'K'
    Lys,
    /// Methionine 	    'M'
    Met,
    /// Phenylalanine 	'F'
    Phe,
    /// Proline 	    'P'
    Pro,
    /// Serine 	        'S'
    Ser,
    /// Threonine 	    'T'
    Thr,
    /// Tryptophan 	    'W'
    Trp,
    /// Tyrosine 	    'Y'
    Tyr,
    /// Valine 	        'V'
    Val,

    /// Stop codon
    Stop,
}

impl From<&str> for AminoAcid {
    fn from(value: &str) -> Self {
        use AminoAcid::*;
        match value {
            "A" => Ala,
            "R" => Arg,
            "N" => Asn,
            "D" => Asp,
            "C" => Cys,
            "Q" => Gln,
            "E" => Glu,
            "G" => Gly,
            "H" => His,
            "I" => Ile,
            "L" => Leu,
            "K" => Lys,
            "M" => Met,
            "F" => Phe,
            "P" => Pro,
            "S" => Ser,
            "T" => Thr,
            "W" => Trp,
            "Y" => Tyr,
            "V" => Val,
            "STOP" => Stop,
            unknown => panic!("unknown amino acid code '{unknown}'"),
        }
    }
}

pub(crate) fn dict(codon: &[RNA]) -> AminoAcid {
    use AminoAcid::*;
    use RNA::*;
    match codon {
        [U, U, U] => Phe,
        [U, U, C] => Phe,
        [U, U, A] => Leu,
        [U, U, G] => Leu,
        [C, U, U] => Leu,
        [C, U, C] => Leu,
        [C, U, A] => Leu,
        [C, U, G] => Leu,
        [A, U, U] => Ile,
        [A, U, C] => Ile,
        [A, U, A] => Ile,
        [A, U, G] => Met,
        [G, U, U] => Val,
        [G, U, C] => Val,
        [G, U, A] => Val,
        [G, U, G] => Val,
        [U, C, U] => Ser,
        [U, C, C] => Ser,
        [U, C, A] => Ser,
        [U, C, G] => Ser,
        [C, C, U] => Pro,
        [C, C, C] => Pro,
        [C, C, A] => Pro,
        [C, C, G] => Pro,
        [A, C, U] => Thr,
        [A, C, C] => Thr,
        [A, C, A] => Thr,
        [A, C, G] => Thr,
        [G, C, U] => Ala,
        [G, C, C] => Ala,
        [G, C, A] => Ala,
        [G, C, G] => Ala,
        [U, A, U] => Tyr,
        [U, A, C] => Tyr,
        [U, A, A] => Stop,
        [U, A, G] => Stop,
        [C, A, U] => His,
        [C, A, C] => His,
        [C, A, A] => Gln,
        [C, A, G] => Gln,
        [A, A, U] => Asn,
        [A, A, C] => Asn,
        [A, A, A] => Lys,
        [A, A, G] => Lys,
        [G, A, U] => Asp,
        [G, A, C] => Asp,
        [G, A, A] => Glu,
        [G, A, G] => Glu,
        [U, G, U] => Cys,
        [U, G, C] => Cys,
        [U, G, A] => Stop,
        [U, G, G] => Trp,
        [C, G, U] => Arg,
        [C, G, C] => Arg,
        [C, G, A] => Arg,
        [C, G, G] => Arg,
        [A, G, U] => Ser,
        [A, G, C] => Ser,
        [A, G, A] => Arg,
        [A, G, G] => Arg,
        [G, G, U] => Gly,
        [G, G, C] => Gly,
        [G, G, A] => Gly,
        [G, G, G] => Gly,
        _ => panic!("unknown codon"),
    }
}
