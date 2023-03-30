use std::collections::HashMap;
use std::fs::read_to_string;
use std::io::{self, stdout, BufWriter, Write};

use sequence::{dict, AminoAcid, Codon, Seq, Sequence, DNA, RNA};

mod sequence;

type CodonToAminoAcidDict = HashMap<Codon, AminoAcid>;

#[allow(dead_code)]
fn read_dictionary(path: &str) -> io::Result<CodonToAminoAcidDict> {
    let contents = read_to_string(path)?;

    let mut dict = CodonToAminoAcidDict::new();
    for line in contents.lines() {
        let mut split = line.split_ascii_whitespace();
        let mut codon = split.next().unwrap().chars();
        let codon = [codon.next(), codon.next(), codon.next()].map(|nt| nt.unwrap().into());
        let amino_acid = split.next().unwrap().into();
        dict.insert(codon, amino_acid);
    }

    Ok(dict)
}

fn read_dna(path: &str, lines: Option<usize>) -> io::Result<Seq<DNA>> {
    let contents = read_to_string(path)?;
    let dna_lines = contents.lines().skip(1); // skip(1) to discard the faste header
    let dna_string: String = if let Some(lines) = lines {
        dna_lines.take(lines).collect()
    } else {
        dna_lines.collect()
    };

    let mut dna = Seq::with_capacity(dna_string.len());
    for nt in dna_string.chars() {
        dna.push(DNA::from(nt))
    }

    Ok(dna)
}

/// Transcribes DNA to RNA.
///
/// DNA must be valid. That is, must only contain 'A', 'C', 'T', or 'G' characters.
///
/// # Panics
///
/// If DNA is invalid---i.e., an invalid nucleotide is encountered---this function panics.
fn transcribe_dna(dna: Seq<DNA>) -> Seq<RNA> {
    let coding_strand = dna;

    let mut rna = Seq::new();
    for nt in coding_strand {
        rna.push(nt.into());
    }

    rna
}

fn find_start_position<T>(start_codon: [T; 3], seq: &Seq<T>, start: usize) -> Option<usize>
where
    T: PartialEq,
{
    let mut n = start;
    while n + 3 <= seq.len() {
        if seq[n..n + 3] == start_codon {
            return Some(n);
        }
        n += 1;
    }

    None
}

/// Translates RNA to an amino acid sequence (protein).
///
/// # Panics
///
/// Length of `rna` must be greater than or equal to the value of `start`.
fn translate(rna: &Seq<RNA>, start: usize) -> Vec<AminoAcid> {
    assert!(start <= rna.len());

    let mut start = start;
    let mut prot = Vec::new();
    while start + 3 <= rna.len() {
        // if let Some(aa) = amino_acid_dict.get(&rna[start..start + 3]) {
        let aa = dict(&rna[start..start + 3]);
        if aa == AminoAcid::Stop {
            break;
        }
        prot.push(aa);
        start += 3
    }

    prot
}

fn translate_all(start_codon: [RNA; 3], rna: &Seq<RNA>) -> Vec<Vec<AminoAcid>> {
    // Find all start positions in the RNA sequence.
    let mut starts = Vec::new();
    let mut pos = 0;
    while pos <= rna.len() {
        if let Some(start) = find_start_position(start_codon, &rna, pos + 1) {
            starts.push(start);
            pos = start
        } else {
            break;
        }
    }

    // Translate all proteins following the start positions.
    let mut proteins = Vec::new();
    for start in starts {
        proteins.push(translate(rna, start));
    }

    proteins
}

fn main() -> io::Result<()> {
    // const AUG: Codon = [RNA::A, RNA::U, RNA::G];
    //
    // let dna = read_dna("E_Coli.txt", None)?;
    // let rna = transcribe_dna(dna);
    //
    // let first_start = find_start_position(AUG, &rna, 0).expect("should contain 'AUG'");
    // let first_prot = translate(&rna, first_start);
    // println!("First protein: {first_prot:?}");
    //
    // let all_prots = translate_all(AUG, &rna);
    // println!("All protein:");
    // let mut stdout = BufWriter::new(stdout());
    // for (n, prot) in all_prots.iter().enumerate() {
    //     writeln!(stdout, "({n}) => {prot:?}", n = n + 1)?;
    // }

    let dna = read_dna("E_Coli.txt", None)?;
    dbg!(dna.len());
    let seq = Sequence::from_vec(dna.clone());
    dbg!(seq.chunks.len());

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use DNA::*;
    static TEST_DNA: [sequence::DNA; 700] = [
        A, G, C, T, T, T, T, C, A, T, T, C, T, G, A, C, T, G, C, A, A, C, G, G, G, C, A, A, T, A,
        T, G, T, C, T, C, T, G, T, G, T, G, G, A, T, T, A, A, A, A, A, A, A, G, A, G, T, C, T, C,
        T, G, A, C, A, G, C, A, G, C, T, T, C, T, G, A, A, C, T, G, G, T, T, A, C, C, T, G, C, C,
        G, T, G, A, G, T, A, A, A, T, T, A, A, A, A, T, T, T, T, A, T, T, G, A, C, T, T, A, G, G,
        T, C, A, C, T, A, A, A, T, A, C, T, T, T, A, A, C, C, A, A, T, A, T, A, G, G, C, A, T, A,
        G, C, G, C, A, C, A, G, A, C, A, G, A, T, A, A, A, A, A, T, T, A, C, A, G, A, G, T, A, C,
        A, C, A, A, C, A, T, C, C, A, T, G, A, A, A, C, G, C, A, T, T, A, G, C, A, C, C, A, C, C,
        A, T, T, A, C, C, A, C, C, A, C, C, A, T, C, A, C, C, A, C, C, A, C, C, A, T, C, A, C, C,
        A, T, T, A, C, C, A, T, T, A, C, C, A, C, A, G, G, T, A, A, C, G, G, T, G, C, G, G, G, C,
        T, G, A, C, G, C, G, T, A, C, A, G, G, A, A, A, C, A, C, A, G, A, A, A, A, A, A, G, C, C,
        C, G, C, A, C, C, T, G, A, C, A, G, T, G, C, G, G, G, C, T, T, T, T, T, T, T, T, C, G, A,
        C, C, A, A, A, G, G, T, A, A, C, G, A, G, G, T, A, A, C, A, A, C, C, A, T, G, C, G, A, G,
        T, G, T, T, G, A, A, G, T, T, C, G, G, C, G, G, T, A, C, A, T, C, A, G, T, G, G, C, A, A,
        A, T, G, C, A, G, A, A, C, G, T, T, T, T, C, T, G, C, G, G, G, T, T, G, C, C, G, A, T, A,
        T, T, C, T, G, G, A, A, A, G, C, A, A, T, G, C, C, A, G, G, C, A, G, G, G, G, C, A, G, G,
        T, G, G, C, C, A, C, C, G, T, C, C, T, C, T, C, T, G, C, C, C, C, C, G, C, C, A, A, A, A,
        T, C, A, C, C, A, A, C, C, A, C, C, T, G, G, T, G, G, C, G, A, T, G, A, T, T, G, A, A, A,
        A, A, A, C, C, A, T, T, A, G, C, G, G, C, C, A, G, G, A, T, G, C, T, T, T, A, C, C, C, A,
        A, T, A, T, C, A, G, C, G, A, T, G, C, C, G, A, A, C, G, T, A, T, T, T, T, T, G, C, C, G,
        A, A, C, T, T, C, T, G, A, C, G, G, G, A, C, T, C, G, C, C, G, C, C, G, C, C, C, A, G, C,
        C, G, G, G, A, T, T, C, C, C, G, C, T, G, G, C, G, C, A, A, T, T, G, A, A, A, A, C, T, T,
        T, C, G, T, C, G, A, C, C, A, G, G, A, A, T, T, T, G, C, C, C, A, A, A, T, A, A, A, A, C,
        A, T, G, T, C, C, T, G, C, A, T, G, G, C, A, T, T, A, G, T, T, T, G, T, T, A, G, G, G, C,
        A, G, T, G, C, C, C, G, G, A,
    ];

    #[test]
    fn sequence_from_vec() {
        let dna_vec = TEST_DNA.to_vec();
        let seq = Sequence::from_vec(dna_vec.clone());

        assert_eq!(dna_vec, seq.to_dna_vec());
    }

    #[test]
    fn get_from_sequence() {
        let dna_vec = TEST_DNA.to_vec();
        let seq = Sequence::from_vec(dna_vec.clone());

        let mut seq_get_vec = Vec::new();
        for i in 0..seq.len() {
            seq_get_vec.push(seq.get(i))
        }

        assert_eq!(dna_vec, seq_get_vec);
    }
}
