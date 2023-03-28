use std::collections::HashMap;
use std::fs::read_to_string;
use std::io;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum DNA {
    A,
    T,
    C,
    G,
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
enum RNA {
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

type Seq<T> = Vec<T>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum AminoAcid {
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

type Codon = [RNA; 3];
type CodonToAminoAcidDict = HashMap<Codon, AminoAcid>;

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

    let mut dna = Seq::new();
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
fn translate(
    amino_acid_dict: &CodonToAminoAcidDict,
    rna: &Seq<RNA>,
    start: usize,
) -> Vec<AminoAcid> {
    assert!(start <= rna.len());

    let mut start = start;
    let mut prot = Vec::new();
    while start + 3 <= rna.len() {
        if let Some(aa) = amino_acid_dict.get(&rna[start..start + 3]) {
            if aa == &AminoAcid::Stop {
                break;
            }

            prot.push(*aa)
        } else {
            return prot;
        }

        start += 3
    }

    prot
}

fn translate_all(
    amino_acid_dict: &CodonToAminoAcidDict,
    start_codon: [RNA; 3],
    rna: &Seq<RNA>,
) -> Vec<Vec<AminoAcid>> {
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
        proteins.push(translate(amino_acid_dict, rna, start));
    }

    proteins
}

fn main() -> io::Result<()> {
    const AUG: Codon = [RNA::A, RNA::U, RNA::G];

    let dict = read_dictionary("Genetic_Code.txt")?;
    let dna = read_dna("E_Coli.txt", None)?;
    let rna = transcribe_dna(dna);

    let first_start = find_start_position(AUG, &rna, 0).expect("should contain 'AUG'");
    let first_prot = translate(&dict, &rna, first_start);
    println!("First protein: {first_prot:?}");

    let all_prots = translate_all(&dict, AUG, &rna);
    println!("All protein:");
    for (n, prot) in all_prots.iter().enumerate() {
        //println!("({n}) => {prot}", n = n + 1);
    }

    Ok(())
}
