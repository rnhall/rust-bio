//Declares a structure with a variable length encoding scheme. Memory is allocated
//for each 4 nucleotides as a single u8 in a vector of u8 values.

use std::fmt;
use std::ops::BitXor;
use std::ops::Not;

//Should I include mutable kmers and immutable kmers?

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Kmer {
    pub k: usize,
    pub sequence: Vec<u8>,
}

pub struct KmerIter {
    pub kmer: Kmer,
    pub position: usize,
    pub nucleotide: u8,
}

pub struct Kmerizer<'a> {
    pub k: usize,
    pub position: usize,
    pub sequence: &'a[u8],
    pub current_kmer: Kmer,
}

impl Kmer {
    pub fn new(len: usize, byte_seq: &[u8]) -> Self {
        let mut kmer = Kmer
        {   k: len, 
            sequence: Vec::new()
        };
        kmer.encode(byte_seq);
        kmer
    }

    pub fn from_literal(str_literal: &str) -> Self {
        Kmer::new(str_literal.len(), str_literal.as_bytes())
    }

    pub fn encode(&mut self, byte_seq: &[u8]) {

        for chunk in byte_seq.chunks(4) {
            let mut bit_seq: u8 = 0;
            for (i, nucleotide) in chunk.iter().enumerate() {
                match nucleotide {
                    //Apparently pow() wants a u32... though none of the values there should
                    //ever be larger than u8...
                    b'T' => {
                        bit_seq += 2u8.pow(((i*2)+1) as u32) + 2u8.pow((i*2) as u32);
                    }
                    b'A' => {
                        bit_seq += 0;
                    }
                    b'G' => {
                        bit_seq += 2u8.pow((i*2) as u32);
                    }
                    b'C' => {
                        bit_seq += 2u8.pow(((i*2)+1) as u32);
                    }
                    _ => {
                        panic!("Non-valid nucleotide detected!");
                    }
                }
            }
            self.sequence.push(bit_seq)
        }
    }

    pub fn decode(&self) -> String {
        let mut counter = 0;
        let mut byte_seq = String::new();
        for mer in self.sequence.iter(){
            let mut div = *mer;
            for _j in 0..4 {
                let rem = (div % 4) as u8;
                div = div / 4;
                match rem {
                    3 => {
                        byte_seq.push('T');
                    }
                    0 => {
                        byte_seq.push('A');
                    }
                    1 => {
                        byte_seq.push('G');
                    }
                    2 => {
                        byte_seq.push('C');
                    }
                    _ => {
                        panic!("Non-valid nucleotide detected!");
                    }
                }
                counter += 1;
                if counter >= self.k {
                    break
                }
            }
        }
        byte_seq
    }

    //Does not consume the Kmer and returns a new Kmer
    pub fn make_complement(&self) -> Kmer {
        let complement = Kmer::new(self.k, self.decode().as_bytes());
        !complement
    }

    //Modifies the existing Kmer
    pub fn complement(&mut self) {
        self.sequence = self.sequence.iter().map(|x| !x).collect();
    }

    pub fn make_reverse_complement(&self) -> Kmer {
        let reverse: String = self.decode().chars().rev().collect();
        let kmer = Kmer::from_literal(reverse.as_str());
        !kmer
    }

    pub fn reverse_complement(&mut self) {
        let reverse: String = self.decode().chars().rev().collect();
        self.encode(reverse.as_bytes());
        self.complement();
    }

    pub fn index(&self, position: usize) -> u8 {
        if self.k < position {
            panic!("Index is greater than kmer length!");
        }
        let bit_mask: u8 = 0b00000011;
        let shift = 2 * (position % 4);
        (self.sequence[position / 4] & (bit_mask << shift)) >> (shift)
    }
}

impl fmt::Display for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        //Create our output string
        let mut sequence = String::new();
        for mer in self.sequence.iter() {
            let mut div = *mer;
            for _i in 0..4 {
                let rem = (div % 4) as u8;
                div = div / 4;
                match rem {
                    3 => {
                        sequence.push('T');
                    }
                    0 => {
                        sequence.push('A');
                    }
                    1 => {
                        sequence.push('G');
                    }
                    2 => {
                        sequence.push('C');
                    }
                    _ => {
                        panic!("Oh no! Detected a non valid nucleotide!");
                    }
                }
            }
        }
        sequence.truncate(4 * (self.k / 4) + (self.k % 4));
        write!(f, "Vmer[{}]: {}", self.k, sequence)
    }
}

//Iterator implementations

impl Iterator for KmerIter {
    type Item = u8;
    fn next(&mut self) -> Option<u8> {
        if self.position == self.kmer.k {
            return None
        }
        self.nucleotide = self.kmer.index(self.position);
        self.position += 1;
        Some(self.nucleotide)
    }
}

impl IntoIterator for Kmer {
    type Item = u8;
    type IntoIter = KmerIter;
    fn into_iter(self) -> Self::IntoIter {
        KmerIter {
            kmer: self,
            position: 0,
            nucleotide: 0,
        }
    }
}

//BITWISE IMPLEMENTATIONS

impl BitXor for Kmer {
    type Output = Self;
    fn bitxor(self, rhs: Self) -> Self::Output {
        assert_eq!(self.k, rhs.k);
        let mut xor_sequence: Vec<u8> = Vec::new();
        for (i, mer) in self.sequence.iter().enumerate() {
            println!("{}", mer ^ rhs.sequence[i]);
            xor_sequence.push(mer ^ rhs.sequence[i]);
        }
        Kmer {
            k: self.k,
            sequence: xor_sequence,
        }
    }
}

impl Not for Kmer {
    type Output = Self;
    fn not(self) -> Self::Output {
        let mut not_sequence: Vec<u8> = Vec::new();
        for mer in self.sequence.iter() {
            not_sequence.push(!mer);
        }
        Kmer {
            k: self.k,
            sequence: not_sequence,
        }
    }
}

//General utility function
pub fn byte_to_nuc(byte: u8) -> char {
    match byte {
        0 => {'A'}
        1 => {'G'}
        2 => {'C'}
        3 => {'T'}
        _ => {panic!("Non-valid nucleotide detected!")}
    }
}

pub fn nuc_to_byte(nuc: char) -> u8 {
    match nuc {
        'A' => {0}
        'G' => {1}
        'C' => {2}
        'T' => {3}
        _ => {panic!("Non-valid nucleotide detected!")}
    }
}

//TESTS

#[cfg(test)]
mod tests {
    use super::Kmer;
    use crate::data_structures::kmer::byte_to_nuc;


    //TODO: WRITE SOME ACTUALLY GOOD UNIT TESTS
    #[test]
    fn test_vmer_instantiations() {
        let _kmer_literal = Kmer::from_literal("ATGCATGCATGCATGCATGCATGC");
        let sequence = String::from("AAAAATTTTTGGGGGCCCCC");
        let k = 5;
        for kmer in sequence.as_bytes().windows(k) {
            let mut kmer1 = Kmer::new(k, kmer);
            //let vmer2 = Vmer::new(k, kmer);
            println!("Kmer:    {}", kmer1);
            //println!("Rev:     {}", kmer1.retain_reverse());
            kmer1.reverse_complement();
            println!("RevComp: {}", kmer1);
        }
    }

    #[test]
    fn test_nucleotide_iterator() {
        let kmer_literal = Kmer::from_literal("ATGCATGCATGCATGCATGCATGC");
        println!("{}", kmer_literal);
        for i in kmer_literal {
            println!("{}", byte_to_nuc(i));
        }
    }
}