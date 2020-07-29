//Declares a structure with a variable length encoding scheme. Memory is allocated
//for each 4 nucleotides as a single u8 in a vector of u8 values.

use std::fmt;
use std::ops::BitXor;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Vmer {
    pub k: usize,
    pub sequence: Vec<u8>,
}

impl Vmer {
    pub fn new(len: usize, byte_seq: &[u8]) -> Self {
        let mut vmer = Vmer
        {   k: len, 
            sequence: Vec::new()
        };
        vmer.encode(byte_seq);
        vmer
    }

    pub fn from_literal(str_literal: &str) -> Self {
        Vmer::new(str_literal.len(), str_literal.as_bytes())
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

    pub fn decode(self) -> String {
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
                        panic!("Oh no! Detected a non valid nucleotide!");
                    }
                }
            }
        }
        byte_seq
    }
}

impl fmt::Display for Vmer {
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

impl BitXor for Vmer {
    type Output = Self;
    fn bitxor(self, rhs: Self) -> Self::Output {
        assert_eq!(self.k, rhs.k);
        let mut xor_sequence: Vec<u8> = Vec::new();
        for (i, mer) in self.sequence.iter().enumerate() {
            xor_sequence.push(mer ^ rhs.sequence[i]);
        }
        Vmer::new(self.k, xor_sequence.as_slice())
    }
}
