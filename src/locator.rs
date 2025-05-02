/// The `Locator` struct and its associated methods are used to locate and align a query sequence
/// against a reference sequence. It provides functionality to calculate alignment details such as
/// percent identity, indels, and aligned strings.
///
/// # Structs
/// - `Locator`: Represents the alignment result, including reference start and end positions,
///   percent identity, indel presence, and aligned strings.
///
/// # Methods
/// - `Locator::new`: Constructs a new `Locator` instance with the given alignment details.
/// - `Locator::build`: Builds a `Locator` instance by aligning a query sequence against a reference
///   sequence using the specified algorithm.
/// - `get_aln`: Performs a semi-global alignment between a query and reference sequence using a
///   scoring function and gap penalties.
/// - `pattern_match`: Uses the Myers bit-parallel algorithm to find approximate matches of a
///   pattern in a text with a maximum allowed distance.
/// - `from_path`: Converts an alignment path into aligned strings, calculates percent identity,
///   and determines the presence of indels.
/// - `algorithm1`: Implements a specific alignment algorithm to align a query sequence against a
///   reference sequence.
///
/// # Modules
/// - `test`: Contains unit tests for the `Locator` struct and its associated methods.
///
/// # Dependencies
/// - `bio`: Provides bioinformatics utilities for alignment and pattern matching.
/// - `std::error::Error`: Used for error handling.
/// - `std::fmt::Display`: Implements display formatting for the `Locator` struct.
///
/// # Usage
/// The `Locator` struct is designed to be used in bioinformatics applications where sequence
/// alignment is required. It supports two algorithms for alignment:
/// - Algorithm 1: A semi-global alignment approach.
/// - Algorithm 2: A combination of pattern matching and refinement.
///
/// The `Locator::build` method determines which algorithm to use based on the query length and
/// user-specified parameters.
///
/// # Example
/// ```rust
/// use virust_locator::locator::Locator;
/// use virust_locator::config::Args;
/// let args = Args {
///     query: vec!["ATTAACAGAGATTTGTGAAGAAATGGAAAAGGAAGGAAAAATTACAAAAATTGGGCCTGAAAATCCATATAACACTCCAATATTTGCCATAAAAAAGAAGGACAGTACTAAGTGGAGAAAATTAGTAGATTTCAGAGAGCTCAATAAAAGAACTCAAGACTTTTGGGAGGTTCAATTAGGAATACCACACCCAGCAGGGTTAAAAAAGAAAAAATCAGTGACAGTACTGGATGTGGGGGATGCATATTTTTCTGTTCCTTTAGATG".to_string()],
///     reference: "HXB2".to_string(),
///     type_query: "nt".to_string(),
///     algorithm: 1,
/// };
///
/// let locator = Locator::build(&args).unwrap().pop().unwrap().unwrap();
/// println!("{}", locator);
/// ```
use crate::config::Args;
use crate::reference::retrieve_reference_sequence;
use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation;
use bio::alignment::pairwise::*;
use bio::pattern_matching::myers::long;
use std::error::Error;
use std::fmt::Display;

#[derive(Debug)]
pub struct Locator {
    pub ref_start: usize, // starting from 1 on reference
    pub ref_end: usize,   // inclusive
    pub percent_identity: f64,
    pub indel: bool,
    pub query_aligned_string: String,
    pub ref_aligned_string: String,
}

impl Display for Locator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.ref_start,
            self.ref_end,
            self.percent_identity,
            self.indel,
            self.query_aligned_string,
            self.ref_aligned_string
        )
    }
}

impl Locator {
    pub fn new(
        ref_start: usize,
        ref_end: usize,
        percent_identity: f64,
        indel: bool,
        query_aligned_string: String,
        ref_aligned_string: String,
    ) -> Self {
        Locator {
            ref_start,
            ref_end,
            percent_identity,
            indel,
            query_aligned_string,
            ref_aligned_string,
        }
    }

    pub fn build(args: &Args) -> Result<Vec<Option<Locator>>, Box<dyn Error>> {
        let query_vec = args
            .query
            .iter()
            .map(|x| x.as_bytes())
            .collect::<Vec<&[u8]>>();

        let ref_seq = retrieve_reference_sequence(&args.reference, &args.type_query)?.sequence;

        let algorithm = args.algorithm;

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

        let mut result_vec = Vec::new();
        for query in query_vec {
            if query.len() < 300 || algorithm == 1 {
                let loc = algorithm1(query, ref_seq, score)?;
                result_vec.push(loc);
            } else {
                let s1 = &query[..100];
                let s2 = &query[query.len() - 100..];

                let aln1 = pattern_match(s1, ref_seq, 30);

                if aln1.is_none() {
                    let loc = algorithm1(query, ref_seq, score)?;
                    result_vec.push(loc);
                    continue;
                }
                let pos_start = aln1.unwrap().ystart as usize;

                let aln2 = pattern_match(s2, ref_seq, 30);

                if aln2.is_none() {
                    let loc = algorithm1(query, ref_seq, score)?;
                    result_vec.push(loc);
                    continue;
                }
                let pos_end = aln2.unwrap().yend as usize;

                let refined_ref = &ref_seq[pos_start..pos_end];

                let mut loc = algorithm1(query, refined_ref, score)?.unwrap();
                loc.ref_start = pos_start + 1;
                loc.ref_end = pos_end;
                result_vec.push(Some(loc));
            }
        }
        return Ok(result_vec);
    }
}

fn get_aln(
    query: &[u8],
    ref_seq: &[u8],
    score: fn(u8, u8) -> i32,
    gap_open: i32,
    gap_extend: i32,
) -> Result<Alignment, Box<dyn Error>> {
    let mut aligner =
        Aligner::with_capacity(query.len(), ref_seq.len(), gap_open, gap_extend, &score);

    Ok(aligner.semiglobal(query, ref_seq))
}

fn pattern_match(pattern: &[u8], text: &[u8], max_dist: usize) -> Option<Alignment> {
    let mut myers = long::Myers::<u64>::new(pattern);
    let mut lazy_matches = myers.find_all_lazy(text, max_dist);
    let mut aln = Alignment::default();
    match lazy_matches.by_ref().min_by_key(|&(_, dist)| dist) {
        Some((best_end, _)) => {
            lazy_matches.alignment_at(best_end, &mut aln);
            return Some(aln);
        }
        None => {
            return None;
        }
    }
}

fn from_path(aln: Alignment, query: &[u8], ref_seq: &[u8]) -> (String, String, f64, bool) {
    let mut ref_string = String::new();
    let mut query_string = String::new();
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut matches = 0;
    for p in aln.path().iter() {
        let (query_pos, ref_pos, state) = p;

        if *state == AlignmentOperation::Match {
            ref_string.push(ref_seq[*ref_pos - 1] as char);
            query_string.push(query[*query_pos - 1] as char);
            matches += 1;
        } else if *state == AlignmentOperation::Subst {
            ref_string.push(ref_seq[*ref_pos - 1] as char);
            query_string.push(query[*query_pos - 1] as char);
            mismatches += 1;
        } else if *state == AlignmentOperation::Ins {
            query_string.push(query[*query_pos - 1] as char);
            ref_string.push('-');
            gaps += 1;
        } else if *state == AlignmentOperation::Del {
            ref_string.push(ref_seq[*ref_pos - 1] as char);
            query_string.push('-');
            gaps += 1;
        }
    }
    let percent_identity = (matches as f64 / (matches + mismatches + gaps) as f64) * 100.0;

    let indel = if gaps > 0 { true } else { false };

    (ref_string, query_string, percent_identity, indel)
}

fn algorithm1(
    query: &[u8],
    ref_seq: &[u8],
    score: fn(u8, u8) -> i32,
) -> Result<Option<Locator>, Box<dyn Error>> {
    let aln = get_aln(query, ref_seq, score, -5, -1)?;
    let ref_start = aln.ystart as usize;
    let ref_end = aln.yend as usize;
    let (ref_aligned_string, query_aligned_string, percent_identity, indel) =
        from_path(aln, query, ref_seq);

    let loc = Locator {
        ref_start: ref_start + 1,
        ref_end,
        percent_identity,
        indel,
        query_aligned_string,
        ref_aligned_string,
    };
    Ok(Some(loc))
}

#[cfg(test)]
mod test {
    use super::*;

    static ONE_LOC: (i32, i32, f64, bool, &'static str, &'static str) = (
        2648,
        3209,
        83.98576512455516,
        true,
        "ATTAACAGAGATTTGTGAAGAAATGGAAAAGGAAGGAAAAATTACAAAAATTGGGCCTGAAAATCCATATAACACTCCAATATTTGCCATAAAAAAGAAGGACAGTACTAAGTGGAGAAAATTAGTAGATTTCAGAGAGCTCAATAAAAGAACTCAAGACTTTTGGGAGGTTCAATTAGGAATACCACACCCAGCAGGGTTAAAAAAGAAAAAATCAGTGACAGTACTGGATGTGGGGGATGCATATTTTTCTGTTCCTTTAGATG-----------------------------------AGTGTAAACAATGAAACACCAGGGATTAGATATCAATATAATGTGCTACCACAGGGGTGGAAAGGATCACCATCAATATTCCAGAGTAGCATGACAAAAATCTTAGAGCCCTTTAGAGCAAAAAACCCAGAAATAGTCATCTATCAATATATGGATGACTTATGTGTAGGATCTGACTTAGAAATAGGGCAACATAGAGCAAAAATAGAGGAGTTAAGAGAACATCTATTGAAGTGGGGATTGACCACACCAGACAAGAAA",
        "ATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAA",
    );

    static TWO_LOC: (i32, i32, f64, bool, &'static str, &'static str) = (
        6585,
        7208,
        83.98576512455516,
        true,
        "ATTAACAGAGATTTGTGAAGAAATGGAAAAGGAAGGAAAAATTACAAAAATTGGGCCTGAAAATCCATATAACACTCCAATATTTGCCATAAAAAAGAAGGACAGTACTAAGTGGAGAAAATTAGTAGATTTCAGAGAGCTCAATAAAAGAACTCAAGACTTTTGGGAGGTTCAATTAGGAATACCACACCCAGCAGGGTTAAAAAAGAAAAAATCAGTGACAGTACTGGATGTGGGGGATGCATATTTTTCTGTTCCTTTAGATG-----------------------------------AGTGTAAACAATGAAACACCAGGGATTAGATATCAATATAATGTGCTACCACAGGGGTGGAAAGGATCACCATCAATATTCCAGAGTAGCATGACAAAAATCTTAGAGCCCTTTAGAGCAAAAAACCCAGAAATAGTCATCTATCAATATATGGATGACTTATGTGTAGGATCTGACTTAGAAATAGGGCAACATAGAGCAAAAATAGAGGAGTTAAGAGAACATCTATTGAAGTGGGGATTGACCACACCAGACAAGAAA",
        "ATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAGACTTCAGGAAGTATACTGCATTTACCATACCTAAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAA",
    );

    static MY_ARGS: (&'static str, &'static str, &'static str, u8) = (
        "ATTAACAGAGATTTGTGAAGAAATGGAAAAGGAAGGAAAAATTACAAAAATTGGGCCTGAAAATCCATATAACACTCCAATATTTGCCATAAAAAAGAAGGACAGTACTAAGTGGAGAAAATTAGTAGATTTCAGAGAGCTCAATAAAAGAACTCAAGACTTTTGGGAGGTTCAATTAGGAATACCACACCCAGCAGGGTTAAAAAAGAAAAAATCAGTGACAGTACTGGATGTGGGGGATGCATATTTTTCTGTTCCTTTAGATGAGTGTAAACAATGAAACACCAGGGATTAGATATCAATATAATGTGCTACCACAGGGGTGGAAAGGATCACCATCAATATTCCAGAGTAGCATGACAAAAATCTTAGAGCCCTTTAGAGCAAAAAACCCAGAAATAGTCATCTATCAATATATGGATGACTTATGTGTAGGATCTGACTTAGAAATAGGGCAACATAGAGCAAAAATAGAGGAGTTAAGAGAACATCTATTGAAGTGGGGATTGACCACACCAGACAAGAAA",
        "HXB2",
        "nt",
        1,
    );

    static MY_ARGS2: (&'static str, &'static str, &'static str, u8) = (
        "AAATTAACCCCACTCTGTGTTGAATTAAATTGTACTAAGTATGAGGGTAATAGTACTACTACCACGAATAGTACTACTGCCACTACGAATAGTACTGCTGCCCCTAACGGGACGGAGACGGGAATGAAAAATTGCTCTTTCTATGTTAACACGGTCACAAACTATAAGGTGCAGAAGAAATATGCACTTTTCTATGATCTTGATATAGTACAAATAGAAGGTAGTAATACTAGCTATAGGATAACAAAGTGTAACACCTCAATCAGCACAGTACAATGCACACATGGTATTAAACCAGTAGTATCAACTCAATTATTGTTAAATGGCAGCTTAGCAGAAGAAAAGATAGTCATCAGATCTAGCAACTTCTCTAGCAACACTGAAAGCATAATAGTACAGCTGAAAAACCCTGTAGAAATTAACTGTACAAGACCCAACAACAATAGAAGACAGAGTATCCATATTGGACCAGGGAGAGCGTTTTTTACAACAGGAGAAATAATAGGAGATATAAGACAA",
        "HXB2",
        "nt",
        1,
    );

    #[test]
    fn test_get_aln() {
        let search_string = b"AAATTAACCCCACTCTGTGTTGAATTAAATTGTACTAAGTATGAGGGTAATAGTACTACTACCACGAATAGTACTACTGCCACTACGAATAGTACTGCTGCCCCTAACGGGACGGAGACGGGAATGAAAAATTGCTCTTTCTATGTTAACACGGTCACAAACTATAAGGTGCAGAAGAAATATGCACTTTTCTATGATCTTGATATAGTACAAATAGAAGGTAGTAATACTAGCTATAGGATAACAAAGTGTAACACCTCAATCAGCACAGTACAATGCACACATGGTATTAAACCAGTAGTATCAACTCAATTATTGTTAAATGGCAGCTTAGCAGAAGAAAAGATAGTCATCAGATCTAGCAACTTCTCTAGCAACACTGAAAGCATAATAGTACAGCTGAAAAACCCTGTAGAAATTAACTGTACAAGACCCAACAACAATAGAAGACAGAGTATCCATATTGGACCAGGGAGAGCGTTTTTTACAACAGGAGAAATAATAGGAGATATAAGACAA";
        let ref_seq = retrieve_reference_sequence("HXB2", "nt").unwrap().sequence;
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let gap_open = -5;
        let gap_extend = -1;

        let aln = get_aln(search_string, ref_seq, score, gap_open, gap_extend).unwrap();
        assert_eq!(aln.ystart, 6584);
        assert_eq!(aln.yend, 7208);
    }

    #[test]
    fn test_locator_1() {
        let targe_loc = Locator::new(
            ONE_LOC.0 as usize,
            ONE_LOC.1 as usize,
            ONE_LOC.2,
            ONE_LOC.3,
            ONE_LOC.4.to_string(),
            ONE_LOC.5.to_string(),
        );

        let my_arg = Args {
            query: vec![MY_ARGS.0.to_string()],
            reference: MY_ARGS.1.to_string(),
            type_query: MY_ARGS.2.to_string(),
            algorithm: MY_ARGS.3,
        };

        let loc = Locator::build(&my_arg).unwrap().pop().unwrap().unwrap();

        assert_eq!(loc.ref_start, targe_loc.ref_start);
        assert_eq!(loc.ref_end, targe_loc.ref_end);
        assert_eq!(loc.percent_identity, targe_loc.percent_identity);
        assert_eq!(loc.indel, targe_loc.indel);
        assert_eq!(loc.query_aligned_string, targe_loc.query_aligned_string);
        assert_eq!(loc.ref_aligned_string, targe_loc.ref_aligned_string);
    }

    #[test]
    #[should_panic]
    fn test_locator_2() {
        let targe_loc = Locator::new(
            ONE_LOC.0 as usize,
            ONE_LOC.1 as usize,
            ONE_LOC.2,
            ONE_LOC.3,
            ONE_LOC.4.to_string(),
            ONE_LOC.5.to_string(),
        );

        let my_arg = Args {
            query: vec![MY_ARGS.0.to_string()],
            reference: MY_ARGS.1.to_string(),
            type_query: MY_ARGS.2.to_string(),
            algorithm: 2,
        };

        let loc = Locator::build(&my_arg).unwrap().pop().unwrap().unwrap();

        assert_eq!(loc.ref_start, targe_loc.ref_start);
        assert_eq!(loc.ref_end, targe_loc.ref_end);
        assert_eq!(loc.percent_identity, targe_loc.percent_identity);
        assert_eq!(loc.indel, targe_loc.indel);
        assert_eq!(loc.query_aligned_string, targe_loc.query_aligned_string);
        assert_eq!(loc.ref_aligned_string, targe_loc.ref_aligned_string);
    }

    #[test]
    fn test_locator_3() {
        let targe_loc = Locator::new(
            TWO_LOC.0 as usize,
            TWO_LOC.1 as usize,
            TWO_LOC.2,
            TWO_LOC.3,
            TWO_LOC.4.to_string(),
            TWO_LOC.5.to_string(),
        );

        let my_arg = Args {
            query: vec![MY_ARGS2.0.to_string()],
            reference: MY_ARGS2.1.to_string(),
            type_query: MY_ARGS2.2.to_string(),
            algorithm: 1,
        };

        let loc = Locator::build(&my_arg).unwrap().pop().unwrap().unwrap();

        assert_eq!(loc.ref_start, targe_loc.ref_start);
        assert_eq!(loc.ref_end, targe_loc.ref_end);
    }

    #[test]
    fn test_locator_4() {
        let targe_loc = Locator::new(
            TWO_LOC.0 as usize,
            TWO_LOC.1 as usize,
            TWO_LOC.2,
            TWO_LOC.3,
            TWO_LOC.4.to_string(),
            TWO_LOC.5.to_string(),
        );

        let my_arg = Args {
            query: vec![MY_ARGS2.0.to_string()],
            reference: MY_ARGS2.1.to_string(),
            type_query: MY_ARGS2.2.to_string(),
            algorithm: 2,
        };

        let loc = Locator::build(&my_arg).unwrap().pop().unwrap().unwrap();

        assert_eq!(loc.ref_start, targe_loc.ref_start);
        assert_eq!(loc.ref_end, targe_loc.ref_end);
    }
}
