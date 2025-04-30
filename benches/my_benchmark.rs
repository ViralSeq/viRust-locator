#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]

use criterion::{Criterion, criterion_group, criterion_main};

use bio::alignment::Alignment;
use bio::alignment::pairwise::*;
use bio::pattern_matching::myers::{Myers, long};
use bio::pattern_matching::*;
use virust_locator::*;

fn run_locator() {
    let query = "AAATTAACCCCGCTTTGTGTTAGTTTAAATTGCACAAATTTGAGGAATAGTACTACTACTAATACCACTATAGAGTCTGGTATGGAAAAAGAAATAAAAAATTGCTCTTTCAATATCACCACAAGCATAAAAGATAAGATGCAGAAAGAACATGCATTGTTTTATAACCTTGATATAACACCAATGGATAATAATGATAATAATAATAATACTAATAGTACTTTTTATAGGTTGATAAGTTGTAACACCTCAGTCACTACACAGGCAGTACAGTACAATGCACACATGGAATTAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAGGAGATAGTAATTAGATCTGACAATTTCACAGACAATGCTAAAAGCATAATAGTACACCTGAATGAATCAGTAGTAATTAATTGTACAAGACCCAACAATAATACAAGGAGAAGTATAAATATGGGACCAGGCAGAGCATTTTATACAACAGGAGATATAATAGGAGATATAAGACGA";
    let reference = "HXB2";
    let type_query = "nt".to_string();
    let algorithm = 1;

    let args = virust_locator::config::Args {
        query: query.to_string(),
        reference: reference.to_string(),
        type_query,
        algorithm,
    };

    // Call the locator function with the parsed arguments
    locator::Locator::build(&args).unwrap();
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("run_locator", |b| {
        b.iter(|| {
            run_locator();
        });
    });
}
criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
