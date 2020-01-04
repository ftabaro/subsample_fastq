extern crate rand;

use std::io::prelude::*;
use std::io::{BufReader,BufWriter};
use std::fs::File;
use std::env;

use rand::distributions::{Distribution, Uniform};

use flate2::Compression;
use flate2::write::GzEncoder;
use flate2::read::GzDecoder;

fn extract_reads(p: &str, selected_lines: &Vec<usize>, outpath: &str) -> Result<(), std::io::Error> {
    //let f = File::open(p)?;
    let d = GzDecoder::new(p)?;
    let reader = BufReader::new(f);

    let fo = File::create(outpath)?;
    let mut writer = BufWriter::new(fo);

    let mut e = GzEncoder::new(Vec::new(), Compression::default());

    for tpl in reader.lines().enumerate() {
        let ( line_number, line ) = tpl;
        if selected_lines.contains(&(line_number+1)) {
            let mut line_content = line.unwrap();
            line_content += "\n";
            e.write_all(line_content.as_bytes())?;
        }
    }
    let compressed_bytes = e.finish()?;
    writer.write(&compressed_bytes)?;
    Ok(())
}

fn count_lines(p: &str) -> Result<usize, std::io::Error> {
    let d = GzDecoder::new(p)?;
    //let f = File::open(p)?;
    let reader = BufReader::new(f);
    Ok(reader.lines().count())
}

fn compute_read_lines(n: &usize) -> Vec<usize> {
    let last_line: usize = n * 4;
    let range: Vec<usize> = (last_line-3..=last_line).collect();
    range
}

fn select_reads(nreads: f32, fraction: f32) -> Vec<usize> {
    let mut selected_lines: Vec<usize> = Vec::new();
    let reads_range = 1..*&nreads as usize;

    if (fraction > 0.0) & (fraction < 1.0){

        let nselected: usize = (nreads * fraction).round() as usize;
        eprintln!("{} ({}%) reads will be sampled.", nselected, fraction * 100.0);

        let mut rng = rand::thread_rng();
        let distribution = Uniform::from(reads_range);

        let mut i = 0;
        let selected_reads: Vec<usize> = Vec::new();
        while i < nselected {
            let n = distribution.sample(&mut rng);
            let range = compute_read_lines(&n);
            if !selected_reads.contains(&n) {
                selected_lines.extend(range);
                i = i+1;
            }
        }

        selected_lines.sort();

    } else if fraction == 1.0 {
        for n in reads_range {
            let range: Vec<usize> = compute_read_lines(&n);
            selected_lines.extend(range);
        }
    } else {
        panic!("Invalid fraction.")
    }

    selected_lines
}

fn main() -> std::io::Result<()> {

    let args: Vec<String> = env::args().collect();

    let fraction = &args[1].parse::<f32>().expect("Fraction is not a number!");

    let len_fq1 = count_lines(&args[2])?;
    let len_fq2 = count_lines(&args[3])?;
    assert_eq!(len_fq1, len_fq2);

    let nreads: f32 = len_fq1 as f32 / 4.0;
    let selected_lines: Vec<usize> = select_reads(nreads, *fraction);

    let subsampled_fq1 = &args[4];
    let subsampled_fq2 = &args[5];

    eprintln!("Output mate 1: {}", subsampled_fq1);
    eprintln!("Output mate 2: {}", subsampled_fq2);

    extract_reads(&args[2], &selected_lines, &subsampled_fq1)?;
    extract_reads(&args[3], &selected_lines, &subsampled_fq2)?;

    Ok(())
}