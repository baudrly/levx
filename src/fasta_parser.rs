use flate2::read::GzDecoder;
use std::fs::File;
// Removed unused std::io::Read
use std::io::{BufRead, BufReader, Error, ErrorKind}; 

const APPROX_CHROMOSOME_CAPACITY: usize = 100 * 1024 * 1024;

pub fn load_chromosomes(path: &str) -> Result<Vec<(String, Vec<u8>)>, Error> {
    let file = File::open(path)
        .map_err(|e| Error::new(e.kind(), format!("Failed to open FASTA file '{}': {}", path, e)))?;
    
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        println!("Detected .gz extension for '{}', reading as gzipped FASTA.", path);
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut chromosomes = Vec::new();
    let mut current_sequence = Vec::new();
    let mut current_header_name: Option<String> = None;

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| Error::new(e.kind(), format!("Error reading line from FASTA file '{}': {}", path, e)))?;
        
        if line.starts_with('>') {
            if let Some(header_name) = current_header_name.take() {
                if !current_sequence.is_empty() {
                    chromosomes.push((header_name, std::mem::take(&mut current_sequence)));
                } else {
                    eprintln!("Warning: Chromosome/sequence entry '{}' in file '{}' had no sequence data. Skipping.", header_name, path);
                }
            }
            let header_content = line[1..].trim();
            let new_chrom_name = header_content.split_whitespace().next().unwrap_or(header_content);
            if new_chrom_name.is_empty() {
                 return Err(Error::new(ErrorKind::InvalidData, format!("Encountered an empty chromosome name after '>' in file '{}'. Line: '{}'", path, line)));
            }
            current_header_name = Some(new_chrom_name.to_string());
            current_sequence = Vec::with_capacity(APPROX_CHROMOSOME_CAPACITY);
        } else if current_header_name.is_some() {
            for char_byte in line.trim().bytes() {
                let upper_char = char_byte.to_ascii_uppercase();
                match upper_char {
                    b'A' | b'C' | b'G' | b'T' | b'N' => current_sequence.push(upper_char),
                    _ => { /* Ignore other characters */ }
                }
            }
        }
    }
    if let Some(header_name) = current_header_name.take() {
        if !current_sequence.is_empty() {
            chromosomes.push((header_name, current_sequence));
        } else {
            eprintln!("Warning: Last chromosome/sequence entry '{}' in file '{}' had no sequence data. Skipping.", header_name, path);
        }
    }

    if chromosomes.is_empty() && (path.ends_with(".fa") || path.ends_with(".fasta") || path.ends_with(".fa.gz") || path.ends_with(".fasta.gz")) {
         println!("Warning: No valid chromosome sequences found in '{}'. Output will be empty if this was the only input.", path);
    }
    Ok(chromosomes)
}