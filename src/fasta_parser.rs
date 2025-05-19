use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Error, Cursor}; // Removed ErrorKind

const APPROX_CHROMOSOME_CAPACITY: usize = 100 * 1024 * 1024;

pub fn load_chromosomes_from_bytes<'a>(
    fasta_data: &'a [u8],
    is_gzipped: bool,
) -> Result<Vec<(String, Vec<u8>)>, Error> {
    let cursor = Cursor::new(fasta_data);
    let reader: Box<dyn BufRead + 'a> = if is_gzipped {
        Box::new(BufReader::new(GzDecoder::new(cursor)))
    } else {
        Box::new(BufReader::new(cursor))
    };
    parse_fasta_from_bufread(reader)
}

#[allow(dead_code)]
pub fn load_chromosomes_from_path(path: &str) -> Result<Vec<(String, Vec<u8>)>, Error> {
    let file = File::open(path)
        .map_err(|e| Error::new(e.kind(), format!("Failed to open FASTA file '{}': {}", path, e)))?;
    
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    parse_fasta_from_bufread(reader)
}

fn parse_fasta_from_bufread<'a>(
    mut reader: Box<dyn BufRead + 'a>
) -> Result<Vec<(String, Vec<u8>)>, Error> {
    let mut chromosomes = Vec::new();
    let mut current_sequence = Vec::with_capacity(APPROX_CHROMOSOME_CAPACITY);
    let mut current_header_name: Option<String> = None;
    let mut line_buffer = String::new();

    loop {
        line_buffer.clear();
        match reader.read_line(&mut line_buffer) {
            Ok(0) => break, // EOF
            Ok(_) => {
                let line_trimmed = line_buffer.trim_end_matches(|c| c == '\r' || c == '\n');
                if line_trimmed.starts_with('>') {
                    if let Some(header_name) = current_header_name.take() {
                        if !current_sequence.is_empty() {
                            chromosomes.push((header_name, std::mem::replace(&mut current_sequence, Vec::with_capacity(APPROX_CHROMOSOME_CAPACITY))));
                        }
                    }
                    let header_content = line_trimmed[1..].trim();
                    let new_chrom_name = header_content.split_whitespace().next().unwrap_or("").trim();
                    if new_chrom_name.is_empty() {
                        current_header_name = Some(format!("UnnamedSequence{}", chromosomes.len() + 1));
                    } else {
                        current_header_name = Some(new_chrom_name.to_string());
                    }
                } else if current_header_name.is_some() {
                    for char_byte in line_trimmed.bytes() {
                        let upper_char = char_byte.to_ascii_uppercase();
                        match upper_char {
                            b'A' | b'C' | b'G' | b'T' | b'N' => current_sequence.push(upper_char),
                            _ => { /* Ignore */ }
                        }
                    }
                }
            }
            Err(e) => return Err(Error::new(e.kind(), format!("Error reading line from FASTA data: {}", e))),
        }
    }

    if let Some(header_name) = current_header_name.take() {
        if !current_sequence.is_empty() {
            chromosomes.push((header_name, current_sequence));
        }
    }
    Ok(chromosomes)
}