use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind};

// A typical human chromosome is < 300 Mbp. X chromosome is ~156 Mbp.
// Pre-allocating avoids multiple reallocations during file reading.
const APPROX_CHROMOSOME_CAPACITY: usize = 200 * 1024 * 1024; // 200 MiB for base pairs

/// Loads the first sequence from a FASTA file (can be gzipped if path ends with .gz).
/// Converts characters to uppercase and keeps only A, C, G, T, N.
pub fn load_chromosome_sequence(path: &str) -> Result<Vec<u8>, Error> {
    let file = File::open(path)
        .map_err(|e| Error::new(e.kind(), format!("Failed to open FASTA file '{}': {}", path, e)))?;
    
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        println!("Detected .gz extension for '{}', attempting to read as gzipped FASTA.", path);
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut sequence = Vec::with_capacity(APPROX_CHROMOSOME_CAPACITY);
    let mut found_sequence_header = false;

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| Error::new(e.kind(), format!("Error reading line from FASTA file '{}': {}", path, e)))?;
        if line.starts_with('>') {
            if found_sequence_header {
                // Found another sequence header, stop after the first one.
                break;
            }
            found_sequence_header = true;
        } else if found_sequence_header {
            for char_byte in line.trim().bytes() {
                let upper_char = char_byte.to_ascii_uppercase();
                match upper_char {
                    b'A' | b'C' | b'G' | b'T' | b'N' => sequence.push(upper_char),
                    _ => { /* Ignore other characters */ }
                }
            }
        }
    }

    if !found_sequence_header && sequence.is_empty() {
         Err(Error::new(ErrorKind::InvalidData, format!("No FASTA sequence header '>' found in the file '{}'.", path)))
    } else if sequence.is_empty() {
        Err(Error::new(ErrorKind::InvalidData, format!("No valid sequence data found after FASTA header in file '{}'.", path)))
    } else {
        // Optional: If memory is extremely tight and exact sizing is critical.
        // sequence.shrink_to_fit(); 
        Ok(sequence)
    }
}