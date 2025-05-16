mod fasta_parser;
mod levenshtein;

// Use specific sub-crates if needed, but top-level arrow should re-export
use arrow::array::{ArrayRef, PrimitiveArray};
use arrow::datatypes::{DataType, Field, Schema, UInt16Type, UInt32Type, UInt8Type};
use arrow::error::{ArrowError, Result as ArrowResult}; // Make sure ArrowError is in scope
use arrow::ipc::writer::FileWriter;
use arrow::record_batch::RecordBatch;

use crossbeam_channel::{bounded, Sender};
// polars::prelude::PolarsError is no longer needed for direct error conversion from ArrowError
use rayon::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use std::sync::Arc;
use std::thread;

const GRID_SPACING: usize = 1_000;
const CHUNK_SIZE_1: usize = 10;
const CHUNK_SIZE_2: usize = 100;

const DIST_THRESHOLD_1: usize = 100_000;
const DIST_THRESHOLD_2: usize = 1_000_000;

const ARROW_BATCH_SIZE: usize = 1 << 16;

struct ResultBatch {
    idx1: Vec<u32>,
    idx2: Vec<u32>,
    dist_val: Vec<u16>,
    dist_type: Vec<u8>,
}

impl ResultBatch {
    fn new() -> Self {
        ResultBatch {
            idx1: Vec::with_capacity(ARROW_BATCH_SIZE),
            idx2: Vec::with_capacity(ARROW_BATCH_SIZE),
            dist_val: Vec::with_capacity(ARROW_BATCH_SIZE),
            dist_type: Vec::with_capacity(ARROW_BATCH_SIZE),
        }
    }

    fn add(&mut self, i1: u32, i2: u32, dv: u16, dt: u8) {
        self.idx1.push(i1);
        self.idx2.push(i2);
        self.dist_val.push(dv);
        self.dist_type.push(dt);
    }

    fn is_full(&self) -> bool {
        self.idx1.len() >= ARROW_BATCH_SIZE
    }

    fn is_empty(&self) -> bool {
        self.idx1.is_empty()
    }
}

#[derive(Debug)]
struct ChannelSendError;
impl std::fmt::Display for ChannelSendError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Worker failed to send batch to writer thread")
    }
}
impl std::error::Error for ChannelSendError {}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args();
    args.next();
    let fasta_path = args.next().ok_or_else(|| "Usage: program <fasta_file> <output_ipc_file> (missing fasta_file)".to_string())?;
    let output_path = args.next().ok_or_else(|| "Usage: program <fasta_file> <output_ipc_file> (missing output_ipc_file)".to_string())?;

    println!("Loading chromosome sequence from: {}", fasta_path);
    let chromosome_sequence = fasta_parser::load_chromosome_sequence(&fasta_path)
        .map_err(|e| format!("Failed to load FASTA file '{}': {}", fasta_path, e))?;
    println!("Chromosome sequence length: {} base pairs.", chromosome_sequence.len());

    let chrom_len = chromosome_sequence.len();
    if chrom_len == 0 {
        eprintln!("Chromosome sequence is empty. Exiting.");
        return Ok(());
    }

    let num_grid_points = chrom_len / GRID_SPACING;
    if num_grid_points < 2 {
        eprintln!("Chromosome too short (length: {} bp) for any pairs on the {} bp grid. Needs at least {} bp. Exiting.",
                  chrom_len, GRID_SPACING, 2 * GRID_SPACING);
        return Ok(());
    }
    let total_pairs = if num_grid_points > 1 { num_grid_points * (num_grid_points - 1) / 2 } else {0};
    println!("Number of {} bp grid points: {}. This will result in {} pairwise comparisons.",
             GRID_SPACING, num_grid_points, total_pairs);

    let chrom_arc = Arc::new(chromosome_sequence);

    let schema = Arc::new(Schema::new(vec![
        Field::new("idx1", DataType::UInt32, false),
        Field::new("idx2", DataType::UInt32, false),
        Field::new("distance", DataType::UInt16, false),
        Field::new("type", DataType::UInt8, false),
    ]));

    let file = File::create(&output_path)
        .map_err(|e| format!("Failed to create output file '{}': {}", output_path, e))?;
    let buf_writer = BufWriter::with_capacity(128 * 1024, file);
    // API Change: FileWriter::new -> FileWriter::try_new
    let mut arrow_writer = FileWriter::try_new(buf_writer, &schema)?;


    let num_threads_for_pool = num_cpus::get();
    println!("Using Rayon thread pool with up to {} threads for computation.", num_threads_for_pool);

    let (tx, rx) = bounded::<ResultBatch>(num_threads_for_pool.max(1) * 2);

    let writer_thread_schema = schema.clone();
    let writer_thread = thread::spawn(move || -> ArrowResult<()> { // ArrowResult is arrow::error::Result
        let mut batches_written = 0;
        let mut total_rows_written = 0;
        loop {
            match rx.recv() {
                Ok(mut result_batch) => {
                    let num_rows_in_batch = result_batch.idx1.len();
                    // API Change: PrimitiveArray::from_vec -> PrimitiveArray::from
                    let col_idx1: ArrayRef = Arc::new(PrimitiveArray::<UInt32Type>::from(std::mem::take(&mut result_batch.idx1)));
                    let col_idx2: ArrayRef = Arc::new(PrimitiveArray::<UInt32Type>::from(std::mem::take(&mut result_batch.idx2)));
                    let col_dist_val: ArrayRef = Arc::new(PrimitiveArray::<UInt16Type>::from(std::mem::take(&mut result_batch.dist_val)));
                    let col_dist_type: ArrayRef = Arc::new(PrimitiveArray::<UInt8Type>::from(std::mem::take(&mut result_batch.dist_type)));

                    let record_batch = RecordBatch::try_new(
                        writer_thread_schema.clone(),
                        vec![col_idx1, col_idx2, col_dist_val, col_dist_type],
                    )?;
                    arrow_writer.write(&record_batch)?;
                    batches_written += 1;
                    total_rows_written += num_rows_in_batch;
                    if batches_written % 100 == 0 {
                        println!("Writer thread: Written {} batches ({} rows total) to IPC file.", batches_written, total_rows_written);
                    }
                }
                Err(_) => {
                    println!("Writer thread: All data received. Finalizing Arrow file.");
                    arrow_writer.finish()?;
                    println!("Writer thread: Arrow file finished. Total batches written: {}, total rows: {}.", batches_written, total_rows_written);
                    break Ok(());
                }
            }
        }
    });

    let computation_result = (0..num_grid_points.saturating_sub(1))
        .into_par_iter()
        .try_for_each_with(tx, |thread_tx, idx1| -> Result<(), ChannelSendError> {
            let mut current_batch = ResultBatch::new();
            let local_chrom_arc = Arc::clone(&chrom_arc);

            let pos1 = idx1 * GRID_SPACING;

            for idx2 in (idx1 + 1)..num_grid_points {
                let pos2 = idx2 * GRID_SPACING;
                let genome_dist = pos2 - pos1;

                let (len_to_compare, dist_type_val) = if genome_dist <= DIST_THRESHOLD_1 {
                    (CHUNK_SIZE_1, 0u8)
                } else if genome_dist <= DIST_THRESHOLD_2 {
                    (CHUNK_SIZE_2, 1u8)
                } else {
                    (GRID_SPACING, 2u8)
                };

                let seq1 = &local_chrom_arc[pos1 .. pos1 + len_to_compare];
                let seq2 = &local_chrom_arc[pos2 .. pos2 + len_to_compare];

                let dist = levenshtein::levenshtein_distance(seq1, seq2);
                current_batch.add(idx1 as u32, idx2 as u32, dist, dist_type_val);

                if current_batch.is_full() {
                    if thread_tx.send(current_batch).is_err() {
                        eprintln!("Error: Worker (processing grid index idx1={}) failed to send batch. Writer thread might be down.", idx1);
                        return Err(ChannelSendError);
                    }
                    current_batch = ResultBatch::new();
                }
            }

            if !current_batch.is_empty() {
                if thread_tx.send(current_batch).is_err() {
                     eprintln!("Error: Worker (processing grid index idx1={}) failed to send final batch. Writer thread might be down.", idx1);
                    return Err(ChannelSendError);
                }
            }
            Ok(())
        });

    if let Err(e) = computation_result {
        eprintln!("A worker thread encountered an error: {}. Computation has been aborted.", e);
    } else {
        println!("All computational tasks completed successfully and sent data to writer thread.");
    }

    match writer_thread.join() {
        Ok(Ok(_)) => println!("Writer thread finished successfully and Arrow IPC file finalized."),
        Ok(Err(arrow_err)) => {
            eprintln!("Writer thread failed with Arrow error: {:?}", arrow_err);
            // API Change: Directly box ArrowError as it implements std::error::Error
            return Err(Box::new(arrow_err));
        }
        Err(panic_payload) => {
            eprintln!("Writer thread panicked: {:?}", panic_payload);
            // Use a simple string error for panics
            let panic_msg = if let Some(s) = panic_payload.downcast_ref::<String>() {
                s.clone()
            } else if let Some(s) = panic_payload.downcast_ref::<&str>() {
                s.to_string()
            } else {
                "Writer thread panicked with an unknown type".to_string()
            };
            return Err(Box::new(std::io::Error::new(std::io::ErrorKind::Other, panic_msg)));
        }
    }

    println!("Program finished. Output written to specified IPC file.");
    Ok(())
}