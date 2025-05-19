mod fasta_parser;
mod levenshtein;

use wasm_bindgen::prelude::*;
use std::sync::Arc;

use arrow::array::{ArrayRef, PrimitiveArray, StringBuilder as ArrowStringBuilder};
use arrow::datatypes::{DataType, Field, Schema, UInt16Type, UInt32Type, UInt8Type};
use arrow::error::Result as ArrowResult;
use arrow::ipc::writer::FileWriter;
use arrow::record_batch::RecordBatch;

use serde::Serialize;

const GRID_SPACING: usize = 1_000;
const CHUNK_SIZE_1: usize = 10;
const CHUNK_SIZE_2: usize = 100;
const DIST_THRESHOLD_1: usize = 100_000;
const DIST_THRESHOLD_2: usize = 1_000_000;
const IPC_ARROW_BATCH_SIZE: usize = 1 << 12; 
const JSON_PLOT_POINT_SAMPLE_RATE: usize = 1;

struct ComputedDistancePointInternal {
    idx1: u32, idx2: u32, distance: u16, resolution_type: u8,
}

#[derive(Serialize)]
struct JsonPlotPoint { x: u32, y: u32, d: u16 }

#[derive(Serialize)]
struct ChromosomeJsonData { name: String, max_idx: u32, points: Vec<JsonPlotPoint> }

// Renamed type parameters to UpperCamelCase
fn process_chromosomes_core<IpcH, JsonH>(
    all_chromosomes_data: Vec<(String, Vec<u8>)>,
    mut ipc_point_handler: IpcH,
    mut json_chrom_handler: JsonH,
    is_ipc_path: bool,
    progress_callback: &js_sys::Function,
) -> Result<(), String>
where
    IpcH: FnMut(&str, ComputedDistancePointInternal) -> Result<(), String>,
    JsonH: FnMut(ChromosomeJsonData) -> Result<(), String>,
{
    let total_chroms = all_chromosomes_data.len();
    for (chrom_idx, (chrom_name, chrom_sequence_data)) in all_chromosomes_data.into_iter().enumerate() {
        let _ = progress_callback.call1(&JsValue::NULL, &JsValue::from_str(&format!(
            "Rust Wasm Core: Chromosome {}/{} '{}' ({} bp)...",
            chrom_idx + 1, total_chroms, chrom_name, chrom_sequence_data.len()
        )));

        let current_chrom_arc: Arc<Vec<u8>> = Arc::new(chrom_sequence_data);
        let chrom_len = current_chrom_arc.len();

        if chrom_len == 0 { continue; }
        let num_grid_points = chrom_len / GRID_SPACING;
        if num_grid_points < 2 { continue; }

        let mut points_for_this_chrom_json: Vec<JsonPlotPoint> = if !is_ipc_path { Vec::new() } else { Vec::with_capacity(0) };
        let mut current_max_idx_for_chrom_json: u32 = 0;
        let mut point_counter_for_json_sampling = 0;
        let mut pairs_processed_for_chrom_progress = 0u64;
        let mut last_progress_update_idx1 = 0;

        for idx1_usize in 0..num_grid_points.saturating_sub(1) {
            let idx1 = idx1_usize as u32;
            let pos1 = idx1_usize * GRID_SPACING;
            for idx2_usize in (idx1_usize + 1)..num_grid_points {
                let idx2 = idx2_usize as u32;
                let pos2 = idx2_usize * GRID_SPACING;
                let genome_dist = pos2 - pos1;
                let (len_to_compare, dist_type_val) = if genome_dist <= DIST_THRESHOLD_1 {
                    (CHUNK_SIZE_1, 0u8)
                } else if genome_dist <= DIST_THRESHOLD_2 {
                    (CHUNK_SIZE_2, 1u8)
                } else {
                    (GRID_SPACING, 2u8)
                };
                
                let seq1 = &current_chrom_arc[pos1 .. pos1 + len_to_compare];
                let seq2 = &current_chrom_arc[pos2 .. pos2 + len_to_compare];
                let dist = levenshtein::levenshtein_distance(seq1, seq2);

                let point_internal = ComputedDistancePointInternal {
                    idx1, idx2, distance: dist, resolution_type: dist_type_val,
                };

                if is_ipc_path {
                    ipc_point_handler(&chrom_name, point_internal)?;
                } else { 
                    point_counter_for_json_sampling += 1;
                    if point_counter_for_json_sampling % JSON_PLOT_POINT_SAMPLE_RATE == 0 {
                        points_for_this_chrom_json.push(JsonPlotPoint { x: idx1, y: idx2, d: dist });
                    }
                    if idx1 > current_max_idx_for_chrom_json { current_max_idx_for_chrom_json = idx1; }
                    if idx2 > current_max_idx_for_chrom_json { current_max_idx_for_chrom_json = idx2; }
                }
                pairs_processed_for_chrom_progress += 1;
            }
            if idx1_usize > last_progress_update_idx1 + (num_grid_points / 20).max(10) || idx1_usize == num_grid_points.saturating_sub(2) {
                let _ = progress_callback.call1(&JsValue::NULL, &JsValue::from_str(&format!(
                    "Rust Wasm Core: Chromosome '{}', grid point {}/{}. ({} pairs processed)", 
                    chrom_name, idx1_usize + 1, num_grid_points.saturating_sub(1), pairs_processed_for_chrom_progress
                )));
                last_progress_update_idx1 = idx1_usize;
            }
        }
        if !is_ipc_path {
             json_chrom_handler(ChromosomeJsonData {
                name: chrom_name.clone(),
                max_idx: current_max_idx_for_chrom_json,
                points: points_for_this_chrom_json,
            })?;
        }
        let _ = progress_callback.call1(&JsValue::NULL, &JsValue::from_str(&format!("Rust Wasm Core: Finished chromosome '{}'. Total pairs: {}", chrom_name, pairs_processed_for_chrom_progress)));
    }
    Ok(())
}

#[wasm_bindgen]
pub fn process_fasta_to_plot_json(
    fasta_content: &[u8], is_gzipped: bool, progress_callback: &js_sys::Function,
) -> Result<String, JsValue> {
    let _ = progress_callback.call1(&JsValue::NULL, &JsValue::from_str("Rust Wasm (JSON): Starting..."));
    let all_chromosomes_data = fasta_parser::load_chromosomes_from_bytes(fasta_content, is_gzipped)
        .map_err(|e: std::io::Error| JsValue::from_str(&format!("FASTA parsing error: {}", e)))?; // Specify Error type for clarity
    if all_chromosomes_data.is_empty() { return Err(JsValue::from_str("No chromosomes in FASTA.")); }
    let _ = progress_callback.call1(&JsValue::NULL, &JsValue::from_str(&format!("Rust Wasm (JSON): Loaded {} chromosomes.", all_chromosomes_data.len())));

    let mut all_plot_data_for_json: Vec<ChromosomeJsonData> = Vec::new();
    
    process_chromosomes_core(
        all_chromosomes_data,
        |_chrom_name, _point_internal| { Ok(()) }, 
        |chrom_json_data| { 
            all_plot_data_for_json.push(chrom_json_data);
            Ok(())
        },
        false, 
        progress_callback
    ).map_err(|e: String| JsValue::from_str(&e))?; // Corrected map_err
    
    let _ = progress_callback.call1(&JsValue::NULL, &JsValue::from_str("Rust Wasm (JSON): Serializing data..."));
    serde_json::to_string(&all_plot_data_for_json)
        .map_err(|e| JsValue::from_str(&format!("JSON serialization error: {}", e)))
}

struct IpcDataBatchInternal {
    chrom_names: ArrowStringBuilder, idx1: Vec<u32>, idx2: Vec<u32>,
    dist_val: Vec<u16>, dist_type: Vec<u8>,
}
impl IpcDataBatchInternal {
    fn new() -> Self {
        IpcDataBatchInternal {
            chrom_names: ArrowStringBuilder::with_capacity(IPC_ARROW_BATCH_SIZE, IPC_ARROW_BATCH_SIZE * 10),
            idx1: Vec::with_capacity(IPC_ARROW_BATCH_SIZE), idx2: Vec::with_capacity(IPC_ARROW_BATCH_SIZE),
            dist_val: Vec::with_capacity(IPC_ARROW_BATCH_SIZE), dist_type: Vec::with_capacity(IPC_ARROW_BATCH_SIZE),
        }
    }
    fn add(&mut self, chrom_name: &str, point: &ComputedDistancePointInternal) {
        self.chrom_names.append_value(chrom_name); self.idx1.push(point.idx1); self.idx2.push(point.idx2);
        self.dist_val.push(point.distance); self.dist_type.push(point.resolution_type);
    }
    fn is_full(&self) -> bool { self.idx1.len() >= IPC_ARROW_BATCH_SIZE }
    fn len(&self) -> usize { self.idx1.len() }
    fn write_to_arrow_and_clear(&mut self, writer: &mut FileWriter<Vec<u8>>, schema: &Arc<Schema>) -> ArrowResult<()> {
        if self.len() == 0 { return Ok(()); }
        let col_chrom_name: ArrayRef = Arc::new(self.chrom_names.finish_cloned());
        let col_idx1: ArrayRef = Arc::new(PrimitiveArray::<UInt32Type>::from(std::mem::take(&mut self.idx1)));
        let col_idx2: ArrayRef = Arc::new(PrimitiveArray::<UInt32Type>::from(std::mem::take(&mut self.idx2)));
        let col_dist_val: ArrayRef = Arc::new(PrimitiveArray::<UInt16Type>::from(std::mem::take(&mut self.dist_val)));
        let col_dist_type: ArrayRef = Arc::new(PrimitiveArray::<UInt8Type>::from(std::mem::take(&mut self.dist_type)));
        let batch = RecordBatch::try_new(schema.clone(), vec![col_chrom_name, col_idx1, col_idx2, col_dist_val, col_dist_type])?;
        writer.write(&batch)?;
        self.idx1.clear(); self.idx2.clear(); self.dist_val.clear(); self.dist_type.clear();
        Ok(())
    }
}

#[wasm_bindgen]
pub fn process_fasta_to_ipc_bytes(
    fasta_content: &[u8], is_gzipped: bool, progress_callback: &js_sys::Function,
) -> Result<Vec<u8>, JsValue> {
    let _ = progress_callback.call1(&JsValue::NULL, &JsValue::from_str("Rust Wasm (IPC): Starting..."));
    let all_chromosomes_data = fasta_parser::load_chromosomes_from_bytes(fasta_content, is_gzipped)
        .map_err(|e: std::io::Error| JsValue::from_str(&format!("FASTA parsing error: {}", e)))?; // Specify Error type
    if all_chromosomes_data.is_empty() { return Err(JsValue::from_str("No chromosomes in FASTA.")); }
    let _ = progress_callback.call1(&JsValue::NULL, &JsValue::from_str(&format!("Rust Wasm (IPC): Loaded {} chromosomes.", all_chromosomes_data.len())));

    let schema = Arc::new(Schema::new(vec![
        Field::new("chromosome", DataType::Utf8, false), Field::new("idx1", DataType::UInt32, false),
        Field::new("idx2", DataType::UInt32, false), Field::new("distance", DataType::UInt16, false),
        Field::new("type", DataType::UInt8, false),
    ]));
    let ipc_buffer_owned: Vec<u8> = Vec::new();
    let mut arrow_writer = FileWriter::try_new(ipc_buffer_owned, &schema)
        .map_err(|e: arrow::error::ArrowError| JsValue::from_str(&format!("Arrow writer init error: {}", e)))?;
    let mut current_ipc_batch = IpcDataBatchInternal::new();

    process_chromosomes_core(
        all_chromosomes_data,
        |chrom_name_for_point, point_internal| {
            current_ipc_batch.add(chrom_name_for_point, &point_internal);
            if current_ipc_batch.is_full() {
                current_ipc_batch.write_to_arrow_and_clear(&mut arrow_writer, &schema)
                    .map_err(|e| format!("IPC batch write error: {}", e))?;
            }
            Ok(())
        },
        |_chrom_json_data| { Ok(()) },
        true, 
        progress_callback
    ).map_err(|e: String| JsValue::from_str(&e))?; // Corrected map_err

    if current_ipc_batch.len() > 0 {
        current_ipc_batch.write_to_arrow_and_clear(&mut arrow_writer, &schema)
            .map_err(|e| JsValue::from_str(&format!("IPC final batch write error: {}", e)))?;
    }
    
    let _ = progress_callback.call1(&JsValue::NULL, &JsValue::from_str("Rust Wasm (IPC): Finalizing IPC stream..."));
    let final_ipc_buffer = arrow_writer.into_inner()
        .map_err(|e| JsValue::from_str(&format!("Arrow IPC finalization error: {}", e)))?;
    let _ = progress_callback.call1(&JsValue::NULL, &JsValue::from_str("Rust Wasm (IPC): Processing complete."));
    Ok(final_ipc_buffer)
}