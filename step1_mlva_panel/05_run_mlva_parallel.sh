#!/bin/bash
# ============================================================================
# MLVA PARALLEL ANALYSIS - COMPLETE WORKFLOW
# ============================================================================
# This is the MAIN script that does EVERYTHING:
# 1. Splits primers into chunks
# 2. Runs MLVA analysis in parallel on all chunks
# 3. Merges all results into final CSV files
#
# JUST RUN THIS ONE SCRIPT - IT HANDLES EVERYTHING!
# ============================================================================

set -e  # Exit on any error

echo ""
echo "============================================================================"
echo "🚀 MLVA PARALLEL ANALYSIS - COMPLETE WORKFLOW"
echo "============================================================================"
echo ""

# ============================================================================
# CONFIGURATION - EDIT THESE VARIABLES
# ============================================================================

PRIMERS_FILE="primers_GCA_045016705.1_TRF_fixed.txt"
GENOMES_DIR="genomes_Dietzia/"
OUTPUT_DIR="splitted_jobs_MLVA"
NUM_CHUNKS=8                    # Split primers into 8 chunks for parallelization
THREADS_PER_CHUNK=2             # Use 2 threads per chunk (8 chunks × 2 = 16 total threads)

# ============================================================================
# VALIDATION
# ============================================================================

echo "📋 CONFIGURATION:"
echo "   Primers file:     $PRIMERS_FILE"
echo "   Genomes directory: $GENOMES_DIR"
echo "   Output directory:  $OUTPUT_DIR"
echo "   Number of chunks:  $NUM_CHUNKS"
echo "   Threads per chunk: $THREADS_PER_CHUNK"
echo ""

# Check if required files exist
if [ ! -f "$PRIMERS_FILE" ]; then
    echo "❌ ERROR: Primers file not found: $PRIMERS_FILE"
    exit 1
fi

if [ ! -d "$GENOMES_DIR" ]; then
    echo "❌ ERROR: Genomes directory not found: $GENOMES_DIR"
    exit 1
fi

echo "✅ All input files found"
echo ""

# ============================================================================
# STEP 1: CREATE OUTPUT DIRECTORY
# ============================================================================

echo "============================================================================"
echo "STEP 1: Creating output directory"
echo "============================================================================"
echo ""

mkdir -p "$OUTPUT_DIR"
echo "✅ Output directory ready: $OUTPUT_DIR"
echo ""

# ============================================================================
# STEP 2: SPLIT PRIMERS INTO CHUNKS
# ============================================================================

echo "============================================================================"
echo "STEP 2: Splitting primers file into $NUM_CHUNKS chunks"
echo "============================================================================"
echo ""

PRIMERS_CHUNKS_DIR="$OUTPUT_DIR/primers_chunks"

python3 split_primers.py \
    -i "$PRIMERS_FILE" \
    -n $NUM_CHUNKS \
    -o "$PRIMERS_CHUNKS_DIR"

echo ""

# ============================================================================
# STEP 3: RUN MLVA ANALYSIS ON ALL CHUNKS IN PARALLEL
# ============================================================================

echo "============================================================================"
echo "STEP 3: Running MLVA analysis on all $NUM_CHUNKS chunks (PARALLEL)"
echo "============================================================================"
echo ""

START_TIME=$(date +%s)
PIDS=()
FAILED_CHUNKS=()

for chunk in $(seq 1 $NUM_CHUNKS); do
    CHUNK_PRIMERS="$PRIMERS_CHUNKS_DIR/primers_chunk_${chunk}.txt"
    CHUNK_OUTPUT="$OUTPUT_DIR/MLVA_chunk_${chunk}"
    
    echo "📊 [$chunk/$NUM_CHUNKS] Launching chunk analysis..."
    
    # Run in background
    ./MLVA_finder_optimized.py \
      -i "$GENOMES_DIR" \
      -o "$CHUNK_OUTPUT" \
      -p "$CHUNK_PRIMERS" \
      --full-locus-name \
      --predicted-PCR-size-table \
      --flanking-seq 200 \
      --threads $THREADS_PER_CHUNK > "$OUTPUT_DIR/chunk_${chunk}.log" 2>&1 &
    
    PIDS+=($!)
done

echo ""
echo "⏳ Waiting for all $NUM_CHUNKS chunks to complete..."
echo "   Monitor progress with: tail -f $OUTPUT_DIR/chunk_*.log"
echo ""

# Wait for all background jobs and check for failures
job_count=0
for pid in "${PIDS[@]}"; do
    job_count=$((job_count + 1))
    if wait $pid; then
        echo "✅ Chunk $job_count completed successfully"
    else
        echo "❌ Chunk $job_count FAILED (exit code: $?)"
        FAILED_CHUNKS+=($job_count)
    fi
done

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo ""
echo "============================================================================"
echo "CHUNK ANALYSIS COMPLETE"
echo "============================================================================"
echo "⏱️  Total time: ${DURATION}s ($((DURATION / 60))m $((DURATION % 60))s)"
echo ""

# Check for failed chunks
if [ ${#FAILED_CHUNKS[@]} -gt 0 ]; then
    echo "⚠️  WARNING: The following chunks failed:"
    printf '   - Chunk %s\n' "${FAILED_CHUNKS[@]}"
    echo ""
    echo "You can rerun failed chunks manually or check the log files:"
    for chunk in "${FAILED_CHUNKS[@]}"; do
        echo "   tail -f $OUTPUT_DIR/chunk_${chunk}.log"
    done
else
    echo "✅ All chunks completed successfully!"
fi

echo ""

# ============================================================================
# STEP 4: MERGE ALL RESULTS INTO FINAL CSV FILES
# ============================================================================

echo "============================================================================"
echo "STEP 4: Merging results from all chunks"
echo "============================================================================"
echo ""

FINAL_RESULTS_DIR="$OUTPUT_DIR/Final_Results"

python3 merge_mlva_results.py \
    -i "$OUTPUT_DIR" \
    -o "$FINAL_RESULTS_DIR"

echo ""

# ============================================================================
# SUMMARY
# ============================================================================

echo "============================================================================"
echo "✅ MLVA ANALYSIS COMPLETE!"
echo "============================================================================"
echo ""

echo "📁 Output directory structure:"
echo "   $OUTPUT_DIR/"
echo "   ├── primers_chunks/              (Split primer files)"
echo "   │   ├── primers_chunk_1.txt"
echo "   │   ├── primers_chunk_2.txt"
echo "   │   └── ... (up to chunk $NUM_CHUNKS)"
echo "   │"
echo "   ├── MLVA_chunk_1/                (Individual chunk results)"
echo "   ├── MLVA_chunk_2/"
echo "   ├── ..."
echo "   │"
echo "   ├── chunk_*.log                  (Log files for each chunk)"
echo "   │"
echo "   └── Final_Results/               ⭐ YOUR FINAL RESULTS HERE"
echo "       ├── Final_MLVA_analysis.csv          (Main results)"
echo "       ├── Final_predicted_PCR_size_table.csv"
echo "       ├── Final_output.csv"
echo "       └── Final_mismatchs.txt"
echo ""

echo "📊 FINAL RESULT FILES:"
if [ -d "$FINAL_RESULTS_DIR" ]; then
    echo ""
    for file in "$FINAL_RESULTS_DIR"/*; do
        if [ -f "$file" ]; then
            filename=$(basename "$file")
            filesize=$(du -h "$file" | cut -f1)
            lines=$(wc -l < "$file" 2>/dev/null || echo "?")
            echo "   ✅ $filename ($filesize, $lines lines)"
        fi
    done
    echo ""
else
    echo "   ⚠️  Final results directory not found"
fi

echo "============================================================================"
echo ""
echo "🎉 SUCCESS! Your MLVA analysis is ready to use!"
echo ""
echo "📌 IMPORTANT:"
echo "   - Main results file: $FINAL_RESULTS_DIR/Final_MLVA_analysis.csv"
echo "   - All individual chunk results are in: $OUTPUT_DIR/MLVA_chunk_*/"
echo "   - You can delete chunk directories to save space once verified"
echo ""
echo "============================================================================"
echo ""
