#!/bin/bash

#SBATCH --job-name=caliby_af3
#SBATCH --output=logs/caliby_af3_%A.out
#SBATCH --error=logs/caliby_af3_%A.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16

#
# Pipeline: Caliby Sequence Design → AlphaFold3 Structure Prediction
# This script runs the Caliby-AF3 pipeline on Sherlock
#
# Usage:
#   sbatch run_caliby_af3_sherlock.sbatch
# ==================== Configuration ====================

# Your workspace directory on Sherlock
WORKDIR="/scratch/users/$USER/protein_design"

# Path to Singularity image (with Caliby) - same image as Caliby-Boltz pipeline
IMG="/scratch/users/$USER/containers/rfdiff_proteinmpnn_boltz_caliby.sif"

# AlphaFold3 paths (from shared group location)
AF3_SCRIPT="/home/groups/possu/alphafold3/run_af3.sh"
AF3_IMAGE="/home/groups/possu/alphafold3/alphafold3_latest.sif"
AF3_MODEL_DIR="/oak/stanford/groups/possu/alphafold3"
AF3_DB_DIR="/oak/stanford/groups/possu/alphafold3/databases"

# Model weights directories (for Caliby)
RFDIFF_WEIGHTS_DIR="/scratch/users/$USER/rfdiff_models/rfdiffusion"
MPNN_WEIGHTS_DIR="/scratch/users/$USER/rfdiff_models/proteinmpnn"
CALIBY_WEIGHTS_DIR="/scratch/users/$USER/rfdiff_models/caliby/model_params/caliby"

# Output logs directory
LOGDIR="$WORKDIR/logs"

# ==================== Pipeline Arguments ====================

# Modify these based on your needs
# Input directory containing scaffold PDB files
SCAFFOLDS_DIR="$WORKDIR/prame_scaffolds_v1"  # Adjust path to your scaffolds directory

# Chain IDs
BINDER_CHAIN="A"
MHC_CHAIN="F"
PEPTIDE_CHAIN="H"

# Number of sequences to design per scaffold
NUM_CALIBY_SEQS=4
NUM_AF3_SEQS=4  # Add this if missing

# Template PDB (you may need to set this)
TEMPLATE_PDB=""  # Set this to your template PDB path if needed

ONTARGET_PEPTIDE_SEQ=""  # Optional: on-target peptide sequence
OFFTARGET_TEMPLATE_PDB=""  # Optional: off-target template PDB
OFFTARGET_PEPTIDE_SEQ="SLLMWITQV"                    # Off-target peptide sequence (optional)

# Caliby checkpoint (will auto-detect if not specified)
CALIBY_CKPT="$CALIBY_WEIGHTS_DIR/caliby.ckpt"

# Output directory
OUT_DIR="$WORKDIR/caliby_af3_outputs"

# ==================== Setup ====================

mkdir -p "$LOGDIR"

# Check if Singularity image exists (for Caliby)
if [ ! -f "$IMG" ]; then
    echo "ERROR: Singularity image not found at $IMG"
    echo "Please build it first using:"
    echo "  cd $WORKDIR/rfdiff_proteinmpnn_image"
    echo "  bash build_singularity_sherlock.sh"
    exit 1
fi

# Check if Caliby weights exist
if [ ! -f "$CALIBY_CKPT" ]; then
    echo "ERROR: Caliby checkpoint not found at $CALIBY_CKPT"
    echo "Please download Caliby weights first:"
    echo "  cd $WORKDIR/rfdiff_proteinmpnn_image"
    echo "  bash download_caliby_weights.sh /scratch/users/$USER/rfdiff_models"
    exit 1
fi

# Check if AlphaFold3 is available
if [ ! -f "$AF3_IMAGE" ]; then
    echo "ERROR: AlphaFold3 Singularity image not found at $AF3_IMAGE"
    echo "AlphaFold3 should be available at: /home/groups/possu/alphafold3/"
    exit 1
fi

# Check if scaffolds directory exists
if [ ! -d "$SCAFFOLDS_DIR" ]; then
    echo "ERROR: Scaffolds directory not found: $SCAFFOLDS_DIR"
    exit 1
fi

# ==================== Fix: Install hydra-core ====================
# The Caliby environment is missing hydra-core. Install it before running.

echo "========================================="
echo "Installing hydra-core in Caliby environment..."
echo "========================================="

singularity exec --nv \
    --bind "$WORKDIR:/workspace" \
    --bind "$RFDIFF_WEIGHTS_DIR:/models/rfdiffusion" \
    --bind "$MPNN_WEIGHTS_DIR:/models/proteinmpnn" \
    --bind "$CALIBY_WEIGHTS_DIR:/opt/caliby/model_params/caliby" \
    "$IMG" \
    bash -c "/opt/caliby/envs/caliby/bin/pip install hydra-core --quiet"

if [ $? -ne 0 ]; then
    echo "WARNING: Failed to install hydra-core. Trying with conda instead..."
    singularity exec --nv \
        --bind "$WORKDIR:/workspace" \
        --bind "$RFDIFF_WEIGHTS_DIR:/models/rfdiffusion" \
        --bind "$MPNN_WEIGHTS_DIR:/models/proteinmpnn" \
        --bind "$CALIBY_WEIGHTS_DIR:/opt/caliby/model_params/caliby" \
        "$IMG" \
        bash -c "source /opt/caliby/envs/caliby/etc/profile.d/conda.sh && conda activate caliby && pip install hydra-core --quiet"
fi

echo "Hydra-core installation complete."
echo ""

# ==================== Environment Variables ====================

export CUDA_VISIBLE_DEVICES=0
export CALIBY_DIR=/opt/caliby
export CALIBY_ENV=/opt/caliby/envs/caliby
export CALIBY_MODEL_DIR="/scratch/users/$USER/rfdiff_models/caliby"

# ==================== Run Pipeline ====================

echo "========================================="
echo "Caliby-AlphaFold3 Pipeline"
echo "========================================="
echo "Scaffolds directory: $SCAFFOLDS_DIR"
echo "Output directory: $OUT_DIR"
echo "Template PDB: $TEMPLATE_PDB"
echo "Caliby checkpoint: $CALIBY_CKPT"
echo "AlphaFold3 image: $AF3_IMAGE"
echo "========================================="
echo ""

# Build command
CMD="python /workspace/run_caliby_af3_pipeline.py \
    --scaffolds_dir $SCAFFOLDS_DIR \
    --binder_chain $BINDER_CHAIN \
    --mhc_chain $MHC_CHAIN \
    --peptide_chain $PEPTIDE_CHAIN \
    --num_caliby_seqs $NUM_CALIBY_SEQS \
    --num_af3_seqs $NUM_AF3_SEQS \
    --template_pdb $TEMPLATE_PDB \
    --caliby_ckpt $CALIBY_CKPT \
    --af3_script $AF3_SCRIPT \
    --af3_image $AF3_IMAGE \
    --af3_model_dir $AF3_MODEL_DIR \
    --af3_db_dir $AF3_DB_DIR \
    --out_dir $OUT_DIR"

# Add optional on-target peptide sequence argument if provided
if [ -n "$ONTARGET_PEPTIDE_SEQ" ]; then
    CMD="$CMD --ontarget_peptide_sequence \"$ONTARGET_PEPTIDE_SEQ\""
fi

# Add optional off-target arguments if provided
if [ -n "$OFFTARGET_TEMPLATE_PDB" ] && [ -f "$OFFTARGET_TEMPLATE_PDB" ]; then
    CMD="$CMD --offtarget_template_pdb $OFFTARGET_TEMPLATE_PDB"
fi

if [ -n "$OFFTARGET_PEPTIDE_SEQ" ]; then
    CMD="$CMD --offtarget_peptide_sequence \"$OFFTARGET_PEPTIDE_SEQ\""
fi

# Run pipeline script inside Caliby container (for Caliby dependencies)
# The script will handle calling AF3 Singularity from within
# Note: Nested Singularity calls may require proper setup
singularity exec --nv \
    --bind "$WORKDIR:/workspace" \
    --bind "$RFDIFF_WEIGHTS_DIR:/models/rfdiffusion" \
    --bind "$MPNN_WEIGHTS_DIR:/models/proteinmpnn" \
    --bind "$CALIBY_WEIGHTS_DIR:/opt/caliby/model_params/caliby" \
    --bind "$AF3_SCRIPT:/opt/af3_script.sh:ro" \
    --bind "$AF3_IMAGE:/opt/af3_image.sif:ro" \
    --bind "$AF3_MODEL_DIR:/opt/af3_models:ro" \
    --bind "$AF3_DB_DIR:/opt/af3_databases:ro" \
    "$IMG" \
    bash -c "$CMD"

EXIT_CODE=$?

echo ""
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "Pipeline completed successfully!"
else
    echo "Pipeline failed with exit code $EXIT_CODE"
fi
echo "Output directory: $OUT_DIR"
echo "========================================="

exit $EXIT_CODE

