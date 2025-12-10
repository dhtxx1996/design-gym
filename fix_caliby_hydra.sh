#!/bin/bash
# Fix script to install hydra-core in Caliby environment
# Add this to your sbatch script before running the pipeline

# Inside the Singularity container, install hydra-core in the Caliby environment
singularity exec --nv \
    --bind "$WORKDIR:/workspace" \
    --bind "$RFDIFF_WEIGHTS_DIR:/models/rfdiffusion" \
    --bind "$MPNN_WEIGHTS_DIR:/models/proteinmpnn" \
    --bind "$CALIBY_WEIGHTS_DIR:/opt/caliby/model_params/caliby" \
    "$IMG" \
    bash -c "/opt/caliby/envs/caliby/bin/pip install hydra-core"

