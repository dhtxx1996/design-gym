#!/bin/bash
# Benchmark script for pH-sensitive protein design task
# Usage:
#   ./bench.sh [difficulty] [model] [run_name]   - Run single benchmark
#   ./bench.sh all [max_parallel]                - Run all combinations in parallel

cd "$(dirname "$0")"
TASK="ph_sensitive_design"

run_single() {
    local diff="$1" model="$2" run="$3"
    case $diff in
        easy|medium|hard) question="question_${diff}.md" ;;
        *) question="$diff" ;;
    esac
    python agent.py --task "$TASK" --question "$question" --model "$model" --output "$run" --overwrite
}

if [ "$1" = "all" ]; then
    MODELS=("gpt-4o" "o3" "gpt-5")
    DIFFS=("easy" "medium" "hard")
    MAX_JOBS="${2:-4}"
    SUBMISSION_ID="$(date +%m%d%H%M%S)"
    
    echo "Running $((${#MODELS[@]} * ${#DIFFS[@]})) benchmarks (max $MAX_JOBS parallel)"
    echo "Submission ID: $SUBMISSION_ID"
    
    for model in "${MODELS[@]}"; do
        for diff in "${DIFFS[@]}"; do
            while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do sleep 1; done
            
            run="${model}_${diff}_${SUBMISSION_ID}"
            logdir="$TASK/agent_output/$run"
            mkdir -p "$logdir"
            
            echo "[START] $run"
            (run_single "$diff" "$model" "$run" > "$logdir/bench.log" 2>&1 && \
                echo "[DONE]  $run" || echo "[FAIL]  $run (see $logdir/bench.log)") &
            jobs=$((jobs + 1))
        done
    done
    wait
    echo "All benchmarks complete!"
else
    DIFF="${1:-medium}"
    MODEL="${2:-gpt-4o}"
    RUN="${3:-${MODEL}_${DIFF}_$(date +%m%d%H%M)}"
    
    echo "=== $TASK: $DIFF / $MODEL -> $RUN ==="
    run_single "$DIFF" "$MODEL" "$RUN"
    echo "Done! Results in: $TASK/agent_output/$RUN/"
fi
