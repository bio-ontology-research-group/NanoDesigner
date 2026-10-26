#!/bin/bash
#
# NanoDesigner - de novo design for a user-supplied target.
#
# Unlike the scripts under NanoDesigner_assessment_experiment_2/, which read
# their input from the benchmark test.json splits, this one takes the JSON
# produced by prepare_NanoDesigner_inputs_Denovo.ipynb.
#
# Usage:
#   bash NanoDesigner_pipeline.sh PATH_TO_JSON [STARTING_ITERATION]
#   bash NanoDesigner_pipeline.sh CASE_NAME   [STARTING_ITERATION]
#
# Pass the JSON written by prepare_NanoDesigner_inputs_Denovo.ipynb directly.
# Its filename (without .json) becomes the case name, and the notebook already
# makes that match the "pdb" field inside the file, which is what the pipeline
# keys its outputs on. Do not rename the file.
#
# Example:
#   bash NanoDesigner_test_Cases/NanoDesigner_pipeline.sh \
#        NanoDesigner_test_Cases/denovo_epitope_info/6lr7_5lrt/6lr7_5lrt_ep_1.json
#
# Giving a bare case name instead looks for
# NanoDesigner_test_Cases/CASE_NAME/{dataset.json,CASE_NAME.json}.

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 PATH_TO_JSON [STARTING_ITERATION]"
    echo "       $0 CASE_NAME    [STARTING_ITERATION]"
    echo "If STARTING_ITERATION is not provided, defaults to 1"
    exit 1
fi

INPUT=$1
START_ITER=${2:-1}

# ============================================================================
#  EDIT THIS BLOCK
# ============================================================================

# Path to your NanoDesigner clone.
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# CDR generation parameters. Pick a de novo config:
#   NanoDesigner_diffab_denovo_1CDR.yml    (CDRH3 only)
#   NanoDesigner_diffab_denovo_3CDRs.yml   (all three CDRs)
#   NanoDesigner_adesigner_denovo_1CDR.yml (ADesigner instead of DiffAb)
# Remember to set the checkpoint path inside the chosen file.
CONFIG="${PROJECT_DIR}/config_files/NanoDesigner_diffab_denovo_1CDR.yml"

# Design parameters (paper settings).
R=50    # randomized nanobodies at initialization
N=15    # top mutants carried to the next iteration
d=100   # docked models generated per nanobody
n=5     # top docked models fed to inference
max_iter=10

# ============================================================================

DYMEAN_CODE_DIR="${PROJECT_DIR}/dyMEAN"
ADESIGNER_CODE_DIR="${PROJECT_DIR}/ADesigner"
DIFFAB_CODE_DIR="${PROJECT_DIR}/diffab"
MAIN_FOLDER="${PROJECT_DIR}/NanoDesigner_test_Cases"

# Resolve the input. A path to an existing .json is used as-is; anything else
# is treated as a case name and looked up under NanoDesigner_test_Cases/.
if [ -f "$INPUT" ]; then
    BASENAME=$(basename "$INPUT" .json)
    DATA_DIR="${MAIN_FOLDER}/${BASENAME}"
    mkdir -p "$DATA_DIR"
    # Keep a copy beside the results so each run records the input it used.
    if [ "$(cd "$(dirname "$INPUT")" && pwd)" != "$(cd "$DATA_DIR" && pwd)" ]; then
        cp "$INPUT" "${DATA_DIR}/${BASENAME}.json"
    fi
    DATASET="${DATA_DIR}/${BASENAME}.json"
else
    BASENAME="$INPUT"
    DATA_DIR="${MAIN_FOLDER}/${BASENAME}"
    if [ -f "${DATA_DIR}/dataset.json" ]; then
        DATASET="${DATA_DIR}/dataset.json"
    elif [ -f "${DATA_DIR}/${BASENAME}.json" ]; then
        DATASET="${DATA_DIR}/${BASENAME}.json"
    else
        echo "Error: '${INPUT}' is not a file, and no input JSON was found for it as a case name."
        echo "Looked for:"
        echo "  ${DATA_DIR}/dataset.json"
        echo "  ${DATA_DIR}/${BASENAME}.json"
        echo
        echo "Generate the JSON with NanoDesigner_test_Cases/prepare_NanoDesigner_inputs_Denovo.ipynb,"
        echo "then pass its path to this script."
        exit 1
    fi
fi

if [ ! -f "$CONFIG" ]; then
    echo "Error: config file not found at $CONFIG"
    exit 1
fi

echo "Case:    ${BASENAME}"
echo "Dataset: ${DATASET}"
echo "Config:  ${CONFIG}"

# Parameters read from the config file.
MODEL=$(sed -n 's/^model_name: //p' "$CONFIG")                      # DiffAb or ADesigner
CKPT=$(sed -n 's/^[ ]*checkpoint: //p' "$CONFIG")
MAX_OBJECTIVE=$(sed -n 's/^[ ]*maximization_objective: //p' "$CONFIG")  # dg for de novo
CDRS=$(sed -n 's/^[ ]*CDRS: //p' "$CONFIG")
INITIAL_CDR=$(sed -n 's/^[ ]*initial_cdr: //p' "$CONFIG")

SAVE_DIR="${DATA_DIR}/${MODEL}_${INITIAL_CDR}_cdr"
HDOCK_DIR="${SAVE_DIR}/HDOCK"
RESULT_DIR="${SAVE_DIR}/results"

eval "$(conda shell.bash hook)"
conda activate nanodesigner1

export DIFFAB_CODE_DIR
export DYMEAN_CODE_DIR
export ADESIGNER_CODE_DIR
export PYTHONPATH="${PROJECT_DIR}:${PYTHONPATH}"
export NANOBODIES_PROJECT_DIR="${PROJECT_DIR}"   # dr_sasa is resolved from this

for ((i=START_ITER; i<=max_iter; i++)); do

    echo "Iteration $i"
    echo "------Running structure prediction and docking (HDOCK) --------"

    SUMMARY_FILE_INFERENCE=${RESULT_DIR}_iteration_$i/summary_iter_${i}.json
    SUMMARY_FILE_PACKED=${RESULT_DIR}_iteration_$i/summary_packed_iter_${i}.json
    SUMMARY_FILE_PACKED_REFINED=${RESULT_DIR}_iteration_$i/summary_packed_refined_iter_${i}.json

    cd ${DYMEAN_CODE_DIR}
    start_time=$SECONDS
    python -m models.pipeline.NanoDesigner_structure_pred_and_docking \
        --dataset_json ${DATASET} \
        --randomized $R \
        --best_mutants $N \
        --cdr_type ${CDRS} \
        --cdr_model ${MODEL} \
        --hdock_models ${HDOCK_DIR}_iter_$i \
        --n_docked_models $d \
        --iteration $i \
        --initial_cdr ${INITIAL_CDR} \
        --csv_dir ${DATA_DIR}/csv_iter_$i \
        --csv_dir_ ${DATA_DIR}/csv_iter_

    wait
    end_time=$SECONDS
    echo "Structure prediction and docking took approximately: $((end_time - start_time)) seconds"

    echo "------ Processing and selection of top docked models --------"
    start_time=$SECONDS

    python -m models.pipeline.NanoDesigner_select_top_docked_models \
        --test_set ${DATASET} \
        --hdock_models ${HDOCK_DIR}_iter_$i \
        --iteration $i \
        --top_n $n

    end_time=$SECONDS
    echo "Selection of top docked models took approximately: $((end_time - start_time)) seconds"

    echo "------Refinement and filtering of docked models --------"
    start_time=$SECONDS

    check_and_retry_stuck_docking() {
        echo "Checking for stuck refinement items..."
        find ${HDOCK_DIR}_iter_$i -name ".*.lock" -mmin +10 -delete 2>/dev/null || true
        python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_refinement_after_docking.py \
            --hdock_models ${HDOCK_DIR}_iter_$i \
            --iteration $i \
            --inference_summary ${SUMMARY_FILE_INFERENCE} \
            --dataset_json ${DATASET} \
            --top_n $n
    }

    find ${HDOCK_DIR}_iter_$i -name ".*.lock" -delete 2>/dev/null || true

    for instance in {1..5}; do
        echo "Starting refinement instance $instance"
        python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_refinement_after_docking.py \
            --hdock_models ${HDOCK_DIR}_iter_$i \
            --iteration $i \
            --inference_summary ${SUMMARY_FILE_INFERENCE} \
            --dataset_json ${DATASET} \
            --top_n $n &
        sleep 3
    done

    wait
    check_and_retry_stuck_docking
    find ${HDOCK_DIR}_iter_$i -name ".*.lock" -delete 2>/dev/null || true

    echo "Refinement completed in $((SECONDS - start_time)) seconds"

    echo "-----------Inference-----------"
    start_time=$SECONDS

    mkdir -p ${RESULT_DIR}_iteration_$i

    if [ "$MODEL" = "DiffAb" ]; then

        conda activate nanodesigner2

        python ${DYMEAN_CODE_DIR}/models/pipeline/cdr_models/NanoDesigner_diffab_inference_original.py \
            --dataset ${DATASET} \
            --config ${CONFIG} \
            --out_dir ${RESULT_DIR}_iteration_$i \
            --hdock_models ${HDOCK_DIR}_iter_$i \
            --diffab_code_dir ${DIFFAB_CODE_DIR} \
            --dymean_code_dir ${DYMEAN_CODE_DIR} \
            --summary_dir ${SUMMARY_FILE_INFERENCE} \
            --gpu 0 \
            --iteration $i

        conda deactivate

    elif [ "$MODEL" = "ADesigner" ]; then

        conda activate nanodesigner1

        cd $ADESIGNER_CODE_DIR
        python ${ADESIGNER_CODE_DIR}/generate_pipeline.py \
            --ckpt ${CKPT} \
            --test_set ${DATASET} \
            --out_dir ${RESULT_DIR}_iteration_$i \
            --hdock_models ${HDOCK_DIR}_iter_$i \
            --rabd_topk 5 \
            --mode "1*1" \
            --rabd_sample 20 \
            --config ${CONFIG} \
            --iteration $i

    else
        echo "Error: Unknown model type '${MODEL}'. Please choose DiffAb or ADesigner. Revise the config file."
        exit 1
    fi

    wait
    end_time=$SECONDS
    echo "Inference took approximately: $((end_time - start_time)) seconds"

    cd ${DYMEAN_CODE_DIR}

    echo "--------- Conducting side chain packing ----------------"
    start_time=$SECONDS

    python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_side_chain_packing.py \
        --summary_json ${SUMMARY_FILE_INFERENCE} \
        --out_file ${SUMMARY_FILE_PACKED} \
        --test_set ${DATASET} \
        --cdr_model ${MODEL}
    wait

    end_time=$SECONDS
    echo "Side chain packing took approximately: $((end_time - start_time)) seconds"

    echo "------Refinement of packed complexes --------"
    start_time=$SECONDS

    # Lock files live beside the summary being written.
    LOCK_DIR=$(dirname "${SUMMARY_FILE_PACKED_REFINED}")

    check_and_retry_stuck_packing() {
        echo "Checking for stuck items..."
        find "${LOCK_DIR}/.locks" -name "*.lock" -mmin +10 -delete 2>/dev/null || true
        python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_refinement_after_side_chain_packing_parallel.py \
            --in_file ${SUMMARY_FILE_PACKED} \
            --out_file ${SUMMARY_FILE_PACKED_REFINED} \
            --cdr_model ${MODEL}
    }

    find "${LOCK_DIR}" -name ".*.lock" -delete 2>/dev/null || true
    find "${LOCK_DIR}/.locks" -name "*.lock" -delete 2>/dev/null || true

    for instance in {1..5}; do
        echo "Starting side chain packing instance $instance"
        python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_refinement_after_side_chain_packing_parallel.py \
            --in_file ${SUMMARY_FILE_PACKED} \
            --out_file ${SUMMARY_FILE_PACKED_REFINED} \
            --cdr_model ${MODEL} &
        sleep 3
    done

    wait
    check_and_retry_stuck_packing

    find "${LOCK_DIR}" -name ".*.lock" -delete 2>/dev/null || true
    find "${LOCK_DIR}/.locks" -name "*.lock" -delete 2>/dev/null || true

    echo "Side chain packing refinement completed in $((SECONDS - start_time)) seconds"

    echo "----------Iteration evaluation-----------"
    start_time=$SECONDS

    mkdir -p ${DATA_DIR}/csv_iter_$i

    python ${DYMEAN_CODE_DIR}/NanoDesigner_metrics_computations.py \
        --test_set ${DATASET} \
        --summary_json ${SUMMARY_FILE_PACKED_REFINED} \
        --hdock_models ${HDOCK_DIR}_iter_$i \
        --cdr_type ${CDRS} \
        --iteration $i \
        --cdr_model ${MODEL} \
        --csv_dir ${DATA_DIR}/csv_iter_$i

    wait
    end_time=$SECONDS
    echo "Evaluation took approximately: $((end_time - start_time)) seconds"

    echo "----------Best mutants selection-----------"
    python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_best_mutant_selection.py \
        --dataset_json ${DATASET} \
        --top_n $N \
        --hdock_models ${HDOCK_DIR}_iter_$i \
        --iteration $i \
        --csv_dir ${DATA_DIR}/csv_iter_$i \
        --objective ${MAX_OBJECTIVE}

    wait

done

echo "Done. Results are under ${SAVE_DIR}"