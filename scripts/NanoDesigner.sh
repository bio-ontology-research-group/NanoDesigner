#!/bin/zsh
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 path_to_json_file"
    echo "Make sure is named as XXXX.json (if nanobody-antigen complex) or XXXX_XXXX.json (if nanobody and antigen comes from different PDB sources (ID: XXXX)." 
    exit 1
fi

JSON_FILE=$1
# Extract the filename without directory and extension
BASENAME=$(basename "$JSON_FILE" .json)


# Define variables
DYMEAN_CODE_DIR=./NanoDesigner/dyMEAN
cd ${DYMEAN_CODE_DIR}
DIFFAB_CODE_DIR=./NanoDesigner/diffab
ADESIGNER_CODE_DIR=./NanoDesigner/ADesigner
CONFIG=./NanoDesigner/config_files/NanoDesigner_ADesigner_codesign_single.yml
MAIN_FOLDER=./NanoDesigner/your_working_directory/ADesigner


#NanoDesigner variables
R=50 # Number of randomized nanobodies (Initialization step)
N=15 # Top best mutants to proceed with to subsequent iterations
d=100 # docked models to generate
n=5 # top docked models to feed to inference stage
max_iter=10


# Extract variables from config file
MODEL=$(sed -n 's/^model_name: //p' "$CONFIG") # options: ADesigner or DiffAb
CKPT=$(sed -n 's/^[ ]*checkpoint: //p' "$CONFIG")
MAX_OBJECTIVE=$(sed -n 's/^[ ]*maximization_objective: //p' "$CONFIG") # options: dg for de novo; ddg for optimization design
CDRS=$(sed -n 's/^[ ]*CDRS: //p' "$CONFIG")
INITIAL_CDR=$(sed -n 's/^[ ]*initial_cdr: //p' "$CONFIG")
DATA_DIR=${MAIN_FOLDER}/${BASENAME}
SAVE_DIR=${DATA_DIR}/${MODEL}_${INITIAL_CDR}_cdr
HDOCK_DIR=${SAVE_DIR}/HDOCK
RESULT_DIR=${SAVE_DIR}/results

mkdir -p $DATA_DIR


#Copy the json file in it (to conserve the structure)
# Update the DATASET variable with the new location
cp "$JSON_FILE" "$DATA_DIR"
DATASET="$DATA_DIR/$(basename "$JSON_FILE")"

# Source .bashrc to ensure Conda is initialized
source ~/miniconda3/etc/profile.d/conda.sh
conda activate nanodesigner1

# source ~/.bashrc
# conda activate nanodesigner2

for ((i=1; i<=max_iter; i++)); do

    echo "Iteration $i"
    echo "------Running docking simulation (Hdock) --------"

    SUMMARY_FILE_INFERENCE=${RESULT_DIR}_iteration_$i/summary_iter_${i}.json
    SUMMARY_FILE_PACKED=${RESULT_DIR}_iteration_$i/summary_packed_iter_${i}.json
    SUMMARY_FILE_PACKED_REFINED=${RESULT_DIR}_iteration_$i/summary_packed_refined_iter_${i}.json

    start_time=$SECONDS
    python -m models.pipeline.NanoDesigner_docking_simulations \
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
    execution_time=$(($end_time - $start_time))
    echo "Docking Simulation took approximately: $execution_time seconds."



    echo "------ Processing and selection of top docked models --------"
    start_time=$SECONDS

    if [ "$i" -eq 1 ]; then
        F=$R
    else
        F=$N
    fi

    k=$(($F+3))  # Added 3 extra rounds to ensure correct processing (in case dSASA computation fails)
    for ((a=1; a<=$k; a++)); do
        python -m models.pipeline.NanoDesigner_select_top_hdock_models \
            --test_set ${DATASET} \
            --hdock_models ${HDOCK_DIR}_iter_$i \
            --iteration $i \
            --top_n $n

    done

    
    wait
    end_time=$SECONDS
    execution_time=$((end_time - start_time))
    echo "Selection of top docked models took approximately: $execution_time seconds"


    echo "------Refinement and filtering of docked models --------"
    start_time=$SECONDS

    # Create file with tasks to conduct with GREASY
    python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_refinement_hdock_models.py \
        --hdock_models ${HDOCK_DIR}_iter_$i \
        --iteration $i \
        --inference_summary ${SUMMARY_FILE_INFERENCE} \
        --dataset_json ${DATASET} \
        --top_n $n

    wait
    end_time=$SECONDS
    execution_time=$((end_time - start_time))
    echo "Refinement and filtering of docked models took approximately: $execution_time seconds"


    echo "-----------Inference-----------"
    start_time=$SECONDS


    if [ "$MODEL" = "DiffAb" ]; then

        conda activate nanodesigner2

        echo "dataset ${DATASET}"
        echo "config ${CONFIG}"
        echo "config ${RESULT_DIR}_iteration_$i"

        mkdir -p ${RESULT_DIR}_iteration_$i

        python ${DYMEAN_CODE_DIR}/models/pipeline/cdr_models/Diffab_for_NanoDesigner.py \
            --dataset ${DATASET} \
            --config ${CONFIG} \
            --out_dir ${RESULT_DIR}_iteration_$i \
            --hdock_models ${HDOCK_DIR}_iter_$i \
            --diffab_code_dir ${DIFFAB_CODE_DIR}  \
            --dymean_code_dir ${DYMEAN_CODE_DIR}  \
            --iteration $i 

        conda deactivate
        wait
        end_time=$SECONDS
        execution_time=$((end_time - start_time))
        echo "Inference took approximately: $execution_time seconds"

    elif [ "$MODEL" = "ADesigner" ]; then

        conda activate nanodesigner1

        cd $ADESIGNER_CODE_DIR
        python ${ADESIGNER_CODE_DIR}/generate_pipeline.py \
            --ckpt ${CKPT} \
            --test_set ${DATASET} \
            --out_dir ${RESULT_DIR}_iteration_$i \
            --hdock_models ${HDOCK_DIR}_iter_$i \
            --rabd_topk 5  \
            --mode "1*1"  \
            --rabd_sample 20 \
            --config ${CONFIG} \
            --iteration $i 


    else
        echo "Error: Unknown model type '${MODEL}'. Please choose from DiffAb or ADesigner. Revise config file"
        exit 1
    fi


    echo --------- Conducting Side chain Packing ----------------

    cd ${DYMEAN_CODE_DIR}

    start_time=$SECONDS
    python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_sidechainpacking.py \
        --summary_json ${SUMMARY_FILE_INFERENCE} \
        --out_file ${SUMMARY_FILE_PACKED} \
        --test_set ${DATASET} \
        --cdr_model ${MODEL}
    wait

    end_time=$SECONDS
    execution_time=$((end_time - start_time))
    echo "Side chain Packing after Side Chain Packing took approximately: $execution_time seconds"
    

    echo --------- Refinement of Packed complexes ---------------

    start_time=$SECONDS

    python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_refinement_spacked_complexes.py \
        --cdr_model ${MODEL} \
        --in_file "${SUMMARY_FILE_PACKED}" \
        --out_file "${SUMMARY_FILE_PACKED_REFINED}"

    wait

    end_time=$SECONDS
    execution_time=$((end_time - start_time))
    echo "Refinement after Side Chain Packing took approximately: $execution_time seconds"


    echo "----------Iter Evaluation-----------"
    python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_evaluation.py \
        --test_set ${DATASET} \
        --summary_json ${SUMMARY_FILE_PACKED_REFINED} \
        --hdock_models ${HDOCK_DIR}_iter_$i \
        --cdr_type ${CDRS} \
        --iteration $i \
        --cdr_model ${MODEL} \
        --csv_dir ${DATA_DIR}/csv_iter_$i
    
    wait
    execution_time=$((end_time - start_time))
    end_time=$SECONDS
    echo "Evaluation took approximately: $execution_time seconds"

    echo "----------Best Mutants Selection-----------"
    python ${DYMEAN_CODE_DIR}/models/pipeline/NanoDesigner_best_mutant_selection.py \
        --dataset_json ${DATASET} \
        --top_n $N \
        --hdock_models ${HDOCK_DIR}_iter_$i \
        --iteration $i \
        --csv_dir ${DATA_DIR}/csv_iter_$i \
        --objective $MAX_OBJECTIVE 

    wait  

done
