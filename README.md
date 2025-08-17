# NanoDesigner: Resolving the Complex-CDR Interdependency with Iterative Refinement

![Alt text](https://github.com/Melissaurious/NanoDesigner/blob/main/NanoDesigner.png)

NanoDesigner is an end-to-end workflow for the design and optimization of nanobodies. It integrates key stages—Structure Prediction, Docking, CDR Generation,and Side-Chain Packing—into an iterative framework based on an Expectation Maximization algorithm. Our method effectively tackles an often overlooked interdependency
challenge where accurate docking presupposes a priori knowledge of the CDR conformation, while effective CDR generation relies on accurate docking outputs to guide its design.


**Journal of Cheminformatics** 17, 120 (2025).  
📄 [Read the paper here](https://doi.org/10.1186/s13321-025-01069-2) | [Journal link](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-025-01069-2)


## Table of Contents
- [Installation](#installation)
- [External Tool Installation](#external-tool-installation)
- [Data Download and Preprocess](#data-download-and-preprocess)
- [Inference-tool assessment](#inference)
- [NanoDesigner](#nanodesigner)
- [Citation](#citation)
- [License](#license)
- [Credits](#credits)



## Installation


```bash
git clone https://github.com/Melissaurious/NanoDesigner.git
cd NanoDesigner
```

### Create and activate the Conda environment for graph-based models
```bash
conda env create -f nanodesigner_1.yml -n nanodesigner1
conda activate nanodesigner1
```

### Create and activate the Conda environment for difussion-based model
```bash
conda env create -f nanodesigner_2.yml -n nanodesigner2
conda activate nanodesigner2
```

## External Tool Installation

The following repositories and software tools are required for NanoDesigner. Clone them into the `NanoDesigner` directory and follow the installation steps provided in their respective webpages:

- [IgFold](https://github.com/Graylab/IgFold) - *Trained models already included in nanodesigner1 conda environment.*
- [DockQ](https://github.com/bjornwallner/DockQ)
- [Rosetta](https://docs.rosettacommons.org/demos/latest/tutorials/install_build/install_build)
- [FoldX](https://foldxsuite.crg.eu/products#foldx) - *Academic license employed.*
- [HDOCK](http://huanglab.phys.hust.edu.cn/software/hdocklite/)
- [dr_sasa_n](https://github.com/nioroso-x3/dr_sasa_n) - *Follow the instructions in the repository to compile this tool.*

After installing the tools, ensure to update the `dyMEAN/configs.py` file with the full paths to the installed tools.

Source code for TMscore evatuation is at `dyMEAN/evaluation/`, please compile as:
```bash
cd dyMEAN
g++ -static -O3 -ffast-math -lm -o evaluation/TMscore evaluation/TMscore.cpp
```


## Data Download and Preprocess

The data download and preparation steps are necessary to replicate our data processing, filtering, and preparation for training and inference. All required instructions are included in the following scripts.

### 1. Preprocess the Data
Use the `download_and_process_data_reduced.py` script to download SAbDab structures, extract CDR information, and analyze antibody-antigen interactions based on solvent-accessible surface area (dSASA).
The preprocessing pipeline performs structure validation, CDR extraction, interaction analysis, and epitope mapping. It automatically handles data download, quality filtering, and generates the final dataset with mapped epitopes from CDR-antigen contacts.

#### Usage
```bash
python download_and_process_data_reduced.py \
    --output_folder <output_directory> \
    --type <Antibody|Nanobody> \
    --numbering <imgt|chothia> \
```

#### Parameters
* `--output_folder`: Directory for all output files (required). Please provide as input ./NanoDesigner/Data_download_and_processing
* `--type`: Structure type to analyze - `Antibody` or `Nanobody` (default: `Nanobody`)
* `--numbering`: Numbering scheme - `imgt` or `chothia` (default: `imgt`)
* `--max_resolution`: Maximum resolution threshold in Å (default: `4.0`)
* `--tsv_file`: Custom SAbDab TSV file (auto-downloaded if not provided). We provide: ./NanoDesigner/Data_download_and_processing/sabdab_summary_all-5.tsv
* `--raw_structures_dir`: Custom PDB structures directory (auto-downloaded if not provided)

#### Output
The main result is `CDRH3_interacting_[type]_[numbering]_unique.json` containing the final filtered dataset with CDR sequences, interaction data, and epitope mappings (input for the next step).

_Nanobody Dataset Analisis: Length and binding involvement of nanobody CDRHs_
<div style="text-align: center;">
  <img src="https://github.com/Melissaurious/NanoDesigner/blob/main/combined_CDR_length_binding_paper2.png" alt="Alt Text" width="700">
</div>


### 2. Split the Data
After downloading and preprocessing the data, use the split_data.py script to create train/validation/test splits with sequence-based clustering.
The splitting script performs sequence similarity clustering on either CDRH3 or antigen sequences, with [40%,30%,20%]  and [95%,80%,60%] threshoolds, respectively. Then creates k-fold cross-validation splits ensuring that similar sequences are kept within the same split. This prevents data leakage where the model might see similar sequences during training and testing.
#### Usage
```bash
python ./Data_download_and_processing/scripts/split_data.py \
    --data_files <path_to_json_file(s)> \
    --immuno_molecule <Antibody|Nanobody|Antibody Nanobody> \
    --cluster_targets <CDRH3|Ag> \
    --out_dir <output_directory> \
```


#### Parameters
* `--data_files`: Nanobody or/and Antibody CDRH3_interacting_[type]_imgt_unique.json (required)
* `--immuno_molecule`: Molecule file to use to construct the clusters, if both files are provided the content will be pooled - `Antibody` or `Nanobody` or both.
* `--cluster_targets`: Type of sequence to use for clustering based on sequence similarity
* `--out_dir`: Directory to save subfodlers with generated cluster per training configuration.


## Inference-tool assessment

We provide trained models for dyMEAN, ADesigner and DiffAb, one or multiple folds for evaluating different training configurations in this link [here](https://drive.google.com/drive/folders/1SZUP4ovqYtHjxIQSJ4-UoO60wW3YN-lj?usp=share_link). Please create a directory called all_checkpoints inside ./NanoDesigner and unzip each file (e.g tar -xzvf ./NanoDesigner/all_checkpoints/DIFFAB.tar.gz). The tested training configurations used in this study are:

**Training Data Variants:**
- Nanobodies-only
- Nanobodies + Antibodies
- Antibodies fine-tuned on Nanobodies

**Clustering Thresholds and sequences:**
- Antigen (Ag) sequence clustering (95%, 80%, 60%)
- CDRH3 sequence clustering (40%, 30%, 20%)

The tested configurations systematically assess model performance across different training data compositions, structural diversity levels (via CDRH3 loop clustering thresholds), and antigenic similarity (through varying antigen sequence clustering thresholds).

### Running Assessments
To replicate our generate evaluation results run:
```bash
./NanoDesigner/scripts/CDRH3_model_assessment.sh
```

#### Variables to update
* `MODEL`: Select the model to run the inference and evaluation for. Options: "dyMEAN" or "ADESIGN" or "DIFFAB".
* `BASE_DIR`: Your working directory (`./NanoDesigner/Tool_assesment_experiment_1`).
* `FINETUNE`: Set to False if you want to evalaute either Nanobody; Nanobody_Antibody; set to TRUE to automatically search for the fine-tuned models.
* `TOTAL_FOLDS`: For Supplementary material tables S2 and S2 set to 0, for Main Text Table 1, set to 9 (thresholds 60% for Ag, 40% for CDRH#) .
* `CLUSTER_TYPE`: Select "Ag" or "CDRH3" depending on the sequence to cluster from.
* `CLUSTER_THRESHOLD`: Set 95, 80 or 60 for Ag and 40, 30 or 20 for CDRH3.


### Checkpoints
The `./NanoDesigner/scripts/CDRH3_model_assessment.sh` will autoamtically locate the trained checkpoints for each tool and training configuration as long as the directory structure is respected as in [here](https://drive.google.com/drive/folders/1SZUP4ovqYtHjxIQSJ4-UoO60wW3YN-lj?usp=share_link)

```
all_checkpoints/
├── ADESIGN/
│   ├── Antibody_clustered_Ag_60/
│   │   └── fold_0/
│   │       └── best.ckpt
│   ├── Antibody_clustered_Ag_80/
│   ├── Antibody_clustered_Ag_95/
│   ├── Antibody_clustered_CDRH3_20/
│   ├── Antibody_clustered_CDRH3_30/
│   ├── Antibody_clustered_CDRH3_40/
│   ├── Nanobody_Antibody_clustered_Ag_60/
│   ├── Nanobody_Antibody_clustered_Ag_60_3CDRs/
│   ├── Nanobody_Antibody_clustered_Ag_80/
│   ├── Nanobody_Antibody_clustered_Ag_95/
│   ├── Nanobody_Antibody_clustered_CDRH3_20/
│   ├── Nanobody_Antibody_clustered_CDRH3_30/
│   ├── Nanobody_Antibody_clustered_CDRH3_40/
│   ├── Nanobody_clustered_Ag_60/
│   ├── Nanobody_clustered_Ag_80/
│   ├── Nanobody_clustered_Ag_95/
│   ├── Nanobody_clustered_CDRH3_20/
│   ├── Nanobody_clustered_CDRH3_30/
│   └── Nanobody_clustered_CDRH3_40/
├── DIFFAB/
│   ├── Antibody_clustered_Ag_60/
│   │   └── fold_0/
│   │       └── 200000.pt
│   └── ... (same subfolders as ADESIGN)
└── dyMEAN/
    ├── Antibody_clustered_Ag_60/
    │   └── fold_0/
    │       └── version_0/
    │           └── checkpoint/
    │               └── epoch144_step49590.ckpt
    └── ... (same subfolders as ADESIGN)

```

### Results Analysis
The performance measures presented in:
- Supplementary Material Tables 2 & 3
- Main Text Table 1

were generated using `NanoDesigner/Tool_assesment_results/generate_tables_experiment_1.ipynb`. This notebook processes the assessment outputs to produce the comparative tables in our publication.



## NanoDesigner

NanoDesigner is an end-to-end workflow designed for both **de novo** and **optimization** cases in nanobody-antigen complex design. The workflow script is located in each of the folders located at `./NanoDesigner/NanoDesigner_assessment_experiment_2`. 

NanoDesigner supports two distinct design approaches: **de novo** and **optimization**. All required input data is generated from the [Data Download and Preprocess](#data-download-and-preprocess)) stage:


- **De Novo Design**: This mode is used when no pre-existing nanobody-antigen complex 3D structure  is available. In this scenario, the design process assumes there is no reference complex, and the design process is guided by the maximization objective ΔG (binding free energy). To assess the effectiveness of NanoDesigner in this scenario we make use of original complexes to evaluate improvement in binding energy compared to a reference complex (ΔΔG) and compute success rate defined as the ratio of designs showing improved affinity values.


- **Optimization Cases**: Optimization mode is applied when you want to improve the binding affinity of an existing nanobody based on a reference nanobody-antigen complex. This approach uses a known complex structure as the starting point, where the design and selection process is based on the maximization objective ΔΔG.

To run evaluate NanoDesigner in any of these modes, please refer to `/NanoDesigner/NanoDesigner_assessment_experiment_2/` folder and update inside the bash script the following variables. Example:

```bash
bash ./NanoDesigner/NanoDesigner_assessment_experiment_2/NanoDesigner_DiffAb_optimization/NanoDesigner_pipeline.sh PDB_ID
```
#### Variables to update
* `BASE_DIR`: Your working directory
#### Files to update
* `CONFIG`: Files containing working parameters for the CDR generation at `./NanoDesigner/config_files`. Please also update the trained model absolute path.


By given the PDB_ID as input the script will automatically look for the the input information required for our method (a json file containing all produced information during [Data Download and Preprocess](#data-download-and-preprocess)). The script also accepts a second input in case the process stops and user desire to start from the last reached iteration. 


We highly encourage to keep a constant number of total number of designs across iterations for simplicity:

```python
R = 50  # Number of randomized nanobodies (Initialization step)
N = 15  # Top best mutants to proceed with to subsequent iterations
d = 100 # Docked models to generate with Hdock
n = 5   # Top docked models to feed to inference stage
k_iteration_1 = 3   # Inference designs at iteration 1 (pre-defined inside config files)
k_iteration_x = 10  # Inference designs at iteration x (pre-defined inside config files)

Rxnxk = 750 (Iteration 1)
Nxnxk = 750 (Iteration X)
```

#### Checkpoints
For proof of concepts of NanoDesigner, please download and employ DiffAb or ADesigner trained models found [here](https://drive.google.com/drive/folders/1kGK3rV138lG8vQpGAtHv5oNP_a11Gr01?usp=share_link). These comes from fold 0 of the best perfoming training configuration found in experiment 1 results: Comnbined dataset Nanobodies + Antibodies, clustered based on antigen sequence similarity with a threshold of 60%.


### Results Analysis
The performance measures presented in:
- Main Text Table 2

were generated using `./NanoDesigner/NanoDesigner_assessment_experiment_2/generate_tables_experiment_2.ipynb`. This notebook processes the assessment outputs to produce the comparative tables in our publication.

### Test Cases


In some cases researchers may want to design a nanobody targeting an antigen of therapeutic interest for which a complex may not be available. An input preparation step is required prior running NanoDesigner. In these cases, the user can provide the desired VHH scaffold and target antigen crystal stucture sources, as well as the desired target epitope (as a sequence; the program also accept non-continues sequences e.g. "YKLV;CLL"). Chain isolation is handled for cases in which the VHH scaffold or antigen are in bound state with other proteins.

In our manuscript we ran Nanodesigner for the de novo design of VHH naobodies targeting different antigens of interest:

**For 6LR7-4OBE and 6LR7-5LRT**, please follow the instructions and run `./NanoDesigner/NanoDesigner_test_Cases/NanoDesigner_test_cases_input_prep.ipynb` notebook, it will ouput one json file per provided epitope sequence for which we recommend to place inside `./NanoDesigner/NanoDesigner_test_Cases/6lr7_5lrt` renamed as dataset.json. You can now proceed to run the pipeline as: 

```bash
bash ./NanoDesigner/NanoDesigner_test_Cases/NanoDesigner_pipeline.sh 6lr7_5lrt
```

**For 7EOW-8PW**H, please follow the instructions and run `./NanoDesigner/NanoDesigner_test_Cases/NanoDesigner_test_cases_input_prep_HER2.ipynb` notebook, it will ouput one json file per provided epitope sequence for which we recommend to place inside `./NanoDesigner/NanoDesigner_test_Cases/7eow_8pwh` renamed as dataset.json. In this special case, the epitope is extracted analyzing the binding interface of HER2 protein bound to two monoclonal therapeutic antibodies prior chain isolation: PERTUZUMAB and TRASTUZUMAB. You can now proceed to run the pipeline as: 

```bash
bash ./NanoDesigner/NanoDesigner_test_Cases/NanoDesigner_pipeline.sh 7eow_8pwh
```


#### Variables to update
* `BASE_DIR`: Your working directory
#### Files to update
* `CONFIG`: Fiiles containing working parameters for the CDR generation at `./NanoDesigner/config_files`. Choose from the de novo alternatives.


### NanoDesigner test cases:
*De novo design escenario:
CDRH3 or 3CDRs design with ΔG optimization objective.

<div style="text-align: center;">
  <img src="https://github.com/Melissaurious/NanoDesigner/blob/main/nanodesigner_test_cases_2025.png" alt="Alt Text" width="700">
</div>



## Citation

Rios Zertuche, M.M., Kafkas, Ş., Renn, D. et al. Nanodesigner: resolving the complex-CDR interdependency with iterative refinement. J Cheminform 17, 120 (2025). https://doi.org/10.1186/s13321-025-01069-2

## License
TODO

## Credits

This codebase is primarily based on the following deep learning tools. We thank the authors for their contributions:

- [Diffab](https://github.com/luost26/diffab)  
- [dyMEAN](https://github.com/THUNLP-MT/dyMEAN)  
- [ADesigner](https://github.com/A4Bio/ADesigner) 

We also acknowledge the rest of tools and software that played a crucial role in the workflow employed in this study.
We sincerely thank the authors of these tools for their invaluable work, which made this project possible.





