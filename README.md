# NanoDesigner: Resolving the Complex-CDR Interdependency with Iterative Refinement

![Alt text](https://github.com/bio-ontology-research-group/NanoDesigner/blob/main/NanoDesigner_.png)

NanoDesigner is an end-to-end workflow for the design and optimization of nanobodies. It integrates key stages—Structure Prediction, Docking, CDR Generation,and Side-Chain Packing—into an iterative framework based on an Expectation Maximization algorithm. Our method effectively tackles an often overlooked interdependency
challenge where accurate docking presupposes a priori knowledge of the CDR conformation, while effective CDR generation relies on accurate docking outputs to guide its design.


## Table of Contents
- [Installation](#installation)
- [External Tool Installation](#external-tool-installation)
- [Data Download and Preprocess](#data-download-and-preprocess)
- [Training and Inference](#training-and-inference)
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
- [FoldX](https://foldxsuite.crg.eu/products#foldx)
- [HDOCK](http://huanglab.phys.hust.edu.cn/software/hdocklite/)
- [dr_sasa_n](https://github.com/nioroso-x3/dr_sasa_n) - *Follow the instructions in the repository to compile this tool.*

After installing the tools, ensure to update the `dyMEAN/configs.py` file with the full paths to the installed tools.

Source code for TMscore evatuation is at `dyMEAN/evaluation/`, please compile as:
```bash
g++ -static -O3 -ffast-math -lm -o evaluation/TMscore evaluation/TMscore.cpp
```


## Data Download and Preprocess

The data download and preparation steps are necessary to replicate our data processing, filtering, and preparation for training. All required instructions are included in the provided Jupyter notebooks.

### 1. Preprocess the Data
- Open the notebook located at `jupyter_notebooks/process_datasets.ipynb`.
- Follow the instructions in the notebook to download and preprocess the datasets.

### 2. Split the Data
- Once preprocessing is complete, open the notebook at `jupyter_notebooks/split_data.ipynb`.
- Use this notebook to split the processed data into training and testing sets.

### Notes:
- Ensure [Jupyter Notebook](https://jupyter.org/install) is installed. To check, run:
  ```bash
  jupyter --version


## Training and Inference

The `scripts` folder contains `.sh` scripts for both training and inference workflows used in the study. These scripts are configured for each of the tools employed in this study and are designed to facilitate a 10-fold cross-validation setup.

### 1. Training
- **Location**: Training scripts for each tool are located in the `scripts` directory.
- **Configuration**: Update the file paths and any necessary parameters inside the scripts. This includes specifying paths for datasets, output directories and additional variables.
- **10-Fold Cross-Validation**: The scripts are pre-configured to implement a 10-fold cross-validation strategy. Refer to *-Data Download and Preprocess*.

### 2. Inference
- **Location**: Inference scripts for each tool are also available in the `scripts` directory.
- **Configuration**: Make sure that the paths across the training and inference scripts match. The folder specified in the training script dictates the location of the generated checkpoints, which will be used during inference.
- **Manual Checkpoint Selection**: For GNN-based tools, selection of the best checkpoint must be done manually. Refer to the instructions provided in the script files for guidance.



To run the training or inference for a specific tool, execute the corresponding script as in the example:
```bash
bash scripts/train_tool.sh
```


**Note:** The 10-fold generated datasets, used to conduct the proof of concept for this project, can be found in [this Google Drive folder](https://drive.google.com/drive/folders/1CzBCQGvpHiBCufGCLoa15-fe9c0Mg1Xq?usp=share_link). For details on how these datasets were generated, please refer to the [Data Download and Preprocess](#data-download-and-preprocess) section.



## NanoDesigner

NanoDesigner is an end-to-end workflow designed for both **de novo** and **optimization** cases in nanobody-antigen complex design. The workflow script is located in the `scripts` folder and can be executed as follows:

```bash
bash scripts/NanoDesigner.sh your_working_directory/denovo_epitope_info/7eow_8pwh_example/7eow_8pwh_ep_1.json
```

The workflow requires a script and a JSON file containing the necessary information for each entry (a nanobody-antigen complex or nanobody scaffold and antigen structure). 

- **De Novo Design**: In cases where the 3D structure of a nanobody-antigen complex is absent (referred to as "de novo" design), the input JSON file can be generated using the notebook `jupyter_notebooks/prepare_NanoDesigner_inputs_Denovo.ipynb`. This notebook guides you through creating a properly formatted JSON file.

- **Optimization Cases**: For existing complexes, simply select a relevant line from the dataset-generated JSON files (prepared during the [Data Download and Preprocess](#data-download-and-preprocess) stage) and use it to create an input JSON file.

All required information for both cases should be obtained during the data download and preprocessing stage. Ensure the configuration files (`config_files`) are updated as needed to reflect your setup.

For proof of concepts of NanoDesigner, please download and employ DiffAb or ADesigner trained models found [here](https://drive.google.com/drive/folders/1kGK3rV138lG8vQpGAtHv5oNP_a11Gr01?usp=share_link).


We highly encourage to keep a constant number of total number of designs across iterations for simplicity:

```python
R = 50  # Number of randomized nanobodies (Initialization step)
N = 15  # Top best mutants to proceed with to subsequent iterations
d = 100 # Docked models to generate with Hdock
n = 5   # Top docked models to feed to inference stage
k_iteration_1 = 3   # Number of predictions obtained from CDR Generation stage at iteration 1
k_iteration_x = 10  # Number of predictions obtained from CDR Generation stage at iteration x

Rxnxk = 750 (Iteration 1)
Nxnxk = 750 (Iteration X)
```

### NanoDesigner test cases:
*De novo design escenario;
*CDRH3 or 3CDRs design with ΔG optimization objective.*

<div style="text-align: center;">
  <img src="https://github.com/Melissaurious/NanoDesigner/blob/main/test_cases.png" alt="Alt Text" width="500">
</div>



## Citation
TODO

## License
TODO

## Credits

This codebase is primarily based on the following deep learning tools. We thank the authors for their contributions:

- [Diffab](https://github.com/luost26/diffab)  
- [dyMEAN](https://github.com/THUNLP-MT/dyMEAN)  
- [ADesigner](https://github.com/A4Bio/ADesigner) 

We also acknowledge the rest of tools and software that played a crucial role in the workflow employed in this study.
We sincerely thank the authors of these tools for their invaluable work, which made this project possible.




