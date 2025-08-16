# APIARC: A Workflow for Integrated Analysis of RNA-seq and ChIP-seq Data  

## Overview  
APIARC is a computational workflow designed to analyze common genes from RNA-seq and ChIP-seq data, enabling the discovery of novel biological insights.  

## Installation  

### 1. Install Mamba and Snakemake  
Follow these steps to set up the environment:  

1. **Access your Linux terminal.**  
2. **Initialize the base environment:**  
   ```bash
   source ~/.bashrc
   ```  
3. **Install Mamba:**  
   ```bash
   conda install -n base -c conda-forge mamba
   ```  
4. **Create a new environment (e.g., `APIARC`) and install Snakemake:**  
   ```bash
   mamba create -c conda-forge -c bioconda -n APIARC snakemake
   ```  
5. **Activate the environment:**  
   ```bash
   conda activate APIARC
   ```  
6. **Install required dependencies using the environment file:**  
   ```bash
   mamba env update -n APIARC -f environment.yml
   ```  

**Note:** After completing the installation, you can activate the environment in future sessions by running:  
```bash
source ~/.bashrc; conda activate APIARC
```  

---

### 2. Deploy and Configure Workflows  

1. **Clone the APIARC repository:**  
   ```bash
   git clone https://github.com/yxiaobo/APIARC
   ```  

APIARC consists of the following workflows:  

#### Reference Workflow  
[`reference`](./reference): Downloads reference data required for all workflows.  

#### RNA-seq Workflow  
[`RNAProj`](./RNAProj): Processes RNA sequencing data.  

#### ChIP-seq Workflow  
[`CHIPProj`](./CHIPProj): Processes Chromatin Immunoprecipitation Sequencing data.  

#### Integrated Analysis Workflow  
[`integratedProj`](./integratedProj): Analyzes common genes identified by `RNAProj` and `CHIPProj`.  

**Important:**  
- Complete `RNAProj` and `CHIPProj` before running `integratedProj`.  
- The workflow concludes after `integratedProj` finishes execution.  

---

## System Requirements  

### Recommended Specifications  
- **Operating System:** Linux (tested on Ubuntu 20.04)  
- **Memory:** ≥ 16 GB RAM  
- **CPU:** ≥ 8 cores (10 cores tested)  
- **Storage:** ≥ 500 GB available space  

### Software Versions (Tested)  
- `mamba=2.0.8`  
- `snakemake=7.32.4`  

---

## Reproducibility  
APIARC is designed for Linux systems. Ensure all dependencies are installed as specified to guarantee reproducibility.  

For optimal performance, adhere to the recommended system specifications.