# Bisulfite Sequencing (Biseq) Processing Pipeline

This pipeline is designed to process bisulfite sequencing (biseq) data through several stages, including genome preparation, read mapping, deduplication, and filtering based on read depth. Below is an overview of the steps involved in this pipeline, along with a description of each directory.

---

## Directory Structure
03_biseq-processing/
├── 03.0_genome-preparation/
├── 03.1_mapping/
├── 03.2_deduplication/
└── 03.3_filtering/

---

### 1. Genome Preparation (`03.0_genome-preparation/`)

**Purpose**: This step involves preparing the genome for bisulfite sequencing analysis. The genome is indexed for alignment in subsequent steps.

- **Files**:
  - `Makefile`: Contains commands for preparing the genome.
  - `aa_targets`: Holds genome-specific files and configurations.


---

### 2. Read Mapping (`03.1_mapping/`)

**Purpose**: In this step, the bisulfite-treated reads are aligned to the prepared reference genome. The mapping step is critical to ensure that the reads are correctly aligned for further analysis.

- **Files**:
  - `Makefile`: Automates the mapping process.
  - `aa_targets`: Contains the list of input files and parameters for the mapping step.

---

### 3. Deduplication (`03.2_deduplication/`)

**Purpose**: After mapping, duplicate reads that likely originated from PCR amplification are removed. Deduplication reduces redundancy and improves the accuracy of downstream analyses.

- **Files**:
  - `Makefile`: Automates the deduplication process.
  - `aa_targets`: Lists the input files and settings for deduplication.

---

### 4. Read Filtering (`03.3_filtering/`)

**Purpose**: This step filters the reads based on certain criteria, such as read depth. Filtering ensures that only high-quality data is used in the final analysis.

- **Subdirectories**:
  - `aa_ctxt-split`: Contains data split by context.
  - `ab_chr-split`: Contains data split by chromosome.
  - `ac_filter`: Filtering by number of C per 100bp-window, 5X read depth.

---

## Usage

Each step of the pipeline can be executed by running the corresponding `Makefile` in each directory. Before running `make`, you may need to generate the necessary targets using the `create-targets.sh` script.

```bash
# Step 1: Generate the necessary targets
~/bin/create-targets.sh ../samples2target

# Step 2: Run the Makefile with multiple threads
make -j<number of threads>
```
---

## Contact

For any issues or questions about this pipeline, please contact Veronica N. Ibanez at [veronicanoeibanez@gmail.com].

---
