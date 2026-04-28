<div align="center">

<h1>CellSavior</h1>

<img src="https://github.com/user-attachments/assets/e090d84c-e869-4226-acf3-f123ddb6321f" width="450"/>

<br>

<h3><b>A tool for cell calling from targeted single-cell DNA sequencing data</b></h3>

</div>
<br>
CellSavior is an open-source algorithm for cell calling in targeted single-cell DNA sequencing data. It identifies true cellular barcodes from noise using ensemble clustering across multiple distance metrics. The method integrates clustering results with confidence-weighted scoring to improve robustness. CellSavior enhances reproducibility and may recover biologically relevant cell populations missed by existing pipelines.
<br><br>



## Install

### Cloning GitHub repository

```
> git clone https://github.com/dauer-tgilab/CellSavior.git
> cd CellSavior
```



## Usage

### Basic Usage

```
> python CellSavior.py
> --input_barcode_distribution /path/to/${sample}.all.barcode.distribution.tsv \
> --input_bam /path/to/${sample}.mapped.bam \
> --output_barcode_distribution /path/to/${sample}.cellsavior.barcode.distribution.tsv \
> --output_bam /path/to/${sample}.cellsavior.bam \
> --output_prefix_for_plots /path/to/${sample}_prefix \
> --thread ${thread}
```



### h5 generation
If single-cell variant calling has been performed, HDF5 file generation is possible. It requires the barcode × variant AF and DP matrices, along with the previously mentioned input files.

```

```



## Data Format

### Input
1. mapped.bam file<br>
   A mapped BAM file contains sequencing reads aligned to a reference genome, with barcode information stored in the RG tag for each read.
   
2. all.barcode.distribution.tsv file
```
barcode  target1  target2  ...  target_n
barcode_1  depth[1,1]  depth[1,2]  ...  depth[1,n]
barcode_2  depth[2,1]  depth[2,2]  ...  depth[2,n] 
...
barcode_m  depth[m,1]  depth[m,2]  ...  depth[m,n]
```

### Output
1. cellsavior.bam file
   The Cellsavior BAM file contains reads corresponding to barcodes identified as valid cells, preserving alignment information and RG-tagged barcode annotations for downstream analysis.
   
2. cellsavior.barcode.distribution.tsv file
```
cell_barcode  target1  target2  ...  target_n
cell_barcode_1  depth[1,1]  depth[1,2]  ...  depth[1,n]
cell_barcode_2  depth[2,1]  depth[2,2]  ...  depth[2,n] 
...
cell_barcode_m  depth[m,1]  depth[m,2]  ...  depth[m,n]
```



## Requirements
CellSavior was developed using Python 3.11.
```
  hdbscan>=0.8.41
  matplotlib-base>=3.10.8
  matplotlib-inline>=0.2.1
  numba>=0.64.0
  numpy>=2.4.2
  pandas>=3.0.1
  scikit-learn>=1.8.0
  scipy>=1.17.1
  seaborn>=0.13.2
  statsmodels>=0.14.6
```
