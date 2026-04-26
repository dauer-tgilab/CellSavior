# CellSavior

<p align="center">

  <img src="https://github.com/user-attachments/assets/e090d84c-e869-4226-acf3-f123ddb6321f" width="400"/>

</p>

<p align="center">

  <b>A tool for cell calling from target-based single-cell DNA sequencing data</b>

</p>

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

