# CellSavior
 A tool for cell calling from target-based single-cell DNA sequencing data
 


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

