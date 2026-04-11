# CellSavior
A tool for cell calling from target-based single-cell DNA sequencing data




=> Basic Usage
  python CellSavior.py \
    --input_barcode_distribution ${input_barcode_distribution} \
    --input_bam ${input_mapped_bam} \
    --output_barcode_distribution ${output_barcode_distribution} \
    --output_bam ${output_cellsavior_bam} \
    --output_prefix_for_plots ${path-to-output_prefix_for_plots}
    --thread ${thread}

