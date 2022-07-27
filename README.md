# Nanopore demultiplexing and variant calling pipeline
Nextflow pipeline designed for Will Hamilton (wh2) for quick onsite demultiplexing and variant calling of Oxford Nanopore data in Ghana
## Usage
```
nextflow run main.nf
  --sequencing_manifest        Manifest containing paths to sequencing directories and sequencing summary files (mandatory)
  --results_dir                Specify results directory - default nextflow_results (optional)
  --help                       Print this help message (optional)
```

## Sequencing manifest format
The sequencing manifest is in a csv format and contains two columns: `sequencing_dir` and `sequence_summary_file`

The sequencing directory is the folder which contains all the Oxford Nanopore sequencing data.

The sequence summary file is required for QC and will be in the sequencing directory.


Example manifest:
```
sequencing_dir,sequence_summary_file
20220721_multiplex_lab_mixes/no_sample/20220721_1106_MN39679_FAT74380_02d931ba,20220721_multiplex_lab_mixes/no_sample/20220721_1106_MN39679_FAT74380_02d931ba/sequencing_summary_FAT74380_e59d6c08.txt
```

Always provide the full paths to files or folders when creating a manifest.