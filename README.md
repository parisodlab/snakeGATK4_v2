# How to launch jobs
## 1. Prepare the input files
The input files are the following:
- Read files in FASTQ format
- Reference genome in FASTA format
- Reference genome index files `samtools faidx ref.fasta` and `bwa index ref.fasta`
- Reference genome dict file `gatk --java-options '-Xmx8g -Xms8g' CreateSequenceDictionary -R ref.fasta -O ref.dict`


## 2. Prepare the configuration file
The configuration file `config_main.yaml` contains the following fields:
- `PROJECT`: the name of the project
- `QC`: the quality control parameters, if required or not
- `TRIMMING`: the trimming parameters, if required or not
- `GENOME`: the genome file path
- `SCAFFOLDS`: whether to process all scaffolds or only a few scaffolds
- `METAFILE`: the metafile path	with the following fields:
    - `sample`: the sample name
    - `prefix`: the prefix of the sample
    - `sample_barcode`: the sample barcode #TODO: check if this is required
    - `flowcell_barcode`: the flowcell barcode
    - `lane`: the lane number
    - `pop`: the population
    - `ploidy`: the ploidy
    - `fq1`: the first read file path
    - `fq2`: the second read file path
    - `flow_cell_type`: the flow cell type #TODO: check if this is required
    - `pcr_indel_model`: the PCR indel model #TODO: check if this is required
- `FINALOUTPUT`: the final output path
- `tmpdir`: the temporary directory path



## 3. Launch the job
The job can be launched using the following command:
```bash
python main_pipeline.py
```

