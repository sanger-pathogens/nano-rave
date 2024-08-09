#!/usr/bin/env nextflow
// Copyright (C) 2022,2023 Genome Research Ltd.

nextflow.enable.dsl = 2

def printHelp() {
    log.info """
    Usage:
        nextflow run main.nf

    Options:
        --sequencing_manifest        Manifest containing paths to sequencing directories and sequencing summary files (mandatory)
        --reference_manifest         Manifest containing reference identifiers and paths to fastq reference files (mandatory)
        --results_dir                Specify results directory [default: ./nextflow_results] (optional)
        --variant_caller             Specify a variant caller to use [medaka (default), medaka_haploid, freebayes, clair3] (optional)
        --clair3_args                Specify clair3 variant calling parameters - must include model e.g. --clair3_args "--model_path /opt/models/r941_prom_sup_g5014" (optional)
        --min_barcode_dir_size       Specify the expected minimum size of the barcode directories, in MB. Must be > 0. [default: 10] (optional)
        --keep_bam_files             Save BAM files in results directory [default: false] (optional)
        --help                       Print this help message (optional)
    """.stripIndent()
}

def validate_path_param(
    param_option, 
    param, 
    type="file", 
    mandatory=true) 
{
    valid_types=["file", "directory"]
    if (!valid_types.any { it == type }) {
            log.error("Invalid type '${type}'. Possibilities are ${valid_types}.")
            return 1
    }
    param_name = (param_option - "--").replaceAll("_", " ")
    if (param) {
        def file_param = file(param)
        if (!file_param.exists()) {
            log.error("The given ${param_name} '${param}' does not exist.")
            return 1
        } else if (
              (type == "file" && !file_param.isFile())
              ||
              (type == "directory" && !file_param.isDirectory())
          ) {
            log.error("The given ${param_name} '${param}' is not a ${type}.")
            return 1
        }
    } else if (mandatory) {
        log.error("No ${param_name} specified. Please specify one using the ${param_option} option.")
        return 1
    }
    return 0
}

def validate_choice_param(param_option, param, choices) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (param) {
        if (!choices.any { it.contains(param.toString()) }) {
            log.error("Please specify the ${param_name} using the ${param_option} option. Possibilities are ${choices}.")
            return 1
        }
    } else {
        log.error("Please specify the ${param_name} using the ${param_option} option")
        return 1
    }
    return 0
}

def validate_number_param(param_option, param) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (param != null) /* Explicit comparison with null, because 0 is an acceptable value */ {
        if (!(param instanceof Number)) {
            log.error("The ${param_name} specified with the ${param_option} option must be a valid number")
            return 1
        }
    } else {
        log.error("Please specify the ${param_name} using the ${param_option} option")
        return 1
    }
    return 0
}

def validate_min_barcode_dir_size(param_option, param) {
    if (validate_number_param(param_option, param) == 1) {
        return 1
    }
    param_name = (param_option - "--").replaceAll("_", " ")
    if (!(param > 0)) {
        log.error("The ${param_name} specified with the ${param_option} option must have a positive value")
        return 1
    }
    return 0
}

def validate_results_dir(results_dir) {
    results_dir = file(results_dir)
    if (results_dir.exists() && !results_dir.isDirectory()) {
        log.error("The given results_dir '${results_dir}' is not a directory.")
        return 1
    }
    return 0
}

def validate_clair3_args(clair3_args) {
    def invalid_options = ['--bam_fn', '--ref_fn', '--threads', '--platform', '--output']
    def required_options = ['--model_path']
    def errors = 0
    if (params.variant_caller != "clair3" && clair3_args) {
        log.error("Clair3 arguments were provided but clair3 was not set as the --variant_caller!")
        errors += 1
    } else if (params.variant_caller == "clair3" && clair3_args) {
        def invalid_options_found = invalid_options.findAll { clair3_args.contains(it) }
        if (invalid_options_found) {
            log.error("The following clair3 options were provided in --clair3_args, but are reserved for use by this pipeline: ${invalid_options_found}")
            errors += 1
        }
        def required_options_found = required_options.findAll { clair3_args.contains(it) }
        def required_options_not_found = required_options.minus(required_options_found)
        if (required_options_not_found) {
            log.error("The following clair3 options were not provided in --clair3_args, but are required for use by this pipeline: ${required_options_not_found}")
            errors += 1
        }
    }
    return errors
}

def validate_parameters() {
    def errors = 0

    errors += validate_path_param("--reference_manifest", params.reference_manifest)
    errors += validate_path_param("--sequencing_manifest", params.sequencing_manifest)
    errors += validate_choice_param("--variant_caller", params.variant_caller, ["medaka", "medaka_haploid", "freebayes", "clair3"])
    errors += validate_min_barcode_dir_size("--min_barcode_dir_size", params.min_barcode_dir_size)
    errors += validate_results_dir(params.results_dir)
    errors += validate_clair3_args(params.clair3_args)

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}

process SORT_FASTQS {
    /* Checks each barcode directory contains sufficient reads and concatenates fastq.gz files for further processing */
    input:
        tuple val(sequencing_dir), val(sequence_summary_file)
    output:
        path("*.fastq.gz"), emit: full_fastq_files
    script:
        """
        threshold=${params.min_barcode_dir_size}
        sample_name=\$(echo $sequencing_dir | awk -F "/" '{ print \$(NF-1) }')
        echo \$sample_name
        for dir in ${sequencing_dir}/fastq_pass/barcode*
        do
          barcode=\$(echo \${dir} | awk -F "/" '{ print \$NF }')
          disk_usage=\$(du --apparent-size -shm \$dir | awk '{ print \$1 }')
          if [ "\${disk_usage}" -ge "\${threshold}" ]
          then
            zcat \${dir}/*.fastq.gz > \${sample_name}_\${barcode}.fastq
          else
            echo "WARN: Skipping '\${barcode}' directory '\${dir}' as it contains <\${threshold}MB of fastq.gz files." >&2
          fi
        done
        gzip *.fastq
        """
}

process NANOPLOT_QC {
    conda "bioconda::nanoplot=1.38.0"
    container "quay.io/biocontainers/nanoplot:1.38.0--pyhdfd78af_0"
    publishDir "${params.results_dir}/qc/nanoplot", mode: 'copy', overwrite: true
    input:
        path(fastq_file)
    output:
        path("*_nanoplot_qc/*")
    script:
        """
        sample=\$(basename $fastq_file | awk -F "." '{ print \$1}')
        NanoPlot -t 2 --fastq $fastq_file -o \${sample}_nanoplot_qc
        """
}

process PYCOQC {
    conda "bioconda::pycoqc=2.5.2"
    container "quay.io/biocontainers/pycoqc:2.5.2--py_0"
    publishDir "${params.results_dir}/qc/pycoqc", mode: 'copy', overwrite: true
    input:
        tuple val(sequencing_dir), val(sequence_summary_file)
    output:
        path("*.html")
        path("*.json")
    script:
        """
        sample=\$(echo $sequencing_dir | awk -F "/" '{ print \$(NF-1) }')
        pycoQC -f $sequence_summary_file -o \${sample}_pycoqc.html -j \${sample}_pycoqc.json
        """
}

process NORMALISE_FASTAS {
    conda "conda-forge::biopython=1.78"
    container "quay.io/biocontainers/biopython:1.78"
    input:
        tuple val(ref_id), path(reference)
    output:
        tuple val(ref_id), path(reference), emit: normalised_ref_ch
    script:
        """
        #!/usr/bin/env python3
        from Bio import SeqIO
        import shutil

        records = SeqIO.parse("${reference}", "fasta")
        SeqIO.write(records, "${reference}.normalised", "fasta")
        shutil.move("${reference}.normalised", "${reference}")
        """
}

process GET_CHROM_SIZES_AND_INDEX {
    conda "bioconda::samtools=1.15.1"
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    input:
        tuple val(ref_id), path(reference)
    output:
        tuple val(ref_id), path("*.fai"), emit: fasta_index_ch
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        samtools faidx ${reference}
        """
}

process MINIMAP2_INDEX {
    conda "bioconda::minimap2=2.17"
    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"
    input:
        tuple val(ref_id), path(reference)
    output:
        path("*.mmi"), emit: mm2_index_ch
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        minimap2 -ax map-ont -t 2 -d ${reference}.mmi ${reference}
        """
}

process MINIMAP2_ALIGN {
    conda "bioconda::minimap2=2.17"
    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"
    input:
        tuple val(ref_id), path(reference)
        path(mm2_index)
        each path(fastq)
    output:
        path("*.sam"), emit: sam_file
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        ref_id=\$(echo $reference | awk -F "seq_" '{ print \$NF }' | sed 's|.fasta||g')
        fname=\$(basename $fastq | awk -F "." '{ print \$1}')
        minimap2 -ax map-ont --MD -t 2 $mm2_index $fastq > \${fname}_\${ref_id}.sam
        """
}

process SAMTOOLS_VIEW_SAM_TO_BAM {
    conda "bioconda::samtools=1.15.1"
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    input:
        tuple val(ref_id), path(reference)
        path(sam_file)
    output:
        path("*.bam"), emit: bam_file
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        filename=\$(basename $sam_file | awk -F "." '{ print \$1}')
        samtools view -b -h -O BAM -@ 2 -o \${filename}.bam $sam_file
        """
}

process SAMTOOLS_SORT_AND_INDEX {
    conda "bioconda::samtools=1.15.1"
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    if (params.keep_bam_files) {
        publishDir "${params.results_dir}/bams", mode: 'copy', overwrite: true, pattern: "*.bam*"
    }
    input:
        tuple val(ref_id), path(reference)
        path(bam_file)
    output:
        tuple path("*.sorted.bam"), path("*.sorted.bam.bai"),  emit: sorted_bam
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        filename=\$(basename $bam_file | awk -F "." '{ print \$1}')
        samtools sort -@ 2 -o \${filename}.sorted.bam -T \${filename} $bam_file
        samtools index \${filename}.sorted.bam
        """
}

process BEDTOOLS_GENOMECOV {
    conda "bioconda::bedtools=2.29.2"
    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    publishDir "${params.results_dir}/genome_coverage", mode: 'copy', overwrite: true, pattern: "*.bedGraph"
    input:
        tuple path(sorted_bam_file), path(sorted_bam_index)

    output:
        path("*.bedGraph"), emit: genome_cov_ch

    script:
        """
        filename=\$(basename $sorted_bam_file | awk -F "." '{ print \$1}')
        bedtools genomecov -split -ibam $sorted_bam_file -bga | bedtools sort > \${filename}.bedGraph
        """
}

process MEDAKA_VARIANT_CALLING {
    conda "bioconda::medaka=1.4.4"
    container "quay.io/biocontainers/medaka:1.4.4--py38h130def0_0"
    input:
        tuple val(ref_id), path(reference)
        tuple path(sorted_bam_file), path(sorted_bam_index)
    output:
        path("*.vcf"), emit: vcf_ch
    script:
        """
        filename=\$(basename $sorted_bam_file | awk -F "." '{ print \$1}')
        
        medaka_variant -f ${reference} -i ${sorted_bam_file}

        mv medaka_variant/round_1.vcf \${filename}.vcf

        # sort vcf by index to stop tabix crying
        cat \${filename}.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > \${filename}_sorted.vcf
        mv \${filename}_sorted.vcf \${filename}.vcf
        """
}

process MEDAKA_HAPLOID_VARIANT_CALLING {
    conda "bioconda::medaka=1.4.4"
    container "quay.io/biocontainers/medaka:1.4.4--py38h130def0_0"
    input:
        tuple val(ref_id), path(reference)
        each path(fastq)
    output:
        path("*.vcf"), emit: vcf_ch
    script:
        """
        filename=\$(basename ${fastq} | awk -F '.' '{ print \$1}')_${ref_id}

        medaka_haploid_variant -r ${reference} -i ${fastq}

        mv medaka/medaka.annotated.vcf \${filename}.vcf

        # sort vcf by index to stop tabix crying
        cat \${filename}.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > \${filename}_sorted.vcf
        mv \${filename}_sorted.vcf \${filename}.vcf
        """
}

process FREEBAYES_VARIANT_CALLING {
    conda "bioconda::freebayes=1.3.7"
    container "quay.io/biocontainers/freebayes:1.3.7--h1870644_0"
    input:
        tuple val(ref_id), path(reference)
        tuple path(sorted_bam_file), path(sorted_bam_index)
    output:
        path("*.vcf"), emit: vcf_ch
    script:
        """
        filename=\$(basename ${sorted_bam_file} | awk -F "." '{ print \$1}')
        
        freebayes -f ${reference} ${sorted_bam_file} > \${filename}.vcf
        
        # sort vcf by index to stop tabix crying
        cat \${filename}.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > \${filename}_sorted.vcf
        mv \${filename}_sorted.vcf \${filename}.vcf
        """
}

process CLAIR3_VARIANT_CALLING {
    conda "bioconda::clair3=1.0.0"
    container "docker.io/hkubal/clair3@sha256:3c4c6db3bb6118e3156630ee62de8f6afef7f7acc9215199f9b6c1b2e1926cf8"  // Includes models
    publishDir "${params.results_dir}/variant_calling/gvcf", mode: 'copy', overwrite: true, pattern: "*.gvcf.gz*"
    input:
        tuple val(ref_id), path(reference), path(reference_index)
        tuple path(sorted_bam_file), path(sorted_bam_index)
    output:
        path("*.vcf"), emit: vcf_ch
        tuple path("*.gvcf.gz"), path("*.gvcf.gz.tbi"),  optional: true, emit: gvcf_ch
    script:
        """
        filename=\$(basename ${sorted_bam_file} | awk -F "." '{ print \$1}')
        
        run_clair3.sh \
            --bam_fn=${sorted_bam_file} \
            --ref_fn=${reference} \
            --threads=${task.cpus} \
            --platform="ont" \
            --output=. \
            ${params.clair3_args}

        if [[ ! -s "merge_output.vcf.gz" ]]; then
            cp pileup.vcf.gz merge_output.vcf.gz
        fi
        # sort vcf by index to stop tabix crying
        zcat merge_output.vcf.gz | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > \${filename}.vcf
        if [[ -f "merge_output.gvcf.gz" ]]; then
            # rename for publishing
            mv merge_output.gvcf.gz \${filename}.gvcf.gz
            mv merge_output.gvcf.gz.tbi \${filename}.gvcf.gz.tbi
        fi
        """
}

process BGZIP_AND_INDEX_VCF {
    conda "bioconda::tabix=1.11"
    container "quay.io/biocontainers/tabix:1.11--hdfd78af_0"
    publishDir "${params.results_dir}/variant_calling/vcf", mode: 'copy', overwrite: true
    input:
        path(vcf_file)

    output:
        path("*.vcf.gz"), emit: bgzip_vcf_file_ch
        path("*.vcf.gz.tbi"), emit: vcf_index_ch

    script:
        """
        filename=\$(basename $vcf_file| awk -F "." '{ print \$1}')
        bgzip -c $vcf_file > \${filename}.vcf.gz
        tabix \${filename}.vcf.gz
        """
}

workflow NANOSEQ {
    take:
        fastq_ch

    main:
        ref_manifest_ch = Channel.fromPath(params.reference_manifest)
        ref_path_ch = ref_manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.reference_id, file(row.reference_path)) }

        NANOPLOT_QC(fastq_ch)

        NORMALISE_FASTAS(ref_path_ch)

        GET_CHROM_SIZES_AND_INDEX(
            NORMALISE_FASTAS.out.normalised_ref_ch
        )

        MINIMAP2_INDEX(
            GET_CHROM_SIZES_AND_INDEX.out.ref_ch
        )

        MINIMAP2_ALIGN(
            MINIMAP2_INDEX.out.ref_ch,
            MINIMAP2_INDEX.out.mm2_index_ch,
            fastq_ch
        )

        SAMTOOLS_VIEW_SAM_TO_BAM(
            MINIMAP2_ALIGN.out.ref_ch,
            MINIMAP2_ALIGN.out.sam_file
        )

        SAMTOOLS_SORT_AND_INDEX(
            SAMTOOLS_VIEW_SAM_TO_BAM.out.ref_ch,
            SAMTOOLS_VIEW_SAM_TO_BAM.out.bam_file,
        )

        BEDTOOLS_GENOMECOV(
            SAMTOOLS_SORT_AND_INDEX.out.sorted_bam
        )

        if (params.variant_caller == "medaka_haploid") {
            MEDAKA_HAPLOID_VARIANT_CALLING(
                NORMALISE_FASTAS.out.normalised_ref_ch,
                fastq_ch
            )
            BGZIP_AND_INDEX_VCF(MEDAKA_HAPLOID_VARIANT_CALLING.out.vcf_ch)
        }

        if (params.variant_caller == "medaka") {
            MEDAKA_VARIANT_CALLING(
                SAMTOOLS_SORT_AND_INDEX.out.ref_ch,
                SAMTOOLS_SORT_AND_INDEX.out.sorted_bam
            )
            BGZIP_AND_INDEX_VCF(MEDAKA_VARIANT_CALLING.out.vcf_ch)
        }

        if (params.variant_caller == "freebayes") {
            FREEBAYES_VARIANT_CALLING(
                SAMTOOLS_SORT_AND_INDEX.out.ref_ch,
                SAMTOOLS_SORT_AND_INDEX.out.sorted_bam
            )
            BGZIP_AND_INDEX_VCF(FREEBAYES_VARIANT_CALLING.out.vcf_ch)
        }

        if (params.variant_caller == "clair3") {
            SAMTOOLS_SORT_AND_INDEX.out.ref_ch
                .combine(GET_CHROM_SIZES_AND_INDEX.out.fasta_index_ch, by: 0)
                .set { clair3_ref_input }
            CLAIR3_VARIANT_CALLING(
                clair3_ref_input,
                SAMTOOLS_SORT_AND_INDEX.out.sorted_bam
            )
            BGZIP_AND_INDEX_VCF(CLAIR3_VARIANT_CALLING.out.vcf_ch)
        }
}

workflow {
    if (params.help) {
        printHelp()
        exit 0
    }

    validate_parameters()

    sequencing_manifest = Channel.fromPath(params.sequencing_manifest)
    sequence_ch = sequencing_manifest.splitCsv(header: true, sep: ',')
                    .map{ row -> tuple(file(row.sequencing_dir), file(row.sequence_summary_file)) }

    SORT_FASTQS(sequence_ch)

    PYCOQC(sequence_ch)

    NANOSEQ(SORT_FASTQS.out.full_fastq_files.flatten())
}
