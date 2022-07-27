#!/usr/bin/env nextflow

def printHelp() {
    log.info """
    Usage:
        nextflow run main.nf

    Options:
        --sequencing_manifest        Manifest containing paths to sequencing directories and sequencing summary files (mandatory)
        --results_dir                Specify results directory - default nextflow_results (optional)
        --help                       Print this help message (optional)
    """.stripIndent()
}

def validate_parameters() {
    def errors = 0

    if (params.reference_manifest) {
        reference_manifest=file(params.reference_manifest)
        if (!reference_manifest.exists()) {
            log.error("The manifest file specified does not exist.")
            errors += 1
        }
    }
    else {
        log.error("No reference manifest file specified. Please specify one using the --reference_manifest option.")
        errors += 1
    }

    if (params.sequencing_manifest) {
        sequencing_manifest=file(params.sequencing_manifest)
        if (!sequencing_manifest.exists()) {
            log.error("The sequencing manifest file specified does not exist.")
            errors += 1
        }
    }
    else {
        log.error("No sequencing manifest file specified. Please specify one using the --sequencing_manifest option.")
        errors += 1
    }


    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}

process SORT_FASTQS {
    input:
        tuple val(sequencing_dir), val(sequence_summary_file)
    output:
        path("*.fastq.gz"), emit: full_fastq_files
    script:
        """
        threshold=10
        sample_name=\$(echo $sequencing_dir | awk -F "/" '{ print \$(NF-1) }')
        echo \$sample_name
        for dir in ${sequencing_dir}/fastq_pass/barcode*
        do
          disk_usage=\$(du -shm \$dir | awk '{ print \$1 }')
          if [ "\${disk_usage}" -gt "\${threshold}" ]
          then
            barcode=\$(echo \${dir} | awk -F "/" '{ print \$NF }')
            zcat \${dir}/*.fastq.gz > \${sample_name}_\${barcode}.fastq
          fi
        done
        gzip *.fastq
        """
}

process NANOPLOT_QC {
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

process GET_CHROM_SIZES_AND_INDEX {
    container "quay.io/biocontainers/samtools:1.13--h8c37831_0"
    input:
        tuple val(ref_id), path(reference)
    output:
        tuple val(ref_id), path("*.sizes"), emit: sizes_ch
        tuple val(ref_id), path("*.fai"), emit: fasta_index_ch
        tuple val(ref_id), path(reference), emit: ref_ch
    script:
        """
        samtools faidx ${reference}
        cut -f 1,2 ${reference}.fai > ${reference}.sizes
        """
}

process MINIMAP2_INDEX {
    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"
    input:
        tuple val(ref_id), path(reference)
        tuple val(ref_id), path(sizes)
    output:
        path("*.mmi"), emit: mm2_index_ch
        path(reference), emit: ref_ch
        path(sizes), emit: sizes_ch
    script:
        """
        minimap2 -ax map-ont -t 2 -d ${reference}.mmi ${reference}
        """
}

process MINIMAP2_ALIGN {
    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"
    input:
        path(mm2_index)
        path(reference)
        path(fastq_file_list)
        path(sizes)
    output:
        path("*.sam"), emit: sam_file
    script:
        """
        ref_id=\$(echo $reference | awk -F "seq_" '{ print \$NF }' | sed 's|.fasta||g')
        for f in ${fastq_file_list}
        do
          fname=\$(basename \$f | awk -F "." '{ print \$1}')
          minimap2 -ax map-ont --MD -t 2 $mm2_index \$f > \${fname}_\${ref_id}.sam
        done
        """

}

process SAMTOOLS_VIEW_SAM_TO_BAM {
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    input:
        path(sam_file)
    output:
        path("*.bam"), emit: bam_file
    script:
        """
        filename=\$(basename $sam_file | awk -F "." '{ print \$1}')
        samtools view -b -h -O BAM -@ 2 -o \${filename}.bam $sam_file
        """
}

process SAMTOOLS_SORT_AND_INDEX {
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    input:
        path(bam_file)
    output:
        path("*.sorted.bam"), emit: sorted_bam
        path("*.sorted.bam.bai"), emit: sorted_bam_index
    script:
        """
        filename=\$(basename $bam_file | awk -F "." '{ print \$1}')
        samtools sort -@ 2 -o \${filename}.sorted.bam -T \${filename} $bam_file
        samtools index \${filename}.sorted.bam
        """
}

process BEDTOOLS_GENOMECOV {
    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    publishDir "${params.results_dir}/genome_coverage", mode: 'copy', overwrite: true, pattern: "*.bedGraph"
    input:
        path(sorted_bam_file)
        path(sorted_bam_index)

    output:
        path("*.bedGraph"), emit: genome_cov_ch

    script:
        """
        filename=\$(basename $sorted_bam_file | awk -F "." '{ print \$1}')
        bedtools genomecov -split -ibam $sorted_bam_file -bg | bedtools sort > \${filename}.bedGraph
        """
}

process SNIFFLES_VARIANT_CALLING {
    container "quay.io/biocontainers/sniffles:1.0.12--h8b12597_1"
    input:
        path(sorted_bam_file)
        path(sorted_bam_index)

    output:
        path("*.vcf"), emit: vcf_ch

    script:
        """
        filename=\$(basename $sorted_bam_file | awk -F "." '{ print \$1}')
        sniffles -m $sorted_bam_file -v \${filename}.vcf
        # sort vcf by index to stop tabix crying
        cat \${filename}.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > \${filename}_sorted.vcf
        mv \${filename}_sorted.vcf \${filename}.vcf
        """
}

process BGZIP_AND_INDEX_VCF {
    container "quay.io/biocontainers/tabix:1.11--hdfd78af_0"
    publishDir "${params.results_dir}/variant_calling", mode: 'copy', overwrite: true
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

    NANOPLOT_QC(fastq_ch.flatMap())

    GET_CHROM_SIZES_AND_INDEX(ref_path_ch)

    MINIMAP2_INDEX(GET_CHROM_SIZES_AND_INDEX.out.ref_ch,
                   GET_CHROM_SIZES_AND_INDEX.out.sizes_ch)

    MINIMAP2_ALIGN(MINIMAP2_INDEX.out.mm2_index_ch,
                   MINIMAP2_INDEX.out.ref_ch,
                   fastq_ch.collect(),
                   MINIMAP2_INDEX.out.sizes_ch)

    SAMTOOLS_VIEW_SAM_TO_BAM(MINIMAP2_ALIGN.out.sam_file.flatMap())

    SAMTOOLS_SORT_AND_INDEX(SAMTOOLS_VIEW_SAM_TO_BAM.out.bam_file)

    BEDTOOLS_GENOMECOV(SAMTOOLS_SORT_AND_INDEX.out.sorted_bam,
                       SAMTOOLS_SORT_AND_INDEX.out.sorted_bam_index)

    SNIFFLES_VARIANT_CALLING(SAMTOOLS_SORT_AND_INDEX.out.sorted_bam,
                             SAMTOOLS_SORT_AND_INDEX.out.sorted_bam_index)

    BGZIP_AND_INDEX_VCF(SNIFFLES_VARIANT_CALLING.out.vcf_ch)
}

workflow {
    if (params.help) {
        printHelp()
        exit 0
    }

    validate_parameters()

    sequencing_manifest = Channel.fromPath(params.sequencing_manifest)
    sequence_ch = sequencing_manifest.splitCsv(header: true, sep: ',')
                    .map{ row -> tuple(row.sequencing_dir, row.sequence_summary_file) }

    SORT_FASTQS(sequence_ch)

    PYCOQC(sequence_ch)

    NANOSEQ(SORT_FASTQS.out.full_fastq_files)
}