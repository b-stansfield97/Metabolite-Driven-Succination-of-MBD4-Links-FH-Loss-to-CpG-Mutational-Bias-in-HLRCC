#!/bin/bash
# =============================================================================
#  HLRCC Whole-Exome Somatic Variant Calling Pipeline
# =============================================================================
#
#  Purpose
#  -------
#  End-to-end tumour-normal somatic variant calling for a single HLRCC
#  whole-exome pair. Aligns paired-end FASTQs to GRCh38 with BWA-MEM,
#  marks duplicates with Picard, calls somatic SNVs/indels with GATK
#  Mutect2 in paired (tumour + matched normal) mode, filters with
#  FilterMutectCalls and annotates to MAF format with Funcotator.
#
#  The resulting MAF (with a `ref_context` field) is the input to
#  HLRCC_Exome_Analysis.r, which performs SBS mutational signature
#  extraction for the manuscript.
#
#  Usage
#  -----
#      ./HLRCC_Exome_Variant_Calling.sh
#
#  All paths / sample names are read from `config.yaml` in the current
#  working directory (see config.yaml.example for a template). Run once
#  per tumour-normal pair, swapping in a different config.yaml each time.
#
#  Inputs (resolved from config.yaml)
#  ----------------------------------
#    reference.fasta              GRCh38 reference FASTA (must be indexed
#                                 with `bwa index`, `samtools faidx`, and
#                                 `gatk CreateSequenceDictionary`).
#    reference.funcotator_sources Funcotator data-sources directory.
#    reference.ref_version        Funcotator --ref-version (e.g. hg38).
#    resources.threads            Threads for BWA / samtools sort.
#    resources.outdir             Output directory (will be created).
#    samples.tumor.{name,fastq1,fastq2}
#    samples.normal.{name,fastq1,fastq2}
#
#  Outputs
#  -------
#    $OUTDIR/bam/<sample>_sorted.bam        BWA-MEM aligned, sorted BAM
#    $OUTDIR/bam/<sample>_rg.bam            With read-group header
#    $OUTDIR/bam/<sample>_dedup.bam         Duplicates marked
#    $OUTDIR/bam/<sample>_dedup_metrics.txt MarkDuplicates metrics
#    $OUTDIR/vcf/<tumor>_unfiltered.vcf.gz  Raw Mutect2 calls
#    $OUTDIR/vcf/<tumor>_filtered.vcf.gz    FilterMutectCalls output
#    $OUTDIR/maf/<tumor>_annotated.maf      Funcotator MAF (input to R)
#
#  Dependencies (must be on PATH)
#  ------------------------------
#    yq      (mikefarah/yq, v4+; YAML parsing)
#    bwa     (>= 0.7.17)
#    samtools (>= 1.10)
#    picard  (CLI wrapper, e.g. via conda)
#    gatk    (>= 4.2.0.0)
#
#  Paper
#  -----
#  Metabolite-Driven Succination of MBD4 Links FH Loss to CpG Mutational
#  Bias in HLRCC.
#
#  Author: Ben Stansfield
# =============================================================================
set -euo pipefail

# -----------------------------------------------------------------------------
# 1. Load configuration
# -----------------------------------------------------------------------------
CONFIG="config.yaml"

REF_GENOME=$(yq e '.reference.fasta' "$CONFIG")
FUNCOTATOR_DATA_SOURCES=$(yq e '.reference.funcotator_sources' "$CONFIG")
FUNCOTATOR_REF_VERSION=$(yq e '.reference.ref_version' "$CONFIG")
THREADS=$(yq e '.resources.threads' "$CONFIG")
OUTDIR=$(yq e '.resources.outdir' "$CONFIG")

TUMOR_SAMPLE=$(yq e '.samples.tumor.name'   "$CONFIG")
TUMOR_FASTQ1=$(yq e '.samples.tumor.fastq1' "$CONFIG")
TUMOR_FASTQ2=$(yq e '.samples.tumor.fastq2' "$CONFIG")

NORMAL_SAMPLE=$(yq e '.samples.normal.name'   "$CONFIG")
NORMAL_FASTQ1=$(yq e '.samples.normal.fastq1' "$CONFIG")
NORMAL_FASTQ2=$(yq e '.samples.normal.fastq2' "$CONFIG")

mkdir -p "$OUTDIR"/{bam,vcf,maf}

# -----------------------------------------------------------------------------
# 2. Tumour: align -> add read groups -> mark duplicates
# -----------------------------------------------------------------------------
# BWA-MEM alignment piped into samtools for BAM conversion and coordinate sort.
bwa mem -M -t "$THREADS" "$REF_GENOME" "$TUMOR_FASTQ1" "$TUMOR_FASTQ2" | \
  samtools view -Sb - | \
  samtools sort -@ "$THREADS" -o "$OUTDIR/bam/${TUMOR_SAMPLE}_sorted.bam"

# Add a read-group header (required by GATK).
picard AddOrReplaceReadGroups \
  I="$OUTDIR/bam/${TUMOR_SAMPLE}_sorted.bam" \
  O="$OUTDIR/bam/${TUMOR_SAMPLE}_rg.bam" \
  RGID=${TUMOR_SAMPLE}_rg RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=$TUMOR_SAMPLE \
  SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

# Mark optical / PCR duplicates.
picard MarkDuplicates \
  I="$OUTDIR/bam/${TUMOR_SAMPLE}_rg.bam" \
  O="$OUTDIR/bam/${TUMOR_SAMPLE}_dedup.bam" \
  M="$OUTDIR/bam/${TUMOR_SAMPLE}_dedup_metrics.txt" \
  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

# -----------------------------------------------------------------------------
# 3. Normal: align -> add read groups -> mark duplicates
# -----------------------------------------------------------------------------
# Same three-step preprocessing applied to the matched normal sample.
bwa mem -M -t "$THREADS" "$REF_GENOME" "$NORMAL_FASTQ1" "$NORMAL_FASTQ2" | \
  samtools view -Sb - | \
  samtools sort -@ "$THREADS" -o "$OUTDIR/bam/${NORMAL_SAMPLE}_sorted.bam"

picard AddOrReplaceReadGroups \
  I="$OUTDIR/bam/${NORMAL_SAMPLE}_sorted.bam" \
  O="$OUTDIR/bam/${NORMAL_SAMPLE}_rg.bam" \
  RGID=${NORMAL_SAMPLE}_rg RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=$NORMAL_SAMPLE \
  SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

picard MarkDuplicates \
  I="$OUTDIR/bam/${NORMAL_SAMPLE}_rg.bam" \
  O="$OUTDIR/bam/${NORMAL_SAMPLE}_dedup.bam" \
  M="$OUTDIR/bam/${NORMAL_SAMPLE}_dedup_metrics.txt" \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT

# -----------------------------------------------------------------------------
# 4. Somatic variant calling with Mutect2 (paired tumour + normal)
# -----------------------------------------------------------------------------
gatk Mutect2 \
  -R "$REF_GENOME" \
  -I "$OUTDIR/bam/${TUMOR_SAMPLE}_dedup.bam"  \
  -I "$OUTDIR/bam/${NORMAL_SAMPLE}_dedup.bam"  \
  -normal "$NORMAL_SAMPLE" \
  -tumor  "$TUMOR_SAMPLE" \
  -O "$OUTDIR/vcf/${TUMOR_SAMPLE}_unfiltered.vcf.gz"

# -----------------------------------------------------------------------------
# 5. Filter raw Mutect2 calls
# -----------------------------------------------------------------------------
gatk FilterMutectCalls \
  -V "$OUTDIR/vcf/${TUMOR_SAMPLE}_unfiltered.vcf.gz" \
  -R "$REF_GENOME" \
  -O "$OUTDIR/vcf/${TUMOR_SAMPLE}_filtered.vcf.gz"

# -----------------------------------------------------------------------------
# 6. Annotate filtered VCF to MAF (Funcotator)
# -----------------------------------------------------------------------------
# The output MAF includes the `ref_context` field (reference allele +/- 10 nt)
# that HLRCC_Exome_Analysis.r uses to build the trinucleotide context for
# mutational signature extraction.
gatk Funcotator \
  -R "$REF_GENOME" \
  -V "$OUTDIR/vcf/${TUMOR_SAMPLE}_filtered.vcf.gz" \
  -O "$OUTDIR/maf/${TUMOR_SAMPLE}_annotated.maf" \
  --output-file-format MAF \
  --data-sources-path "$FUNCOTATOR_DATA_SOURCES" \
  --ref-version "$FUNCOTATOR_REF_VERSION"