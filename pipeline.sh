#!/bin/bash# pipeline.sh
# Author: Nitin Sharma 
# Affiliation: BTech Biotechnology, NIT Allahabad
# Date: 2025-11-16
# Purpose: Full pipeline to align reads, call variants and produce round_1.vcf
# Usage: bash pipeline.sh

set -euo pipefail

# -----------------------
# Config / overrides
# -----------------------
REF="Pfalciparum_3D7.fasta"
THREADS=${THREADS:-8}
PREFIX=${PREFIX:-alignment}
OUTDIR=${OUTDIR:-medaka_out}
VCF=${VCF:-round_1.vcf}

# -----------------------
# Autodetect reads
# -----------------------
# If READ1/READ2 pre-set, use them. Otherwise try common names.
if [ -n "${READ1:-}" ] && [ -n "${READ2:-}" ]; then
  READ1="$READ1"
  READ2="$READ2"
elif [ -f "SRR26194785_1.fastq" ] && [ -f "SRR26194785_2.fastq" ]; then
  READ1="SRR26194785_1.fastq"
  READ2="SRR26194785_2.fastq"
elif [ -f "SRR26194785.fastq" ]; then
  READ1="SRR26194785.fastq"
  READ2=""
else
  # pick any single fastq in directory as fallback
  candidate=$(ls *.fastq 2>/dev/null | head -n1 || true)
  if [ -n "$candidate" ]; then
    READ1="$candidate"
    READ2=""
  else
    echo "No fastq files found in the current folder. Place SRR26194785.fastq or SRR26194785_1/2.fastq here." >&2
    exit 2
  fi
fi

# -----------------------
# Decide mapping preset
# -----------------------
if [ -n "${PRESET:-}" ]; then
  PRESET="$PRESET"
else
  # If paired reads detected, assume short reads (sr); otherwise assume nanopore (map-ont)
  if [ -n "${READ2:-}" ]; then
    PRESET="sr"
  else
    PRESET="map-ont"
  fi
fi

echo "Reference:    $REF"
echo "Reads:        $READ1 ${READ2:+$READ2}"
echo "Threads:      $THREADS"
echo "Preset:       $PRESET"

# -----------------------
# 1) Align reads to reference (minimap2)
# -----------------------
echo "Running minimap2 alignment..."
if [ -n "${READ2:-}" ]; then
  # paired
  minimap2 -ax ${PRESET} "$REF" "$READ1" "$READ2" > ${PREFIX}.sam
else
  # single
  minimap2 -ax ${PRESET} "$REF" "$READ1" > ${PREFIX}.sam
fi

# -----------------------
# 2) Convert, sort, index (samtools)
# -----------------------
echo "Converting/sorting/indexing with samtools..."
samtools view -bS ${PREFIX}.sam > ${PREFIX}.bam
samtools sort -o ${PREFIX}.sorted.bam ${PREFIX}.bam
samtools index ${PREFIX}.sorted.bam

# remove intermediates to save space
rm -f ${PREFIX}.sam ${PREFIX}.bam || true

# -----------------------
# 3) Variant calling
#    - If nanopore single-end (map-ont) and medaka available -> medaka_variant
#    - Else use samtools mpileup + bcftools call
# -----------------------
# Helper: check if medaka_variant exists on PATH
medaka_bin="$(command -v medaka_variant || true)"
bcftools_bin="$(command -v bcftools || true)"

if [ "$PRESET" = "map-ont" ] && [ -n "$medaka_bin" ] && { [ -z "${READ2:-}" ] || [ "${FORCE_MEDAKA:-0}" = "1" ]; }; then
  echo "Running medaka_variant (Nanopore variant caller)..."
  # medaka_variant expects reads (fastq), reference, and output dir
  mkdir -p "$OUTDIR"
  # medaka_variant's CLI varies by version; common invocation:
  medaka_variant -i "$READ1" -f "$REF" -o "$OUTDIR" -t $THREADS || {
    echo "medaka_variant failed. See $OUTDIR for logs. Trying to find a VCF inside $OUTDIR..."
  }

  # find medaka VCF
  if [ -f "${OUTDIR}/variants.vcf" ]; then
    mv "${OUTDIR}/variants.vcf" "$VCF"
  elif [ -f "${OUTDIR}/round_1.vcf" ]; then
    mv "${OUTDIR}/round_1.vcf" "$VCF"
  else
    # try any vcf in output dir
    vcf_candidate=$(ls -1 ${OUTDIR}/*.vcf 2>/dev/null | head -n1 || true)
    if [ -n "$vcf_candidate" ]; then
      cp "$vcf_candidate" "$VCF"
    else
      echo "Could not find a VCF in $OUTDIR after medaka. Listing $OUTDIR:" >&2
      ls -l "$OUTDIR" >&2
      exit 3
    fi
  fi

elif [ -n "$bcftools_bin" ]; then
  echo "Using samtools + bcftools for variant calling (short-read / fallback)."
  # Create raw BCF via mpileup then call SNPs/indels
  # -Ou produces uncompressed BCF in memory
 # use bcftools mpileup (modern workflow)
 # adjust -q / -Q / -C to filter mapping/base quality if you want
    bcftools mpileup -f "$REF" ${PREFIX}.sorted.bam -Ou \
  | bcftools call -mv -Oz -o "${VCF}.gz"

# index and write plain VCF
    bcftools index "${VCF}.gz"
    bcftools view "${VCF}.gz" -Ov -o "$VCF"

  # cleanup
  rm -f "${VCF}.gz" "${VCF}.gz.csi" || true

else
  echo "ERROR: No suitable variant caller found. medaka_variant not found and bcftools not installed." >&2
  echo "Install medaka (for nanopore) or bcftools (for short reads). Example:" >&2
  echo "  conda install -c bioconda medaka" >&2
  echo "  sudo apt install bcftools" >&2
  exit 4
fi

echo "VCF produced: $VCF"
