#!/bin/bash
# ============================================================
#  nanovir — Nanopore Virus Analysis Pipeline
#  Version: 1.0.0  |  License: MIT
#  Author:  Thibaut Armel Cherif Gnimadi
#  GitHub:  https://github.com/armel001/nanovir
# ============================================================

set -euo pipefail

VERSION="1.0.0"

# Colors
GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BOLD='\033[1m'; NC='\033[0m'

log()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] $1${NC}"; }
warn()  { echo -e "${YELLOW}[$(date '+%H:%M:%S')] WARNING: $1${NC}"; }
error() { echo -e "${RED}[$(date '+%H:%M:%S')] ERROR: $1${NC}" >&2; exit 1; }

# ── DEFAULTS ─────────────────────────────────────────────────
SAMPLE=""; READS=""; REF=""; MODEL_MEDAKA=""
THREADS=8; OUTDIR=""; CONFIG_FILE=""
MIN_DEPTH=10; MIN_QUAL=20; MIN_LEN=200; MAX_LEN=0
PRIMERS=""     # BED file of primer positions (amplicon mode)
SKIP_QC=false; SKIP_MEDAKA=false

# ── USAGE ────────────────────────────────────────────────────
usage() {
cat << EOF
${BOLD}nanovir v${VERSION}${NC} — Nanopore Virus Analysis Pipeline

${BOLD}Usage:${NC}
  $(basename "$0") -s SAMPLE -r READS -R REFERENCE -m MODEL [options]
  $(basename "$0") -c config.sh

${BOLD}Required:${NC}
  -s, --sample    Sample name
  -r, --reads     Input FASTQ / FASTQ.GZ
  -R, --ref       Reference genome FASTA
  -m, --model     Medaka model (ignored if --skip-medaka)

${BOLD}Options:${NC}
  -t, --threads   CPU threads [default: 8]
  -o, --outdir    Output directory [default: results_SAMPLE]
  -c, --config    Config file (see config/config.example.sh)
  --min-depth     Min depth for variant calling & masking [default: 10]
  --min-qual      Min variant quality score [default: 20]
  --min-len       Min read length for Chopper [default: 200]
  --max-len       Max read length for Chopper [default: 0 = disabled]
                  ⚠ Recommended for amplicon sequencing (set close to amplicon size)
  --primers       BED file of primer positions for soft-clipping (amplicon mode)
                  Format: chrom  start  end  name  score  strand
  --skip-qc       Skip NanoStat QC
  --skip-medaka   Skip Medaka polishing
  -h, --help      Show this help
  -v, --version   Show version

${BOLD}Amplicon mode:${NC}
  For amplicon sequencing, always set --max-len near your amplicon size and
  provide --primers to remove primer sequences from the alignments.
  Example (ARTIC ~400 bp amplicons):
    $(basename "$0") -s SARS2 -r reads.fastq.gz -R refs/sars2.fasta \\
      -m r1041_e82_400bps_hac_variant_g632 \\
      --min-len 300 --max-len 500 --primers refs/artic_primers.bed

${BOLD}Medaka models:${NC}
  R9.4.1  Fast/HAC/SUP : r941_min_fast_g507 / r941_min_hac_g507 / r941_min_sup_g507
  R10.4.1 Fast/HAC/SUP : r1041_e82_400bps_fast_variant_g632
                          r1041_e82_400bps_hac_variant_g632
                          r1041_e82_400bps_sup_variant_g632
  Run: medaka tools list_models

${BOLD}Examples:${NC}
  $(basename "$0") -s MPXV_001 -r reads.fastq.gz -R refs/mpxv.fasta -m r1041_e82_400bps_hac_variant_g632
  $(basename "$0") -c config/mpxv_sample.sh
  $(basename "$0") -s SARS2 -r reads.fastq.gz -R refs/sars2.fasta -m r1041_e82_400bps_sup_variant_g632 -t 24

EOF
  exit 0
}

# ── ARGUMENT PARSING ─────────────────────────────────────────
[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--sample)      SAMPLE="$2";       shift 2 ;;
    -r|--reads)       READS="$2";        shift 2 ;;
    -R|--ref)         REF="$2";          shift 2 ;;
    -m|--model)       MODEL_MEDAKA="$2"; shift 2 ;;
    -t|--threads)     THREADS="$2";      shift 2 ;;
    -o|--outdir)      OUTDIR="$2";       shift 2 ;;
    -c|--config)      CONFIG_FILE="$2";  shift 2 ;;
    --min-depth)      MIN_DEPTH="$2";    shift 2 ;;
    --min-qual)       MIN_QUAL="$2";     shift 2 ;;
    --min-len)        MIN_LEN="$2";      shift 2 ;;
    --max-len)        MAX_LEN="$2";      shift 2 ;;
    --primers)        PRIMERS="$2";      shift 2 ;;
    --skip-qc)        SKIP_QC=true;      shift ;;
    --skip-medaka)    SKIP_MEDAKA=true;  shift ;;
    -h|--help)        usage ;;
    -v|--version)     echo "nanovir v${VERSION}"; exit 0 ;;
    *) error "Unknown option: $1 (use -h for help)" ;;
  esac
done

# ── CONFIG FILE ──────────────────────────────────────────────
if [[ -n "$CONFIG_FILE" ]]; then
  [[ -f "$CONFIG_FILE" ]] || error "Config file not found: $CONFIG_FILE"
  # shellcheck source=/dev/null
  source "$CONFIG_FILE"
fi

# ── VALIDATION ───────────────────────────────────────────────
[[ -z "$SAMPLE" ]]      && error "Sample name required (-s)"
[[ -z "$READS" ]]       && error "Reads file required (-r)"
[[ -z "$REF" ]]         && error "Reference FASTA required (-R)"
[[ -f "$READS" ]]       || error "Reads file not found: $READS"
[[ -f "$REF" ]]         || error "Reference not found: $REF"
[[ "$SKIP_MEDAKA" == false && -z "$MODEL_MEDAKA" ]] \
  && error "Medaka model required (-m). Use --skip-medaka to bypass."

[[ -n "$PRIMERS" && ! -f "$PRIMERS" ]] && error "Primers BED file not found: $PRIMERS"

[[ -z "$OUTDIR" ]] && OUTDIR="results_${SAMPLE}"

# ── DEPENDENCY CHECK ─────────────────────────────────────────
check_deps() {
  local tools=("minimap2" "samtools" "bcftools" "mosdepth")
  [[ "$SKIP_QC"     == false ]] && tools+=("NanoStat")
  [[ "$SKIP_MEDAKA" == false ]] && tools+=("medaka_consensus")
  tools+=("chopper")

  local missing=0
  for t in "${tools[@]}"; do
    command -v "$t" &>/dev/null || { warn "Missing tool: $t"; ((missing++)); }
  done
  [[ $missing -gt 0 ]] && error "$missing tool(s) not found. See docs/installation.md"
}

# ── BANNER ───────────────────────────────────────────────────
print_banner() {
  echo ""
  echo -e "${BOLD}╔══════════════════════════════════╗${NC}"
  echo -e "${BOLD}║        nanovir  v${VERSION}           ║${NC}"
  echo -e "${BOLD}╚══════════════════════════════════╝${NC}"
  printf "  %-16s %s\n" "Sample:"    "$SAMPLE"
  printf "  %-16s %s\n" "Reads:"     "$READS"
  printf "  %-16s %s\n" "Reference:" "$REF"
  printf "  %-16s %s\n" "Model:"     "${MODEL_MEDAKA:-N/A}"
  printf "  %-16s %s\n" "Threads:"   "$THREADS"
  printf "  %-16s %s\n" "Output:"    "$OUTDIR"
  local len_range="≥${MIN_LEN} bp"
  [[ "$MAX_LEN" -gt 0 ]] && len_range="${MIN_LEN}–${MAX_LEN} bp"
  printf "  %-16s depth≥%s  QUAL≥%s  len %s\n" "Filters:" "$MIN_DEPTH" "$MIN_QUAL" "$len_range"
  [[ -n "$PRIMERS" ]] && printf "  %-16s %s\n" "Primers BED:" "$PRIMERS"
  echo ""
}

# ═══════════════════════════════════════════════════════════════
# PIPELINE STEPS
# ═══════════════════════════════════════════════════════════════

# ── STEP 1: QC ───────────────────────────────────────────────
step_qc() {
  if [[ "$SKIP_QC" == true ]]; then
    log "[1/8] QC skipped"
    return
  fi
  log "[1/8] Quality control (NanoStat)..."
  NanoStat --fastq "$READS" \
           --outdir "${OUTDIR}/qc" \
           --name "${SAMPLE}.txt" \
           --threads "$THREADS" \
           2>>"${OUTDIR}/logs/nanostat.log"
}

# ── STEP 2: TRIMMING ─────────────────────────────────────────
step_trim() {
  local len_info="min ${MIN_LEN} bp"
  [[ "$MAX_LEN" -gt 0 ]] && len_info="${MIN_LEN}–${MAX_LEN} bp"
  log "[2/8] Read filtering (Chopper, length: ${len_info})..."
  TRIMMED_READS="${OUTDIR}/trimmed/${SAMPLE}.fastq.gz"

  local chopper_args="--minlength $MIN_LEN --quality 8 --threads $THREADS"
  [[ "$MAX_LEN" -gt 0 ]] && chopper_args="$chopper_args --maxlength $MAX_LEN"

  if [[ "$READS" == *.gz ]]; then gunzip -c "$READS"
  else cat "$READS"; fi \
  | chopper $chopper_args \
            2>>"${OUTDIR}/logs/chopper.log" \
  | gzip > "$TRIMMED_READS"
}

# ── STEP 3: ALIGNMENT ────────────────────────────────────────
step_align() {
  log "[3/8] Alignment (minimap2)..."
  ALIGNED_BAM="${OUTDIR}/alignment/${SAMPLE}.bam"

  minimap2 -ax map-ont -t "$THREADS" --secondary=no \
           -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ONT" \
           "$REF" "$TRIMMED_READS" \
           2>>"${OUTDIR}/logs/minimap2.log" \
  | samtools view -F 4 -b -@ "$THREADS" \
  | samtools sort -@ "$THREADS" -o "$ALIGNED_BAM"

  samtools index -@ "$THREADS" "$ALIGNED_BAM"
  samtools flagstat "$ALIGNED_BAM" > "${OUTDIR}/alignment/${SAMPLE}_flagstat.txt"

  local mapped total
  mapped=$(grep "reads mapped:" "${OUTDIR}/alignment/${SAMPLE}_flagstat.txt" | awk '{print $1}')
  total=$(grep  "in total"      "${OUTDIR}/alignment/${SAMPLE}_flagstat.txt" | awk '{print $1}')
  log "    → ${mapped} / ${total} reads mapped"
}

# ── STEP 4: PRIMER CLIPPING (amplicon mode only) ──────────────
step_clip_primers() {
  if [[ -z "$PRIMERS" ]]; then
    log "[4/8] Primer clipping skipped (no --primers BED provided)"
    return
  fi

  log "[4/8] Clipping primer sequences (samtools ampliconclip)..."
  local clipped="${OUTDIR}/alignment/${SAMPLE}_clipped.bam"

  samtools ampliconclip \
    -b "$PRIMERS" \
    --soft-clip \
    --both-ends \
    -@ "$THREADS" \
    -o - \
    "$ALIGNED_BAM" \
    2>>"${OUTDIR}/logs/ampliconclip.log" \
  | samtools sort -@ "$THREADS" -o "$clipped"

  samtools index -@ "$THREADS" "$clipped"

  # Replace the main BAM with the clipped version for all downstream steps
  ALIGNED_BAM="$clipped"
  log "    → Primers soft-clipped. BAM updated: $clipped"
}

# ── STEP 5: COVERAGE ─────────────────────────────────────────
step_coverage() {
  log "[5/8] Coverage analysis (mosdepth)..."

  mosdepth --threads "$THREADS" --no-abbrev \
           "${OUTDIR}/coverage/${SAMPLE}" \
           "$ALIGNED_BAM" \
           2>>"${OUTDIR}/logs/mosdepth.log"

  samtools depth -a "$ALIGNED_BAM" > "${OUTDIR}/coverage/${SAMPLE}_depth.txt"
  samtools coverage "$ALIGNED_BAM"  > "${OUTDIR}/coverage/${SAMPLE}_coverage.txt"

  # BED of positions below MIN_DEPTH (used for masking)
  awk -v d="$MIN_DEPTH" '$3 < d {print $1"\t"($2-1)"\t"$2}' \
      "${OUTDIR}/coverage/${SAMPLE}_depth.txt" \
      > "${OUTDIR}/coverage/${SAMPLE}_low_cov.bed"

  local mean pct masked
  mean=$(awk 'NR>1{s+=$7;n++} END{printf "%.1f",s/n}' "${OUTDIR}/coverage/${SAMPLE}_coverage.txt")
  pct=$(awk  'NR>1{s+=$6;n++} END{printf "%.1f",s/n}' "${OUTDIR}/coverage/${SAMPLE}_coverage.txt")
  masked=$(wc -l < "${OUTDIR}/coverage/${SAMPLE}_low_cov.bed")
  log "    → mean depth: ${mean}x | breadth: ${pct}% | bases to mask: ${masked}"
}

# ── STEP 6: MEDAKA ───────────────────────────────────────────
step_medaka() {
  if [[ "$SKIP_MEDAKA" == true ]]; then
    log "[6/8] Medaka skipped — using alignment BAM for variant calling"
    VARIANT_BAM="$ALIGNED_BAM"
    return
  fi
  log "[6/8] Consensus polishing (Medaka)..."

  medaka_consensus -i "$TRIMMED_READS" \
                   -d "$REF" \
                   -o "${OUTDIR}/medaka" \
                   -t "$THREADS" \
                   -m "$MODEL_MEDAKA" \
                   2>>"${OUTDIR}/logs/medaka.log"

  VARIANT_BAM="${OUTDIR}/medaka/calls_to_draft.bam"
  [[ -f "$VARIANT_BAM" ]] || error "Medaka BAM not found — check ${OUTDIR}/logs/medaka.log"
}

# ── STEP 7: VARIANT CALLING ──────────────────────────────────
step_variants() {
  log "[7/8] Variant calling (bcftools)..."

  # Raw calls (ploidy 1 = haploid/viral)
  bcftools mpileup -f "$REF" --max-depth 50000 --min-BQ 20 \
                   -a FORMAT/AD,FORMAT/DP -Ou "$VARIANT_BAM" \
                   2>>"${OUTDIR}/logs/bcftools.log" \
  | bcftools call -mv --ploidy 1 -Oz \
                  -o "${OUTDIR}/variants/${SAMPLE}_raw.vcf.gz" \
                  2>>"${OUTDIR}/logs/bcftools.log"
  bcftools index "${OUTDIR}/variants/${SAMPLE}_raw.vcf.gz"

  # Normalize (left-align indels)
  bcftools norm -f "$REF" -Oz \
                "${OUTDIR}/variants/${SAMPLE}_raw.vcf.gz" \
                -o "${OUTDIR}/variants/${SAMPLE}_norm.vcf.gz" \
                2>>"${OUTDIR}/logs/bcftools.log"
  bcftools index "${OUTDIR}/variants/${SAMPLE}_norm.vcf.gz"

  # Filter
  bcftools filter -e "QUAL<${MIN_QUAL} || DP<${MIN_DEPTH}" -Oz \
                  "${OUTDIR}/variants/${SAMPLE}_norm.vcf.gz" \
                  -o "${OUTDIR}/variants/${SAMPLE}_filtered.vcf.gz" \
                  2>>"${OUTDIR}/logs/bcftools.log"
  bcftools index "${OUTDIR}/variants/${SAMPLE}_filtered.vcf.gz"

  local snps indels
  snps=$(bcftools   view -v snps   -H "${OUTDIR}/variants/${SAMPLE}_filtered.vcf.gz" 2>/dev/null | wc -l)
  indels=$(bcftools view -v indels -H "${OUTDIR}/variants/${SAMPLE}_filtered.vcf.gz" 2>/dev/null | wc -l)
  log "    → SNPs: ${snps} | Indels: ${indels} (after filtering)"
}

# ── STEP 8: CONSENSUS ────────────────────────────────────────
step_consensus() {
  log "[8/8] Generating consensus FASTA..."

  local out_fa="${OUTDIR}/consensus/${SAMPLE}_consensus.fa"
  local low_cov="${OUTDIR}/coverage/${SAMPLE}_low_cov.bed"
  local mask_args=""

  if [[ -s "$low_cov" ]]; then
    mask_args="--mask ${low_cov} --mask-with N"
  fi

  # shellcheck disable=SC2086
  bcftools consensus -f "$REF" --mark-del '-' $mask_args \
                     "${OUTDIR}/variants/${SAMPLE}_filtered.vcf.gz" \
                     > "$out_fa" 2>>"${OUTDIR}/logs/bcftools.log"

  # Rename FASTA header
  sed -i "1s/.*/>$SAMPLE/" "$out_fa"
  log "    → ${OUTDIR}/consensus/${SAMPLE}_consensus.fa"
}

# ── SUMMARY REPORT ───────────────────────────────────────────
write_report() {
  local report="${OUTDIR}/summary_report.txt"
  {
    echo "nanovir v${VERSION} — Summary Report"
    echo "Date   : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Sample : $SAMPLE"
    echo ""
    echo "── Inputs ────────────────────────────────────"
    echo "  Reads     : $READS"
    echo "  Reference : $REF"
    echo "  Model     : ${MODEL_MEDAKA:-N/A (skipped)}"
    echo ""
    echo "── Parameters ────────────────────────────────"
    echo "  Threads   : $THREADS"
    echo "  Min depth : $MIN_DEPTH"
    echo "  Min QUAL  : $MIN_QUAL"
    echo "  Min length: $MIN_LEN bp"
    echo ""
    echo "── Alignment ─────────────────────────────────"
    [[ -f "${OUTDIR}/alignment/${SAMPLE}_flagstat.txt" ]] \
      && cat "${OUTDIR}/alignment/${SAMPLE}_flagstat.txt"
    echo ""
    echo "── Coverage ──────────────────────────────────"
    [[ -f "${OUTDIR}/coverage/${SAMPLE}_coverage.txt" ]] \
      && cat "${OUTDIR}/coverage/${SAMPLE}_coverage.txt"
    echo ""
    echo "── Variants (filtered) ───────────────────────"
    local snps indels
    snps=$(bcftools   view -v snps   -H "${OUTDIR}/variants/${SAMPLE}_filtered.vcf.gz" 2>/dev/null | wc -l)
    indels=$(bcftools view -v indels -H "${OUTDIR}/variants/${SAMPLE}_filtered.vcf.gz" 2>/dev/null | wc -l)
    echo "  SNPs   : $snps"
    echo "  Indels : $indels"
    echo ""
    echo "── Output files ──────────────────────────────"
    echo "  Consensus : ${OUTDIR}/consensus/${SAMPLE}_consensus.fa"
    echo "  VCF       : ${OUTDIR}/variants/${SAMPLE}_filtered.vcf.gz"
    echo "  BAM       : ${OUTDIR}/alignment/${SAMPLE}.bam"
    echo "  Log       : $LOG_FILE"
    echo "──────────────────────────────────────────────"
  } > "$report"

  echo ""
  cat "$report"
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
main() {
  local t0=$SECONDS

  print_banner
  check_deps

  # Create output dirs
  mkdir -p "${OUTDIR}"/{logs,qc,trimmed,alignment,coverage,medaka,variants,consensus}
  LOG_FILE="${OUTDIR}/logs/nanovir_$(date '+%Y%m%d_%H%M%S').log"
  exec > >(tee -a "$LOG_FILE") 2>&1

  # Declare shared state variables (set inside step functions)
  TRIMMED_READS=""
  ALIGNED_BAM=""
  VARIANT_BAM=""

  step_qc
  step_trim
  step_align
  step_clip_primers
  step_coverage
  step_medaka
  step_variants
  step_consensus
  write_report

  local elapsed=$(( SECONDS - t0 ))
  echo ""
  echo -e "${BOLD}${GREEN}Done in $(( elapsed/60 ))m $(( elapsed%60 ))s — results in ${OUTDIR}/${NC}"
}

main "$@"
