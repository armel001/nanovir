#!/bin/bash
# nanovir — fichier de configuration exemple
# Copier ce fichier et l'adapter pour chaque sample :
#   cp config/config.example.sh config/mon_sample.sh
#   ./nanovir.sh -c config/mon_sample.sh

# ── Requis ───────────────────────────────────────────────────
SAMPLE="MPXV_UK_001"
READS="data/Ech001.fastq.gz"
REF="refs/Monkeypox_MT903345.fasta"
MODEL_MEDAKA="r1041_e82_400bps_hac_variant_g632"

# ── Performance ──────────────────────────────────────────────
THREADS=16

# ── Sortie ───────────────────────────────────────────────────
OUTDIR="results/${SAMPLE}"

# ── Filtres ──────────────────────────────────────────────────
MIN_DEPTH=10     # Bases < cette profondeur → masquées en N dans le consensus
MIN_QUAL=20      # Variants < ce score QUAL → supprimés
MIN_LEN=200      # Reads < cette longueur (bp) → supprimés par Chopper

# ── Étapes optionnelles ──────────────────────────────────────
SKIP_QC=false       # true = passer NanoStat
SKIP_MEDAKA=false   # true = passer Medaka (plus rapide, moins précis)
