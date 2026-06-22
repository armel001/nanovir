# nanovir 🧬

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bash](https://img.shields.io/badge/Shell-Bash-blue.svg)](https://www.gnu.org/software/bash/)
[![ONT](https://img.shields.io/badge/Sequencing-Oxford%20Nanopore-important.svg)](https://nanoporetech.com/)

**Auteur :** Thibaut Armel Cherif Gnimadi — [@armel001](https://github.com/armel001)

Pipeline Bash pour l'analyse de génomes viraux depuis des données Oxford Nanopore.

---

## Ce que fait nanovir

```
Reads FASTQ
    │
    ├─ [1] QC           NanoStat   — statistiques qualité des reads
    ├─ [2] Filtrage     Chopper    — filtre par longueur / qualité
    ├─ [3] Alignement   minimap2   — mapping sur le génome de référence
    ├─ [4] Couverture   mosdepth   — profondeur par base, masque les zones froides
    ├─ [5] Polissage    Medaka     — correction des erreurs Nanopore
    ├─ [6] Variants     bcftools   — appel SNPs / indels
    └─ [7] Consensus    bcftools   — séquence consensus finale (N = zones masquées)
```

Compatible avec tous les flowcells (R9.4.1, R10.4.1…) et tous les virus ayant un génome de référence disponible.

---

## Installation

### 1. Cloner le dépôt

```bash
git clone https://github.com/armel001/nanovir.git
cd nanovir
```

### 2. Créer l'environnement conda

```bash
conda env create -f environment.yml
conda activate nanovir
```

### 3. Installer Medaka (via pip)

```bash
conda activate nanovir

pip install medaka               # CPU uniquement
pip install "medaka[cuda11]"     # GPU avec CUDA 11
pip install "medaka[cuda12]"     # GPU avec CUDA 12
```

Vérifier :

```bash
medaka --version
medaka tools list_models
```

### 4. Rendre le script exécutable

```bash
chmod +x nanovir.sh
```

---

## Utilisation

```bash
./nanovir.sh -s SAMPLE -r READS -R REFERENCE -m MODEL [options]
```

### Paramètres

| Paramètre | Description | Défaut |
|-----------|-------------|--------|
| `-s` | Nom du sample | *(requis)* |
| `-r` | Fichier FASTQ ou FASTQ.GZ | *(requis)* |
| `-R` | Génome de référence FASTA | *(requis)* |
| `-m` | Modèle Medaka | *(requis sauf si `--skip-medaka`)* |
| `-t` | Nombre de threads | `8` |
| `-o` | Répertoire de sortie | `results_SAMPLE` |
| `-c` | Fichier de configuration | — |
| `--min-depth` | Profondeur min pour les variants et le masquage | `10` |
| `--min-qual` | Score QUAL min pour les variants | `20` |
| `--min-len` | Longueur min des reads (Chopper) | `200 bp` |
| `--skip-qc` | Passer l'étape NanoStat | — |
| `--skip-medaka` | Passer Medaka (plus rapide, moins précis) | — |
| `-h` | Aide | — |
| `-v` | Version | — |

---

## Modèles Medaka

| Flowcell | Mode | Modèle |
|----------|------|--------|
| R9.4.1 | Fast | `r941_min_fast_g507` |
| R9.4.1 | HAC | `r941_min_hac_g507` |
| R9.4.1 | SUP | `r941_min_sup_g507` |
| R10.4.1 | Fast | `r1041_e82_400bps_fast_variant_g632` |
| R10.4.1 | HAC | `r1041_e82_400bps_hac_variant_g632` |
| R10.4.1 | SUP | `r1041_e82_400bps_sup_variant_g632` |

> Utilisez le même mode que votre basecaller (Dorado/Guppy). SUP = meilleure précision.

```bash
medaka tools list_models   # lister tous les modèles installés
```

---

## Exemples

### Monkeypox (MPXV)

```bash
./nanovir.sh \
  -s MPXV_UK_001 \
  -r data/Ech001.fastq.gz \
  -R refs/Monkeypox_MT903345.fasta \
  -m r1041_e82_400bps_hac_variant_g632 \
  -t 16
```

### SARS-CoV-2 — paramètres stricts

```bash
./nanovir.sh \
  -s SARS2_001 \
  -r data/sample.fastq.gz \
  -R refs/SARS-CoV-2_MN908947.fasta \
  -m r1041_e82_400bps_sup_variant_g632 \
  -t 24 \
  --min-depth 30 \
  --min-qual 30
```

### Sans Medaka (mode rapide)

```bash
./nanovir.sh -s Flu_H3N2 -r data/flu.fastq.gz -R refs/h3n2.fasta \
  -m r1041_e82_400bps_hac_variant_g632 --skip-medaka -t 8
```

### Avec un fichier de configuration

```bash
cp config/config.example.sh config/mon_sample.sh
# Éditer config/mon_sample.sh
./nanovir.sh -c config/mon_sample.sh
```

### Boucle multi-samples

```bash
for fq in data/*.fastq.gz; do
  sample=$(basename "$fq" .fastq.gz)
  ./nanovir.sh -s "$sample" -r "$fq" -R refs/ref.fasta \
    -m r1041_e82_400bps_hac_variant_g632 -t 16
done
```

---

## Résultats

```
results_SAMPLE/
├── logs/                    Log complet horodaté
├── qc/                      Rapport NanoStat
├── trimmed/                 Reads filtrés
├── alignment/               BAM + statistiques d'alignement
├── coverage/                Profondeur par base, zones masquées
├── medaka/                  Consensus Medaka
├── variants/                VCF brut, normalisé, filtré
├── consensus/
│   └── SAMPLE_consensus.fa  ← Consensus final (N = faible couverture)
└── summary_report.txt       ← Rapport de synthèse
```

---

## Dépendances

| Outil | Rôle | Version min |
|-------|------|-------------|
| minimap2 | Alignement | ≥ 2.26 |
| samtools | Traitement BAM | ≥ 1.17 |
| bcftools | Variants & consensus | ≥ 1.17 |
| mosdepth | Couverture | ≥ 0.3.8 |
| NanoStat | QC reads | ≥ 1.6 |
| chopper | Filtrage reads | ≥ 0.7 |
| medaka | Polissage consensus | ≥ 1.11 |

---

## Citation

Si vous utilisez nanovir dans une publication, citez les outils sous-jacents :

- **minimap2** : Li H. (2018) *Bioinformatics* 34:3094
- **samtools/bcftools** : Danecek P. et al. (2021) *GigaScience* 10:giab008
- **Medaka** : Oxford Nanopore Technologies — github.com/nanoporetech/medaka
- **mosdepth** : Pedersen & Quinlan (2018) *Bioinformatics* 34:867
- **NanoStat** : De Coster et al. (2018) *Bioinformatics* 34:2666
- **chopper** : De Coster & Rademakers (2023) *Bioinformatics* 39:btad311

---

## Licence

MIT — voir [LICENSE](LICENSE)

---

## Contact

**Thibaut Armel Cherif Gnimadi**
GitHub : [github.com/armel001](https://github.com/armel001)
# nanovir
