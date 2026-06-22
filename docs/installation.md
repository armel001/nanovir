# Guide d'installation — nanovir

## 1. Prérequis

- Linux (Ubuntu 20.04+)
- [Conda ou Mamba](https://github.com/conda-forge/miniforge)
- GPU NVIDIA optionnel (pour Medaka GPU)

Installer Miniforge si nécessaire :

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh
source ~/.bashrc
```

---

## 2. Cloner et installer

```bash
git clone https://github.com/your-org/nanovir.git
cd nanovir

# Créer l'environnement (mamba est plus rapide que conda)
mamba env create -f environment.yml
conda activate nanovir
```

---

## 3. Installer Medaka

Medaka s'installe séparément via pip dans l'environnement actif :

```bash
conda activate nanovir

# CPU uniquement
pip install medaka

# GPU CUDA 11 (V100, A100...)
pip install "medaka[cuda11]"

# GPU CUDA 12
pip install "medaka[cuda12]"
```

Vérification :

```bash
medaka --version
medaka tools list_models
```

---

## 4. Rendre le script exécutable

```bash
chmod +x nanovir.sh
./nanovir.sh --help
```

---

## Mise à jour

```bash
conda activate nanovir
conda env update -f environment.yml --prune
pip install --upgrade medaka
```

---

## Cluster HPC (SLURM)

```bash
#!/bin/bash
#SBATCH --job-name=nanovir
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --gres=gpu:1     # supprimer pour CPU uniquement

module load conda
conda activate nanovir

./nanovir.sh \
  -s "$SAMPLE" \
  -r "$READS" \
  -R "$REF" \
  -m r1041_e82_400bps_hac_variant_g632 \
  -t "$SLURM_CPUS_PER_TASK"
```

---

## Problèmes courants

**`medaka: command not found`**
```bash
conda activate nanovir && pip install medaka
```

**`Error: Medaka model not found`**
```bash
medaka tools list_models   # choisir le modèle correspondant
```

**Très peu de reads mappés (< 50%)**
→ Vérifier que la référence correspond bien au virus/souche séquencé.

**Toutes les positions masquées en N**
→ La couverture est inférieure à `--min-depth`. Baisser ce seuil ou vérifier la référence.

**Erreur Medaka BAM introuvable**
```bash
cat results_SAMPLE/logs/medaka.log
```
