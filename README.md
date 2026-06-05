<p align="center">
  <p align="center">
    <img width=200 align="center" src="./assets/figures/hillerlab.png" >
  </p>

  <span>
    <h1 align="center">
        xorf
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/xorf" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/xorf?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> end-to-end robust and comprehensive ORF prediction pipeline  </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://www.genome.gov/genetics-glossary/Open-Reading-Frame">orf</a> .
        <a href="https://github.com/hillerlab/xorf/blob/main/assets/pipeline/xorf.mermaid">pipeline</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

---

<div align="center">

<pre>
genome       ─────────────────────────────────────────────
ORFs           ATG════TAA     ATG══════════TGA   ATG══TAG

translationAi  █████░░░░      ████████░░░░       ██░░░░░
RNASamba       ████░░░░░      ███████░░░░░       █░░░░░░
TRANSAID       ██████░░░      █████████░░░       ██░░░░░
Netstart2      ███░░░░░░      ██████░░░░░        █░░░░░░
BLAST          ███████░░      ████████░░░        ░░░░░░░

features       └──────┘       └──────────┘       └─────┘
                   │                │                │
                   ▼                ▼                ▼
              GBoost: +        GBoost: +        GBoost: −
</pre>

</div>

---

> [!IMPORTANT]
> - **Softmask both genomes** (lowercase, do NOT hardmask). RepeatModeler 2 per genome is recommended; add WindowMasker if you see runaway LASTZ runtimes. We also provide a soft-masking solution in [softmask](https://github.com/hillerlab/softmask).
> - **Scaffold names**: no spaces; avoid dots (rename `NC_00000.1` → `NC_00000`)
> - Inputs accepted: `.fasta` or `.2bit`.
> - **Container image**: We offer a pre-built container image for the whole pipeline as well as individual modules. By default the pipeline runs with [ghcr.io/hillerlab/xorf:latest](https://github.com/hillerlab/containers/pkgs/container/xorf). Additional images can be found at [containers](https://github.com/hillerlab/containers) and nextflow modules at [core](https://github.com/hillerlab/core).
> - **UCSC replacement**: As of >=3.1.0, the pipeline uses [chaintools](https://github.com/alejandrogzi/chaintools), a Rust library to work with .chain files.

---

## Usage

> [!NOTE]
> Requirements: Nextflow ≥ 25.04.6, Docker or Apptainer, Java.

```bash
git clone https://github.com/hillerlab/softmask.git
cd softmask
```

Edit `params.json` (set `genome`, `assembly_prefix`), then:

```bash
# Docker
nextflow run main.nf -params-file params.json -profile docker

# Apptainer / Singularity
nextflow run main.nf -params-file params.json -profile apptainer
```

Smoke test:
```bash
nextflow run main.nf -profile test,apptainer
```

> [!NOTE]
> You can also specify these options directly in `params.json`.

A helper sh script is provided to run the pipeline on a SLURM cluster. See details below.

<details>
<summary>Click to expand</summary>


Edit the path variables at the top of `assets/hpc/xorf.sh` (cache dir, container image, manifest path), then submit:

```bash
sbatch --array=1-<N> xorf.sh
```

Each array task spawns one Nextflow head job that submits all compute as child SLURM jobs.

REPEATMASKER run as SLURM job arrays. Partition routing, array sizes, and resource tiers are documented inline in `nextflow.config` — edit there to match your cluster.

</details>

---

## Output

```
results/
├── 00_concat/       *bed
├──── 00_concat/raw/ *bed
├── 01_duplicates/   *bed
├── 02_results/      *bed
└── pipeline_info/    timeline, trace, DAG
```

---

## Where to edit

| File | What |
|------|------|
| `params.json` | Genome paths, alignment settings, checkpoints — per run |
| `nextflow.config` | Compute resources, profiles, container, SLURM — rarely |

