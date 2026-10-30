# cannabinoid-synthase-workflow

[![DOI](https://zenodo.org/badge/1439704713.svg)](https://doi.org/10.5281/zenodo.23806790)

A reproducible, multi-level (L1–L7) in silico workflow for identifying functional
copies of the cannabinoid acid synthases THCAS, CBDAS and CBCAS in
haplotype-resolved *Cannabis sativa* genome assemblies.

The synthase loci on chromosome 7 contain many partial, pseudogenized and
divergent copies, and annotation-based and genome-based evidence often disagree.
This workflow evaluates each synthase-like copy against seven criteria, resolves
conflicting evidence by rules fixed before the analysis, and assigns every copy
to one of five classes, recording why each copy was or was not selected as a
functional candidate.

> **Scope.** The workflow is a sequence- and genome-based screen. It does not
> demonstrate biochemical activity. "Functional candidate" means only that a copy
> satisfies all of L1–L7 and is prioritised for experimental validation.

This code accompanies:

> Akahoshi Y, et al. Manuscript submitted to *GigaByte* (2026). The full citation will be added on publication.

## Requirements

- Linux (tested on Ubuntu 26.04 LTS under WSL2)
- Conda
- Python 3.10.14, BLAST+ 2.16.0, MAFFT 7.526, IQ-TREE 3.1.2, Biopython 1.84,
  Matplotlib 3.10.9 (all pinned in `environment.yml`)

```bash
conda env create -f environment.yml
conda activate plant_genomics
```

SignalP 6.0, used in the paper for signal-peptide prediction, is an external web
service and is not part of the automated workflow.

## Quick start: example dataset

```bash
bash scripts/run_example.sh
```

The example dataset (`example_data/`) is a 1.2 Mb window of chromosome 7 of
haplotype SAN2b (54.0–55.2 Mb), with every coordinate re-based to 1
(see `example_data/README_offset.txt`). It contains the whole THCAS-series
cluster of this haplotype and runs through the same code path as the full
analysis in a few minutes.

**Expected result:** `example_out/SAN2b/06_classification.tsv` with 11 copies.
SAN2b g313950 is classified as a *Pseudogenized THCAS-like copy* on the evidence
of a premature stop codon in the genomic hit. The 11 classifications are
identical to those obtained with the full SAN2b assembly as input.

To check the result automatically, run `bash tests/test_example.sh`; it reruns
the example and compares every classification with
`tests/expected_SAN2b_classification.tsv`.

The example dataset was built from the full assembly with
`scripts/make_example_dataset.py`. It is an excerpt of the SAN2b assembly,
annotation and predicted proteins of Lynch et al. (2025), which are released
under CC0 1.0 on Figshare+; please cite the original data when using it.

## Workflow

| Step | Script | Evaluates | Main output |
|---|---|---|---|
| 01 | `01_search_synthase_candidates.sh` | BLASTP / tBLASTn search (E ≤ 1e-10) | `blastp_hits.tsv`, `tblastn_hits.tsv` |
| 02 | `02_extract_gene_models.py` | Genomic hits → overlapping GFF3 gene models | `02_gene_models.tsv` |
| 03 | `03_evaluate_orf.py` | L2, L6 from the annotated protein; MAFFT alignment | `03_orf.tsv`, `03_candidates.aligned.fasta` |
| 03g | `03g_genomic_orf_check.py` | L2, L6 from the genome (evidence-priority rule) | `03g_genomic_orf.tsv` |
| 04 | `04_check_functional_residues.py` | L3 (FAD-related), L4 (active-site) residues | `04_functional_residues.tsv` |
| 05 | `05_run_phylogeny.sh` | ML tree (IQ-TREE, rooted on BBE) | `*.contree` |
| 05b | `05b_assign_clade.py` | L5 (clade membership, UFBoot ≥ 95) | `05b_l5.tsv` |
| 05c | `05c_build_evidence.py` | Joins all levels (no decisions) | `06_evidence.tsv` |
| 06 | `06_assign_L1_L7.py` | Five-class classification | `06_classification.tsv` |

`scripts/run_haplotype.sh` runs steps 03–06 for one haplotype.

### Evaluation levels

| Level | Criterion |
|---|---|
| L1 | BLAST homology to THCAS/CBDAS/CBCAS (identity ≥ 30 %, alignment ≥ 200 aa) |
| L2 | Full-length ORF (≥ 500 aa, no internal gaps) |
| L3 | FAD-related residues conserved (His114, Cys176, His184; THCAS precursor numbering) |
| L4 | Active-site residues conserved (His292, Tyr417, Tyr484) |
| L5 | Smallest clade containing the copy and exactly one reference is supported by UFBoot ≥ 95 |
| L6 | No premature stop codon, frameshift or large deletion (> 20 aa) |
| L7 | Haplotype origin confirmed and locus structure interpretable |

### Classes

| Class | Definition |
|---|---|
| Functional THCAS/CBDAS/CBCAS candidate | Satisfies all of L1–L7 |
| THCAS/CBDAS/CBCAS-like copy | High homology but not all of L1–L7 |
| Divergent …-like copy | Satisfies L2 and L6 but not L4 or L5 |
| Pseudogenized …-like copy | Premature stop or frameshift, or lacks L2, L3 and L4 |
| No detectable … sequence | No significant hit (haplotype-level) |

### Evidence-priority rules (fixed before the analysis)

- L2 and L6 are judged from the annotated protein, except when a single tBLASTn
  HSP covers the reference end to end (≤ 3 residues missing at each terminus),
  is gapless and in frame, and translates without an internal stop; the genomic
  evidence is then used. Multiple HSPs are not joined, start codons are not
  searched for and introns are not predicted.
- The series name (THCAS/CBDAS/CBCAS) comes from the nearest BLAST reference and
  is kept separate from the phylogenetic position (L5).
- Pseudogenization requires direct evidence (a premature stop or frameshift).
  Split CDS or in-frame deletions alone are not used.
- For truncated copies, L2 takes precedence over L5.
- Levels that cannot be evaluated are not forced; such copies do not satisfy all
  of L1–L7.

## Running on the full assemblies

The assemblies, GFF3 annotations and protein FASTA files are those of the
*Cannabis sativa* pangenome (Lynch et al. 2025, *Nature* 643:1001–1010,
doi:10.1038/s41586-025-09065-0), available from the Salk Institute Michael Lab
resources (https://resources.michael.salk.edu/) and Figshare+
(doi:10.25452/figshare.plus.c.7248427.v1). They are not redistributed here.

The scripts work in a project directory `ROOT` (default `~/cbdas_project`)
that holds the input data, with this repository cloned inside it:

```
$ROOT/
  cannabinoid-synthase-workflow/   (this repository)
  00_genomes/<group>/<H>/<H>.softmasked.fasta.gz   (or <H>.unmasked.fasta.gz)
                        <H>.all.proteins.fasta.gz
                        <H>.high_confidence.gff3.gz
  refs/synthase_query.fasta        (copy from this repository)
  refs/BBE_EcaBBE_P30986.fasta     (copy from this repository)
  work/<H>/                        (steps 01-02, created by the commands below)
  results/<H>/                     (steps 03-06, created by run_haplotype.sh)
```

For each haplotype `H` (SODLa, SODLb, EH23a, EH23b, KOMPa, KOMPb, SAN2a, SAN2b),
from `$ROOT`:

```bash
S=cannabinoid-synthase-workflow/scripts
G=$(dirname "$(find 00_genomes -name "$H.all.proteins.fasta.gz")")
bash $S/01_search_synthase_candidates.sh --refs refs/synthase_query.fasta \
  --proteins "$G/$H.all.proteins.fasta.gz" --genome "$G/$H.softmasked.fasta.gz" \
  --outdir "work/$H/01_search" --threads 16
python3 $S/02_extract_gene_models.py --hits "work/$H/01_search/tblastn_hits.tsv" \
  --gff "$G/$H.high_confidence.gff3.gz" --out "work/$H/02_gene_models.tsv"
ROOT="$PWD" bash $S/run_haplotype.sh "$H"
```

Each haplotype takes about 3 minutes for steps 03-06 on a desktop PC.
After all eight haplotypes are done, from `$ROOT`:

```bash
python3 $S/summarize_no_model_hits.py > no_gene_model_regions.tsv
python3 $S/build_additional_file_1.py > Additional_file_1_synthase_copy_detail.tsv
bash $S/07_cross_haplotype_tree.sh      # Fig. 3 tree (sequences gathered by 07a)
```

`no_gene_model_regions.tsv` lists tBLASTn hits that overlap no GFF3 gene model;
these are reported in Additional file 1 but not classified.

The classifications in the paper (Table 4, Fig. 2, Additional file 1) were
regenerated with this version of the code for all eight haplotypes and are
identical to the published results.

## Figures

`figures/` contains the scripts that draw Figs. 1–4 of the paper, together with
their PDF/PNG output and the Fig. 3 tree (`Fig3_cross_haplotype.contree`).

| Figure | Script |
|---|---|
| Fig. 1 Cannabinoid acid biosynthesis and chemotypes | `Fig1_pathway_v11.py` |
| Fig. 2 Chromosome 7 synthase clusters | `Fig2_chr7_locus_map_v25.py` |
| Fig. 3 Maximum-likelihood tree | `Fig3_phylogenetic_tree_v6.py` |
| Fig. 4 Santhica 27 copy architecture | `Fig4_SAN2_disruption_v5.py` |

## Citation

If you use this workflow, please cite the paper above and the archived version
of this repository: doi:10.5281/zenodo.23806790 (all versions; each release also has its own DOI on Zenodo).

## License

MIT License. See `LICENSE`.
