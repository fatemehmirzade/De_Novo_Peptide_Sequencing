# Characterizing Unexpected Peptides in Human Proteomics Data Using De Novo Peptide Sequencing

our pipeline take the unidentified medoid spectra of MassIVE dataset MSV000088598, run Casanovo on it, filter the predictions, splits them per dataset, removes human peptides and assigns the rest to taxa with Unipept.

| Step | Spectra or peptides |
|------|---------------------|
| medoid spectra in MSV000088598 | 44,697,225 |
| after removing ANN-SoLo annotated spectra | 39,413,436 |
| casanovo predictions (charge 1 to 4) | 35,461,875 |
| after score and length filtering | 2,069,180 |
| after removing human peptides | 814,117 |

all results cover 216 public human proteomics datasets.

## repository layout

| Folder | Content |
|--------|---------|
| `codes/` | main pipeline scripts (steps 1 to 13 below) |
| `text_mining_pipeline/` | metadata extraction from the papers |
| `papers_text/` | text of the papers used as input for text mining |
| `text_mining_results/` | extracted metadata per dataset |
| `human_match/` | BLAST hits against the human proteome, per dataset |
| `filtered_tsv(No_Human)/` | peptides left after removing the human ones |
| `unipept_results/` | taxonomic assignments |
| `terminus/` | N- and C-terminal residue analysis |

## pipeline (Codes/)

run the scripts in this order -> each script has its input and output paths as constants at the top of the file so edit those before running.

1. **[filter_mgf.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/filter_mgf.py)**
   remove the spectra that ANN-SoLo already identified & reads the `spectra_ref` index of every PSM in the mzTab and streams the MGF, keeping only the spectra with no PSM.
   input: `cluster_ident_2.mgf`, `cluster_ident_n.mgf` and the matching mzTab files output: `cluster_unannotated_2.mgf`, `cluster_unannotated_n.mgf`
   
3. **[run_casanovo.sh](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/Run_casanovo.sh)**
   slurm script that runs Casanovo 5.1.2 on both unannotated MGF files with default settings (`casanovo sequence <file>.mgf --model casanovo_v5_0_0.ckpt`). Casanovo only supports precursor charges 1 to 4 and skips the rest -> Output: one mzTab per MGF.

4. **[convert_mztab_tsv.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/convert_mztab_tsv.py)**
   converts the PSM section of each Casanovo mzTab into a TSV file.

5. **[filter_by_casanovo_score.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/filter_by_casanovo_score.py)** keeps a PSM if the Casanovo score is at least 0.6 or between -0.4 and 0.0.

6. **[remove_modifications.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/remove_modifications.py)** removes modifications from the sequences, so `LFM[Oxidation]GK` becomes `LFMGK`.

7. **[filter_by_length.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/filter_by_length.py)** keeps peptides of 10 residues or more. Shorter peptides match too many unrelated proteins.

8. **[split_datasets.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/split_datasets.py)** splits the PSMs per MassIVE dataset. The MSV ID is read from the spectrum title in the MGF, Output: `<MSV ID>.tsv` and `<MSV ID>_matched.mgf`.

9. **[merge_datasets.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/merge_datasets.py)** merges the results of the two cluster size groups (`_2` and `_n`) for each dataset step 7 output has to be placed in the `cluster_ident_2_unannotated_data` and `cluster_ident_n_unannotated_data` folders first.

10. **[extract_sequence.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/extract_sequence.py)** write a sequence only TSV (input for Unipept) and a FASTA file (input for BLAST) per dataset.

11. **[tsv_to_fasta.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/tsv_to_fasta.py)** converts a sequence TSV to FASTA. Useful to rebuild the FASTA files without rerunning step 9.

12. **[terminus_plot.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/terminus_plot.py)** plots the first and last residue of all peptides, as a quality check (Supplementary Figure S1) -> in the paper, K and R make up 91.72% of the C-termini, as expected for trypsin. it uses the step 6 output, so it can run right after step 6.

13. **[filter_homo_by_sequences.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/filter_homo_by_sequences.py)** removes the peptides that matched a human protein, leaving the foreign peptides & it needs `Human_match/<MSV ID>_human_matches.txt`, which comes from a BLAST run of the step 9 FASTA files against the human proteome (`blastp-short`, 100% query coverage, at least 90% identity).

14. **[Unipept_taxon.sh](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Codes/Unipept_taxon.sh)** slurm script running `unipept pept2lca --equate --all` on each dataset. It gives the lowest common ancestor, its rank and the full lineage for every peptide. Failed files are retried three times and finished files are skipped, so the job can be resubmitted after an interruption.

## text mining (Text_mining_pipeline/)

extract metadata from the papers of the 216 datasets, using DocETL agents& it is used to check whether the taxonomic results agree with what the papers report.

1. **[prepare_data.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Text_mining_pipeline/prepare_data.py)** read the paper text, splits it into sections (abstract, methods, supplementary and others), clean it and write `papers_dataset.json`.

2. **DocETL** one YAML per metadata category. Each one runs a map operation over the papers and extracts a set of fields:
   
   example for running -> docetl run 01_biological_info.yaml
   
   - `01_biological_info.yaml`: organism, strain, age, sex, organism part, specimen, treatment
   - `02_ms_instruments.yaml`: instrument, acquisition method, fragmentation, mass tolerances
   - `03_sample_prep.yaml`: labeling, enzyme, sample preparation
   - `04_separation.yaml`: chromatography and separation
   - `05_data_analysis.yaml`: search engine, database, analysis parameters
   - `06_clinical_experimental.yaml`: clinical and experimental design
   - `07_factor_values.yaml`: experimental factors

4. **[merge_and_generate_ann.py](https://github.com/fatemehmirzade/De_Novo_Peptide_Sequencing/blob/main/Text_mining_pipeline/merge_and_generate_ann.py)** merge the seven outputs per paper and writes the final annotation table.

input text goes in `Papers_text/`, output goes in `Text_mining_results/`.

## data

input spectra and previous annotations: MassIVE `MSV000088598` (doi:10.25345/C52K34, `ftp://massive-ftp.ucsd.edu/v04/MSV000088598/`) and files are `cluster_ident_2.mgf` and `cluster_ident_n.mgf` and the ANN-SoLo results are `cluster_ident_2.mztab` and `cluster_ident_n.mztab`-> together about 168 GB.

Casanovo output for the 39,413,436 unidentified medoid spectra, plus supplementary tables and files: Zenodo, doi:10.5281/zenodo.22305568.

## requirements

- python 3.10 or newer. `pandas`, `matplotlib` and `seaborn` for steps 8 and 11, the other python scripts use the standard library only.
- Casanovo 5.1.2 with the `casanovo_v5_0_0.ckpt` model and a GPU.
- unipept CLI (needs Node.js).
- BLAST+ for the human filtering step.
- docETL and an API key for the text mining pipeline.
