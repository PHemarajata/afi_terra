# AFI Terra Pipeline

Metagenomic Rickettsiales detection workflow.

Main steps:

1 Human read removal  
2 fastp read cleaning  
3 Taxonomic classification  
4 minimap2 alignment to Rickettsiales 16S panel  
5 Interpretation using NTC-aware rules  

Single-sample workflow:

`wdl/AFI_Rickettsiales_Main.wdl`

- `mode`: `validation` or `routine`
- `classifier_mode`: `single` (Centrifuger) or `double` (Kraken2 16G + Kraken2 Rick)
- `use_human_scrub`: `true` or `false` (defaults to `true`)

Batch workflow (mixed validation + routine in one submission):

`wdl/AFI_Rickettsiales_Batch.wdl`

- Input: `Array[SampleSpec] samples`
- Scatters each sample through `AFI_Rickettsiales_Main`
- Supports per-sample `mode`, `classifier_mode`, `expected_taxon`, and optional scrub override

Example inputs:

- `wdl/inputs/validation_single.example.json`
- `wdl/inputs/validation_double.example.json`
- `wdl/inputs/routine_single.example.json`
- `wdl/inputs/routine_double.example.json`
- `wdl/inputs/batch_mixed.example.json`

TSV to batch JSON helper:

- Template sample sheet: `wdl/inputs/batch_samples.template.tsv`
- Builder script: `scripts/build_batch_inputs_json.py`
- Example command:

`python3 scripts/build_batch_inputs_json.py \
	--sample-sheet wdl/inputs/batch_samples.template.tsv \
	--out-json wdl/inputs/batch_generated.example.json \
	--rickettsiales-panel gs://YOUR_BUCKET/ref/rickettsiales_16S_panel.fasta \
	--ntc-background gs://YOUR_BUCKET/ref/ntc_background.tsv \
	--kraken-db-16g gs://YOUR_BUCKET/db/kraken2_16g \
	--kraken-db-rick gs://YOUR_BUCKET/db/kraken2_rickettsiales \
	--centrifuger-db gs://YOUR_BUCKET/db/centrifuger_bact_arch_plus_rickettsiales`

Direct from AFI mapping table (auto mode from expected results):

`python3 scripts/build_batch_inputs_json.py \
	--mapping-tsv AFI_optimizeProtocol.tsv \
	--out-json wdl/inputs/batch_from_AFI_optimizeProtocol.auto.json \
	--rickettsiales-panel gs://YOUR_BUCKET/ref/rickettsiales_16S_panel.fasta \
	--ntc-background gs://YOUR_BUCKET/ref/ntc_background.tsv \
	--kraken-db-16g gs://YOUR_BUCKET/db/kraken2_16g \
	--kraken-db-rick gs://YOUR_BUCKET/db/kraken2_rickettsiales \
	--centrifuger-db gs://YOUR_BUCKET/db/centrifuger_bact_arch_plus_rickettsiales \
	--fastq-uri-prefix gs://YOUR_BUCKET/fastq \
	--default-classifier-mode double \
	--mode-policy auto \
	--per-sample-use-human-scrub true`

Compare single vs double outputs after Terra runs:

- Script: `scripts/compare_single_double_outputs.py`
- Accepts either specific summary TSV files or parent directories (recursive search).

Example (compare both validation and routine outputs):

`python3 scripts/compare_single_double_outputs.py \
	--single-validation /path/to/single_run_outputs \
	--double-validation /path/to/double_run_outputs \
	--single-routine /path/to/single_run_outputs \
	--double-routine /path/to/double_run_outputs \
	--out-prefix comparison/single_vs_double`

Outputs:

- `comparison/single_vs_double.summary.tsv`
- `comparison/single_vs_double.validation_differences.tsv`
- `comparison/single_vs_double.routine_differences.tsv`
