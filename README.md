# AFI Terra Pipeline

Metagenomic Rickettsiales detection workflow.

Terra import options (to avoid missing-import errors):

1) Recommended: Dockstore + GitHub

- This repo includes .dockstore.yml with both workflows:
	- /wdl/AFI_16S_Main.wdl
	- /wdl/AFI_16S_Batch.wdl
- In Dockstore, link GitHub repository PHemarajata/afi_terra.
- Create or refresh a workflow version from the Git tag/branch.
- In Terra, import from Dockstore (TRS) instead of uploading a single WDL file.
- Terra will pull the full descriptor set and resolve imports automatically.

2) Direct Terra upload (fallback)

- Upload a zip that contains all imported files with relative paths preserved.
- At minimum include:
	- wdl/AFI_16S_Batch.wdl
	- wdl/AFI_16S_Main.wdl
	- wdl/tasks/*.wdl
	- NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl
- Uploading only AFI_16S_Batch.wdl will fail because it imports AFI_16S_Main.wdl.

Docker images used by WDL tasks:

- AFI core image (fastp, minimap2, samtools, kraken2, python, pandas, scripts): `phemarajata614/afi-terra:0.1.0`
- Centrifuger image (single-classifier path): `phemarajata614/centrifuger:1.1`

Build and push AFI core image:

`docker build --platform linux/amd64 -t phemarajata614/afi-terra:0.1.0 .`

`docker push phemarajata614/afi-terra:0.1.0`

Optional latest tag:

`docker tag phemarajata614/afi-terra:0.1.0 phemarajata614/afi-terra:latest`

`docker push phemarajata614/afi-terra:latest`

Quick local tool check:

`docker run --rm phemarajata614/afi-terra:0.1.0 bash -lc "fastp --version && minimap2 --version && samtools --version | head -n 1 && kraken2 --version && python3 -c 'import pandas; print(pandas.__version__)'"`

Terra setup note:

- You can keep WDL default image tags as above, or override task-level docker inputs where exposed (e.g., in classify/validate tasks) during method configuration.

Main steps:

1 Human read removal  
2 fastp read cleaning  
3 Taxonomic classification  
4 minimap2 alignment to Rickettsiales 16S panel  
5 Interpretation using NTC-aware rules  

Single-sample workflow:

`wdl/AFI_16S_Main.wdl`

- `mode`: `validation` or `routine`
- `use_human_scrub`: `true` or `false` (defaults to `true`)
- `centrifuger_resource_profile`: `balanced` by default; use `low_cost`, `balanced`, `high_sensitivity`, or `custom`
- Other performance knobs: `fastp_threads`, `minimap_threads`, and `minimap_sort_memory_per_thread`

Batch workflow (mixed validation + routine in one submission):

`wdl/AFI_16S_Batch.wdl`

- Inputs: parallel arrays from the Terra sample table (`run_ids`, `sample_ids`, FASTQs, `sample_types`, `modes`, and `expected_taxa`)
- Scatters each sample through preprocessing, Centrifuger, 16S alignment, interpretation, and summary tasks
- Supports per-sample `run_id`, `mode`, and `expected_taxa`, plus run-wide scrub and performance settings
- Top-level outputs include `sample_result_manifest.tsv` and `sample_results_bundle.tar.gz`; use these in Terra set tables to find and download sample-specific results without expanding long `Array[File]` cells.

Centrifuger resource profiles:

- `low_cost`: 4 threads, `48G`, `local-disk 200 HDD`
- `balanced` default: 8 threads, `64G`, `local-disk 250 HDD`
- `high_sensitivity`: 16 threads, `128G`, `local-disk 500 HDD`
- `custom`: honors manual `classify_threads`, `centrifuger_memory`, and `centrifuger_disks` inputs

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
	--centrifuger-db centrifuger_bact_arch_plus_rickettsiales \
	--centrifuger-resource-profile balanced`

Direct from AFI mapping table (auto mode from expected results):

`python3 scripts/build_batch_inputs_json.py \
	--mapping-tsv AFI_optimizeProtocol.tsv \
	--out-json wdl/inputs/batch_from_AFI_optimizeProtocol.auto.json \
	--rickettsiales-panel gs://YOUR_BUCKET/ref/rickettsiales_16S_panel.fasta \
	--centrifuger-db centrifuger_bact_arch_plus_rickettsiales \
	--centrifuger-resource-profile balanced \
	--fastq-uri-prefix gs://YOUR_BUCKET/fastq \
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

Run-by-run operation with NTC-dependent interpretation:

Because ntc_background is run-dependent, execute one sequencing run at a time.

Suggested 2-pass pattern per run:

1) Build run-specific batch JSON (filter by run_id):

`python3 scripts/build_batch_inputs_json.py \
	--mapping-tsv AFI_optimizeProtocol.tsv \
	--run-id 3 \
	--out-json wdl/inputs/run3.batch.json \
	--fastq-uri-prefix gs://YOUR_BUCKET/fastq \
	--centrifuger-resource-profile balanced \
	--mode-policy auto \
	--per-sample-use-human-scrub true`

2) First Terra pass for that run (include PC8 + NTC + clinical/validation samples) to produce NTC metrics outputs.

3) Build run-specific ntc_background.tsv from NTC metrics files:

`python3 scripts/build_ntc_background_from_metrics.py \
	--ntc-metrics path/to/NTC1.metrics.tsv \
	--ntc-metrics path/to/NTC2.metrics.tsv \
	--out wdl/inputs/run3.ntc_background.tsv`

4) Upload run3.ntc_background.tsv to GCS and re-run the same run3 batch config using that file for AFI_16S_Batch.ntc_background.

Notes:

- Keep NTC and PC8 in each run submission so controls are processed with the same run.
- In Terra set tables, prefer the `sample_result_manifest` and `sample_results_bundle` outputs for sample-level review and download.
- In mode auto, rows with expected_results become validation and rows without expected_results become routine.
- Control type handling:
	- Validation mode supports PC_MIX8, MIXED4, and PC_SINGLE as positive control classes.
	- Routine mode treats only PC_MIX8 as routine positive control class.
- Multi-expected validation controls:
	- For MIXED4/PC_MIX8/PC_SINGLE controls with multiple expected taxa, pass a delimited list in `expected_taxa` (sample-sheet) or `expected_taxon` (JSON/WDL).
	- Supported delimiters: `;` or `,` or `|`.
	- Concordance is `Concordant` only when all expected taxa are detected as Confirmed/Probable.
