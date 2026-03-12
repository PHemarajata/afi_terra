# AFI Terra Pipeline

Metagenomic Rickettsiales detection workflow.

Terra import options (to avoid missing-import errors):

1) Recommended: Dockstore + GitHub

- This repo includes .dockstore.yml with both workflows:
	- /wdl/AFI_Rickettsiales_Main.wdl
	- /wdl/AFI_Rickettsiales_Batch.wdl
- In Dockstore, link GitHub repository PHemarajata/afi_terra.
- Create or refresh a workflow version from the Git tag/branch.
- In Terra, import from Dockstore (TRS) instead of uploading a single WDL file.
- Terra will pull the full descriptor set and resolve imports automatically.

2) Direct Terra upload (fallback)

- Upload a zip that contains all imported files with relative paths preserved.
- At minimum include:
	- wdl/AFI_Rickettsiales_Batch.wdl
	- wdl/AFI_Rickettsiales_Main.wdl
	- wdl/tasks/*.wdl
	- NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl
- Uploading only AFI_Rickettsiales_Batch.wdl will fail because it imports AFI_Rickettsiales_Main.wdl.

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

Run-by-run operation with NTC-dependent interpretation:

Because ntc_background is run-dependent, execute one sequencing run at a time.

Suggested 2-pass pattern per run:

1) Build run-specific batch JSON (filter by run_id):

`python3 scripts/build_batch_inputs_json.py \
	--mapping-tsv AFI_optimizeProtocol.tsv \
	--run-id 3 \
	--out-json wdl/inputs/run3.batch.json \
	--fastq-uri-prefix gs://YOUR_BUCKET/fastq \
	--default-classifier-mode single \
	--mode-policy auto \
	--per-sample-use-human-scrub true`

2) First Terra pass for that run (include PC8 + NTC + clinical/validation samples) to produce NTC metrics outputs.

3) Build run-specific ntc_background.tsv from NTC metrics files:

`python3 scripts/build_ntc_background_from_metrics.py \
	--ntc-metrics path/to/NTC1.metrics.tsv \
	--ntc-metrics path/to/NTC2.metrics.tsv \
	--out wdl/inputs/run3.ntc_background.tsv`

4) Upload run3.ntc_background.tsv to GCS and re-run the same run3 batch config using that file for AFI_Rickettsiales_Batch.ntc_background.

Notes:

- Keep NTC and PC8 in each run submission so controls are processed with the same run.
- In mode auto, rows with expected_results become validation and rows without expected_results become routine.
- Control type handling:
	- Validation mode supports PC_MIX8, MIXED4, and PC_SINGLE as positive control classes.
	- Routine mode treats only PC_MIX8 as routine positive control class.
- Multi-expected validation controls:
	- For MIXED4/PC_MIX8/PC_SINGLE controls with multiple expected taxa, pass a delimited list in `expected_taxa` (sample-sheet) or `expected_taxon` (JSON/WDL).
	- Supported delimiters: `;` or `,` or `|`.
	- Concordance is `Concordant` only when all expected taxa are detected as Confirmed/Probable.
