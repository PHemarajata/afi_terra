# 16S AFI Unknown 1 Terra Notes

Use `AFI_16S_Batch` in Terra with `Run workflow with inputs defined by file paths`, then upload one of the JSON files in this folder.

Files prepared for this run:

- `16s_afi_unknown_1.batch_samples.tsv`: editable sample sheet with all routine samples, `NTC_S11`, and `PC_S12`
- `16s_afi_unknown_1.batch.first_pass.single.json`: first-pass Terra inputs using the single-classifier path
- `16s_afi_unknown_1.batch.second_pass.single.template.json`: second-pass Terra inputs after building a run-specific NTC background TSV
- `ntc_background.placeholder.tsv`: minimal placeholder TSV to upload for the first pass

Docker override fields are now exposed at the workflow level:

- `AFI_16S_Batch.afi_core_docker`: image used by `fastp`, minimap2, metrics extraction, interpretation, validation, and Kraken tasks
- `AFI_16S_Batch.centrifuger_docker`: image used only by the Centrifuger classification task

Rebuild and push a replacement AFI core image from this repo:

1. `docker login`
2. `bash scripts/build_push_afi_core_image.sh phemarajata614 0.2.1 linux/amd64`
3. The JSON files are already set to `phemarajata614/afi-terra:0.2.1`

Equivalent manual commands:

1. `docker build --platform linux/amd64 -t phemarajata614/afi-terra:0.2.1 .`
2. `docker push phemarajata614/afi-terra:0.2.1`

Sample handling encoded in the JSON files:

- routine unknowns: `sample_type = clinical`, `mode = routine`
- positive control: `PC_S12` with `sample_type = PC_MIX8`, `mode = routine`
- negative control: `NTC_S11` with `sample_type = NTC`, `mode = routine`
- excluded: `Undetermined*`

Before first pass:

1. Upload `ntc_background.placeholder.tsv` to a GCS location.
2. The panel FASTA path is set to `gs://fc-36ebdf16-bf31-4fef-9963-fc780c8f7367/uploads/centrifuger_db/rickettsiales_panel_16S.clean.fa` in both JSON templates.
3. The first-pass placeholder NTC path is set to `gs://fc-80712e02-4823-47c2-bdea-80127f018355/uploads/16S_afi_unknown_1/ntc_background.placeholder.tsv` in the first-pass JSON.
4. The AFI core image override is set to `phemarajata614/afi-terra:0.2.1` in both JSON templates.
5. The Centrifuger image is set to `phemarajata614/centrifuger:1.1.0` in both JSON templates.
6. The bucket currently contains `centrifuger_index.tar.gz`, not an extracted index prefix, so the JSONs now provide both:
	- `AFI_16S_Batch.centrifuger_db = centrifuger_bact_arch_plus_rickettsiales`
	- `AFI_16S_Batch.centrifuger_db_archives = [gs://fc-36ebdf16-bf31-4fef-9963-fc780c8f7367/uploads/centrifuger_db/centrifuger_index.tar.gz]`
7. The Centrifuger task now requests `128G` RAM and `local-disk 500 HDD`, because the archive in GCS is about `67 GiB` compressed and Terra was previously launching the task on a `1 CPU / 2 GB` VM.
8. Upload `16s_afi_unknown_1.batch.first_pass.single.json` into Terra and launch the batch run.

After first pass:

1. Collect the `metrics` outputs for the NTC samples.
2. Build a run-specific background TSV with `scripts/build_ntc_background_from_metrics.py`.
3. Upload that TSV to GCS.
4. Replace `AFI_16S_Batch.ntc_background` in `16s_afi_unknown_1.batch.second_pass.single.template.json` with the real uploaded TSV path.
5. Keep the same working `AFI_16S_Batch.afi_core_docker` override in the second-pass JSON.
6. Re-run the same batch with the second-pass JSON.

Note: `RunCentrifuger` now supports archive-backed Terra runs by extracting `AFI_16S_Batch.centrifuger_db_archives` and resolving the local prefix from `AFI_16S_Batch.centrifuger_db` before invoking `centrifuger -x`.

Reason for the new Centrifuger changes: the attached Terra log showed `centrifuger` segfaulting while pointed at `gs://.../centrifuger_bact_arch_plus_rickettsiales`, but `gsutil ls` against that bucket showed only `centrifuger_index.tar.gz` and `rickettsiales_panel_16S.clean.fa`. The same log also showed Terra running the classifier on `custom-1-2048`, which is far too small for this database.

Reason for the new docker override: the Terra run failed in `FastpClean` with `fastp: command not found`, which indicates the default `phemarajata614/afi-terra:0.1.0` image currently pulled by Terra does not match the Dockerfile in this repository.

If you want to use the double-classifier path instead, change each sample's `classifier_mode` to `double`, then add both `AFI_16S_Batch.kraken_db_16g` and `AFI_16S_Batch.kraken_db_rick` to the JSON.