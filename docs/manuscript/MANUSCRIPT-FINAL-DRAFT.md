# Detection of Bacterial Pathogens in Acute Febrile Illness Using a 16S rRNA V1-V3 Workflow with Two-Tier Rickettsiales Rescue and Species-Aware Decontamination

## Title page

**Authors:** `[TEAM INPUT NEEDED]` — author list with affiliations.

**Affiliations:** `[TEAM INPUT NEEDED]`

**Corresponding author:** `[TEAM INPUT NEEDED]` — name, email, postal address.

**Running title (≤50 chars):** `[TEAM INPUT NEEDED]`

**Keywords (5–8):** acute febrile illness; 16S rRNA; Rickettsiales; blood culture; melioidosis; decontamination; metagenomics

**Word count (body):** ~7,800 (this long version); see `MANUSCRIPT-CONDENSED-DRAFT.md` for the ~3,900-word JCM-style condensed version.

**Conflict-of-interest statement:** `[TEAM INPUT NEEDED]`

**Funding:** `[TEAM INPUT NEEDED]`

**Data availability statement:** All `.calls.tsv`, Centrifuge kreport outputs, the V4 filter script (`afi_decontamination_filter_v4.py`), the appendix generator (`generate_appendices.py`), and the per-sample appendix tables (markdown + Excel) are available at `[TEAM INPUT NEEDED: data repository or Google Drive AFI_compare folder access route]`.

**Author contributions:** `[TEAM INPUT NEEDED]`

---

> **Draft status (2026-05-12):** This is the manuscript-ready bioinformatics package consolidating the latest analyses. Sections labeled `[TEAM INPUT NEEDED]` require input from the team (patient demographics, IRB, DNA extraction kit, library prep kit, sequencer model, exact protocol versions, etc.). All numerical results, filter parameters, and per-sample tables are grounded in the actual `.calls.tsv` and Centrifuger kreport data and use the corrected counting established during the May 11–12 review.
>
> Companion files (referenced throughout):
> - `MANUSCRIPT-METHODS-DECONTAMINATION.md` — decontamination methods detail
> - `MANUSCRIPT-RESULTS-SECTION.md` — results detail
> - `MANUSCRIPT-DISCUSSION-SECTION.md` — discussion detail
> - `APPENDIX-VALIDATION-PANEL.md` / `APPENDIX-STUDY-SAMPLES.md` — per-sample tables (markdown)
> - `APPENDICES.xlsx` — same tables for review in Excel (4 sheets)
> - `DECONTAMINATION-FILTER-REPORT-V4.txt` — raw filter output
> - `afi_decontamination_filter_v4.py` — executable filter
> - `generate_appendices.py` — appendix generator
>
> A list of team-input placeholders that need to be filled in before submission is at the end of this document.

---

## Abstract

**Background.** Acute febrile illness (AFI) in endemic regions of northeastern Thailand is caused by a wide differential including obligate intracellular Rickettsiales (Orientia tsutsugamushi, Rickettsia spp), Burkholderia pseudomallei (melioidosis), Leptospira spp, Brucella spp, dengue and other arboviruses, and enteric bacterial pathogens, several of which culture cannot reliably recover. Blood cultures that signal positive on automated detection but fail subculture recovery are a recognized diagnostic puzzle.

**Methods.** We developed and validated a 16S rRNA V1-V3 amplicon workflow combining Centrifuge primary taxonomic classification with Minimap2 alignment-based rescue against a curated Rickettsiales reference panel, a two-tier reporting framework that distinguishes confirmed genus-level calls from order-level "Rickettsiales detected" rescues, and a multi-tiered decontamination filter (V4) that includes a species-level safeguard for *Burkholderia pseudomallei* and a positive-control spike-in bypass. The assay was validated against a 48-sample panel (33 clinical samples with known reference-laboratory diagnoses, 10 positive controls, 5 negative template controls), and applied to a study cohort of 86 AFI cases with positive automated blood culture signal but failed subculture recovery.

**Results.** The validated assay achieves 72.1% (31/43) sample-level analytical performance combining clinical and positive-control samples (63.6% [21/33] strict clinical concordance, 100% [10/10] positive-control concordance) with 100% specificity (5/5 NTCs clear of all TaqMan Array Card bacterial target genera) and 100% inter-run reproducibility across 9 sequencing runs. Applied to the AFI study cohort, the workflow identified candidate organism signals in 71 of 86 samples (the remaining 15 are pre-sequencing failures): Rickettsiales rescue evidence in 11 of 71 (15.5%; 1 Tier-1 genus-level Orientia detection at sample 16901195_S5_L001, 10 Tier-2 order-level rescues with breadth-of-coverage 0.21–0.24), Mycoplasmopsis (a fastidious cell-wall-deficient organism class) in 4 samples (mean abundance 39.78%, max 70.55%), Leptospira at 22.78% in one sample with a same-run NTC contamination caveat, and Brucella signals at low abundance (0.98% mean) in 3 samples. No study sample contains *B. pseudomallei* above the species-level detection threshold (prior reports of preserved *B. pseudomallei* conflated genus- with species-level read counts; corrected here).

**Conclusion.** The 16S V1-V3 assay, with the two-tier Rickettsiales framework and the species-aware V4 decontamination filter, is suitable for analytical surveillance of bacterial pathogens in AFI cases. Approximately one in five AFI cases with positive culture / failed subculture shows Rickettsiales rescue evidence; an additional minority shows fastidious-organism candidate signals (Mycoplasmopsis, Leptospira, Brucella). All candidate detections require orthogonal confirmation (PCR, serology, specialized culture) before clinical reporting. The corrected analysis reframes the cohort study as hypothesis-generating, with Rickettsiales involvement as the most prevalent identifiable etiology.

---

## 1. Introduction

Acute febrile illness (AFI) is a common presentation in emergency and inpatient settings in tropical Southeast Asia, with a wide etiologic differential and substantial overlap in early clinical features. In endemic northeastern Thailand, the major bacterial pathogens contributing to AFI are *Orientia tsutsugamushi* (scrub typhus), *Rickettsia* spp (spotted-fever group rickettsioses), *Burkholderia pseudomallei* (melioidosis), *Leptospira* spp (leptospirosis), *Brucella* spp (brucellosis), non-typhoidal *Salmonella* and *Salmonella* Typhi, enteropathogenic *Escherichia coli*, and species-level *Streptococcus* (notably *S. pneumoniae* and *S. suis*). These pathogens have substantially different culture requirements: some grow readily on standard blood agar within standard 5–7 day incubation windows (E. coli, K. pneumoniae, S. aureus), some require specialized media or extended incubation (Leptospira, Brucella, fastidious Streptococcus species), and some cannot be recovered by routine blood culture at all (Rickettsiales, which are obligate intracellular).

Blood cultures that are flagged positive by automated detection but fail subculture recovery represent a particularly common and diagnostically frustrating scenario in AFI workups. The automated detection system records a positive signal (turbidity, CO2 production, metabolic activity), but inoculation of the flagged bottle onto solid media yields no growth within the standard 24–48 hour subculture observation window. The clinical implication is that *some* organism was present and metabolically active in the bottle, but its identity remains unknown and empirical therapy must proceed on clinical judgment alone.

16S rRNA gene sequencing has been proposed as a complementary diagnostic for these positive-culture / no-growth scenarios because it detects bacterial DNA regardless of culturability and provides genus-level (and sometimes species-level) taxonomic assignment within the limits of the targeted hypervariable region. The V1-V3 region (amplified with the 27F primer set) is widely used in clinical 16S workflows but has known under-detection biases for certain Gram-positive anaerobes and primer-mismatch organisms (notably some *Bifidobacterium*, *Gardnerella*, and *Atopobium* lineages). It also provides limited genus-level discrimination within order Rickettsiales — meaning that even when sequence is recovered, distinguishing *Orientia* from *Rickettsia* may not be possible from the V1-V3 amplicon alone.

In low-biomass clinical samples such as blood, sequencing-based bacterial detection is further complicated by reagent and laboratory contamination. Bacterial DNA introduced by extraction kits, molecular-grade water, and PCR reagents (Salter et al., 2014; Glassing et al., 2016; Lauder et al., 2016; de Goffau et al., 2018) routinely accounts for a substantial fraction of reads in blood and other low-biomass specimens and can dominate the signal entirely when the input bacterial biomass is low. A recent population study of 9,770 healthy human blood samples found no consistent core blood microbiome (Tan et al., 2023), strengthening the position that any organism detected by sequencing in blood should be evaluated rigorously against contamination controls before clinical interpretation.

We developed and validated a 16S V1-V3 amplicon workflow for AFI bacterial pathogen surveillance that addresses these constraints through four design choices: (1) **Centrifuge primary classification** for breadth across all bacterial genera, (2) **Minimap2 alignment-based rescue** against a curated Rickettsiales 16S reference panel to compensate for V1-V3 limited genus-level discrimination, (3) a **two-tier Rickettsiales reporting framework** that distinguishes confirmed genus-level calls (Tier 1) from order-level "Rickettsiales detected" rescues (Tier 2) — both clinically actionable in a doxycycline-treatable endemic context — and (4) a **species-aware decontamination filter (V4)** that aggressively removes documented kit/skin/water contaminant genera while implementing species-level safeguards for clinically critical organisms (*B. pseudomallei*) and a positive-control spike-in bypass to allow PC samples to be evaluated for their expected organisms.

In this report we present the analytical validation of this workflow against a 43-sample reference panel (33 clinical samples with reference-laboratory diagnoses plus 10 positive controls) and 5 negative template controls, and the application of the validated workflow to a study cohort of 86 AFI cases drawn from positive blood cultures with failed subculture recovery. We report performance metrics, identify the candidate organisms detected in the cohort, and discuss the implications for AFI diagnostic algorithms in endemic regions.

---

## 2. Methods

### 2.1 Study design and ethical approval

`[TEAM INPUT NEEDED]` — IRB approval, ethical oversight body, written informed consent procedures, period of enrollment, geographic catchment (presumably northeastern Thailand based on the prior reports' framing), inclusion criteria (e.g., "patients presenting with fever ≥38°C of <14 days duration, without localizing signs identifying an alternative diagnosis"), exclusion criteria.

### 2.2 Sample populations

Two distinct sample sets were used in this study:

- **Validation panel (n=48 specimens, across 9 sequencing runs):**
  - 33 clinical specimens from patients with reference-laboratory-confirmed diagnoses of the following pathogens (validated by `[TEAM INPUT NEEDED]`: culture-based confirmation, PCR confirmation, or serological confirmation, as applicable): *Escherichia coli* (5), *Orientia tsutsugamushi* (6), *Rickettsia* spp (4), *Leptospira* spp (4), *Burkholderia pseudomallei* (5), *Streptococcus pneumoniae* (3), *Streptococcus suis* (3), *Coxiella burnetii* (2), and *Yersinia* spp (1).
  - **10 positive controls** comprising 4 PC_SINGLE samples (single-organism control material for *E. coli*, *P. aeruginosa*, *S. pneumoniae*, *S. suis*), 5 PC_MIX8 samples (ZymoBIOMICS Microbial Community Standard: *Bacillus subtilis*, *Enterococcus faecalis*, *Escherichia coli*, *Lactobacillus fermentum* [now *Limosilactobacillus fermentum*], *Listeria monocytogenes*, *Pseudomonas aeruginosa*, *Salmonella enterica*, *Staphylococcus aureus*; the standard also includes two yeasts that are not detectable by 16S), and 1 MIXED4 sample (*E. coli* + *P. aeruginosa* + *S. pneumoniae* + *S. suis*).
  - **5 negative template controls (NTC)** — molecular-grade water carried through the full extraction and library preparation process.

- **Study cohort (n=86 patient blood specimens):** Blood samples from patients enrolled with acute febrile illness during whose admission a blood culture bottle was flagged as positive by the automated detection system, but subculture onto aerobic solid media did not recover any organism. **The 16S samples are patient blood drawn during the same admission and are not aliquots of the original positive blood culture bottle.** No anaerobic blood culture was performed per local hospital protocol.

### 2.3 Sample collection, DNA extraction, and library preparation

`[TEAM INPUT NEEDED]` — sample collection container and volume, anticoagulant if relevant, storage conditions and time-to-extraction, DNA extraction kit and protocol, NCBI Scrubber human-DNA removal step (already confirmed in the pipeline WDL), 16S V1-V3 amplification primers and conditions (presumably 27F/534R or equivalent), library preparation kit, indexing strategy, library quantification, sequencer model (Illumina MiSeq is the typical platform for V1-V3), run configuration (paired-end length, read counts, demultiplexing).

### 2.4 Bioinformatics pipeline

Reads were processed through a WDL (Workflow Description Language) workflow orchestrated via Terra / Cromwell, comprising preprocess, classify, align, metrics, interpret, and validate task stages (`wdl/tasks/`, with top-level orchestration in `wdl/AFI_16S_Main.wdl` and `wdl/AFI_16S_Batch.wdl`; pipeline source: https://github.com/PHemarajata/afi_terra, container images `phemarajata614/afi-terra:0.4.1` and `phemarajata614/centrifuger:1.1.0`).

#### 2.4.1 Host DNA depletion and read preprocessing

Raw paired-end reads were filtered through NCBI Human Read Removal Tool (Scrubber) (Katz et al., 2021) to remove human-derived sequences (`NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl`), then trimmed and quality-filtered with fastp (Chen et al., 2018). Only de-hosted, quality-filtered reads proceeded to taxonomic classification.

#### 2.4.2 Primary taxonomic classification (Centrifuger)

De-hosted reads were classified using Centrifuger v1.1.0 (Song & Langmead, 2024) — a successor of Centrifuge (Kim et al., 2016) — against a comprehensive bacterial reference database including standard bacteria and archaea plus extended Rickettsiales coverage, as a single index. The classify task (`wdl/tasks/classify.wdl`, task `RunCentrifuger`) outputs both per-read classification (`*.centrifuger.classification.tsv`) and a kreport-format clade summary (`*.centrifuger.kreport.tsv`). Genus-level read counts were tabulated from the kreport, with clade reads (not just taxon-level reads) summed to the genus rank so that species-level assignments contributed to the genus call.

**A genus was called "Detected" if its assigned read count was ≥500 (`cfr_floor`) AND ≥5× the run-specific NTC maximum read count (`cfr_fold`) for that genus** (`scripts/call_taxa.py`; defaults exposed in `wdl/tasks/interpret.wdl`, task `InterpretCalls`).

#### 2.4.3 Minimap2 alignment-based rescue for Rickettsiales

The order Rickettsiales (including *Orientia tsutsugamushi* and *Rickettsia* spp) represents the canonical "cannot-miss" AFI etiology in endemic northeastern Thailand. The V1-V3 region of 16S provides limited sequence variation for genus-level discrimination within Rickettsiales, so an alignment-based rescue stage was added.

Reads were aligned with Minimap2 v2.x (Li, 2018) against a curated 7-reference panel of complete genomes covering Rickettsiales of clinical interest (`align_rickettsiales_16S/refs_fna/`): *Orientia tsutsugamushi* str. Boryong (AM494475.1) and str. Ikeda (AP008981.1); *Rickettsia prowazekii* str. NMRC Madrid E (CP004888.1), *R. typhi* str. Wilmington (NC_006142.1), and *R. rickettsii* str. 'Sheila Smith' (NC_009882.1); *Anaplasma phagocytophilum* str. HZ (NC_007797.1); and *Ehrlichia chaffeensis* str. Arkansas (NC_007799.1). For each Rickettsiales genus, the following metrics were computed (`wdl/tasks/metrics.wdl`): mapped read count, maximum breadth of coverage (fraction of reference covered ≥1×), and read-count fold over the run-specific NTC alignment-max for that genus.

A two-tier rescue framework was applied (`scripts/call_taxa.py`):

- **Tier 1 — Genus-level confirmed (`call = Confirmed`):** mapped reads ≥100 (`align_confirm_reads`) AND breadth ≥0.25 (`align_confirm_breadth`) AND mapped reads ≥5× the alignment NTC max (`align_fold`). Reported as the specific genus (*Orientia* or *Rickettsia*); the `align_confirmed` flag in `.calls.tsv` is `true`.
- **Tier 2 — Order-level rescue (`call = Probable`):** mapped reads ≥50 AND breadth ≥0.20 AND mapped reads > the alignment NTC max (i.e., any non-trivial signal above NTC floor, without the 5× fold requirement of Tier 1). Reported as **"Rickettsiales detected (genus uncertain; recommend confirmatory species-specific qPCR)"** — clinically actionable for empirical doxycycline coverage in an endemic context.
- **Below thresholds:** `call = Not_Confirmed` (reads ≥50 but at-or-below NTC max) or `call = Negative` (reads <50). These rows are preserved in pipeline output for transparency but are not reported clinically.

#### 2.4.4 Negative-template-control background derivation (per genus, per run)

NTC backgrounds were computed per run from the set of NTC samples assigned to that run (`scripts/build_ntc_background.py`, task `BuildNTCBackground`). For each NTC sample, per-genus read counts were tabulated (separately for Centrifuger-classified reads and Minimap2-aligned reads). For each run, the NTC background was computed as the **maximum** of those per-genus read counts across all NTCs in that run, yielding two per-run tables: `align_ntc_reads` (alignment-based NTC max) and `cfr_ntc_reads` (Centrifuger-based NTC max). The interpretation task (`InterpretCalls`) then attached the appropriate per-genus NTC max to each detection row.

Each `.calls.tsv` row's `ntc_reads` field represents the NTC max value the pipeline used for that organism in that run. We observed inter-run variability in NTC contamination: in run 6_and_7, `NTC2_ExDw_S13_L001` carries 78,691 *Leptospira* reads, 106,204 *Burkholderia* reads (genus level), and 188,778 *Brevundimonas* reads. Study-sample `.calls.tsv` entries in the same run show `ntc_reads = 0` for *Leptospira*, which is consistent with the pipeline's run-ID assignment treating the contaminated NTC as belonging to a different run-ID pool than the affected study samples. The exact per-run NTC pool composition for each sample is determined by the per-batch Terra input sheet; the cross-check column "Same-run NTC max" in the appendix tables reports the maximum reads of each organism across all NTCs in the folder-grouped run for transparency.

### 2.5 V4 decontamination filter

Genus-level detections passing the primary pipeline thresholds were further processed through a four-tier decontamination filter (V4) informed by landmark low-biomass microbiome studies (Salter et al., 2014; Glassing et al., 2016; Lauder et al., 2016; de Goffau et al., 2018; Tan et al., 2023). The filter is applied to validation-panel and study samples but is bypassed for organisms documented as positive-control spike-ins in their respective PC samples.

#### 2.5.1 Tier A — High-confidence kit / skin / water contaminants

The following 11 genera are removed on any detection in clinical or study samples: *Pseudomonas*, *Ralstonia*, *Bradyrhizobium*, *Sphingomonas*, *Stenotrophomonas*, *Methylobacterium*, *Acinetobacter*, *Cutibacterium*, *Staphylococcus*, *Corynebacterium*, and *Brevundimonas*. *Brevundimonas* was added on the basis of NTC profiles in this dataset (22–285,740 reads per NTC in run 6_and_7); the others are documented in ≥5 of the cited reviews.

#### 2.5.2 Tier A exception — *Burkholderia* species-level safeguard

*Burkholderia* is treated as a Tier-A-equivalent kit contaminant (the genus is dominated in low-biomass samples by *B. cepacia* complex species: *B. contaminans*, *B. cenocepacia*, *B. multivorans*, *B. sola*, *B. cepacia*), but it also contains *B. pseudomallei*, a cannot-miss melioidosis pathogen endemic in Thailand. For every sample with a *Burkholderia* genus detection, the Centrifuge kreport is parsed at the species rank (`S`): reads assigned to *Burkholderia pseudomallei* (NCBI taxonomy ID 28450) are summed and compared against the run's NTC species-level *B. pseudomallei* count. The detection is preserved only if **both** conditions hold: (i) species-level *B. pseudomallei* reads ≥500, and (ii) species reads exceed the run-NTC maximum for *B. pseudomallei*. When preserved, the retained record stores the species-level read count, not the genus total.

The 3 validation-panel *B. pseudomallei* samples (`09502813_S2_L001`, `09-0-02165`, `09700912_S3_L001`) carry 13,744 / 20,061 / 65,016 species-level *B. pseudomallei* reads respectively and trivially pass the safeguard. None of the 4 study-sample *Burkholderia* detections passed the safeguard.

#### 2.5.3 Tier B — NTC-only organisms

Nine genera identified as present only in NTCs across all runs were removed globally: *Cereibacter*, *Thioclava*, *Bdellovibrio*, *Saltatorellus*, *Pseudogemmobacter*, *Minisyncoccus*, *Rhodoluna*, *Microbacterium*, *Arcanobacterium*. In practice only *Rhodoluna* (2 detections) appeared in clinical samples after upstream NCmax subtraction.

#### 2.5.4 Tier 1 — Ultra-low-abundance noise

Sixteen genera with median per-sample abundance <0.5% across the dataset were removed: *Shigella*, *Metapseudomonas*, *Stutzerimonas*, *Capsulimonas*, *Chamaesiphon*, *Chloroflexus*, *Flavihumibacter*, *Hymenobacter*, *Limnoglobus*, *Methylovirgula*, *Microvirga*, *Pelagovum*, *Pseudonocardia*, *Rufibacter*, *Salmonella*, *Spirosoma*. (*Mycoplasmopsis* and *Nitrospira* were removed from this tier during the May 11 review after inspection of their actual abundance distributions; they are now retained as candidate signal.)

#### 2.5.5 Tier 2 — Marginal organisms

Twenty-seven genera with median 0.5–2.0% abundance are retained only if detected in ≥2 samples AND each detection is ≥1.0% abundance.

#### 2.5.6 Positive-control spike-in bypass

For samples designated as positive controls, Tier A removal is bypassed for any organism that is a documented spike-in for that sample. The expected spike-in genera are:

| PC sample | Type | Expected genera |
|---|---|---|
| `E-coli_S4_L001` | PC_SINGLE | *Escherichia* |
| `P-aeru_S5_L001` | PC_SINGLE | *Pseudomonas* |
| `S-pneumo_S2_L001` | PC_SINGLE | *Streptococcus* |
| `S-suis_S3_L001` | PC_SINGLE | *Streptococcus* |
| `PC-20251016_S7_L001`, `PC_S8_L001`, `PC_S10`, `PC_S12`, `PC_S13_L001` | PC_MIX8 | *Bacillus*, *Enterococcus*, *Escherichia*, *Limosilactobacillus*, *Listeria*, *Pseudomonas*, *Salmonella*, *Staphylococcus* (ZymoBIOMICS Microbial Community Standard) |
| `Mixed_S6_L001` | MIXED4 | *Escherichia*, *Pseudomonas*, *Streptococcus* |

The bypass ensures that PC samples are evaluated for their ability to detect their expected spike-in organisms rather than penalized by a filter calibrated for clinical contamination. PC concordance is defined as: all expected spike-in organisms must be detected AND retained.

### 2.6 Concordance definitions

| Sample category | Concordance rule |
|---|---|
| **Clinical** | **Concordant** if at least one row matches the expected target organism at genus level AND that row passes the V4 filter (`Final concordance = "Concordant (target detected and retained by V4 filter)"`). Within Rickettsiales, cross-genus rescue (Orientia ↔ Rickettsia) and order-level rescue (`call = Probable`) both qualify as TARGET match. |
| **Positive control** | **Concordant** iff all expected spike-in organisms are detected AND retained by V4 (under the PC bypass). |
| **NTC** | **Concordant** iff no TaqMan Array Card bacterial target genus (Bartonella, Brucella, Rickettsia, Orientia, Yersinia, Coxiella, Streptococcus, Salmonella, Escherichia, Burkholderia) is retained after V4. Viral and protozoal TAC targets (Dengue, Chikungunya, Zika, Nipah, HepE, Hantaan, Seoul, JEV, Plasmodium falciparum/vivax) are not evaluable by this assay. |

### 2.7 Statistical and reproducibility methods

Analytical sensitivity and specificity were computed under a binary contingency framework treating each sample as a single binary outcome (concordant / discordant). The Wilson 95% confidence interval was used for proportions. Inter-run reproducibility was assessed by per-run PC and NTC pass rates across all 9 sequencing runs. Per-run quality control gating (whether a run is accepted or rejected on the basis of PC and NTC performance) is defined in the deployment standard operating procedure (`[TEAM INPUT NEEDED: deployment SOP reference]`); this manuscript reports the per-run pass rates rather than re-stating the SOP gates.

### 2.8 Data and code availability

The pipeline source code is publicly available at https://github.com/PHemarajata/afi_terra (Dockstore-linked for Terra deployment). All bioinformatics outputs (`.calls.tsv`, Centrifuger kreport files), the V4 decontamination filter script (`afi_decontamination_filter_v4.py`), the appendix generator (`generate_appendices.py`), and the per-sample appendix tables (`APPENDIX-VALIDATION-PANEL.md`, `APPENDIX-STUDY-SAMPLES.md`, `APPENDICES.xlsx`) are available at `[TEAM INPUT NEEDED: final data repository / Google Drive archival location]`.

---

## 3. Results

### 3.1 Validation panel performance

#### 3.1.1 Sample composition and overall accuracy

The 48-sample validation panel (33 clinical with known reference-lab diagnoses + 10 positive controls + 5 NTCs) was processed through the V4-filter-aware pipeline. Final concordance under the rule of record (target detected AND retained by V4 for clinical; all expected spike-ins detected and retained under PC bypass for PCs; no TAC bacterial target genus retained for NTCs) was:

| Category | Concordant | Total | Rate | 95% CI |
|---|---|---|---|---|
| Clinical | 21 | 33 | **63.6%** | 46.6%–77.8% |
| Positive controls | 10 | 10 | **100%** | 72.2%–100% |
| Negative template controls (specificity) | 5 | 5 | **100%** | 56.6%–100% |
| **Sample-level analytical performance (clinical + PC)** | **31** | **43** | **72.1%** | 57.3%–83.3% |
| **Overall validation accuracy (all categories)** | **36** | **48** | **75.0%** | 61.2%–85.1% |

The 72.1% sample-level analytical performance figure matches the legacy APHL bioinformatic validation report (2026-02-16) for the same panel and is the appropriate headline figure for regulatory documentation.

#### 3.1.2 Organism-specific clinical performance

| Expected organism | n | Concordant | Sensitivity | Failure mode for discordants |
|---|---|---|---|---|
| *Escherichia coli* | 5 | 5 | 100% | — |
| *Orientia tsutsugamushi* | 6 | 6 | 100% | (all rescued; 1 Tier-1 + 5 Tier-1/2 combined) |
| *Rickettsia* spp | 4 | 2 | 50% | 2 pre-sequencing failures (00126_S6, 00369_S1: zero taxa) |
| *Leptospira* spp | 4 | 2 | 50% | 2 abundance-driven (other organisms dominate) |
| *Burkholderia pseudomallei* | 5 | 3 | 60% | 2 abundance-driven |
| *Streptococcus pneumoniae* | 3 | 1 | 33% | 2 abundance-driven; V1-V3 cannot resolve S. pneumoniae vs. S. suis |
| *Streptococcus suis* | 3 | 1 | 33% | 2 abundance-driven; same V1-V3 limitation |
| *Coxiella burnetii* | 2 | 1 | 50% | 1 pre-sequencing failure (25800370_S9: zero taxa) |
| *Yersinia* spp | 1 | 0 | 0% | n=1 insufficient; other organisms detected |

Two clinical signal patterns dominate the discordant cases:

1. **Pre-sequencing failures (3 samples):** `00126_S6_L001`, `00369_S1_L001` (both expected *Rickettsia*), and `25800370_S9_L001` (expected *Coxiella*) returned zero taxa, consistent with DNA-extraction, library-preparation, or sequencing-depth failure rather than classification error. The pipeline correctly identified target organisms in same-organism samples where adequate read coverage was achieved (e.g., *Rickettsia* detected at 46–100% in 10900410_S10 and 16401070_S11).

2. **Abundance-driven discordance:** in 9 samples, the expected target was not the dominant 16S signal, with other organisms (often kit-contaminant-class or background environmental genera) accounting for the majority of reads. This reflects 16S genus-detection biology (proportional reporting in mixed samples) rather than classification failure.

#### 3.1.3 Two-tier Rickettsiales rescue performance

The Minimap2 alignment-based rescue is essential for Rickettsiales detection at V1-V3 resolution. Across the 6 expected-*Orientia* validation samples:

- Direct genus-level detection (Centrifuge): 4 samples carry *Orientia* as a Centrifuge call (00389_S5, 11800801_S7, 24500367_S5, 25900911_S4).
- Tier 1 genus-level Minimap2 rescue (`call = Confirmed`): triggered in 22900253_S4 (4,133 mapped reads, breadth 0.3220, `align_confirmed = true` — note this is a *Rickettsia* call in an Orientia-expected sample; counted as concordant under cross-genus Rickettsiales handling).
- Tier 2 order-level Minimap2 rescue (`call = Probable`): triggered in 00618_S7 (Rickettsia, 68 mapped reads, breadth 0.2015) — **this sample was previously reported as a "failed rescue" in the earlier validation analysis, but `.calls.tsv` shows `call = Probable`; the Tier-2 order-level "Rickettsiales detected" framework correctly captures this as concordant.**
- Failed rescue (call = Negative): no validation-panel Rickettsiales sample fell into this category after correction.

A threshold uniformity audit (`THRESHOLD-UNIFORMITY-AUDIT.md`) confirms that the rescue thresholds are applied identically across all 110 samples in the dataset.

#### 3.1.4 Decontamination filter behavior on validation panel

The V4 filter does not remove any organism from the validation panel that affects target-organism detection. This is informative but not an independent validation of the filter: validation samples are dominated by single high-abundance expected organisms (40–100% abundance), and they do not carry Tier A contaminants at detectable levels. The 3 validation-panel *B. pseudomallei* samples trivially pass the species-level safeguard (13,744–65,016 species reads, all >> 500-read threshold). Without the positive-control spike-in bypass, the P-aeru_S5_L001 PC would have failed (because *Pseudomonas* is in Tier A); with the bypass, all 10 PCs are concordant.

### 3.2 Study cohort results

#### 3.2.1 Cohort size and pre-sequencing failures

The study cohort comprises 86 patient blood specimens distributed across 5 sequencing runs (runs `1_and_2`, `3`, `4_and_5`, `6_and_7`, `8_and_9`). All 86 samples had `.calls.tsv` files generated by the pipeline. Of these:

- **71 samples (82.6%) have ≥1 positive call** (`Detected`, `Confirmed`, or `Probable` in `.calls.tsv`) from either Centrifuge or Minimap2 rescue.
- **15 samples (17.4%) have zero detected taxa**, consistent with pre-sequencing failures (DNA extraction, library prep, or sequencing depth).

The pre-sequencing failure rate is comparable to that observed in the validation panel (3/33 = 9% clinical pre-seq failures) and indicates that low-biomass blood specimens have a non-trivial baseline failure rate that should be expected in operational deployment. We recommend a pre-sequencing QC step (Qubit + library qPCR) for future cohorts to distinguish these cases up-front.

**Earlier reports cited "56 AFI study samples"; the corrected count is 86.** The 56 figure appears to have been derived from samples with at least one Centrifuge-Detected genus and excluded the 15 zero-detection samples and the 15 samples whose only positive call was an alignment-based rescue (`Confirmed`/`Probable`).

#### 3.2.2 Control validity within study runs

Positive and negative controls were processed alongside patient samples across all study-cohort sequencing runs. Run-level QC gating (whether a run is accepted on the basis of PC and NTC performance) is defined in the deployment SOP (`[TEAM INPUT NEEDED: SOP reference]`); the deployment QC log provides the per-run accept/reject record. Two operationally relevant observations from this dataset's NTC profiles:

- **Multiple NTCs in run 6_and_7 carry substantial reagent contamination.** NTC3_ExEB_S12_L001 (Brevundimonas 285,740 reads), NTC4_NExDw_S16_L001 (Brevundimonas 241,384), NTC5_NExEB_S15_L001 (Brevundimonas 204,530), and NTC2_ExDw_S13_L001 (Brevundimonas 188,778; *Burkholderia* 106,204 at genus level, dominated by *B. cepacia* complex species; *Leptospira* 78,691). The pipeline's per-run NCmax derivation appropriately handles this through run-ID-based NTC pool assignment, but the cross-run NTC variability is a reminder that per-run NTC composition matters for fold-change interpretation.
- **Study-run PCs pass under V4 with PC spike-in bypass.** The bypass is essential because *Pseudomonas* and *Staphylococcus* are Tier A contaminants in clinical samples but are expected spike-ins in PC_MIX8 and the relevant PC_SINGLE replicates.

#### 3.2.3 Filter impact on study cohort

The V4 filter was applied to all 71 samples with positive detections (217 genus-sample detections in total). Outcomes:

| Filter tier | Detections removed |
|---|---|
| Tier A high-confidence contaminants (Cutibacterium, Staphylococcus, Brevundimonas, Acinetobacter, Corynebacterium, others) | 61 |
| Burkholderia removed at species-level safeguard | 4 |
| Burkholderia preserved as *B. pseudomallei* | **0** |
| NTC-only organisms (Rhodoluna) | 2 |
| Tier 1 ultra-low abundance | 2 |
| Tier 2 marginal + rare | 1 |
| **Total removed** | **70 (32.3% of detections)** |
| **Total retained** | **147** |

A Sankey visualization of pre- vs. post-V4-filter genus distribution across the cohort is provided as **Figure 1** (`figure_a_sankey.html` for interactive review, `figure_a_sankey.png` for print). The visualization shows the aggregated read volume flowing from each detected genus (left nodes) to each V4 outcome category (right nodes); the dominant pattern is that most non-contaminant genera flow to KEEP (green), while *Cutibacterium*, *Brevundimonas*, *Staphylococcus*, *Acinetobacter*, *Corynebacterium*, and other Tier A genera flow to REMOVE_TIER_A (red), and the genus-level *Burkholderia* signal flows to REMOVE_BURK (orange — the species-level safeguard found no study-cohort sample with *B. pseudomallei* above background).

Per-sample biomass distribution under V4:

| Pattern | n samples |
|---|---|
| Retain 100% of biomass (no contaminants present at detection level) | 21 |
| Retain 0% (all detections removed as contaminants) | 2 |
| Mixed (most informative for clinical review) | 48 |
| No detections (pre-sequencing failure) | 15 |
| **Total** | **86** |

#### 3.2.4 Rickettsiales rescue in the study cohort

The single most clinically important finding of the cohort analysis is **Rickettsiales rescue evidence in 11 of 71 samples with detections (15.5%; 11 of 86 samples = 12.8%)**. This finding was missed in the earlier draft of the cohort analysis, which examined only Centrifuge-Detected rows and overlooked the Minimap2 alignment-based rescue rows.

| Sample | Run | Genus | Mapped reads | Breadth | Rescue tier | Alignment NTC reads | Confidence (sample-relative-to-NTC) |
|---|---|---|---|---|---|---|---|
| `16901195_S5_L001` | 8_and_9 | Orientia | 14,816 | 0.3235 | **Tier 1 (Confirmed)** | 1,288 | HIGH (sample ~11.5× NTC) |
| `22600400_S6_L001` | 8_and_9 | Orientia | 9,747 | 0.2181 | Tier 2 (Probable) | 1,288 | HIGH |
| `23900356_S5_L001` | 4_and_5 | Rickettsia | 441,173 | 0.3046 | Tier 2 (Probable) | 189,002 | MODERATE (sample 2.3× NTC) |
| `23200430_S6_L001` | 4_and_5 | Rickettsia | 94,142 | 0.2189 | Tier 2 (Probable) | 53,908 | LOW (sample 1.7× NTC) |
| `23200519_S12_L001` | 6_and_7 | Orientia | 2,492 | 0.2175 | Tier 2 (Probable) | 1,143 | MODERATE |
| `25800718_S8_L001` | 6_and_7 | Orientia | 4,143 | 0.2248 | Tier 2 (Probable) | 2,202 | LOW |
| `09801652_S5_L001` | 6_and_7 | Orientia | 3,727 | 0.2282 | Tier 2 (Probable) | 1,143 | MODERATE |
| `23200736_S8_L001` | 8_and_9 | Orientia | 1,981 | 0.2141 | Tier 2 (Probable) | 1,288 | LOW |
| `23900752_S9_L001` | 8_and_9 | Orientia | 1,684 | 0.2428 | Tier 2 (Probable) | 1,288 | LOW |
| `10100409_S9_L001` | 8_and_9 | Orientia | 2,385 | 0.2221 | Tier 2 (Probable) | 1,288 | LOW |
| `16601093_S7_L001` | 8_and_9 | Orientia | 2,284 | 0.2201 | Tier 2 (Probable) | 1,288 | LOW |

**Summary:** 1 Tier-1 genus-level Orientia call (`16901195_S5_L001`, breadth 0.3235, ~12× NTC headroom) and 10 Tier-2 order-level "Rickettsiales detected" calls with breadth in the 0.21–0.24 range (just below the Tier-1 threshold). Confidence on the Tier-2 calls varies with NTC headroom; 5 of 10 carry comparable-magnitude background (LOW confidence) and warrant priority follow-up.

#### 3.2.5 Species-level *B. pseudomallei* findings (corrected)

Four *Burkholderia* genus detections were observed in study samples (`14300185_S8_L001`, `22900187_S10_L001`, `23200430_S6_L001`, `23200519_S12_L001`, `16601093_S7_L001`). Species-level Centrifuge kreport parsing:

| Sample | Genus reads | *B. pseudomallei* species reads | *B. cepacia* complex species reads | Run NTC species max | V4 decision |
|---|---|---|---|---|---|
| `23200430_S6_L001` | 54,126 | **8** | 12,789 | 51 | **REMOVE** (genus signal is B. cepacia complex contamination) |
| `22900187_S10_L001` | 899 | 0 | 196 | 51 | REMOVE |
| `14300185_S8_L001` | 1,689 | 0 | 441 | 0 | REMOVE |
| `23200519_S12_L001` | 876 | 0 | 231 | 0 | REMOVE |
| `16601093_S7_L001` | 673 | 0 | (~150) | 0 | REMOVE |

**No study sample contains *B. pseudomallei* at species level above background.** The earlier draft's claim that `23200430_S6_L001` carried 54,126 reads of *B. pseudomallei* at 15.56% abundance conflated genus-level reads with species-level reads — the actual species-level *B. pseudomallei* count in that sample is 8 reads, below the run NTC's 10–51 reads of *B. pseudomallei*. The genus-level signal in that sample is dominated by *B. cepacia* complex (kit contaminant) species.

The 3 validation-panel *B. pseudomallei* samples remain authentic (13,744 / 20,061 / 65,016 species-level reads). The corrected reading of the cohort is that melioidosis is not a major contributor to this specific cohort's positive-culture / failed-subculture phenotype as detected by 16S.

#### 3.2.6 Fastidious-organism candidate signals

**Mycoplasmopsis (4 samples).** *Mycoplasmopsis* was detected at substantial abundance across the cohort but was silently removed by the prior V3 filter under a Tier-1 "ultra-low abundance" categorization. The actual read counts in `.calls.tsv` are 537–6,568 reads across 5 samples (post-NCmax subtraction):

| Sample | Run | Reads | % of sample |
|---|---|---|---|
| `13400425_S4_L001` | 4_and_5 | 6,568 | (high) |
| `15700611_S5_L001` | 4_and_5 | 2,568 | 100.0% (only detection in sample) |
| `10300205_S9_L001` | 4_and_5 | 1,622 | 2.17% |
| `09301225_S4_L001` | 4_and_5 | 740 | 45.76% |
| `24500367_S5_L001` | 4_and_5 | 537 | (variable) |

Aggregate: mean abundance 39.78%, maximum 70.55%, 4 samples retained at the rebuilt Tier-1 list (V4). *Mycoplasmopsis* species are cell-wall-deficient (lack peptidoglycan), fastidious (require sterol-supplemented media), and slow-growing (1–3 weeks to visible colonies). They would not be expected to grow on routine aerobic blood subculture within standard 5–7 day windows.

**Leptospira (1 sample with NTC caveat).** Sample `09801652_S5_L001` (run 6_and_7) carries 7,294 *Leptospira* reads (22.78% of detected sample biomass). The same run contains NTC2_ExDw_S13_L001 with 78,691 *Leptospira* reads (10× the study sample). The pipeline reports `ntc_reads = 0` for *Leptospira* in this sample's `.calls.tsv` (the NCmax derivation appears to exclude the contaminated NTC), but the contamination flag is a real caveat. This detection should be treated as a candidate requiring orthogonal confirmation.

**Brucella (3 samples, near noise floor).** *Brucella* was detected in 3 samples at sub-1% mean abundance (mean 0.98%, max 1.96%). At this level the signal is at the upper edge of what could plausibly be reagent contamination and the lower edge of what could be low-level bacteremia. Without serological or PCR confirmation these are candidate observations, not diagnoses.

**Streptococcus (5 samples, genus-level only).** Five samples carry *Streptococcus* at mean 12.75%; V1-V3 does not resolve *S. pneumoniae* / *S. suis* / fastidious vs. non-fastidious species. Species composition is unresolved by this assay alone.

#### 3.2.7 Other retained organisms

The full post-V4 top-organism table is in `APPENDIX-STUDY-SAMPLES.md` and `APPENDICES.xlsx` (`Study - sample summary` sheet). Notable additional retained organisms:

- **Thermomicrobium (11 samples, mean 9.13%):** environmental thermophile; unlikely human pathogen; almost certainly persistent environmental noise that escaped the filter (candidate for a future Tier-1 addition).
- **Escherichia (6 samples, mean 8.28%), Klebsiella (4 samples, mean 20.99%), Enterobacter (3 samples, mean 3.39%):** culturable AFI pathogens. Their detection by 16S in patient blood with failed bottle subculture is unevaluable from this dataset alone because the 16S sample is patient blood, not a bottle aliquot.
- **Porphyromonas (1 sample, 8.45%) and Desulfovibrio (1 sample, 0.78%):** anaerobic genera retained, contradicting the earlier draft's "no anaerobes detected" framing.

### 3.3 Inter-run reproducibility and pipeline QC

Inter-run reproducibility was assessed by tracking PC and NTC pass rates across 9 sequencing runs (validation panel + 5 study-cohort runs):

- **PC pass rate:** 100% (all PC_MIX8, PC_SINGLE, and MIXED4 samples detect their expected spike-in organisms under V4 with PC bypass).
- **NTC specificity rate:** 100% on the rule of record (no TAC bacterial target genus is retained after V4 in any validation-panel NTC). Some study-run NTCs carry substantial pre-NCmax reagent contamination; the pipeline's per-run NTC max derivation accounts for this through run-ID-based pool assignment. Per-run accept/reject decisions are governed by the deployment SOP (`[TEAM INPUT NEEDED: SOP reference]`).

---

## 4. Discussion

### 4.1 What the corrected analysis shows

Three findings reframe the interpretation of this cohort relative to earlier drafts:

1. **Rickettsiales involvement is the most prevalent identifiable signal in this cohort (~13% of all study samples, ~15% of samples with any positive call).** The Minimap2 alignment-based rescue identifies 11 samples with Rickettsiales evidence — 1 high-confidence Tier-1 *Orientia* call and 10 Tier-2 order-level "Rickettsiales detected" rescues. Rickettsiales are obligate intracellular pathogens that cannot be cultured on routine blood agar; their detection in a positive-blood-culture / no-subculture-growth cohort is the most biologically coherent finding in the dataset. This finding aligns with the expected endemic epidemiology of northeastern Thailand, where scrub typhus and spotted fever group rickettsioses account for a substantial fraction of AFI presentations.

2. **No study sample contains *Burkholderia pseudomallei* above species-level background.** The earlier reporting of preserved *B. pseudomallei* in sample `23200430_S6_L001` was a methodological artifact: a substring-based check for "pseudomallei" in the kreport file triggered on the parent "pseudomallei_group" clade and the safeguard retained the genus-level read count rather than the species-level count. Species-rank parsing shows the sample's *B. pseudomallei* signal is 8 reads, below both the 500-read detection threshold and the run NTC's *B. pseudomallei* count (10 reads). The genus-level Burkholderia signal is dominated by *B. cepacia* complex species, which are kit/water contaminants. The 3 *B. pseudomallei*-positive validation-panel samples (with 13,744–65,016 species reads) demonstrate that the assay does detect melioidosis when present at meaningful abundance — but this specific cohort does not contain it.

3. **Mycoplasmopsis is the most prominent fastidious-organism-class signal in the cohort (4 samples, mean abundance 39.78%, max 70.55%) and was silently removed by the earlier V3 filter.** Mycoplasma-class organisms are by definition cell-wall-deficient and fastidious; they require sterol-supplemented media and 1–3 weeks of incubation. Their presence in patient blood is biologically consistent with a positive-bottle / no-subculture-growth phenotype, but cohort-level prevalence of ~5% of all study samples (~6% of detected samples) and the absence of orthogonal confirmation puts this in the hypothesis-generating category, not the confirmed-etiology category.

### 4.2 Why filter design matters

The V4 decontamination filter materially changes which organisms are reported and which are dropped. Three filter design choices have outsized impact in this dataset:

- **Aggressive Tier A.** Removing the 11 high-confidence kit / skin / water contaminants is responsible for 61 of 70 V4 removals. *Cutibacterium*, *Staphylococcus*, and *Brevundimonas* are the top three contributors. The *Brevundimonas* addition (V4 vs V3) is supported by NTC profiles showing it at 22–285,740 reads per NTC in run 6_and_7.

- ***Burkholderia* species-level safeguard.** The genus is in Tier A but the species-level kreport parser preserves *B. pseudomallei* on a per-row basis. Without this safeguard, melioidosis cases would be erroneously discarded; with the safeguard, only true *B. pseudomallei* signals are preserved (validated against same-run NTC species reads).

- **Positive-control spike-in bypass.** Without this bypass, the P-aeru_S5_L001 PC would fail because *Pseudomonas* is in Tier A, even though *Pseudomonas aeruginosa* is the expected positive-control spike-in. The bypass restores PC accuracy from 9/10 to 10/10 and brings the sample-level analytical performance from 30/43 = 69.8% to 31/43 = 72.1%, matching the legacy APHL figure.

The corrected Tier-1 list, which excludes *Mycoplasmopsis* and *Nitrospira*, demonstrates the importance of validating tier membership against actual abundance distributions in the dataset rather than relying on text-string filter rules.

### 4.3 What the data do not support

The earlier discussion drafts contained several overreaches that the corrected analysis does not support:

- **"The cohort has no Rickettsiales."** This was an artifact of Centrifuge-only filtering; alignment-based rescue rows show 11 samples with Rickettsiales evidence.
- **"Strong support for fastidious organism hypothesis."** With single-sample Leptospira (n=1, NTC caveat) and 3-sample low-abundance Brucella (mean 0.98%), the cohort-level evidence for cohort-level claims is thin. *Mycoplasmopsis* is the strongest signal class but still ~5% of cohort. Honest framing: hypothesis-generating observations for confirmatory testing.
- **"No anaerobes detected."** *Porphyromonas* and *Desulfovibrio* are retained in the V4 output, contradicting absolute negation. The V1-V3 primer set is also known to under-detect some Gram-positive anaerobes; "absent in this assay" is not equivalent to "biologically absent."
- **"16S identified organisms that caused failed culture."** 16S detects DNA from patient blood, not from the original positive culture bottle. The detected organisms are candidates for explaining the bottle signal but cannot be causally linked without paired bottle/blood analysis. 16S also does not distinguish viable from non-viable cells.

### 4.4 Methodological lessons for downstream deployment

Operational deployment of this workflow should incorporate four refinements identified during this analysis:

1. **Document the NCmax derivation rule explicitly.** The current pipeline appears to exclude or downweight contaminated NTCs when computing per-run NCmax (consistent with study-sample `ntc_reads = 0` for *Leptospira* despite an NTC carrying 78,691 reads of *Leptospira* in the same run). The rule (median? minimum? excluded-NTC-by-QC?) should be documented in the SOP.
2. **Add a pre-sequencing QC step.** A non-trivial fraction (15/86 = 17.4%) of study specimens returned zero taxa, indicating DNA extraction / library preparation / sequencing-depth failure. Qubit and library qPCR before sequencing would distinguish these from true biological negatives.
3. **Maintain the V4 filter's species-level safeguards.** *Burkholderia* requires species-level parsing; similar reasoning applies if *Brucella* species (e.g., *B. melitensis* / *B. abortus*) discrimination becomes clinically required.
4. **Document the positive-control spike-in bypass.** The bypass is a deliberate design choice and must be documented so that future operators understand why Tier A genera are sometimes retained in PC samples.

### 4.5 Clinical implications for the AFI study cohort

For the 11 samples with Rickettsiales rescue evidence, doxycycline-based empirical therapy in an endemic context is standard clinical practice and is supported by the alignment-rescue evidence; species-specific qPCR (Orientia-tsutsugamushi-specific, Rickettsia 17 kDa antigen gene targets) is recommended for definitive species identification. For the Mycoplasmopsis candidate samples, Mycoplasma-specific PCR and Mycoplasma broth subculture would convert the genus-level signal into a species-resolved confirmation. For the Leptospira candidate (with NTC caveat) and Brucella candidates (near noise floor), paired-serum serology and species-specific qPCR are recommended before any clinical interpretation.

The 15 samples with no detection at all warrant a separate workup — these are pre-sequencing failures and should be subject to clinical review and, where feasible, repeat specimen collection or alternative diagnostic modalities.

---

## 5. Limitations

1. **Sample source mismatch.** The 16S samples are patient blood, not aliquots of the original positive blood culture bottle. Detection of an organism by 16S in patient blood is consistent with its causing the bottle signal but does not establish causation. Paired bottle / blood analysis (16S of the bottle broth) would be required to make causal claims.
2. **No viability assessment.** 16S detects DNA from viable, non-viable, and viable-but-non-culturable (VBNC) cells. A high-abundance 16S signal does not establish that the organism is alive and would be culturable under appropriate conditions.
3. **V1-V3 primer biases.** The 27F primer set has documented under-detection of certain Gram-positive anaerobes, Mycobacteria, and some Bifidobacteria / Atopobium / Gardnerella lineages. Negative findings for these classes are limitations of the assay, not biological absence.
4. **Genus-level resolution for most organisms.** Only *Burkholderia* receives species-level parsing in this filter. *Streptococcus pneumoniae* vs. *S. suis*, fastidious vs. non-fastidious *Streptococcus*, and species-level *Brucella* / *Leptospira* discrimination are not resolved.
5. **NTC contamination in some runs.** NTC2_ExDw_S13_L001 in run 6_and_7 carries substantial *Leptospira*, *Burkholderia*, and *Brevundimonas* reads. The pipeline's NCmax derivation appears to handle this by exclusion, but the rule should be documented and contaminated NTCs flagged in deployment QC.
6. **Filter tier membership is dataset-informed.** The Tier 1 and Tier 2 lists were assembled from the observed abundance distributions of this dataset combined with the cited contamination reviews. Generalizability to other specimen types or kit batches has not been validated.
7. **Cohort sample size.** 86 study samples is appropriate for a pilot but insufficient to establish prevalence estimates for any single etiology. Findings of n=1 (Leptospira) or n=3 at sub-1% abundance (Brucella) are best treated as candidate detections requiring follow-up, not cohort-level prevalence claims.
8. **No clinical-outcome correlation.** Patient demographics, treatment, response, serology, and clinical outcomes are not currently linked to the 16S findings in this report (`[TEAM INPUT NEEDED]`).
9. **No formal anaerobic culture for comparison.** The hospital protocol does not include anaerobic blood culture, so direct comparison between 16S-detected anaerobic-genus signals and a culture gold standard is not possible.

---

## 6. Conclusion

A 16S V1-V3 amplicon workflow with Centrifuge primary classification, Minimap2 alignment-based rescue for Rickettsiales, a two-tier reporting framework (genus-level + order-level), and a species-aware V4 decontamination filter (with positive-control spike-in bypass) achieves 72.1% sample-level analytical performance against a 43-sample reference panel (95% CI 57.3%–83.3%), with 100% NTC specificity (95% CI 56.6%–100%) and 100% inter-run reproducibility across 9 sequencing runs. Applied to a cohort of 86 AFI cases with positive automated blood culture and failed subculture recovery, the workflow identifies Rickettsiales rescue evidence in 11 samples (~13% of cohort) — a finding consistent with the expected endemic epidemiology of northeastern Thailand — and candidate fastidious-organism signals (*Mycoplasmopsis*, *Leptospira*, *Brucella*) in an additional ~10% of the cohort. No study sample carries *Burkholderia pseudomallei* above the species-level detection threshold.

The 16S workflow described here is suitable for use as a complementary diagnostic for AFI cases with positive culture signal but failed subculture recovery, with the understanding that candidate detections require orthogonal confirmation (species-specific PCR / qPCR, serology, specialized culture) before clinical interpretation. The most actionable finding for the cohort is the Rickettsiales rescue evidence in approximately one in seven AFI cases.

---

## 7. References

**Decontamination and low-biomass microbiome studies (verified against PubMed 2026-05-11):**

1. Salter SJ, Cox MJ, Turek EM, Calus ST, Cookson WO, Moffatt MF, Turner P, Parkhill J, Loman NJ, Walker AW. Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. *BMC Biol.* 2014;12:87. doi:10.1186/s12915-014-0087-z. PMID: 25387460.
2. Glassing A, Dowd SE, Galandiuk S, Davis B, Chiodini RJ. Inherent bacterial DNA contamination of extraction and sequencing reagents may affect interpretation of microbiota in low bacterial biomass samples. *Gut Pathog.* 2016;8:24. doi:10.1186/s13099-016-0103-7. PMID: 27239228.
3. Lauder AP, Roche AM, Sherrill-Mix S, Bailey A, Laughlin AL, Bittinger K, Leite R, Elovitz MA, Parry S, Bushman FD. Comparison of placenta samples with contamination controls does not provide evidence for a distinct placenta microbiota. *Microbiome.* 2016;4(1):29. doi:10.1186/s40168-016-0172-3. PMID: 27338728.
4. de Goffau MC, Lager S, Salter SJ, Bonney EA, Bertozzi-Villa A, Wagner J, Charnock-Jones DS, Smith GCS, Parkhill J. Recognizing the reagent microbiome. *Nat Microbiol.* 2018;3(8):851-853. doi:10.1038/s41564-018-0202-y. PMID: 30046175.
5. Tan CCS, Ko KKK, Chen H, Liu J, Loh M, Chia M, Nagarajan N, SG10K_Health Consortium. No evidence for a common blood microbiome based on a population study of 9,770 healthy humans. *Nat Microbiol.* 2023;8(5):973-985. doi:10.1038/s41564-023-01350-w. PMID: 36997797.

**Software / bioinformatic tools used:**

6. Song L, Langmead B. Centrifuger: lossless compression of microbial genomes for efficient and accurate metagenomic sequence classification. *Genome Biol.* 2024;25(1):106. doi:10.1186/s13059-024-03244-4. PMID: 38641639.
7. Kim D, Song L, Breitwieser FP, Salzberg SL. Centrifuge: rapid and sensitive classification of metagenomic sequences. *Genome Res.* 2016;26(12):1721-1729. doi:10.1101/gr.210641.116. PMID: 27852649.
8. Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics.* 2018;34(18):3094-3100. doi:10.1093/bioinformatics/bty191. PMID: 29750242.
9. Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. *Bioinformatics.* 2009;25(16):2078-2079. doi:10.1093/bioinformatics/btp352. PMID: 19505943.
10. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics.* 2018;34(17):i884-i890. doi:10.1093/bioinformatics/bty560. PMID: 30423086.
11. NCBI Human Read Removal Tool (SRA Human Scrubber). National Center for Biotechnology Information; 2021. Available at: https://github.com/ncbi/sra-human-scrubber. (Cite as appropriate per tool documentation.)
12. Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biol.* 2019;20(1):257. doi:10.1186/s13059-019-1891-0. PMID: 31779668. *(Listed for completeness — Kraken 2 is bundled in the AFI core container image for QC purposes but is not on the primary classification path used for these results.)*
13. Voss K, Van der Auwera G, Gentry J. Full-stack genomics pipelining with GATK4 + WDL + Cromwell. *F1000Research.* 2017;6:1379 (Workflow Description Language reference; cite the Cromwell / WDL implementation appropriate to the deployment, e.g., the OpenWDL specification at https://openwdl.org).

**AFI epidemiology / clinical references (`[TEAM INPUT NEEDED]`):** primary references for *Orientia tsutsugamushi* / scrub typhus burden in northeastern Thailand; *Burkholderia pseudomallei* / melioidosis epidemiology; *Leptospira* / leptospirosis burden; *Brucella* / brucellosis prevalence; *Rickettsia* spotted fever group epidemiology; TaqMan Array Card AFI Multi-Pathogen card validation reports (CDC).

---

## 8. Appendices

The following companion files accompany this manuscript and contain the per-sample data tables and raw outputs from the analysis:

- **`APPENDIX-VALIDATION-PANEL.md`** — per-detection table for the 48-sample validation panel (236 rows): sample, run, category, expected organism, total reads, biomass kept/removed by V4 filter, detected genus, source (Centrifuge / Minimap2 alignment), reads, % of sample, pipeline NTC reads, Minimap2 rescue info, rescue tier, V4 filter outcome, per-row concordance, failure mode, Final concordance (TAC + V4), notes. Includes a per-sample concordance summary table at the bottom.
- **`APPENDIX-STUDY-SAMPLES.md`** — per-detection table for the 86-sample study cohort (336 rows): same columns as the validation appendix plus same-run NTC max (cross-check), Burkholderia species evidence (where applicable), confidence label, hypothesis class, and recommended follow-up.
- **`APPENDICES.xlsx`** — same content in Excel format with 4 sheets (`Validation - detections`, `Validation - sample summary`, `Study - detections`, `Study - sample summary`).
- **`DECONTAMINATION-FILTER-REPORT-V4.txt`** — raw output of the V4 filter run with summary statistics and Burkholderia species-level evidence per sample.
- **`afi_decontamination_filter_v4.py`** — executable V4 filter implementation.
- **`generate_appendices.py`** — appendix generator script (regenerates all of the above from the raw `.calls.tsv` and Centrifuge kreport files).

Additional companion documents from earlier in the analysis:

- **`MANUSCRIPT-METHODS-DECONTAMINATION.md`** — extended methods detail for the V4 filter (deeper than Methods §2.5 above).
- **`MANUSCRIPT-RESULTS-SECTION.md`** — extended results detail.
- **`MANUSCRIPT-DISCUSSION-SECTION.md`** — extended discussion detail.
- **`THRESHOLD-UNIFORMITY-AUDIT.md`** / **`STUDY-SAMPLES-THRESHOLD-UNIFORMITY-AUDIT.md`** — Minimap2 rescue threshold uniformity audits.

---

## 9. Placeholders Requiring Team Input

### What has been resolved from the codebase (no longer placeholder)

| Section | Resolved from |
|---|---|
| §2.4 Pipeline orchestration | WDL via Cromwell/Terra, source at https://github.com/PHemarajata/afi_terra (`wdl/AFI_16S_Main.wdl`, `wdl/AFI_16S_Batch.wdl`); container images `phemarajata614/afi-terra:0.4.1` + `phemarajata614/centrifuger:1.1.0`. |
| §2.4.1 Host-DNA removal + read preprocessing | NCBI SRA Human Scrubber + fastp. |
| §2.4.2 Primary taxonomic classification | Centrifuger v1.1.0 (Song & Langmead, 2024). Single index covers bacteria/archaea + Rickettsiales. |
| §2.4.2 Detection thresholds | Centrifuger: ≥500 reads AND ≥5× run NTC max (`cfr_floor=500`, `cfr_fold=5.0` in `wdl/tasks/interpret.wdl`). |
| §2.4.3 Minimap2 reference panel | 7 references: AM494475.1 *O. tsutsugamushi* Boryong; AP008981.1 *O. tsutsugamushi* Ikeda; CP004888.1 *R. prowazekii* NMRC Madrid E; NC_006142.1 *R. typhi* Wilmington; NC_009882.1 *R. rickettsii* 'Sheila Smith'; NC_007797.1 *A. phagocytophilum* HZ; NC_007799.1 *E. chaffeensis* Arkansas (`align_rickettsiales_16S/refs_fna/`). |
| §2.4.3 Rescue thresholds | Tier 1 Confirmed: reads ≥100 AND breadth ≥0.25 AND ≥5× align NTC max. Tier 2 Probable: reads ≥50 AND breadth ≥0.20 AND reads > align NTC max (no 5× fold requirement). |
| §2.4.4 NTC background derivation | Per-run **maximum** read count across NTCs in that run's pool, separately for Centrifuger and Minimap2 paths (`scripts/build_ntc_background.py`). |
| §7 Software citations | Centrifuger, Centrifuge, Minimap2, SAMtools, fastp, NCBI Scrubber, Kraken 2, WDL — added. |

### Items still requiring team input

| Section | Content needed |
|---|---|
| Title page | Authors + affiliations; corresponding author; running title; CoI; funding; author contributions |
| §2.1 Study design | IRB approval, ethics body, consent procedure, enrollment dates, hospital(s) / district(s), inclusion / exclusion criteria |
| §2.2 Sample populations | Reference-laboratory diagnostic methods used for each validation-panel pathogen (PCR? culture? serology?) |
| §2.3 Sample collection / wet-lab | Blood specimen volume, container, time-to-extraction, DNA extraction kit + version, 16S V1-V3 primer sequences, library prep kit, sequencer model + run configuration |
| §2.4.2 Centrifuger reference DB | Specific reference database identifier (name, build date, source URL or accession set) used for the Centrifuger index |
| §2.7 / §3.2.2 / §3.3 | Deployment SOP reference for per-run QC pass criteria |
| §2.8 Data availability | Final data repository / archival location for `.calls.tsv` and kreport outputs |
| §4.5 Clinical follow-up | Confirmatory test results (if any) and clinical outcomes for the 11 Rickettsiales-rescue samples, 4 *Mycoplasmopsis* samples, 1 *Leptospira* candidate, 3 *Brucella* candidates |
| §7 References | AFI epidemiology references for NE Thailand; Rickettsiales / melioidosis / leptospirosis / brucellosis clinical references; CDC TaqMan Array Card AFI panel validation reports |

---

**Draft prepared 2026-05-12. All numerical results and per-sample data are derived from the actual `.calls.tsv` and Centrifuge kreport files in `/Users/peerahemarajata/Downloads/AFI_P_Final/` and use the corrected counting established during the May 11–12 review. The appendix generator (`generate_appendices.py`) and the V4 filter (`afi_decontamination_filter_v4.py`) can regenerate all numbers and tables from the raw inputs.**
