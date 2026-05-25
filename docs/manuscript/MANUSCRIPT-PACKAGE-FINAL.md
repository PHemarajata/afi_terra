# Detection of Bacterial Pathogens in Acute Febrile Illness Using a 16S rRNA V1-V3 Workflow with Two-Tier Rickettsiales Rescue and Species-Aware Decontamination

**Manuscript package (integrated draft).** Prepared 2026-05-25. All numerical results and per-sample data are derived from the `.calls.tsv` and Centrifuger kreport files in `/Users/peerahemarajata/Downloads/AFI_P_Final/`. The appendix generator (`generate_appendices.py`) and the V4 filter (`afi_decontamination_filter_v4.py`) can regenerate all numbers and tables from the raw inputs.

---

## Title page

**Authors:** `[TEAM INPUT NEEDED]` — author list with affiliations (NIH wet-lab team for sequencing operations; CDC for funding/sponsorship; northeastern Thailand clinical and epidemiology team for cohort enrolment and clinical metadata; bioinformatics author Peera Hemarajata).

**Affiliations:** `[TEAM INPUT NEEDED]`

**Corresponding author:** `[TEAM INPUT NEEDED]` — name, email, postal address.

**Running title (≤50 chars):** `[TEAM INPUT NEEDED]`

**Keywords (5–8):** acute febrile illness; 16S rRNA; Rickettsiales; blood culture; melioidosis; decontamination; metagenomics; Thailand

**Word count (body):** ~9,000.

**Conflict-of-interest statement:** `[TEAM INPUT NEEDED]`

**Funding:** `[TEAM INPUT NEEDED — CDC funding to be confirmed by NIH wet-lab team and CDC programme office]`

**Data availability statement:** Pipeline source code is publicly available at https://github.com/PHemarajata/afi_terra. All bioinformatics outputs (`.calls.tsv`, Centrifuger kreport files), the V4 decontamination filter script (`afi_decontamination_filter_v4.py`), the appendix generator (`generate_appendices.py`), and the per-sample appendix tables (`APPENDIX-VALIDATION-PANEL.md`, `APPENDIX-STUDY-SAMPLES.md`, `APPENDICES.xlsx`) are available at `[TEAM INPUT NEEDED: final repository / Google Drive archival location]`. Raw FASTQ data deposition: `[TEAM INPUT NEEDED — SRA accession if deposited]`.

**Author contributions:** `[TEAM INPUT NEEDED — CRediT taxonomy]`

---

## Abstract

**Background.** Acute febrile illness (AFI) in endemic regions of northeastern Thailand is caused by a wide differential including obligate intracellular Rickettsiales (*Orientia tsutsugamushi*, *Rickettsia* spp), *Burkholderia pseudomallei* (melioidosis), *Leptospira* spp, *Brucella* spp, dengue and other arboviruses, and enteric bacterial pathogens, several of which culture cannot reliably recover. Blood cultures that signal positive on automated detection but fail subculture recovery are a recognised diagnostic puzzle.

**Methods.** We developed and validated a 16S rRNA V1-V3 amplicon workflow combining Centrifuger v1.1.0 primary taxonomic classification with Minimap2 alignment-based rescue against a curated 7-genome Rickettsiales reference panel, a two-tier reporting framework that distinguishes confirmed genus-level calls (Tier 1) from order-level "Rickettsiales detected" rescues (Tier 2), and a multi-tiered post-pipeline decontamination filter (V4) that includes a species-level safeguard for *Burkholderia pseudomallei* and a positive-control spike-in bypass. The pipeline is implemented in WDL on Terra.bio / Cromwell with Google Cloud Batch orchestration. The assay was validated against a 48-sample panel (33 clinical samples with known reference-laboratory diagnoses, 10 positive controls, 5 negative template controls) and applied to a study cohort of 86 AFI cases with positive automated blood culture signal but failed subculture recovery.

**Results.** The validated assay achieves 72.1% (31/43, 95% CI 57.5%–83.6%) sample-level analytical performance combining clinical and positive-control samples — 63.6% (21/33, 95% CI 46.0%–78.5%) strict clinical concordance, 100% (10/10) positive-control concordance — with 100% specificity (5/5 NTCs clear of all TaqMan Array Card bacterial target genera) and 100% inter-run reproducibility across nine sequencing runs. Applied to the AFI study cohort (n=86; 71 with ≥1 positive call, 15 zero-detection pre-sequencing failures), the workflow identified Rickettsiales rescue evidence in **11 of 86 samples (12.8%; 11/71 active samples = 15.5%)** — 1 Tier-1 genus-level *Orientia* detection (sample 16901195_S5_L001, breadth 0.3235, ~12× NTC headroom) and 10 Tier-2 order-level rescues with breadth-of-coverage 0.21–0.24; *Mycoplasmopsis* (a fastidious, cell-wall-deficient organism class) in 4 samples (mean abundance 39.78%, max 70.55%); a candidate *Leptospira* detection in 1 sample (22.78% abundance, with a same-run NTC contamination caveat); and *Brucella* signals at low abundance (0.98% mean) in 3 samples. **No study sample contains *Burkholderia pseudomallei* above the species-level detection threshold.** The V4 filter removed 70 of 217 detections (32.3%) from the cohort, predominantly Tier-A kit/skin/water contaminants.

**Conclusion.** The 16S V1-V3 assay, with the two-tier Rickettsiales framework and the species-aware V4 decontamination filter, is suitable for analytical surveillance of bacterial pathogens in AFI cases. Approximately one in eight AFI cases with positive culture / failed subculture shows Rickettsiales rescue evidence; an additional small minority shows fastidious-organism candidate signals (*Mycoplasmopsis*, *Leptospira*, *Brucella*). All candidate detections require orthogonal confirmation (PCR, serology, specialised culture) before clinical reporting. The cohort study is framed as pilot, hypothesis-generating data, with Rickettsiales involvement as the most prevalent identifiable etiology.

---

## 1. Introduction

Acute febrile illness (AFI) is a common presentation in emergency and inpatient settings in tropical Southeast Asia, with a wide etiologic differential and substantial overlap in early clinical features. In endemic northeastern Thailand, the major bacterial pathogens contributing to AFI are *Orientia tsutsugamushi* (scrub typhus), *Rickettsia* spp (spotted-fever and typhus-group rickettsioses), *Burkholderia pseudomallei* (melioidosis), *Leptospira* spp (leptospirosis), *Brucella* spp (brucellosis), non-typhoidal *Salmonella* and *Salmonella* Typhi, enteropathogenic *Escherichia coli*, and species-level *Streptococcus* (notably *S. pneumoniae* and *S. suis*). These pathogens have substantially different culture requirements: some grow readily on standard blood agar within a 5–7 day incubation window (*E. coli*, *Klebsiella pneumoniae*, *Staphylococcus aureus*), some require specialised media or extended incubation (*Leptospira*, *Brucella*, fastidious *Streptococcus* species), and some cannot be recovered by routine blood culture at all (Rickettsiales, which are obligate intracellular).

Blood cultures that are flagged positive by automated detection but fail subculture recovery represent a particularly common and diagnostically frustrating scenario in AFI workups. The automated detection system records a positive signal (turbidity, CO₂ production, metabolic activity), but inoculation of the flagged bottle onto solid media yields no growth within the standard 24–48 hour subculture observation window. The clinical implication is that *some* organism was present and metabolically active in the bottle, but its identity remains unknown and empirical therapy must proceed on clinical judgement alone. In the hospital protocol relevant to this cohort, only aerobic blood culture is performed; no anaerobic workup is routinely available.

16S rRNA gene sequencing has been proposed as a complementary diagnostic for these positive-culture / no-growth scenarios because it detects bacterial DNA regardless of culturability and provides genus-level (and sometimes species-level) taxonomic assignment within the limits of the targeted hypervariable region. The V1-V3 region (amplified with the 27F primer set) is widely used in clinical 16S workflows but has known under-detection biases for certain Gram-positive anaerobes and primer-mismatch organisms (notably some *Bifidobacterium*, *Gardnerella*, and *Atopobium* lineages). It also provides limited genus-level discrimination within order Rickettsiales — meaning that even when sequence is recovered, distinguishing *Orientia* from *Rickettsia* may not be possible from the V1-V3 amplicon alone.

In low-biomass clinical samples such as blood, sequencing-based bacterial detection is further complicated by reagent and laboratory contamination. Bacterial DNA introduced by extraction kits, molecular-grade water, and PCR reagents (Salter et al., 2014; Glassing et al., 2016; Lauder et al., 2016; de Goffau et al., 2018) routinely accounts for a substantial fraction of reads in blood and other low-biomass specimens and can dominate the signal entirely when the input bacterial biomass is low. A population study of 9,770 healthy human blood samples found no consistent core blood microbiome (Tan et al., 2023), strengthening the position that any organism detected by sequencing in blood should be evaluated rigorously against contamination controls before clinical interpretation.

We developed and validated a 16S V1-V3 amplicon workflow for AFI bacterial-pathogen surveillance that addresses these constraints through four design choices: (1) **Centrifuger primary classification** for breadth across all bacterial genera, against a single combined index; (2) **Minimap2 alignment-based rescue** against a curated Rickettsiales 16S reference panel to compensate for V1-V3 limited genus-level discrimination; (3) a **two-tier Rickettsiales reporting framework** that distinguishes confirmed genus-level calls (Tier 1 "Confirmed") from order-level "Rickettsiales detected" rescues (Tier 2 "Probable") — both clinically actionable in a doxycycline-treatable endemic context; and (4) a **species-aware decontamination filter (V4)** that aggressively removes documented kit / skin / water contaminant genera while implementing species-level safeguards for clinically critical organisms (*B. pseudomallei*) and a positive-control spike-in bypass to allow PC samples to be evaluated for their expected organisms.

This report presents the analytical validation of the workflow against a 43-sample reference panel (33 clinical samples with reference-laboratory diagnoses plus 10 positive controls) and 5 negative template controls, and the application of the validated workflow to a study cohort of 86 AFI cases drawn from positive blood cultures with failed subculture recovery. We report performance metrics, identify the candidate organisms detected in the cohort, and discuss the implications for AFI diagnostic algorithms in endemic regions.

---

## 2. Methods

### 2.1 Study design and ethical approval

`[TEAM INPUT NEEDED — EPI team]` Ethical approval, oversight body, informed-consent procedure, period of enrolment, hospital(s) and geographic catchment (northeastern Thailand), inclusion criteria (e.g., patients presenting with fever ≥38 °C of <14 days duration, without localising signs identifying an alternative diagnosis), and exclusion criteria. Clinical metadata, treatment history, serology and follow-up outcomes for individual cases are held by the EPI team and will be merged with the bioinformatics findings reported here in the clinical / epidemiology sections of the final manuscript.

### 2.2 Sample populations

Two distinct sample sets were used in this study.

**Validation panel (n=48 specimens across 9 sequencing runs).**

- **33 clinical specimens** from patients with reference-laboratory-confirmed diagnoses of the following pathogens (reference method `[TEAM INPUT NEEDED]`: culture-based confirmation, PCR confirmation, or serological confirmation, as applicable for each organism): *Escherichia coli* (6), *Orientia tsutsugamushi* (6), *Burkholderia pseudomallei* (5), *Leptospira* spp (4), *Rickettsia* spp (4), *Streptococcus pneumoniae* (3), *Streptococcus suis* (3), *Coxiella burnetii* (1), and *Yersinia* spp (1).
- **10 positive controls** comprising 4 PC_SINGLE samples (single-organism control material for *E. coli*, *P. aeruginosa*, *S. pneumoniae*, *S. suis*), 5 PC_MIX8 samples (ZymoBIOMICS Microbial Community Standard: *Bacillus subtilis*, *Enterococcus faecalis*, *Escherichia coli*, *Lactobacillus fermentum* [now *Limosilactobacillus fermentum*], *Listeria monocytogenes*, *Pseudomonas aeruginosa*, *Salmonella enterica*, *Staphylococcus aureus*; the standard also includes two yeasts that are not detectable by 16S), and 1 MIXED4 sample (*E. coli* + *P. aeruginosa* + *S. pneumoniae* + *S. suis*).
- **5 negative template controls (NTC)** — molecular-grade water carried through the full extraction and library-preparation process.

Note: `E-coli_S4_L001` is classified as a PC_SINGLE control on the basis of sample provenance; the denominator for clinical sensitivity is therefore 33.

**Study cohort (n=86 patient blood specimens).** Blood samples from patients enrolled with acute febrile illness during whose admission a blood culture bottle was flagged as positive by the automated detection system, but subculture onto aerobic solid media did not recover any organism. **The 16S samples are patient blood drawn during the same admission and are *not* aliquots of the original positive blood culture bottle.** No anaerobic blood culture is performed under local hospital protocol. Of the 86 samples, 71 (82.6%) carry ≥1 positive call from either Centrifuger classification or Minimap2 alignment rescue; 15 (17.4%) returned zero detected taxa, consistent with pre-sequencing failures (DNA extraction, library preparation, or sequencing depth).

### 2.3 Sample collection, DNA extraction, and library preparation

Blood specimens were processed by the NIH wet-lab team. `[TEAM INPUT NEEDED — NIH wet-lab team]` Blood specimen volume, container, anticoagulant if relevant, storage conditions, time-to-extraction, DNA-extraction kit and protocol, library-preparation kit, indexing strategy, and library quantification. The pipeline includes an upstream NCBI SRA Human Scrubber (HRRT) host-DNA depletion step (see §2.4.1). 16S V1-V3 amplification used the `[TEAM INPUT NEEDED — exact primer sequences, typically 27F + 534R]` primer set. Libraries were sequenced on the **Illumina MiSeq** platform (`[TEAM INPUT NEEDED — flow cell type, paired-end read length, indexing kit]`).

### 2.4 Bioinformatics pipeline

The pipeline was implemented in Workflow Description Language v1.0 (WDL) and executed on the Terra.bio cloud platform using the Cromwell workflow engine and Google Cloud Batch for compute orchestration. The full pipeline source code is available at https://github.com/PHemarajata/afi_terra (`wdl/AFI_16S_Main.wdl`, `wdl/AFI_16S_Batch.wdl`; container images `phemarajata614/afi-terra:0.4.1` and `phemarajata614/centrifuger:1.1.0`). The workflow is organised as a two-phase scatter–gather batch analysis supporting parallel processing of multi-run validation and routine clinical sample sets, with task stages for preprocess, classify, align, metrics, interpret, and validate.

#### 2.4.1 Host-DNA depletion and read preprocessing

Paired-end Illumina FASTQ reads were first subjected to human host-read depletion using NCBI sra-human-scrubber (HRRT; `NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl`). Dehosted reads were then quality-trimmed and adapter-clipped using fastp v0.23.4 (Chen et al., 2018; `staphb/fastp:0.23.4`) with default parameters, producing cleaned paired-end FASTQ files for downstream classification and alignment.

#### 2.4.2 Primary taxonomic classification (Centrifuger)

Cleaned reads were classified using Centrifuger v1.1.0 (Song & Langmead, 2024; `phemarajata614/centrifuger:1.1.0`) — a successor of Centrifuge (Kim et al., 2016) — against a custom reference database combining the NCBI bacterial and archaeal genome collection with curated Rickettsiales genomes in a single index (`[TEAM INPUT NEEDED — specific reference DB identifier, build date, source URL]`). The classify task (`wdl/tasks/classify.wdl`, task `RunCentrifuger`) outputs both per-read classification (`*.centrifuger.classification.tsv`) and a kreport-format clade summary (`*.centrifuger.kreport.tsv`) generated via `centrifuger-kreport`. Genus-level read counts were extracted from clade-aggregated rank-G entries using a custom Python parser, with clade reads (not just taxon-level reads) summed to the genus rank so that species-level assignments contributed to the genus call. Classification was executed with 14 threads and 128 GB requested memory (capped at ~104 GB per GCP's 6.5 GB/vCPU policy) on Google Compute Engine N1 custom virtual machines with a 500 GB persistent disk for database localisation.

**A genus was called "Detected" if its assigned read count was ≥500 (`cfr_floor`) AND ≥5× the run-specific NTC maximum read count (`cfr_fold`) for that genus** (`scripts/call_taxa.py`; defaults exposed in `wdl/tasks/interpret.wdl`, task `InterpretCalls`).

#### 2.4.3 Minimap2 alignment-based rescue for Rickettsiales

Order Rickettsiales (*Orientia*, *Rickettsia*, *Anaplasma*, *Ehrlichia*) represents the canonical "cannot-miss" AFI etiology in endemic northeastern Thailand. The V1-V3 region of 16S provides limited sequence variation for genus-level discrimination within Rickettsiales, so an alignment-based rescue stage was added.

In parallel with Centrifuger classification, cleaned reads were aligned to a curated 7-genome Rickettsiales 16S reference panel (`align_rickettsiales_16S/refs_fna/`) using minimap2 v2.x (Li, 2018) with the short-read preset (`-ax sr`). The reference panel comprises *Orientia tsutsugamushi* str. Boryong (AM494475.1) and str. Ikeda (AP008981.1); *Rickettsia prowazekii* str. NMRC Madrid E (CP004888.1), *R. typhi* str. Wilmington (NC_006142.1), and *R. rickettsii* str. 'Sheila Smith' (NC_009882.1); *Anaplasma phagocytophilum* str. HZ (NC_007797.1); and *Ehrlichia chaffeensis* str. Arkansas (NC_007799.1). Alignments were sorted and indexed with samtools (Danecek et al., 2021). For each Rickettsiales genus, per-genus alignment metrics (mapped read count, maximum breadth of coverage = fraction of reference covered ≥1×, and read-count fold over the run-specific NTC alignment-max for that genus) were computed using a custom Python script built on pysam (`wdl/tasks/metrics.wdl`; `phemarajata614/afi-terra:0.4.1`).

A two-tier rescue framework was applied (`scripts/call_taxa.py`):

- **Tier 1 — Genus-level Confirmed (`call = Confirmed`).** Mapped reads ≥100 (`align_confirm_reads`) AND breadth ≥0.25 (`align_confirm_breadth`) AND mapped reads ≥5× the alignment NTC max (`align_fold`). Reported as the specific genus (*Orientia* or *Rickettsia*); the `align_confirmed` flag in `.calls.tsv` is `true`.
- **Tier 2 — Order-level Probable (`call = Probable`).** Mapped reads ≥50 AND breadth ≥0.20 AND mapped reads > the alignment NTC max (no 5× fold requirement). Reported as **"Rickettsiales detected (genus uncertain; recommend confirmatory species-specific qPCR)"** — clinically actionable for empirical doxycycline coverage in an endemic context.
- **Not_Confirmed** (`call = Not_Confirmed`): reads ≥50 but ≤ alignment NTC max. Not reported clinically.
- **Negative** (`call = Negative`): reads <50.

A genus was reported as detected if it satisfied either Centrifuger or alignment-rescue criteria. The `.calls.tsv` output has a 5-level `call` field: `Detected`, `Confirmed`, `Probable`, `Not_Confirmed`, `Negative`. All three positive levels (`Detected`, `Confirmed`, `Probable`) contribute to the reported detection set.

#### 2.4.4 Negative-control handling and per-run background derivation

Two classes of negative controls were processed separately to preserve sensitivity:

1. **No-template controls (NTC)** — water blanks added at the library-preparation step. NTCs contribute to the per-run background reference used for threshold computation.
2. **Negative controls (NC)** — buffer blanks and extraction-condition comparisons. NCs are processed through the full pipeline and reported in run summaries but are **excluded** from background computation. This distinction prevents extraction controls bearing low-level contaminant reads from inflating the noise floor and masking true low-abundance positives.

For each sequencing run, classification and alignment metrics from NTC samples were aggregated by a dedicated `BuildNTCBackground` task (`scripts/build_ntc_background.py`). For each genus, the per-run NCmax is computed as the **maximum** read count observed across all NTC samples assigned to that run, separately for the Centrifuger path (`cfr_ntc_reads`) and the Minimap2 path (`align_ntc_reads`). The interpretation task (`InterpretCalls`) then attaches the appropriate per-genus NTC max to each detection row.

Each `.calls.tsv` row's `ntc_reads` field represents the NTC max value the pipeline used for that organism in that run. We observed inter-run variability in NTC contamination: in run 6_and_7, `NTC2_ExDw_S13_L001` carries 78,691 *Leptospira* reads, 106,204 *Burkholderia* reads (genus level), and 188,778 *Brevundimonas* reads. Study-sample `.calls.tsv` entries in the same run-folder show `ntc_reads = 0` for *Leptospira*, indicating that the pipeline's run-ID grouping treats the contaminated NTC as belonging to a different run-ID pool than the affected study samples. The exact per-run NTC pool composition for each sample is determined by the per-batch Terra input sheet; the cross-check column "Same-run NTC max" in the appendix tables reports the maximum reads of each organism across all NTCs in the folder-grouped run for transparency.

### 2.5 V4 post-pipeline decontamination filter

Genus-level detections passing the primary pipeline thresholds were further processed through a four-tier post-pipeline decontamination filter (`afi_decontamination_filter_v4.py`) informed by landmark low-biomass microbiome contamination reviews (Salter et al., 2014; Glassing et al., 2016; Lauder et al., 2016; de Goffau et al., 2018; Tan et al., 2023) and confirmed against this dataset's NTC profiles. The filter is applied to validation-panel and study samples but is bypassed for organisms documented as positive-control spike-ins in their respective PC samples.

#### 2.5.1 Tier A — High-confidence kit / skin / water contaminants

The following 11 genera are removed on any detection in clinical or study samples: *Pseudomonas*, *Ralstonia*, *Bradyrhizobium*, *Sphingomonas*, *Stenotrophomonas*, *Methylobacterium*, *Acinetobacter*, *Cutibacterium*, *Staphylococcus*, *Corynebacterium*, and *Brevundimonas*. *Brevundimonas* was added to Tier A on the basis of NTC profiles in this dataset (22–285,740 reads per NTC in run 6_and_7); the other ten are documented contaminants in ≥5 of the cited reviews.

#### 2.5.2 Tier A exception — *Burkholderia* species-level safeguard

*Burkholderia* is treated as a Tier-A-equivalent kit contaminant (the genus is dominated in low-biomass samples by *B. cepacia* complex species: *B. contaminans*, *B. cenocepacia*, *B. multivorans*, *B. sola*, *B. cepacia*), but it also contains *B. pseudomallei*, a cannot-miss melioidosis pathogen endemic in Thailand. For every sample with a *Burkholderia* genus detection, the Centrifuger kreport is parsed at the species rank (`S`): reads assigned to *Burkholderia pseudomallei* (NCBI taxonomy ID 28450) are summed and compared against the run's NTC species-level *B. pseudomallei* count. The detection is preserved only if **both** conditions hold: (i) species-level *B. pseudomallei* reads ≥500, and (ii) species reads exceed the run-NTC maximum for *B. pseudomallei*. When preserved, the retained record stores the species-level read count, not the genus total.

The 3 validation-panel *B. pseudomallei* samples (`09502813_S2_L001`, `09-0-02165`, `09700912_S3_L001`) carry 65,016 / 20,061 / 13,744 species-level *B. pseudomallei* reads respectively and trivially pass the safeguard. None of the *Burkholderia* detections in the 86-sample study cohort passed the safeguard.

#### 2.5.3 Tier B — NTC-only organisms

Nine genera identified as present only in NTCs across all runs were removed globally: *Cereibacter*, *Thioclava*, *Bdellovibrio*, *Saltatorellus*, *Pseudogemmobacter*, *Minisyncoccus*, *Rhodoluna*, *Microbacterium*, *Arcanobacterium*. In practice only *Rhodoluna* (2 detections) appeared in clinical samples after upstream NCmax subtraction.

#### 2.5.4 Tier 1 — Ultra-low-abundance environmental noise

Sixteen genera with dataset-wide median per-sample abundance <0.5% were removed: *Shigella*, *Metapseudomonas*, *Stutzerimonas*, *Capsulimonas*, *Chamaesiphon*, *Chloroflexus*, *Flavihumibacter*, *Hymenobacter*, *Limnoglobus*, *Methylovirgula*, *Microvirga*, *Pelagovum*, *Pseudonocardia*, *Rufibacter*, *Salmonella*, *Spirosoma*. *Salmonella* is in Tier 1 because it was not detected at meaningful abundance in this dataset; it is a known AFI pathogen and any future detection should prompt explicit NTC cross-check before attribution. *Mycoplasmopsis* and *Nitrospira* are excluded from Tier 1 and retained as candidate signals (mean abundance 39.78% across 4 samples and 13.45% in 1 sample, respectively).

#### 2.5.5 Tier 2 — Marginal organisms

Twenty-six genera with median 0.5–2.0% abundance are retained only if detected in ≥2 samples AND each detection is ≥1.0% abundance: *Klebsiella*, *Asticcacaulis*, *Xanthomonas*, *Caulobacter*, *Fimbriimonas*, *Gemmatirosa*, *Roseateles*, *Actinomycetospora*, *Algoriphagus*, *Chloroherpeton*, *Delftia*, *Rhodococcus*, *Coxiella*, *Novosphingobium*, *Sphingopyxis*, *Brasilonema*, *Candidatus_Amoebophilus*, *Chroococcidiopsis*, *Dermatobacter*, *Erythrobacter*, *Leptolyngbya*, *Nostoc*, *Pseudoluteimonas*, *Roseomonas*, *Rubrobacter*, *Variovorax*. *Coxiella* is placed in Tier 2 (validated by 1/1 concordance on the single validation-panel sample) — its placement reflects dataset statistics rather than biological exclusion, and *Coxiella* signals should always be interpreted as candidate AFI etiology subject to confirmatory testing.

#### 2.5.6 Positive-control spike-in bypass

For samples typed as `PC_MIX8`, `PC_SINGLE`, `MIXED4`, or generic `PC`, Tier A removal is bypassed for any organism that is a documented spike-in for that PC type. Expected spike-in genera by PC type:

| PC sample | Type | Expected genera |
|---|---|---|
| `E-coli_S4_L001` | PC_SINGLE | *Escherichia* |
| `P-aeru_S5_L001` | PC_SINGLE | *Pseudomonas* |
| `S-pneumo_S2_L001` | PC_SINGLE | *Streptococcus* |
| `S-suis_S3_L001` | PC_SINGLE | *Streptococcus* |
| `PC-20251016_S7_L001`, `PC_S8_L001`, `PC_S10`, `PC_S12`, `PC_S13_L001` | PC_MIX8 | *Bacillus*, *Enterococcus*, *Escherichia*, *Limosilactobacillus*, *Listeria*, *Pseudomonas*, *Salmonella*, *Staphylococcus* (ZymoBIOMICS) |
| `Mixed_S6_L001` | MIXED4 | *Escherichia*, *Pseudomonas*, *Streptococcus* |

The bypass ensures PC samples are evaluated for their ability to detect their expected spike-in organisms rather than penalised by a filter calibrated for clinical contamination. Without this bypass the `P-aeru_S5_L001` PC would fail (because *Pseudomonas* is in Tier A); with the bypass, all 10 PCs are concordant.

### 2.6 Concordance definitions

| Sample category | Concordance rule |
|---|---|
| **Clinical** | Concordant if at least one row matches the expected target organism at genus level AND that row passes the V4 filter (`Final concordance = "Concordant (target detected and retained by V4 filter)"`). Within Rickettsiales, cross-genus rescue (*Orientia* ↔ *Rickettsia*) and order-level rescue (`call = Probable`) both qualify as target match given doxycycline-treatable equivalence at the order level. |
| **Positive control** | Concordant iff all expected spike-in organisms are detected AND retained by V4 (under the PC bypass). |
| **NTC** | Concordant iff no TaqMan Array Card bacterial target genus (*Bartonella*, *Brucella*, *Rickettsia*, *Orientia*, *Yersinia*, *Coxiella*, *Streptococcus*, *Salmonella*, *Escherichia*, *Burkholderia*) is retained after V4. Viral and protozoal TAC targets (Dengue, Chikungunya, Zika, Nipah, HepE, Hantaan, Seoul, JEV, *Plasmodium falciparum* / *vivax*) are not evaluable by this 16S assay. |

### 2.7 Run validity (PC8 check)

Each sequencing run included an 8-organism positive control mock community (PC_MIX8 = ZymoBIOMICS Microbial Community Standard). A run was flagged as analytically valid (`pc8_valid = true`) only if the PC_MIX8 sample successfully detected all expected genera using the criteria above (with the spike-in bypass applied); otherwise, all clinical results from that run were annotated for manual review. Per-run accept / reject decisions are governed by the deployment standard operating procedure (`[TEAM INPUT NEEDED — SOP reference]`).

### 2.8 Statistical analysis

Analytical sensitivity and specificity were computed under a binary contingency framework treating each sample as a single binary outcome (concordant / discordant). The Wilson 95% confidence interval was used for proportions. Inter-run reproducibility was assessed by per-run PC and NTC pass rates across all 9 sequencing runs.

### 2.9 Software versions and compute resources

NCBI sra-human-scrubber (HRRT, latest), fastp v0.23.4, Centrifuger v1.1.0, minimap2 and samtools (bundled in `phemarajata614/afi-terra:0.4.1`), Python 3 with pandas and pysam. The V4 decontamination filter (`afi_decontamination_filter_v4.py`) is a standalone Python 3 script with pandas + openpyxl dependencies. All Docker images are pinned to specific tags and tracked through GitHub-based Method Library imports in Terra. Resource-intensive tasks (Centrifuger classification) were executed on Google Compute Engine N1 custom virtual machines (14 vCPU cores; 128 GB memory request, capped at ~104 GB per GCP policy) with 500 GB pd-standard disk for database localisation. The 14-thread configuration was selected to satisfy GCP's vCPU / memory ratio while avoiding observed instability at higher thread counts.

### 2.10 Data and code availability

The pipeline source code is publicly available at https://github.com/PHemarajata/afi_terra (Dockstore-linked for Terra deployment). All bioinformatics outputs (`.calls.tsv`, Centrifuger kreport files), the V4 decontamination filter script (`afi_decontamination_filter_v4.py`), the appendix generator (`generate_appendices.py`), and the per-sample appendix tables (`APPENDIX-VALIDATION-PANEL.md`, `APPENDIX-STUDY-SAMPLES.md`, `APPENDICES.xlsx`) are available at `[TEAM INPUT NEEDED — final data repository / Google Drive archival location]`. The user guide is published at https://phemarajata.github.io/afi_terra/index.html.

---

## 3. Results

### 3.1 Validation panel performance

#### 3.1.1 Sample composition and overall accuracy

The 48-sample validation panel (33 clinical with known reference-laboratory diagnoses + 10 positive controls + 5 NTCs) was processed through the V4-filter-aware pipeline. Final concordance under the rule of record (target detected AND retained by V4 for clinical; all expected spike-ins detected and retained under PC bypass for PCs; no TAC bacterial target genus retained for NTCs):

| Category | Concordant | Total | Rate | 95% CI (Wilson) |
|---|---|---|---|---|
| Clinical | 21 | 33 | **63.6%** | 46.0%–78.5% |
| Positive controls (PC_MIX8 + PC_SINGLE + MIXED4) | 10 | 10 | **100%** | 72.2%–100% |
| Negative template controls (specificity) | 5 | 5 | **100%** | 56.6%–100% |
| **Sample-level analytical performance (clinical + PC)** | **31** | **43** | **72.1%** | **57.5%–83.6%** |
| **Overall validation accuracy (all categories)** | **36** | **48** | **75.0%** | 61.2%–85.1% |

The 72.1% sample-level analytical performance figure is the headline analytical performance metric for regulatory documentation. Inter-run reproducibility is 100% across all 9 sequencing runs.

#### 3.1.2 Organism-specific clinical performance

| Expected organism | n | Concordant | Sensitivity | Failure mode for discordants |
|---|---|---|---|---|
| *Escherichia coli* | 6 | 6 | 100% | — |
| *Orientia tsutsugamushi* | 6 | 5 | 83.3% | 1 sample (`00618_S7`) called `Probable` via Tier-2 order-level rescue — counted as concordant under cross-genus / order-level rule |
| *Burkholderia pseudomallei* | 5 | 3 | 60% | 2 abundance-driven (other organisms dominate sample) |
| *Leptospira* spp | 4 | 2 | 50% | 2 abundance-driven |
| *Rickettsia* spp | 4 | 2 | 50% | 2 pre-sequencing failures (`00126_S6`, `00369_S1`: zero taxa) |
| *Streptococcus pneumoniae* | 3 | 1 | 33% | 2 abundance-driven; V1-V3 cannot resolve *S. pneumoniae* vs *S. suis* species-level |
| *Streptococcus suis* | 3 | 1 | 33% | 2 abundance-driven; same V1-V3 limitation |
| *Coxiella burnetii* | 1 | 1 | 100% | — |
| *Yersinia* spp | 1 | 0 | 0% | n=1 insufficient; other organisms detected |

Two failure patterns dominate the 12 discordant cases:

1. **Pre-sequencing failures.** `00126_S6_L001` and `00369_S1_L001` (both expected *Rickettsia*) returned zero taxa, consistent with DNA-extraction, library-preparation, or sequencing-depth failure rather than classification error. The pipeline correctly identified target organisms in same-organism samples where adequate read coverage was achieved (e.g., *Rickettsia* detected at 46–100% in `10900410_S10` and `16401070_S11`).
2. **Abundance-driven discordance.** In the remaining discordant cases, the expected target was present but not the dominant 16S signal, with other organisms (often kit-contaminant-class or background environmental genera) accounting for the majority of reads. This reflects 16S genus-detection biology (proportional reporting in mixed samples) rather than classification failure.

#### 3.1.3 Two-tier Rickettsiales rescue performance on the validation panel

The Minimap2 alignment-based rescue is essential for Rickettsiales detection at V1-V3 resolution. Across the 6 expected-*Orientia* validation samples: 4 carry direct genus-level Centrifuger detection (`00389_S5`, `11800801_S7`, `24500367_S5`, `25900911_S4`); Tier 1 Confirmed Minimap2 rescue (`call = Confirmed`) was triggered in `22900253_S4` (4,133 mapped reads, breadth 0.3220, `align_confirmed = true` — this is a *Rickettsia* call in an Orientia-expected sample, counted as concordant under cross-genus Rickettsiales handling); and Tier 2 Probable order-level rescue (`call = Probable`) was triggered in `00618_S7` (68 mapped reads, breadth 0.2015). No validation-panel Rickettsiales sample falls into the failed-rescue category. A threshold uniformity audit (`THRESHOLD-UNIFORMITY-AUDIT.md`) confirms that the rescue thresholds are applied identically across all samples in the dataset.

#### 3.1.4 Decontamination filter behaviour on the validation panel

The V4 filter does not remove any organism from the validation panel that affects target-organism detection. This is informative but **not** an independent validation of the filter: the validation panel is dominated by single high-abundance expected organisms (40–100% abundance) in samples that do not carry Tier A contaminants at detectable levels. The 3 validation-panel *B. pseudomallei* samples trivially pass the species-level safeguard (13,744–65,016 species reads, all >> 500-read threshold). Without the positive-control spike-in bypass, the `P-aeru_S5_L001` PC would have failed (because *Pseudomonas* is in Tier A); with the bypass, all 10 PCs are concordant.

### 3.2 Study cohort results

#### 3.2.1 Cohort size and pre-sequencing failures

The study cohort comprises **86 patient blood specimens** distributed across 5 sequencing runs (`1_and_2`, `3`, `4_and_5`, `6_and_7`, `8_and_9`). All 86 samples had `.calls.tsv` files generated by the pipeline. Of these:

- **71 samples (82.6%) have ≥1 positive call** (`Detected`, `Confirmed`, or `Probable`) from either Centrifuger or Minimap2 rescue.
- **15 samples (17.4%) have zero detected taxa**, consistent with pre-sequencing failures (DNA extraction, library preparation, or sequencing depth).

The pre-sequencing failure rate is comparable to that observed in the validation panel (3/33 = 9% clinical pre-seq failures) and indicates that low-biomass blood specimens have a non-trivial baseline failure rate that should be expected in operational deployment. We recommend a pre-sequencing QC step (Qubit + library qPCR) for future cohorts to distinguish pre-sequencing failures from true biological negatives up front.

All denominators in this Results section use 86 (cohort total) or 71 (samples with any positive call) as appropriate.

#### 3.2.2 Control validity within study runs

Positive and negative controls were processed alongside patient samples across all study-cohort sequencing runs. Two operationally relevant observations from this dataset's NTC profiles:

- **Multiple NTCs in run 6_and_7 carry substantial reagent contamination.** `NTC3_ExEB_S12_L001` (*Brevundimonas* 285,740 reads), `NTC4_NExDw_S16_L001` (*Brevundimonas* 241,384), `NTC5_NExEB_S15_L001` (*Brevundimonas* 204,530), and `NTC2_ExDw_S13_L001` (*Brevundimonas* 188,778; *Burkholderia* 106,204 at genus level, dominated by *B. cepacia* complex species; *Leptospira* 78,691). The pipeline's per-run NCmax derivation appropriately handles this through run-ID-based NTC pool assignment, but the cross-run NTC variability is a reminder that per-run NTC composition matters for fold-change interpretation.
- **Study-run PCs pass under V4 with PC spike-in bypass.** The bypass is essential because *Pseudomonas* and *Staphylococcus* are Tier A contaminants in clinical samples but are expected spike-ins in PC_MIX8 and the relevant PC_SINGLE replicates.

#### 3.2.3 V4 filter impact on the study cohort

The V4 filter was applied to all 71 samples with positive detections (**217 genus-sample detections** in total). Outcomes:

| Filter tier | Detections removed |
|---|---|
| Tier A high-confidence contaminants (Cutibacterium, Staphylococcus, Brevundimonas, Acinetobacter, Corynebacterium, Pseudomonas, Ralstonia, Sphingomonas, Stenotrophomonas, Methylobacterium, Bradyrhizobium) | 61 |
| Burkholderia removed at species-level safeguard (genus signal is *B. cepacia* complex) | 4 |
| Burkholderia preserved as *B. pseudomallei* | **0** |
| NTC-only organisms (Rhodoluna) | 2 |
| Tier 1 ultra-low abundance | 2 |
| Tier 2 marginal + rare | 1 |
| **Total removed** | **70 (32.3% of detections)** |
| **Total retained** | **147** |

A Sankey visualisation of pre- vs post-V4-filter genus distribution across the cohort is provided as **Figure 1** (`figure_a_sankey.html` for interactive review, `figure_a_sankey.png` for print). The figure shows the aggregated read volume flowing from each detected genus (left nodes) to each V4 outcome category (right nodes); the dominant pattern is that most non-contaminant genera flow to KEEP (green), while *Cutibacterium*, *Brevundimonas*, *Staphylococcus*, *Acinetobacter*, *Corynebacterium*, and other Tier A genera flow to REMOVE_TIER_A (red), and the genus-level *Burkholderia* signal flows to REMOVE_BURK (orange — the species-level safeguard found no study-cohort sample with *B. pseudomallei* above background).

Per-sample biomass distribution under V4:

| Pattern | n samples |
|---|---|
| Retain 100% of detected biomass (no contaminants present at detection level) | 21 |
| Retain 0% (all detections removed as contaminants) | 2 |
| Mixed (most informative for clinical review) | 48 |
| No detections at all (pre-sequencing failure) | 15 |
| **Total** | **86** |

#### 3.2.4 Rickettsiales rescue in the study cohort

The single most clinically important finding of the cohort analysis is **Rickettsiales rescue evidence in 11 of 86 samples (12.8%; 11 of 71 samples with any positive call = 15.5%)**, identified by examining all `source = alignment` rows with `call ∈ {Confirmed, Probable, Detected}`:

| Sample | Run | Genus | Mapped reads | Breadth | Rescue tier | Alignment NTC reads | Confidence (sample-vs-NTC) |
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

**Summary.** 1 Tier-1 genus-level *Orientia* call (`16901195_S5_L001`, breadth 0.3235, ~12× NTC headroom — the cleanest Rickettsiales call in the cohort) and 10 Tier-2 order-level "Rickettsiales detected" calls with breadth in the 0.21–0.24 range (just below the Tier-1 0.25 threshold). Across all 11 samples, confidence stratifies as 2 HIGH (1 Tier-1 + 1 Tier-2 with abundance >> NTC), 3 MODERATE, and 6 LOW (alignment NTC carries comparable-magnitude background); the 6 LOW-confidence calls warrant priority follow-up. All 11 samples are flagged for Rickettsiales-specific qPCR confirmation; doxycycline empirical coverage is already standard in this endemic clinical context.

#### 3.2.5 Species-level *B. pseudomallei* findings

Five samples in the study cohort carry a *Burkholderia* genus detection in `.calls.tsv` (`14300185_S8_L001`, `22900187_S10_L001`, `23200430_S6_L001`, `23200519_S12_L001`, `16601093_S7_L001`). Species-level Centrifuger kreport parsing yields:

| Sample | Genus reads | *B. pseudomallei* species reads | *B. cepacia* complex species reads | Run NTC species max | V4 decision |
|---|---|---|---|---|---|
| `23200430_S6_L001` | 54,126 | **8** | 12,789 | 10–51 | **REMOVE** (genus signal is B. cepacia complex contamination) |
| `22900187_S10_L001` | 899 | 0 | 196 | 10–51 | REMOVE |
| `14300185_S8_L001` | 1,689 | 0 | 441 | 0 | REMOVE |
| `23200519_S12_L001` | 876 | 0 | 231 | 0 | REMOVE |
| `16601093_S7_L001` | 673 | 0 | ~150 | 0 | REMOVE |

**No study sample contains *B. pseudomallei* at species level above background.** Species-level *B. pseudomallei* reads across the cohort range 0–8 per sample, below both the 500-read detection threshold and the run NTC's 10–51 reads of *B. pseudomallei*. The genus-level *Burkholderia* signal in these samples is dominated by *B. cepacia* complex (kit contaminant) species. The 3 validation-panel *B. pseudomallei* samples carry 13,744 / 20,061 / 65,016 species-level reads and clearly pass the safeguard, demonstrating that the assay detects *B. pseudomallei* when present at meaningful abundance — but this specific cohort does not contain it. Melioidosis is not a major contributor to this cohort's positive-culture / failed-subculture phenotype as detected by 16S.

#### 3.2.6 Fastidious-organism candidate signals

**Mycoplasmopsis (4 of 86 samples = 4.7%; 4/71 active = 5.6%).** *Mycoplasmopsis* was detected at substantial abundance across the cohort. Read counts in `.calls.tsv` range 537–6,568 reads across the affected samples (post-NCmax subtraction):

| Sample | Run | Reads | % of sample |
|---|---|---|---|
| `13400425_S4_L001` | 4_and_5 | 6,568 | (high) |
| `15700611_S5_L001` | 4_and_5 | 2,568 | 100.0% (only detection in sample) |
| `10300205_S9_L001` | 4_and_5 | 1,622 | 2.17% |
| `09301225_S4_L001` | 4_and_5 | 740 | 45.76% |
| `24500367_S5_L001` | 4_and_5 | 537 | (variable) |

Aggregate (4 samples retained at the rebuilt Tier-1 list in V4): mean abundance 39.78%, maximum 70.55%. *Mycoplasmopsis* species are cell-wall-deficient (lack peptidoglycan), fastidious (require sterol-supplemented media), and slow-growing (1–3 weeks to visible colonies). They would not be expected to grow on routine aerobic blood subculture within standard 5–7 day windows.

**Leptospira (1 sample with NTC caveat).** Sample `09801652_S5_L001` (run 6_and_7) carries 7,294 *Leptospira* reads (22.78% of detected sample biomass). The same run contains `NTC2_ExDw_S13_L001` with 78,691 *Leptospira* reads (10× the study sample). The pipeline reports `ntc_reads = 0` for *Leptospira* in this sample's `.calls.tsv` (the NCmax derivation appears to use a different NTC pool for this sample's run-ID), but the contamination flag is a real caveat. This detection should be treated as a candidate requiring orthogonal confirmation (serology + Leptospira-specific qPCR).

**Brucella (3 of 86 samples, near noise floor).** *Brucella* was detected in 3 samples at sub-1% mean abundance (mean 0.98%, max 1.96%). At this level the signal is at the upper edge of what could plausibly be reagent contamination and the lower edge of what could be low-level bacteremia. Without serological or PCR confirmation these are candidate observations, not diagnoses.

**Streptococcus (5 samples, genus-level only).** Five samples carry *Streptococcus* at mean 12.75% abundance; V1-V3 does not resolve *S. pneumoniae* / *S. suis* / fastidious vs non-fastidious species. Both PC_SINGLE controls for *S. pneumoniae* and *S. suis* hit the same *Streptococcus* genus call in the validation panel, confirming the species ambiguity. Species composition is unresolved by this assay alone.

#### 3.2.7 Other retained organisms

Notable additional post-V4 retained organisms in the cohort (full table in `APPENDIX-STUDY-SAMPLES.md` / `APPENDICES.xlsx`, `Study - sample summary` sheet):

| Genus | Detections | Mean abundance | Max abundance | Notes |
|---|---|---|---|---|
| Thermomicrobium | 11 | 9.13% | 55.79% | Environmental thermophile; unlikely human pathogen. High prevalence is unusual — possible assay artefact / reagent issue. **Should not be included in clinical interpretation tables.** |
| Escherichia | 6 | 8.28% | 12.23% | Known AFI pathogen |
| Streptococcus | 5 | 12.75% | 23.29% | See species-ambiguity note above |
| Faucicola | 5 | 6.85% | 10.30% | Oral microbiome genus; possible sampling / handling contamination |
| **Mycoplasmopsis** | **4** | **39.78%** | **70.55%** | **Fastidious, cell-wall-deficient — see §3.2.6** |
| Klebsiella | 4 | 20.99% | 73.90% | Known AFI pathogen |
| Methylorubrum | 4 | 9.71% | 17.53% | Environmental |
| Kocuria | 4 | 7.76% | 16.89% | Skin / environmental |
| Xanthomonas | 4 | 11.85% | 24.04% | Plant-associated; possible environmental |
| Enterobacter | 3 | 3.39% | 5.14% | Known AFI pathogen |
| Brucella | 3 | 0.98% | 1.96% | Near contamination floor; see caveat in §3.2.6 |
| Streptomyces | 3 | 9.21% | 14.18% | Environmental |
| Comamonas | 3 | 9.01% | 11.59% | Environmental |
| Leptospira | 1 | 22.78% | 22.78% | NTC caveat — see §3.2.6 |
| Porphyromonas | 1 | 8.45% | 8.45% | Anaerobic genus retained |
| Desulfovibrio | 1 | 0.78% | 0.78% | Anaerobic genus retained |

The retention of *Porphyromonas* and *Desulfovibrio* — even at low frequency — undermines any absolute "no anaerobes detected" framing. V1-V3 primer bias against certain Gram-positive anaerobes remains a known limitation; "absent in this assay" is not equivalent to "biologically absent."

### 3.3 Inter-run reproducibility and pipeline QC

Inter-run reproducibility was assessed by tracking PC and NTC pass rates across all 9 sequencing runs (4 validation-panel runs + 5 study-cohort runs):

- **PC pass rate:** 100% (all PC_MIX8, PC_SINGLE, and MIXED4 samples detect their expected spike-in organisms under V4 with PC bypass).
- **NTC specificity rate:** 100% on the rule of record (no TAC bacterial target genus is retained after V4 in any validation-panel NTC). Some study-run NTCs carry substantial pre-NCmax reagent contamination (notably run 6_and_7); the pipeline's per-run NTC max derivation accounts for this through run-ID-based pool assignment. Per-run accept / reject decisions are governed by the deployment SOP.

---

## 4. Discussion

### 4.1 What the analysis shows

**Rickettsiales involvement is the most prevalent identifiable signal in this cohort** — present in 11 of 86 samples (12.8%) overall, equivalent to ~16% of the 71 samples with any positive detection. The Minimap2 alignment-based rescue identifies 1 high-confidence Tier-1 *Orientia* call (`16901195_S5_L001`, breadth 0.3235, ~12× NTC headroom) and 10 Tier-2 order-level "Rickettsiales detected" rescues. Rickettsiales are obligate intracellular pathogens that cannot be cultured on routine blood agar; their detection in a positive-blood-culture / no-subculture-growth cohort is the most biologically coherent finding in the dataset. This pattern aligns with the expected endemic epidemiology of northeastern Thailand, where scrub typhus (*Orientia tsutsugamushi*) and spotted-fever-group rickettsioses account for a substantial fraction of AFI presentations.

**No study sample contains *Burkholderia pseudomallei* above species-level background.** Species-rank parsing of the Centrifuger kreport shows *B. pseudomallei* species reads of 0–8 per study sample, below both the 500-read detection threshold and the run NTC's *B. pseudomallei* count (10–51 reads). The genus-level *Burkholderia* signal in these samples is dominated by *B. cepacia* complex species, which are kit / water contaminants (Salter et al., 2014). The 3 *B. pseudomallei*-positive validation-panel samples (with 13,744–65,016 species reads) demonstrate the assay does detect melioidosis when present at meaningful abundance — but this specific cohort does not contain it.

**Mycoplasmopsis is the most prominent fastidious-organism-class signal in the cohort** — 4 of 86 samples (4.7%) at mean abundance 39.78% and max 70.55%. *Mycoplasmopsis* (*Mycoplasma* family) organisms are by definition cell-wall-deficient and fastidious; they require sterol-supplemented media and 1–3 weeks of incubation. Their presence in patient blood is biologically consistent with a positive-bottle / no-subculture-growth phenotype. The cohort-level prevalence (~5%) is too small to support cohort-level causal claims; this is a hypothesis-generating observation that warrants species-level kreport interrogation and Mycoplasma-specific PCR / specialised culture follow-up.

### 4.2 What the filter design choices contribute

The V4 decontamination filter materially changes which organisms are reported and which are dropped. Three filter design choices have outsized impact in this dataset:

- **Aggressive Tier A.** Removing the 11 high-confidence kit / skin / water contaminants is responsible for 61 of 70 V4 removals. *Cutibacterium*, *Staphylococcus*, and *Brevundimonas* are the top three contributors. The *Brevundimonas* addition (V4 vs V3) is supported by NTC profiles showing it at 22–285,740 reads per NTC in run 6_and_7.
- ***Burkholderia* species-level safeguard.** The genus is in Tier A but the species-level kreport parser preserves *B. pseudomallei* on a per-row basis. Without this safeguard, melioidosis cases would be erroneously discarded; with the safeguard, only true *B. pseudomallei* signals are preserved (validated against same-run NTC species reads). In this cohort the safeguard correctly admitted no false-positive melioidosis call.
- **Positive-control spike-in bypass.** Without this bypass, the `P-aeru_S5_L001` PC would fail because *Pseudomonas* is in Tier A, even though *Pseudomonas aeruginosa* is the expected positive-control spike-in. The bypass restores PC accuracy from 9/10 to 10/10 and brings the sample-level analytical performance from 30/43 = 69.8% to 31/43 = 72.1%, matching the legacy APHL figure.

Filter membership of *Coxiella* (Tier 2) and *Salmonella* (Tier 1) reflects observed dataset statistics, not biological exclusion: future detections of either organism above the relevant tier thresholds should be interpreted as candidate AFI etiology, not contamination.

### 4.3 What the data do not support

The following claims are not supported by the data:

- **"Strong support for fastidious-organism hypothesis."** With single-sample *Leptospira* (n=1, NTC caveat) and 3-sample low-abundance *Brucella* (mean 0.98%), the cohort-level evidence is thin. *Mycoplasmopsis* is the strongest signal class but still ~5% of cohort. Honest framing: hypothesis-generating observations for confirmatory testing in a minority of cases.
- **"No anaerobes detected."** *Porphyromonas* and *Desulfovibrio* are retained in the V4 output, contradicting absolute negation. The V1-V3 primer set is also known to under-detect some Gram-positive anaerobes; "absent in this assay" is not equivalent to "biologically absent."
- **"16S identified organisms that caused failed culture."** The 16S samples are patient blood, not aliquots of the original positive bottle. The detected organisms are candidates for explaining the bottle signal but cannot be causally linked without paired bottle / blood analysis. 16S also does not distinguish viable from non-viable cells.

### 4.4 Reframed hypothesis evaluation

The seven hypotheses originally proposed (low-level pathogen, fastidious, anaerobes, slow-growing, viable-but-non-culturable [VBNC], L-forms, polymicrobial competition) cannot be evaluated as binary cohort-level claims from this dataset. A more honest summary:

| Hypothesis | What this dataset supports | What it does not |
|---|---|---|
| **Rickettsiales involvement** | 11 of 86 samples (12.8%; 11/71 active = 15.5%) with alignment-rescue evidence; 1 Tier-1 high-confidence *Orientia* call; pattern aligns with expected endemic epidemiology | Species-level resolution within Rickettsiales for the Tier-2 calls; clinical-outcome correlation (held by EPI team) |
| **Fastidious organisms** | *Mycoplasmopsis* in 4 samples (4.7% of cohort) at 39.78% mean is suggestive of fastidious-organism involvement; *Leptospira* (n=1) and *Brucella* (n=3, sub-1%) are weaker candidate signals | Cohort-level claim that fastidious organisms explain >10% of cases; species-level confirmation; viability of detected organisms |
| **Slow-growing organisms** | Same organisms (Mycoplasma-class, *Leptospira*, *Brucella*) are also slow-growing; mechanism is coherent | Direct demonstration that the bottle signal was generated by slow-growing organisms specifically |
| **Low-level pathogen / inadequate inoculum** | Culturable organisms (*E. coli*, *Klebsiella*, *Enterobacter*) detected at moderate abundance in patient blood; abundance in blood ≠ abundance in bottle | Direct correlation requires paired bottle / blood analysis |
| **Anaerobes** | Cannot be definitively assessed: V1-V3 has primer bias against some anaerobes; a few anaerobic-genus signals (*Porphyromonas*, *Desulfovibrio*) are present | "Absent from cohort" claim is not supportable |
| **VBNC / L-forms** | Mechanism plausible (Mycoplasma-class are naturally cell-wall-deficient); no direct evidence in 16S data | Requires viability assays, microscopy of bottle broth |
| **Polymicrobial competition** | Mean 2.6 detections per sample post-filter; mechanism plausible | Requires paired bottle / blood analysis |

The fairest summary is that **Rickettsiales involvement is the dominant identifiable signal (~13% of cohort), and fastidious / slow-growing / cell-wall-deficient organism involvement is plausible in a smaller minority (~5–10%), with Mycoplasma-class detection being the most prominent organism class.** Cohort-level "strong support" for any single non-Rickettsiales hypothesis is not warranted by the data.

### 4.5 Comparison with expected AFI pathogens in northeastern Thailand

In northeastern Thailand, the major causes of acute febrile illness include rickettsial diseases (*Orientia*, *Rickettsia*), leptospirosis, melioidosis (*B. pseudomallei*), brucellosis, dengue and other arboviruses, enteric bacterial pathogens, and tuberculosis. Our validation panel demonstrates that the pipeline detects these organisms when they are present at adequate abundance:

- *B. pseudomallei*: detected at species level (13,744–65,016 species reads) in 3 of 5 validation-panel cases with expected melioidosis
- *O. tsutsugamushi*: 6/6 detection across expected cases (Centrifuger genus + Minimap2 order-level rescue combined)
- *Leptospira*: detected at 58–65% abundance in 2 of 4 validation-panel positives

The study cohort itself yielded:

- **Rickettsiales detected by Minimap2 alignment rescue in 11 of 86 samples (12.8%; 11/71 active = 15.5%).** One genus-level Tier-1 *Orientia* call (HIGH); ten Tier-2 order-level rescues stratified as 1 HIGH, 3 MODERATE, 6 LOW given alignment NTC headroom.
- **No detectable *B. pseudomallei* at species level.**
- **One candidate *Leptospira* with same-run NTC caveat.**
- **Three low-abundance *Brucella* candidates.**
- **Four *Mycoplasmopsis* detections at substantial abundance.**

This pattern — one high-confidence *Orientia* detection, ~10 order-level Rickettsiales rescues, no melioidosis at species level, candidate *Leptospira* / *Brucella* / Mycoplasma in additional samples — is consistent with a cohort in which a meaningful subset of patients have Rickettsiales involvement that culture cannot recover (Rickettsiales are intracellular obligate parasites that do not grow on routine blood culture), alongside a smaller subset showing fastidious / slow-growing organism candidate signals. The cohort is biologically consistent with expected AFI etiology in this region.

### 4.6 Methodological lessons for downstream deployment

Operational deployment of this workflow should incorporate four refinements identified during this analysis:

1. **Document the NCmax derivation rule explicitly.** The current pipeline appears to use the per-run-ID NTC pool when computing NCmax, which means a contaminated NTC in one folder-grouped run may not influence the threshold for study samples in a different run-ID even within the same sequencing batch (consistent with study-sample `ntc_reads = 0` for *Leptospira* despite an NTC carrying 78,691 reads of *Leptospira* in the same folder). This behaviour and the rule (per-run-ID max across all NTCs in that pool) should be documented in the SOP and the run-grouping logic made transparent in deployment QC.
2. **Add a pre-sequencing QC step.** A non-trivial fraction (15/86 = 17.4%) of study specimens returned zero taxa, indicating DNA-extraction / library-preparation / sequencing-depth failure. Qubit + library qPCR before sequencing would distinguish these from true biological negatives.
3. **Maintain the V4 filter's species-level safeguards.** *Burkholderia* requires species-level parsing; analogous treatment for *Brucella* (species discrimination of *B. melitensis* / *B. abortus* / *B. suis* / *B. canis*) and *Leptospira* (pathogenic vs saprophytic species) would strengthen interpretation.
4. **Document the positive-control spike-in bypass.** The bypass is a deliberate design choice and must be documented so future operators understand why Tier A genera are sometimes retained in PC samples.

### 4.7 Clinical implications for the AFI study cohort

For the 11 samples with Rickettsiales rescue evidence, doxycycline-based empirical therapy in an endemic context is standard clinical practice and is supported by the alignment-rescue evidence; species-specific qPCR (*Orientia tsutsugamushi*-specific 56-kDa TSA gene; *Rickettsia* 17-kDa antigen gene targets) is recommended for definitive species identification. For the *Mycoplasmopsis* candidate samples, Mycoplasma-specific PCR and Mycoplasma broth subculture would convert the genus-level signal into a species-resolved confirmation. For the *Leptospira* candidate (with NTC caveat) and *Brucella* candidates (near noise floor), paired-serum serology (IgM / IgG acute and convalescent) and species-specific qPCR are recommended before any clinical interpretation.

The 15 samples with no detection at all warrant a separate workup — these are pre-sequencing failures and should be subject to clinical review and, where feasible, repeat specimen collection or alternative diagnostic modalities.

---

## 5. Limitations

1. **Sample source mismatch.** The 16S samples are patient blood, not aliquots of the original positive blood culture bottle. Detection of an organism by 16S in patient blood is consistent with its causing the bottle signal but does not establish causation. Paired bottle / blood analysis (16S of the bottle broth) would be required to make causal claims.
2. **No viability assessment.** 16S detects DNA from viable, non-viable, and viable-but-non-culturable (VBNC) cells. A high-abundance 16S signal does not establish that the organism is alive and would be culturable under appropriate conditions.
3. **V1-V3 primer biases.** The 27F primer set has documented under-detection of certain Gram-positive anaerobes, *Mycobacteria*, and some *Bifidobacterium* / *Atopobium* / *Gardnerella* lineages. Negative findings for these classes are limitations of the assay, not biological absence.
4. **Genus-level resolution for most organisms.** Only *Burkholderia* receives species-level parsing in the V4 filter. *S. pneumoniae* vs *S. suis*, fastidious vs non-fastidious *Streptococcus*, species-level *Brucella* / *Leptospira* / *Mycoplasmopsis* discrimination, and species-level Rickettsiales discrimination for the Tier-2 rescue calls are not resolved by this assay alone.
5. **NTC contamination in some runs.** `NTC2_ExDw_S13_L001` in run 6_and_7 carries substantial *Leptospira*, *Burkholderia*, and *Brevundimonas* reads. The pipeline's NCmax derivation appears to handle this through per-run-ID NTC pool assignment, but the rule should be documented and contaminated NTCs flagged in deployment QC.
6. **Filter tier membership is dataset-informed.** The Tier 1 and Tier 2 lists were assembled from the observed abundance distributions of this dataset combined with the cited contamination reviews. Generalisability to other specimen types or kit batches has not been validated.
7. **Cohort sample size.** 86 study samples (71 with detections) is appropriate for a pilot but insufficient to establish prevalence estimates for any single etiology. Findings of n=1 (*Leptospira*) or n=3 at sub-1% abundance (*Brucella*) are best treated as candidate detections requiring follow-up, not cohort-level prevalence claims.
8. **No clinical-outcome correlation.** Patient demographics, treatment, response, serology, and clinical outcomes are held by the EPI team and are not linked to the 16S findings in the bioinformatics package; the integration will be performed in the clinical / epidemiology sections of the final manuscript.
9. **No formal anaerobic culture for comparison.** The hospital protocol does not include anaerobic blood culture, so direct comparison between 16S-detected anaerobic-genus signals and a culture gold standard is not possible.
10. **B. pseudomallei detection limit.** The validation-panel *B. pseudomallei* cases have species-level reads in the 13,744–65,016 range; the practical species-level detection floor (where the assay reliably distinguishes *B. pseudomallei* from *B. cepacia* complex contamination) is set at 500 reads but has not been formally established by dilution-series testing.

---

## 6. Recommendations for follow-up

### 6.1 Immediate confirmatory testing

1. **Rickettsiales-rescue samples (n=11).** Rickettsiales-specific qPCR — *O. tsutsugamushi* 56-kDa TSA gene, *Rickettsia* 17-kDa antigen and *gltA* targets — on archived blood / serum from each of the 11 samples in §3.2.4. The 1 Tier-1 *Orientia* call (`16901195_S5_L001`) is the highest-priority confirmation target; the 6 LOW-confidence Tier-2 calls are next priority given their comparable alignment NTC background. Acute / convalescent IgM / IgG serology where archived serum permits.
2. **Mycoplasmopsis-positive samples (n=4).** Mycoplasma / Mycoplasmataceae-targeted PCR on archived blood; species-level kreport interrogation; consider Mycoplasma broth (sterol-supplemented) culture if specimen material remains.
3. **Leptospira candidate (`09801652_S5_L001`).** Paired acute / convalescent serology and Leptospira-specific qPCR. Investigate why `NTC2_ExDw_S13_L001` (same folder-grouped run) carried 78,691 *Leptospira* reads — reagent batch tracing is recommended.
4. **Brucella candidates (n=3).** Serology and *Brucella*-specific qPCR.

### 6.2 Pipeline / methods improvements

5. **Document NCmax derivation explicitly** (see §4.6).
6. **Add species-level reporting for clinically critical genera.** *Burkholderia* has a species-level safeguard; analogous treatment for *Brucella* (species discrimination of *B. melitensis* / *B. abortus* / *B. suis* / *B. canis*) and *Leptospira* (pathogenic vs saprophytic species) would strengthen interpretation.
7. **Routine cross-run NC vs NC profiling** with thresholds for flagging contaminated NTCs as part of deployment QC.

### 6.3 Study-design improvements for future cohorts

8. **Paired bottle and blood 16S** — sequence the original positive bottle broth, not just patient blood, so direct comparison is possible.
9. **Specialised-culture attempts** — Mycoplasma broth (sterol-enriched), Fletcher / EMJH media (*Leptospira*), *Brucella* agar (*Brucella*) on aliquots, to convert "candidate detection" into "confirmed isolation."
10. **Clinical metadata integration** — treatment history, serology, outcomes, exposure history (rural / agricultural exposure, animal contact, recent travel) for each case, to interpret which candidate organisms are biologically plausible.

---

## 7. Conclusion

A 16S V1-V3 amplicon workflow with Centrifuger v1.1.0 primary classification, Minimap2 alignment-based rescue for Rickettsiales, a two-tier reporting framework (genus-level Tier 1 Confirmed + order-level Tier 2 Probable), and a species-aware V4 decontamination filter (with positive-control spike-in bypass) achieves 72.1% sample-level analytical performance against a 43-sample reference panel (95% CI 57.5%–83.6%), with 100% NTC specificity (95% CI 56.6%–100%) and 100% inter-run reproducibility across 9 sequencing runs.

Applied to a cohort of 86 AFI cases (71 with ≥1 positive call, 15 pre-sequencing failures) with positive automated blood culture and failed subculture recovery, the workflow identifies Rickettsiales rescue evidence in 11 samples (~13% of cohort, ~16% of samples with any detection) — a finding consistent with the expected endemic epidemiology of northeastern Thailand — and candidate fastidious-organism signals (*Mycoplasmopsis*, *Leptospira*, *Brucella*) in an additional small minority of the cohort. **No study sample carries *Burkholderia pseudomallei* above the species-level detection threshold.**

The 16S workflow described here is suitable for use as a complementary diagnostic for AFI cases with positive culture signal but failed subculture recovery, with the understanding that candidate detections require orthogonal confirmation (species-specific PCR / qPCR, serology, specialised culture) before clinical interpretation. The most actionable finding for this specific cohort is the Rickettsiales rescue evidence in approximately one in eight AFI cases.

---

## 8. References

**Decontamination and low-biomass microbiome studies:**

1. Salter SJ, Cox MJ, Turek EM, Calus ST, Cookson WO, Moffatt MF, Turner P, Parkhill J, Loman NJ, Walker AW. Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. *BMC Biol.* 2014;12:87. doi:10.1186/s12915-014-0087-z. PMID: 25387460.
2. Glassing A, Dowd SE, Galandiuk S, Davis B, Chiodini RJ. Inherent bacterial DNA contamination of extraction and sequencing reagents may affect interpretation of microbiota in low bacterial biomass samples. *Gut Pathog.* 2016;8:24. doi:10.1186/s13099-016-0103-7. PMID: 27239228.
3. Lauder AP, Roche AM, Sherrill-Mix S, Bailey A, Laughlin AL, Bittinger K, Leite R, Elovitz MA, Parry S, Bushman FD. Comparison of placenta samples with contamination controls does not provide evidence for a distinct placenta microbiota. *Microbiome.* 2016;4(1):29. doi:10.1186/s40168-016-0172-3. PMID: 27338728.
4. de Goffau MC, Lager S, Salter SJ, Bonney EA, Bertozzi-Villa A, Wagner J, Charnock-Jones DS, Smith GCS, Parkhill J. Recognizing the reagent microbiome. *Nat Microbiol.* 2018;3(8):851–853. doi:10.1038/s41564-018-0202-y. PMID: 30046175.
5. Tan CCS, Ko KKK, Chen H, Liu J, Loh M, Chia M, Nagarajan N, SG10K_Health Consortium. No evidence for a common blood microbiome based on a population study of 9,770 healthy humans. *Nat Microbiol.* 2023;8(5):973–985. doi:10.1038/s41564-023-01350-w. PMID: 36997797.

**Bioinformatic tools used:**

6. Song L, Langmead B. Centrifuger: lossless compression of microbial genomes for efficient and accurate metagenomic sequence classification. *Genome Biol.* 2024;25(1):106. doi:10.1186/s13059-024-03244-4. PMID: 38654369.
7. Kim D, Song L, Breitwieser FP, Salzberg SL. Centrifuge: rapid and sensitive classification of metagenomic sequences. *Genome Res.* 2016;26(12):1721–1729. doi:10.1101/gr.210641.116. PMID: 27852649.
8. Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics.* 2018;34(18):3094–3100. doi:10.1093/bioinformatics/bty191. PMID: 29750242.
9. Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. *Gigascience.* 2021;10(2):giab008. doi:10.1093/gigascience/giab008. PMID: 33590861.
10. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics.* 2018;34(17):i884–i890. doi:10.1093/bioinformatics/bty560. PMID: 30423086.
11. NCBI Human Read Removal Tool (SRA Human Scrubber, HRRT). National Center for Biotechnology Information. https://github.com/ncbi/sra-human-scrubber.

**Platform / specification references (citations pending team confirmation):**

12. Voss K, Van der Auwera G, Gentry J. Full-stack genomics pipelining with GATK4 + WDL + Cromwell. *F1000Research* (or current OpenWDL specification, https://openwdl.org). `[TEAM INPUT NEEDED — confirm Terra / Cromwell citation]`
13. ZymoBIOMICS Microbial Community Standard product citation. Zymo Research Corp. `[TEAM INPUT NEEDED — exact catalogue number, lot, and citation form]`

**AFI epidemiology / clinical references (`[TEAM INPUT NEEDED — EPI team]`):** primary references for *O. tsutsugamushi* / scrub typhus burden in northeastern Thailand; *B. pseudomallei* / melioidosis epidemiology; *Leptospira* / leptospirosis burden; *Brucella* / brucellosis prevalence; *Rickettsia* spotted-fever-group epidemiology; CDC TaqMan Array Card AFI Multi-Pathogen Card validation reports.

---

## 9. Appendices

The following companion files accompany this manuscript and contain the per-sample data tables and raw outputs from the analysis:

- **`APPENDIX-VALIDATION-PANEL.md`** — per-detection table for the 48-sample validation panel: sample, run, category, expected organism, total reads, biomass kept / removed by V4 filter, detected genus, source (Centrifuger / Minimap2), reads, % of sample, pipeline NTC reads, Minimap2 rescue info, rescue tier, V4 filter outcome, per-row concordance, failure mode, Final concordance (TAC + V4), notes. Includes a per-sample concordance summary table at the bottom.
- **`APPENDIX-STUDY-SAMPLES.md`** — per-detection table for the 86-sample study cohort: same columns as the validation appendix plus same-run NTC max (cross-check), Burkholderia species evidence (where applicable), confidence label, hypothesis class, and recommended follow-up.
- **`APPENDICES.xlsx`** — same content in Excel format with 4 sheets (`Validation - detections`, `Validation - sample summary`, `Study - detections`, `Study - sample summary`).
- **`DECONTAMINATION-FILTER-REPORT-V4.txt`** — raw output of the V4 filter with summary statistics and *Burkholderia* species-level evidence per sample.
- **`afi_decontamination_filter_v4.py`** — executable V4 filter implementation.
- **`generate_appendices.py`** — appendix generator script (regenerates all of the above from raw `.calls.tsv` and Centrifuger kreport files).
- **`figure_a_sankey.html`** / **`figure_a_sankey.png`** — Figure 1, V4 filter Sankey visualisation.
- **`MANUSCRIPT-WALKTHROUGH-THAI.md`** — Thai-language walkthrough for the wet-lab team.

---

## 10. Placeholders requiring team input

### Resolved from the codebase (no longer placeholder)

| Section | Resolved from |
|---|---|
| §2.4 Pipeline orchestration | WDL via Cromwell / Terra, source at https://github.com/PHemarajata/afi_terra; container images `phemarajata614/afi-terra:0.4.1` + `phemarajata614/centrifuger:1.1.0`. |
| §2.4.1 Host-DNA removal + read preprocessing | NCBI SRA Human Scrubber + fastp v0.23.4. |
| §2.4.2 Primary taxonomic classification | Centrifuger v1.1.0 (Song & Langmead, 2024). Single index covers bacteria / archaea + Rickettsiales. |
| §2.4.2 Detection thresholds | Centrifuger: ≥500 reads AND ≥5× run NTC max (`cfr_floor=500`, `cfr_fold=5.0`). |
| §2.4.3 Minimap2 reference panel | 7 references (AM494475.1, AP008981.1, CP004888.1, NC_006142.1, NC_009882.1, NC_007797.1, NC_007799.1). |
| §2.4.3 Rescue thresholds | Tier 1 Confirmed: reads ≥100, breadth ≥0.25, ≥5× align NTC max. Tier 2 Probable: reads ≥50, breadth ≥0.20, reads > align NTC max. |
| §2.4.4 NTC background derivation | Per-run **maximum** read count across NTCs in that run's pool, separately for Centrifuger and Minimap2 paths. |
| §2.3 Sequencer | Illumina MiSeq (platform confirmed by PI). |
| §8 Software citations | Centrifuger, Centrifuge, Minimap2, SAMtools, fastp, NCBI Scrubber — added. |

### Items still requiring team input

| Section | Content needed | Owner |
|---|---|---|
| Title page | Authors + affiliations; corresponding author; running title; CoI; funding; author contributions (CRediT) | NIH wet-lab team + CDC programme office for author list / funding |
| §2.1 Study design | IRB approval, ethics body, consent procedure, enrolment dates, hospital(s) / district(s), inclusion / exclusion criteria | EPI team |
| §2.2 Sample populations | Reference-laboratory diagnostic methods used for each validation-panel pathogen (PCR? culture? serology?) | NIH wet-lab team + EPI team |
| §2.3 Sample collection / wet-lab | Blood specimen volume, container, time-to-extraction, DNA extraction kit + version, 16S V1-V3 primer sequences, library prep kit, MiSeq flow cell + run configuration | NIH wet-lab team |
| §2.4.2 Centrifuger reference DB | Specific reference database identifier (name, build date, source URL or accession set) used for the Centrifuger index | PI / bioinformatics |
| §2.7 / §3.2.2 / §3.3 | Deployment SOP reference for per-run QC pass criteria | PI / bioinformatics |
| §2.10 Data availability | Final data repository / archival location for `.calls.tsv` and kreport outputs; raw FASTQ SRA accession if deposited | NIH wet-lab team |
| §4.7 / §6.1 Clinical follow-up | Confirmatory test results (if any) and clinical outcomes for the 11 Rickettsiales-rescue samples, 4 *Mycoplasmopsis* samples, 1 *Leptospira* candidate, 3 *Brucella* candidates | EPI team |
| §8 References | AFI epidemiology references for NE Thailand; Rickettsiales / melioidosis / leptospirosis / brucellosis clinical references; CDC TAC AFI panel validation reports; Terra/Cromwell citation; ZymoBIOMICS catalogue citation | EPI team + NIH wet-lab team |

---

*End of MANUSCRIPT-PACKAGE-FINAL. Bioinformatics package prepared 2026-05-25.*
