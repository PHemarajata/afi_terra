# Methods — AFI Terra 16S Metagenomic Pipeline

*Draft methods section for manuscript inclusion. Last revised 2026-05-13 to incorporate the post-pipeline V4 decontamination filter and the corrected analytical validation summary. See notes at bottom for citations and style adjustments.*

---

## Bioinformatic Pipeline

The Acute Febrile Illness (AFI) 16S metagenomic pipeline was implemented in Workflow Description Language v1.0 (WDL) and executed on the Terra.bio cloud platform using the Cromwell workflow engine and Google Cloud Batch for compute orchestration. The full pipeline source code is available at https://github.com/PHemarajata/afi_terra. The workflow is organized as a two-phase scatter–gather batch analysis supporting parallel processing of multi-run validation and routine clinical sample sets.

## Sample Preprocessing

Paired-end Illumina FASTQ reads were first subjected to human host read depletion using NCBI sra-human-scrubber (HRRT). Dehosted reads were then quality-trimmed and adapter-clipped using fastp v0.23.4 (`staphb/fastp:0.23.4`) with default parameters, producing cleaned paired-end FASTQ files for downstream classification and alignment.

## Taxonomic Classification

Cleaned reads were classified using Centrifuger v1.1.0 (`phemarajata614/centrifuger:1.1.0`) against a custom reference database combining the NCBI bacterial and archaeal genome collection with curated Rickettsiales genomes, ensuring comprehensive taxonomic coverage in a single classifier (replacing an earlier dual-database Kraken2 design). Classification outputs were converted to kreport format using `centrifuger-kreport`, and genus-level read counts were extracted from clade-aggregated rank-G entries using a custom Python parser. Classification was executed with 14 threads and 128 GB allocated memory on Google Compute Engine N1 custom virtual machines, with 500 GB persistent disk for database localization.

## Confirmatory 16S Alignment

In parallel with classification, cleaned reads were aligned to a curated Rickettsiales 16S rRNA reference panel (7 reference genomes: *Orientia tsutsugamushi* str. Boryong [AM494475.1] and str. Ikeda [AP008981.1]; *Rickettsia prowazekii* str. NMRC Madrid E [CP004888.1], *R. typhi* str. Wilmington [NC_006142.1], and *R. rickettsii* str. 'Sheila Smith' [NC_009882.1]; *Anaplasma phagocytophilum* str. HZ [NC_007797.1]; *Ehrlichia chaffeensis* str. Arkansas [NC_007799.1]) using minimap2 with the short-read preset (`-ax sr`). Alignments were sorted and indexed with samtools, and per-genus alignment metrics — mapped read count and maximum breadth of coverage — were computed using a custom Python script built on pysam (`phemarajata614/afi-terra:0.4.1`).

## Negative Control Handling and Background Computation

Two classes of negative controls were processed separately to preserve sensitivity:

1. **No-template controls (NTC):** water blanks added at the library preparation step. NTCs contribute to the per-run background reference used for threshold computation.
2. **Negative controls (NC):** buffer blanks and extraction-condition comparisons. NCs are processed through the full pipeline and reported in run summaries but are **excluded** from background computation. This distinction prevents extraction controls bearing low-level contaminant reads from inflating the noise floor and masking true low-abundance positives.

For each sequencing run, classification and alignment metrics from NTC samples were aggregated by a dedicated `BuildNTCBackground` task to establish per-genus background read counts. For each genus, the per-run NCmax (NTC maximum) is computed as the maximum read count observed across all NTC samples assigned to that run, separately for the Centrifuger-classified path (`cfr_ntc_reads`) and the Minimap2-aligned path (`align_ntc_reads`).

## Interpretation Logic

Each sample's classification and alignment results were evaluated against the NTC background using a two-module decision tree with a two-tier Rickettsiales rescue framework:

- **Module 1 — Centrifuger-based broad detection.** A genus was called positive (`Detected`) if classified reads exceeded an absolute floor (≥500 reads, `cfr_floor`) AND surpassed the NTC background by ≥5-fold (`cfr_fold = 5`).
- **Module 2 — Alignment-based confirmation and rescue for *Orientia* and *Rickettsia*.**
  - **Tier 1 — Confirmed (genus-level call).** ≥100 mapped reads (`align_confirm_reads`) AND breadth ≥0.25 (`align_confirm_breadth`) AND mapped reads ≥5× alignment NCmax (`align_fold`). Reported as the specific genus (*Orientia* or *Rickettsia*).
  - **Tier 2 — Probable (order-level "Rickettsiales detected" rescue).** ≥50 mapped reads AND breadth ≥0.20 AND mapped reads > alignment NCmax (no 5× fold requirement). Reported as "Rickettsiales detected (genus uncertain; recommend confirmatory species-specific qPCR)" — clinically actionable for empirical doxycycline coverage in endemic regions where V1–V3 sequence variation does not allow confident genus-level discrimination of Rickettsiales.
  - **Not_Confirmed (equivocal).** ≥50 reads but ≤ alignment NCmax. Not reported clinically.
  - **Negative.** <50 reads.

A genus was reported as detected if it satisfied either module's criteria.

## V4 Post-Pipeline Decontamination Filter

To address cross-study reagent / skin / water contamination beyond per-run NTC subtraction, pipeline outputs were further processed through a four-tier post-pipeline decontamination filter (`afi_decontamination_filter_v4.py`). Tier membership is informed by landmark low-biomass microbiome contamination reviews (Salter 2014, Glassing 2016, Lauder 2016, de Goffau 2018, Tan 2023) and confirmed against this dataset's NTC profiles.

- **Tier A — High-confidence kit / skin / water contaminants.** Eleven genera are removed on any detection in clinical or study samples: *Pseudomonas*, *Ralstonia*, *Bradyrhizobium*, *Sphingomonas*, *Stenotrophomonas*, *Methylobacterium*, *Acinetobacter*, *Cutibacterium*, *Staphylococcus*, *Corynebacterium*, *Brevundimonas*.
- **Tier A exception — *Burkholderia* species-level safeguard.** For every *Burkholderia* genus detection, the Centrifuger kreport is parsed at the species rank for *Burkholderia pseudomallei* (NCBI taxonomy ID 28450). The detection is preserved only if (i) species-level reads ≥ 500 AND (ii) species reads exceed the run-NTC maximum for *B. pseudomallei* at species rank. When preserved, the retained record stores the species-level read count rather than the genus total.
- **Tier B — NTC-only organisms.** Nine genera observed only in NTC samples (*Cereibacter*, *Thioclava*, *Bdellovibrio*, *Saltatorellus*, *Pseudogemmobacter*, *Minisyncoccus*, *Rhodoluna*, *Microbacterium*, *Arcanobacterium*) are removed globally.
- **Tier 1 — Ultra-low-abundance environmental noise.** Sixteen genera with dataset-wide median per-sample abundance < 0.5% are removed.
- **Tier 2 — Marginal organisms.** Twenty-seven genera with median 0.5–2.0% abundance are retained only if detected in ≥ 2 samples AND each detection is ≥ 1.0% abundance.
- **Positive-control spike-in bypass.** For positive-control samples, Tier A removal is bypassed for any organism that is a documented spike-in for that sample's PC type (PC_SINGLE single-organism controls, PC_MIX8 = ZymoBIOMICS Microbial Community Standard [*Bacillus*, *Enterococcus*, *Escherichia*, *Limosilactobacillus*, *Listeria*, *Pseudomonas*, *Salmonella*, *Staphylococcus*], MIXED4 = *E. coli* + *P. aeruginosa* + *S. pneumoniae* + *S. suis*). PC concordance is then defined as all expected spike-in organisms being both detected and retained.

## Run Validity (PC8 Check)

Each sequencing run included an 8-organism positive control mock community (PC_MIX8 = ZymoBIOMICS Microbial Community Standard). A run was flagged as analytically valid (`pc8_valid = true`) only if the PC_MIX8 sample successfully detected all expected genera using the criteria above (with the spike-in bypass applied); otherwise, all clinical results from that run were annotated for manual review.

## Sample Set Composition and Modes

The pipeline accepts a sample table with the following types: NTC, NC, PC_MIX8 (8-organism positive control), PC_SINGLE (single-organism positive), MIXED4 (4-organism mock), PC (generic positive control), and clinical. Each sample is annotated with one of two operating modes: **validation** (with expected-taxa annotation for sensitivity/specificity computation) or **routine** (for clinical reporting). The batch workflow supports samples from multiple sequencing runs in a single submission, with per-run grouping for background computation and run-validity assessment.

## Concordance Definitions

Sample-level concordance was evaluated under the following rules:

- **Clinical samples** were concordant if the expected target organism was detected at the genus level AND retained by the V4 filter. Within Rickettsiales, both order-level rescue (`call = Probable`) and cross-genus rescue (e.g., a *Rickettsia* rescue call in an *Orientia*-expected sample) qualified as target matches given the doxycycline-treatable equivalence at the order level.
- **Positive controls** were concordant if all expected spike-in organisms for that PC type were detected AND retained under the V4 filter with PC bypass.
- **Negative template controls** were concordant if no TaqMan Array Card (TAC) bacterial target genus (*Bartonella*, *Brucella*, *Rickettsia*, *Orientia*, *Yersinia*, *Coxiella*, *Streptococcus*, *Salmonella*, *Escherichia*, *Burkholderia*) was retained after V4. Viral and protozoal TAC targets are not evaluable by this 16S assay.

## Statistical Analysis

Analytical sensitivity and specificity were computed under a binary contingency framework treating each sample as a single binary outcome (concordant / discordant). Wilson 95% confidence intervals were used for proportions. Inter-run reproducibility was assessed by per-run PC and NTC pass rates across all sequencing runs.

## Compute Resource Allocation

Resource-intensive tasks (Centrifuger classification) were executed on Google Compute Engine N1 custom virtual machines (14 vCPU cores; 128 GB memory request, automatically capped at ~104 GB per GCP's 6.5 GB/vCPU policy) with 500 GB pd-standard disk for database localization. The 14-thread configuration was selected to satisfy GCP's vCPU/memory ratio while avoiding observed instability at higher thread counts.

## Software Versions

NCBI sra-human-scrubber (HRRT, latest), fastp v0.23.4, Centrifuger v1.1.0, minimap2 and samtools (bundled in `phemarajata614/afi-terra:0.4.1`), Python 3 with pandas and pysam. The V4 decontamination filter (`afi_decontamination_filter_v4.py`) is a standalone Python 3 script with pandas + openpyxl dependencies. All Docker images are pinned to specific tags and tracked through GitHub-based Method Library imports in Terra.

---

## Notes for Manuscript Adaptation

- **Add citations** for: NCBI HRRT, fastp (Chen 2018, *Bioinformatics* 34:i884), Centrifuger (Song & Langmead 2024, *Genome Biol* 25:106), Centrifuge (Kim et al. 2016, *Genome Res* 26:1721 — for the methodological predecessor), minimap2 (Li 2018, *Bioinformatics* 34:3094), samtools (Li et al. 2009 / Danecek 2021), Terra/Cromwell, and the contamination references for the V4 filter rationale (Salter 2014 *BMC Biol* 12:87 PMID 25387460; Glassing 2016 *Gut Pathog* 8:24 PMID 27239228; Lauder 2016 *Microbiome* 4:29 PMID 27338728; de Goffau 2018 *Nat Microbiol* 3:851 PMID 30046175; Tan 2023 *Nat Microbiol* 8:973 PMID 36997797).
- **Adjust the 14-thread justification** if you prefer to omit the GCP-specific operational detail (useful for reproducibility but may read as too implementation-specific for some journals).
- **Tense:** Currently written in past tense as appropriate for a Methods section. Adjust if your target journal prefers present tense.
- **Length:** ~1,250 words; tighten or expand any section as needed. The V4 filter and concordance-definition sections are new in this revision (2026-05-13).

## See Also

- `MANUSCRIPT-FINAL-DRAFT.md` — long-form manuscript draft (~7,800 words) incorporating these methods
- `MANUSCRIPT-CONDENSED-DRAFT.md` — JCM-style condensed manuscript (~3,900 words)
- `MANUSCRIPT-WALKTHROUGH-THAI.md` — Thai-language walkthrough for the wet-lab team
- `APPENDIX-VALIDATION-PANEL.md` / `APPENDIX-STUDY-SAMPLES.md` / `APPENDICES.xlsx` — per-sample detection tables
- Pipeline source: https://github.com/PHemarajata/afi_terra
- User guide: https://phemarajata.github.io/afi_terra/index.html
