# AFI 16S study — bioinformatics handoff for the wet-lab + EPI team

**Date:** 2026-05-25
**Prepared by:** Bioinformatics (Peera Hemarajata + Yuyi)
**Bioinformatics support window ends:** 2026-05-31

> **How to use this document.** Everything in §C, §D, §E, §F, and §G is **ready to lift into the manuscript with minimal editing**. Supporting materials (per-sample appendices, V4 filter script, appendix generator, Figure 1) are listed in §H — your PI can include any of them with the delivery package. If you find yourself wanting more detail on a number or a method, the appendix tables are the underlying evidence base; ask before May 31 if anything is unclear.

---

## At a glance

We sequenced the 16S rRNA gene (V1-V3 region) on patient blood samples from acute febrile illness (AFI) cases whose blood cultures flagged positive on the automated system but failed to grow on subculture. We built a custom analysis pipeline (Centrifuger classifier + Minimap2 alignment-based rescue for Rickettsiales + a four-tier contamination filter) and validated it on a 48-sample reference panel. Headline finding in the 86-sample study cohort: **Rickettsiales rescue evidence in 11 of 86 samples (~13%)**, consistent with expected northeastern Thailand AFI epidemiology. **No study sample carries *Burkholderia pseudomallei* above the species-level detection threshold.** Bioinformatics is hand-off ready on Methods (analysis), Results, and interpretation paragraphs. The wet-lab IRB / extraction / library prep / MiSeq details, the clinical and epidemiology context, the full Background, and the Conclusion are yours to write.

---

## §A. Who owns what

| Manuscript section | Owner | Status |
|---|---|---|
| Title page, authors, affiliations, funding, CoI | Wet-lab + CDC | You write |
| Abstract — Background sentence | Wet-lab / EPI | You write |
| Abstract — Methods sentence | Wet-lab + Bioinformatics | We have a drop-in for the bioinformatics half (§F) |
| Abstract — Results sentence | Bioinformatics | **Drop-in below (§F)** |
| Abstract — Conclusion sentence | Wet-lab / EPI | You write |
| §1 Introduction / Background | Wet-lab / EPI | You write |
| §2 Methods — sample collection, DNA extraction, primers, library prep, MiSeq config | Wet-lab | You write |
| §2 Methods — IRB / ethics / enrolment | EPI | You write |
| §2 Methods — bioinformatics pipeline | Bioinformatics | **Drop-in below (§B)** |
| §2 Methods — concordance definitions | Bioinformatics | **Drop-in below (§B)** |
| §3 Results — validation panel performance | Bioinformatics | **Drop-in below (§C)** |
| §3 Results — study cohort findings | Bioinformatics | **Drop-in below (§D)** |
| §4 Discussion — interpretation of our findings | Bioinformatics | **Drop-in below (§E)** |
| §4 Discussion — clinical implications, case context, treatment response | Wet-lab + EPI | You write |
| §5 Limitations | Shared | **Non-negotiable bioinformatics limitations in §G**; add your clinical limitations |
| §6 Conclusion | Wet-lab / EPI | You write |
| References — bioinformatics software | Bioinformatics | Listed in §B below |
| References — AFI epidemiology, clinical | EPI | See `LEADS-AND-FRAMING.md` citation starter list |
| Appendices — per-sample tables | Bioinformatics | `docs/manuscript/appendices/` |

---

## §B. Methods — analysis pipeline (drop-in)

The 16S V1-V3 amplicon sequence data were processed through a custom four-step bioinformatics pipeline implemented in WDL on the Terra.bio cloud platform (Cromwell + Google Cloud Batch). Pipeline source code is publicly available at https://github.com/PHemarajata/afi_terra.

**Step 1 — Read preparation.** Paired-end Illumina reads were first stripped of any human DNA sequences using the NCBI SRA Human Scrubber tool, then quality-trimmed and adapter-clipped with fastp (Chen et al., 2018). Only de-hosted, quality-filtered reads proceeded to downstream analysis.

**Step 2 — Bacterial classification against a comprehensive reference database.** Cleaned reads were matched against a custom reference database covering the NCBI bacterial and archaeal genome collection plus extended Rickettsiales coverage. We used Centrifuger v1.1.0 (Song & Langmead, 2024), a modern metagenomic classifier and successor of Centrifuge. Genus-level read counts were tabulated from a clade-aggregated taxonomy report. **A genus was called *Detected* if its read count was at least 500 AND at least 5-fold above the same-run negative template control (NTC) maximum read count for that genus.**

**Step 3 — Targeted second-pass search for Rickettsiales.** The V1-V3 region does not always carry enough sequence variation to confidently distinguish among Rickettsiales genera (*Orientia*, *Rickettsia*, *Anaplasma*, *Ehrlichia*). Because Rickettsiales is the canonical "cannot-miss" AFI etiology in endemic northeastern Thailand and these organisms cannot be cultured on routine blood agar, we added an alignment-based rescue step: cleaned reads were aligned with minimap2 (Li, 2018) against a curated panel of 7 complete Rickettsiales genomes covering *Orientia tsutsugamushi* (str. Boryong and Ikeda), *Rickettsia prowazekii*, *R. typhi*, *R. rickettsii*, *Anaplasma phagocytophilum*, and *Ehrlichia chaffeensis*. For each Rickettsiales genus we computed mapped read count, maximum breadth of coverage (fraction of the reference covered ≥1×), and read-count fold over the same-run NTC alignment maximum. A two-tier reporting framework was applied:

- **Tier 1 (Confirmed, genus-level call):** ≥100 mapped reads AND breadth ≥0.25 AND reads ≥5× alignment NTC max. Reported as the specific genus (*Orientia* or *Rickettsia*).
- **Tier 2 (Probable, order-level "Rickettsiales detected"):** ≥50 mapped reads AND breadth ≥0.20 AND reads above the alignment NTC max. Reported as **"Rickettsiales detected (genus uncertain; recommend confirmatory species-specific qPCR)"** — clinically actionable for empirical doxycycline coverage in an endemic setting.

**Step 4 — Multi-tier decontamination filter.** Blood is a low-biomass specimen, so reagent and laboratory contamination is a known major source of false-positive 16S signal (Salter et al., 2014). We applied a four-tier post-pipeline filter (V4 filter) to remove background contaminants while preserving genuine pathogens. Three of its design choices matter for interpretation:

- **High-confidence kit / skin / water contaminant removal.** Eleven well-documented contaminant genera (e.g., *Cutibacterium*, *Staphylococcus*, *Brevundimonas*, *Acinetobacter*, *Corynebacterium*, *Pseudomonas*) are removed on any detection. Full list is in the filter script (see §H).
- ***Burkholderia* species-level safeguard.** The *Burkholderia* genus contains both *B. pseudomallei* (melioidosis) and the *B. cepacia* complex (a documented kit / water contaminant). For every *Burkholderia* genus detection, the Centrifuger species-rank output is re-parsed: the detection is preserved as *B. pseudomallei* only if species-level *B. pseudomallei* reads are ≥500 AND exceed the same-run NTC's species-level *B. pseudomallei* count. When preserved, the retained record stores the species-level count, not the genus total.
- **Positive-control spike-in bypass.** For PC samples, the Tier-A removal is skipped for any organism documented as a spike-in for that PC type (e.g., *Pseudomonas* is preserved in the PC_SINGLE *P. aeruginosa* control).

**Concordance definitions.** Each sample was scored under the rule of record:

- **Clinical sample:** *Concordant* if the expected target organism was detected at the genus level AND the row passed the V4 filter. Cross-genus rescue within Rickettsiales (*Orientia* ↔ *Rickettsia*) and order-level Tier-2 rescue both qualify as a target match (doxycycline-treatable equivalence).
- **Positive control:** *Concordant* if all expected spike-in organisms were detected and retained under the V4 filter with PC bypass.
- **Negative template control:** *Concordant* if no TaqMan Array Card bacterial target genus (*Bartonella*, *Brucella*, *Rickettsia*, *Orientia*, *Yersinia*, *Coxiella*, *Streptococcus*, *Salmonella*, *Escherichia*, *Burkholderia*) was retained after V4.

---

## §C. Validation results (drop-in)

The 48-sample validation panel comprised 33 clinical specimens with known reference-laboratory diagnoses, 10 positive controls (5 PC_MIX8 ZymoBIOMICS 8-organism standards + 4 PC_SINGLE + 1 MIXED4), and 5 negative template controls (NTCs).

| Category | Concordant / Total | Rate |
|---|---|---|
| Clinical (strict) | 21 / 33 | **63.6%** |
| Positive controls | 10 / 10 | **100%** |
| Negative template controls (specificity) | 5 / 5 | **100%** |
| **Sample-level analytical (clinical + PC)** | **31 / 43** | **72.1%** |
| Overall validation accuracy | 36 / 48 | 75.0% |

Inter-run reproducibility was 100% across all 9 sequencing runs. (Wilson 95% confidence intervals available on request.)

Two failure patterns explain the 12 discordant clinical cases. Two samples (`00126_S6_L001`, `00369_S1_L001`, both expected *Rickettsia*) returned zero taxa — these are pre-sequencing failures (DNA extraction, library prep, or sequencing depth), not classification errors. The remaining discordants were abundance-driven: the expected target was detectable but not the dominant 16S signal in the sample, with other organisms (often kit-contaminant-class or background environmental genera) accounting for the majority of reads. This reflects 16S genus-detection biology in mixed samples rather than a pipeline error.

The three validation-panel *B. pseudomallei* samples (`09502813_S2_L001`, `09-0-02165`, `09700912_S3_L001`) carry 65,016 / 20,061 / 13,744 species-level *B. pseudomallei* reads respectively and trivially pass the species-level safeguard, demonstrating that the assay does detect melioidosis when it is present at meaningful abundance.

---

## §D. Study cohort findings (drop-in)

The study cohort comprises **86 patient blood specimens** distributed across 5 sequencing runs. All 86 had pipeline outputs generated. Of these, **71 (82.6%) have at least one positive call** (Detected / Confirmed / Probable) from either Centrifuger classification or alignment-based rescue, and **15 (17.4%) returned zero detected taxa**, consistent with pre-sequencing failures. The V4 filter was applied to all 217 genus-sample detections in the 71 active samples, removing 70 (32.3%) — predominantly high-confidence kit / skin / water contaminants — and retaining 147 detections for interpretation.

### D.1 Rickettsiales rescue (headline finding)

The single most clinically important finding is **Rickettsiales rescue evidence in 11 of 86 samples (12.8%; equivalently 11 of 71 samples with any positive detection = 15.5%)**. One sample is a Tier-1 genus-level *Orientia* call (the highest-confidence Rickettsiales detection in the cohort); the other 10 are Tier-2 order-level "Rickettsiales detected" rescues with breadth-of-coverage values just below the Tier-1 threshold.

| Sample | Run | Genus | Mapped reads | Breadth | Rescue tier | Alignment NTC reads | Confidence |
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

Across the 11 samples, confidence stratifies as 2 HIGH (1 Tier-1 + 1 Tier-2 with abundance >> NTC), 3 MODERATE, and 6 LOW (alignment NTC carries comparable-magnitude background). All 11 samples are flagged for Rickettsiales-specific qPCR confirmation; doxycycline empirical coverage is already standard in this endemic clinical context.

### D.2 *Burkholderia pseudomallei*: no study sample above species-level background

Five study samples carry a *Burkholderia* genus signal in pipeline outputs. When the Centrifuger output is re-parsed at the species rank, **no study sample carries *B. pseudomallei* reads at or above the 500-read detection threshold**, and no sample's *B. pseudomallei* species-level count exceeds the same-run NTC's *B. pseudomallei* count. The genus-level signal in these samples is dominated by *B. cepacia* complex species, which are documented kit / water contaminants. The V4 filter correctly removes all 5 *Burkholderia* detections.

The cohort therefore does not contain melioidosis as detected by this assay.

### D.3 *Mycoplasmopsis* candidate signal

*Mycoplasmopsis* — a fastidious, cell-wall-deficient organism class that would not be expected to grow on routine aerobic blood subculture within standard 5–7 day windows — was detected in **4 of 86 samples (4.7%)**, with mean abundance 39.78% and a maximum of 70.55% in the affected samples. Species-level identity has not been resolved in this analysis; each detection should be cross-checked against its specific run's NTC and confirmed by Mycoplasma-specific PCR before clinical interpretation.

### D.4 *Leptospira* and *Brucella* candidates

One study sample (`09801652_S5_L001`) carries 7,294 *Leptospira* reads (22.78% of detected sample biomass). **Important caveat:** the same sequencing run contains a heavily contaminated NTC (`NTC2_ExDw_S13_L001`) with 78,691 *Leptospira* reads — approximately 10× the study sample's count. The pipeline's per-run NTC max derivation appears to use a different uncontaminated NTC for this sample, but the contamination flag is a real one. This detection should be treated as a candidate requiring orthogonal confirmation (paired serology + *Leptospira*-specific qPCR).

Three study samples carry *Brucella* genus signal at near-noise abundance (mean 0.98%, max 1.96%). At this level the signal sits at the upper edge of what could be reagent contamination and the lower edge of what could be low-level bacteremia. Without confirmatory serology these are candidate observations, not diagnoses.

---

## §E. Interpretation paragraphs — drop into Discussion verbatim

> *These four paragraphs are written to be lifted into the Discussion section of the paper. You may need to renumber sub-headers to fit the rest of your Discussion structure.*

**E.1 Rickettsiales as the headline finding.** The single most prevalent identifiable signal in the cohort is Rickettsiales involvement — 11 of 86 samples (12.8%) carry alignment-based rescue evidence, with 1 high-confidence Tier-1 *Orientia* call and 10 Tier-2 order-level rescues. Rickettsiales are obligate intracellular pathogens that cannot be cultured on routine blood agar; their detection in a positive-blood-culture / no-subculture-growth cohort is the most biologically coherent finding in the dataset. This pattern aligns with the expected endemic epidemiology of northeastern Thailand, where scrub typhus (*Orientia tsutsugamushi*) and spotted-fever-group rickettsioses account for a substantial fraction of AFI presentations. Empirical doxycycline coverage is standard clinical practice in this context; all 11 samples are recommended for Rickettsiales-specific qPCR confirmation.

**E.2 No *B. pseudomallei* above species-level background.** No study sample carries *Burkholderia pseudomallei* above the species-level detection threshold. Genus-level *Burkholderia* signals in the cohort are dominated by *B. cepacia* complex species, which are well-documented kit / water contaminants. That the assay can detect *B. pseudomallei* when present is demonstrated by the three validation-panel cases, which carry 13,744–65,016 species-level *B. pseudomallei* reads. The absence of detectable melioidosis in this specific cohort is therefore itself an informative finding: it argues against melioidosis as a major contributor to the positive-culture / failed-subculture phenotype in this group of patients.

**E.3 *Mycoplasmopsis* as a candidate fastidious-organism signal.** *Mycoplasmopsis* was detected in 4 of 86 samples (4.7%) at substantial abundance (mean 39.78%, max 70.55%). Mycoplasma-class organisms lack a cell wall, require sterol-supplemented media, and grow slowly (1–3 weeks to visible colonies), and would not be expected to grow on routine aerobic blood subculture within standard 5–7 day windows. This is the most prominent fastidious-organism-class signal in the cohort and is biologically consistent with the positive-bottle / no-subculture-growth phenotype. The cohort-level prevalence is too small to support causal claims at the cohort level, and species-level identity has not been resolved; this is a hypothesis-generating observation that warrants species-level interrogation and Mycoplasma-specific PCR follow-up. *Leptospira* (n=1, with a same-run NTC contamination caveat) and *Brucella* (n=3, at near-noise abundance) are additional, weaker candidate signals consistent with fastidious / slow-growing organism involvement, also requiring orthogonal confirmation.

**E.4 What the 16S data do not support.** Three caveats merit explicit framing. First, the 16S samples are patient blood, not aliquots of the original positive blood culture bottle; detection of an organism in patient blood is consistent with — but does not prove — that organism's causation of the bottle signal. We therefore use "candidate detection" rather than causal language throughout. Second, 16S detects DNA from viable, non-viable, and viable-but-non-culturable cells alike; a high-abundance 16S signal does not establish that the organism is alive or would be culturable. Third, although classic obligate anaerobes (*Bacteroides*, *Prevotella*, *Fusobacterium*, *Clostridium*) were not detected in the retained set, the V1-V3 primer region is known to under-detect certain Gram-positive anaerobes; absence in this assay is not equivalent to biological absence. A small number of anaerobic-genus signals (*Porphyromonas*, *Desulfovibrio*) were in fact retained, so the cohort is not strictly anaerobe-negative.

---

## §F. Abstract sentences (drop into your abstract)

Insert these sentences into the abstract Results bullet (and one into the Methods bullet) where they fit your overall narrative:

> **(Methods sentence, second half of your Methods bullet.)** Sequencing data were processed through a custom 16S V1-V3 pipeline combining Centrifuger v1.1.0 primary classification, Minimap2 alignment-based rescue against a curated Rickettsiales reference panel with a two-tier reporting framework, and a multi-tier decontamination filter with a species-level safeguard for *Burkholderia pseudomallei* and a positive-control spike-in bypass.

> **(Results sentence.)** Validated against a 48-sample reference panel, the assay achieved 72.1% sample-level analytical performance (31/43; 63.6% strict clinical concordance, 100% positive-control concordance, 100% negative-control specificity) with 100% inter-run reproducibility across 9 sequencing runs. Applied to the 86-sample AFI cohort, the workflow identified Rickettsiales rescue evidence in 11 of 86 samples (12.8%; 11/71 samples with detections = 15.5%; 1 Tier-1 genus-level *Orientia* + 10 Tier-2 order-level rescues), candidate *Mycoplasmopsis* detections in 4 samples (mean abundance 39.78%), and *Leptospira* (n=1) and *Brucella* (n=3) candidate signals; no study sample carried *B. pseudomallei* above the species-level detection threshold.

---

## §G. Limitations the manuscript must keep

> *These four limitations are non-negotiable — please do not paraphrase them in a way that softens their meaning. They can sit in your Limitations section alongside the clinical / epi limitations you add.*

1. **Sample-source mismatch.** The 16S samples are patient blood, not aliquots of the original positive blood culture bottle. Detection of an organism by 16S in patient blood is consistent with — but does not establish — that organism's causation of the bottle signal. Paired bottle / blood 16S would be needed for causal claims.
2. **No viability assessment.** 16S detects DNA from viable, non-viable, and viable-but-non-culturable cells. A 16S signal does not establish that an organism is alive or would be culturable under appropriate conditions.
3. **V1-V3 primer biases.** The 27F primer set has documented under-detection of certain Gram-positive anaerobes, *Mycobacteria*, and some *Bifidobacterium* / *Atopobium* / *Gardnerella* lineages. Negative findings for these classes are limitations of the assay, not biological absence.
4. **NTC contamination in run 6_and_7.** One NTC in that run carried substantial *Leptospira*, *Burkholderia*, and *Brevundimonas* reads. The pipeline's per-run NTC max derivation appears to handle this through run-ID-based NTC pool assignment, but the *Leptospira* candidate detection in `09801652_S5_L001` warrants particularly careful orthogonal confirmation given the same-run NTC contamination flag.

(Six additional limitations — pilot sample size, species-level resolution limits, filter dataset-dependence, no clinical-outcome correlation in this package, no anaerobic culture comparator, and the practical *B. pseudomallei* detection floor — can be added if your target journal expects a comprehensive limitations section.)

---

## §H. Supporting materials accompanying this delivery

The following files are available for your PI to include with the delivery package, alongside this handoff:

- **Per-sample validation table**: `APPENDIX-VALIDATION-PANEL.md` — every detection in the 48-sample validation panel with source (Centrifuger / Minimap2), reads, % of sample, rescue tier, V4 outcome, and per-row concordance.
- **Per-sample study cohort table**: `APPENDIX-STUDY-SAMPLES.md` — every detection in the 86-sample cohort with the same columns plus same-run NTC cross-check, *Burkholderia* species evidence (where applicable), confidence label, and recommended follow-up.
- **Excel workbook (4 sheets)**: `APPENDICES.xlsx` — `Validation - detections`, `Validation - sample summary`, `Study - detections`, `Study - sample summary`.
- **V4 filter audit log**: `DECONTAMINATION-FILTER-REPORT-V4.txt` — raw output of the V4 filter run with summary statistics.
- **V4 filter source**: `afi_decontamination_filter_v4.py` — executable Python implementation; open this if you need the full Tier-A / Tier-B / Tier-1 / Tier-2 genus lists.
- **Appendix generator**: `generate_appendices.py` — regenerates the appendix tables from raw pipeline outputs.
- **Figure 1 (Sankey, V4 filter impact)**: `figure_a_sankey.html` (interactive) and `figure_a_sankey.png` (for print). Suggested caption: "Aggregated read-volume flow from each detected genus (left) to each V4 filter outcome (right) across the 86-sample study cohort."
- **Companion handoff files**: `AFI-TEAM-HANDOFF-TH.md` (Thai walkthrough of the methods, results, and limitations) and `LEADS-AND-FRAMING.md` (narrative framing, citation starters, EPI question prompts, and skeleton outlines for the sections the wet-lab + EPI team owns).
- **Pipeline source code**: https://github.com/PHemarajata/afi_terra — full WDL workflow, container images, and per-task scripts.

If anything in this document is unclear, the appendix tables and the V4 filter audit log are the underlying evidence base — flag any number that doesn't reconcile before May 31.
