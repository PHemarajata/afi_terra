# A 16S rRNA V1-V3 Workflow with Two-Tier Rickettsiales Rescue and Species-Aware Decontamination for Acute Febrile Illness Pathogen Surveillance

## Title page

**Authors:** `[PLACEHOLDER - clinical team]`
**Affiliations:** `[PLACEHOLDER - clinical team]`
**Corresponding author:** `[PLACEHOLDER - clinical team]`
**Running title:** `[PLACEHOLDER - ≤50 chars]`
**Keywords:** acute febrile illness; 16S rRNA; Rickettsiales; blood culture; melioidosis; decontamination
**Conflict of interest / Funding:** `[PLACEHOLDER - clinical team]`
**Data availability:** All sequencing outputs, V4 filter script, and per-sample appendices are available at `[PLACEHOLDER - data repository]`.

**Target word count:** ~5,000 words body. Companion long-form manuscript: `MANUSCRIPT-FINAL-DRAFT.md`.

---

## Abstract (~250 words)

**Background.** Acute febrile illness (AFI) in endemic northeastern Thailand has a broad bacterial differential — including obligate intracellular Rickettsiales, *Burkholderia pseudomallei*, *Leptospira*, *Brucella*, and fastidious or slow-growing organisms — many of which routine blood culture cannot recover. Positive automated blood culture signal with failed subculture growth is a recognized diagnostic gap.

**Methods.** We developed a 16S V1-V3 amplicon workflow combining Centrifuge primary classification with Minimap2 alignment-based rescue against a curated Rickettsiales 16S reference panel, a two-tier reporting framework (genus-level vs. order-level "Rickettsiales detected"), and a multi-tiered decontamination filter (V4) implementing a species-level safeguard for *B. pseudomallei* (parsing Centrifuge kreport at the species rank, comparing against same-run NTC species reads, with a 500-read threshold) and a positive-control spike-in bypass. The workflow was validated against a 48-sample panel (33 clinical, 10 positive controls, 5 NTC) and applied to 86 AFI cases with positive blood culture but failed subculture.

**Results.** Sensitivity was 31/43 (72.1%, 95% CI 57.3–83.3%; 63.6% [21/33] clinical-only and 100% [10/10] positive-control subgroups); specificity 5/5 (100%); positive predictive value 31/31 (100%); inter-run reproducibility 100% across 9 runs. In the AFI cohort, 71/86 samples carried ≥1 detection; 15 had zero taxa (pre-sequencing failures). Rickettsiales rescue evidence was identified in 11 samples (1 Tier-1 genus-level *Orientia*, 10 Tier-2 order-level "Rickettsiales detected"). *Mycoplasmopsis* — a fastidious cell-wall-deficient organism class — was retained in 4 samples (mean 39.78%, max 70.55%). *Leptospira* at 22.78% in one sample and *Brucella* at sub-1% abundance in three samples emerged as candidate detections. No study sample carried *B. pseudomallei* above species-level background.

**Conclusion.** The workflow detects clinically relevant AFI bacterial pathogens with good analytical performance and Rickettsiales-optimized rescue. Approximately 13% of AFI cases with positive culture / no growth show Rickettsiales rescue evidence; an additional minority show fastidious-organism candidate signals. Candidate detections require orthogonal confirmation before clinical reporting.

---

## 1. Introduction (~300 words)

Acute febrile illness in tropical Southeast Asia presents a wide etiologic differential that includes *Orientia tsutsugamushi* (scrub typhus), *Rickettsia* spp, *Burkholderia pseudomallei* (melioidosis), *Leptospira*, *Brucella*, and fastidious or slow-growing organisms. These pathogens have substantially different culture requirements: some grow readily on routine blood agar, some require specialized media or extended incubation, and Rickettsiales — obligate intracellular — cannot be cultured by routine methods at all. Blood cultures that are flagged positive by automated detection but fail subculture recovery (positive signal / no growth) are a common diagnostic problem, with clinical management proceeding on empirical therapy alone.

16S rRNA gene sequencing has been proposed as a complementary diagnostic for these scenarios because it detects bacterial DNA regardless of culturability. The V1-V3 region (27F primer set) is the most widely used clinical V1-V3 target but has known under-detection biases for some Gram-positive anaerobes and limited genus-level discrimination within Rickettsiales. In low-biomass clinical samples, reagent and laboratory contamination (1–4) further complicates interpretation; a recent 9,770-sample healthy-human blood study found no consistent core blood microbiome (5), strengthening the case for aggressive contaminant filtering before clinical interpretation.

We developed and validated a 16S V1-V3 amplicon workflow for AFI bacterial surveillance addressing these constraints through (i) Centrifuge primary classification for breadth, (ii) Minimap2 alignment-based rescue against a curated Rickettsiales reference panel, (iii) a two-tier reporting framework (genus-level + order-level "Rickettsiales detected"), and (iv) a species-aware decontamination filter (V4) with a *B. pseudomallei* species-level safeguard and a positive-control spike-in bypass. Here we present the analytical validation against a 43-sample reference panel plus 5 NTCs and the application to a cohort of 86 AFI cases with positive blood culture / failed subculture.

---

## 2. Methods (~1,000 words)

### 2.1 Study design and samples

`[PLACEHOLDER - clinical team]`: IRB approval, ethics body, enrollment period, geographic catchment, inclusion/exclusion criteria.

Two specimen sets were analyzed:

- **Validation panel (n=48 specimens across 9 sequencing runs):** 33 clinical specimens with reference-laboratory-confirmed diagnoses (*E. coli* n=5; *O. tsutsugamushi* n=6; *Rickettsia* spp n=4; *Leptospira* spp n=4; *B. pseudomallei* n=5; *S. pneumoniae* n=3; *S. suis* n=3; *C. burnetii* n=2; *Yersinia* spp n=1); 10 positive controls (4 PC_SINGLE; 5 PC_MIX8 = ZymoBIOMICS Microbial Community Standard [*Bacillus*, *Enterococcus*, *Escherichia*, *Limosilactobacillus*, *Listeria*, *Pseudomonas*, *Salmonella*, *Staphylococcus*]; 1 MIXED4 = *E. coli* + *P. aeruginosa* + *S. pneumoniae* + *S. suis*); 5 NTCs.
- **Study cohort (n=86 patient blood specimens):** patients with AFI for whom automated blood culture flagged positive but subculture onto aerobic solid media did not recover an organism. No anaerobic culture per local protocol. **The 16S samples are patient blood drawn during the same admission — not aliquots of the original positive blood culture bottle.**

### 2.2 Sample processing and sequencing

`[PLACEHOLDER - clinical team]`: blood specimen volume and container, DNA extraction kit and protocol, NCBI Scrubber host-DNA removal, 16S V1-V3 primer sequences (27F-based), library preparation kit, sequencing platform (Illumina MiSeq is the typical V1-V3 platform), and run configuration.

### 2.3 Bioinformatics pipeline

Reads were de-hosted via NCBI Scrubber, then classified by Centrifuge `[PLACEHOLDER - version and reference DB]` against a comprehensive bacterial reference. A genus was called "Detected" if its read count was ≥500 AND ≥5× the run-specific negative-template-control maximum (NCmax). Run-specific NCmax derivation is documented in the deployment SOP `[PLACEHOLDER - clinical team: specify median / max / excluded-by-QC rule]`.

For Rickettsiales, Centrifuge calls were supplemented by Minimap2 alignment-based rescue against a curated Rickettsiales 16S reference panel. For each genus, mapped reads, breadth of coverage, and fold over NCmax were computed. A two-tier rescue framework applies:

- **Tier 1 — Genus-level confirmed (`call = Confirmed`):** ≥100 mapped reads AND breadth ≥0.25 AND ≥5× NCmax AND V1-V3 provides clear genus-level discrimination. Reported as *Orientia* or *Rickettsia*.
- **Tier 2 — Order-level rescue (`call = Probable`):** alignment provides Rickettsiales-order-level evidence but V1-V3 does not allow confident genus assignment. Reported as "Rickettsiales detected (genus uncertain; recommend confirmatory qPCR)" — clinically actionable for empirical doxycycline coverage.
- **Below thresholds (`call = Negative`):** preserved in pipeline output for transparency but not reported clinically.

A threshold uniformity audit confirmed identical rescue threshold application across all 110 samples in the dataset (see Appendix and `THRESHOLD-UNIFORMITY-AUDIT.md`).

### 2.4 V4 decontamination filter

Genus detections passing the primary thresholds were processed through a four-tier decontamination filter informed by landmark low-biomass microbiome studies (1–5):

- **Tier A — High-confidence kit/skin/water contaminants (11 genera removed):** *Pseudomonas*, *Ralstonia*, *Bradyrhizobium*, *Sphingomonas*, *Stenotrophomonas*, *Methylobacterium*, *Acinetobacter*, *Cutibacterium*, *Staphylococcus*, *Corynebacterium*, *Brevundimonas*.
- **Tier A exception — *Burkholderia* species-level safeguard:** For every *Burkholderia* genus detection, the Centrifuge kreport is parsed at species rank for *B. pseudomallei* (NCBI taxonomy 28450). The detection is retained only if species reads ≥500 AND species reads exceed the same-run NTC species maximum. When retained, the record stores species-level reads (not genus total).
- **Tier B — NTC-only organisms (9 genera):** *Cereibacter*, *Thioclava*, *Bdellovibrio*, *Saltatorellus*, *Pseudogemmobacter*, *Minisyncoccus*, *Rhodoluna*, *Microbacterium*, *Arcanobacterium*.
- **Tier 1 — Ultra-low-abundance noise (median <0.5%):** 16 genera removed (full list in `MANUSCRIPT-METHODS-DECONTAMINATION.md`).
- **Tier 2 — Marginal organisms:** 27 genera retained only if detected in ≥2 samples AND each detection is ≥1.0% abundance.
- **Positive-control spike-in bypass:** For PC samples, Tier A removal is bypassed for any organism that is a documented spike-in (e.g., *Pseudomonas* in P-aeru_S5_L001 and PC_MIX8; *Staphylococcus* in PC_MIX8). PC concordance requires all expected spike-ins to be detected AND retained.

### 2.5 Concordance and statistical analysis

| Category | Concordance rule |
|---|---|
| Clinical | Concordant if at least one row matches the expected target organism AND passes V4. Within Rickettsiales, cross-genus rescue and order-level rescue both qualify. |
| Positive control | Concordant if all expected spike-in organisms are detected AND retained under V4 with PC bypass. |
| NTC | Concordant if no TaqMan Array Card bacterial target genus is retained (Bartonella, Brucella, Rickettsia, Orientia, Yersinia, Coxiella, Streptococcus, Salmonella, Escherichia, Burkholderia). |

Sensitivity and specificity were computed as binary per-sample outcomes with Wilson 95% confidence intervals. Inter-run reproducibility was assessed by per-run PC and NTC pass rates across 9 sequencing runs.

---

## 3. Results (~1,400 words)

### 3.1 Validation panel performance

The 48-sample validation panel was processed through the V4-filter-aware pipeline. Each sample was scored under the rule of record and the result was used to populate a standard 2×2 contingency table (clinical + PC = 43 expected-positive; 5 NTCs as expected-negative; TP = 31, FN = 12, TN = 5, FP = 0):

| Parameter | Value | 95% CI (Wilson) |
|---|---|---|
| **Sensitivity** (clinical + PC samples) | **31 / 43 = 72.1%** | 57.3% – 83.3% |
| &nbsp;&nbsp;Clinical-only subgroup | 21 / 33 = 63.6% | 46.0% – 78.5% |
| &nbsp;&nbsp;Positive-control subgroup | 10 / 10 = 100.0% | 72.2% – 100.0% |
| **Specificity** (NTC samples) | **5 / 5 = 100.0%** | 56.6% – 100.0% |
| **Positive predictive value (PPV)** | 31 / 31 = 100.0% | 89.0% – 100.0% |
| **Negative predictive value (NPV)**\* | 5 / 17 = 29.4% | 13.3% – 53.1% |
| Overall analytical accuracy | 36 / 48 = 75.0% | 61.2% – 85.1% |

*\* NPV is computed against the validation-panel composition (43 expected-positive : 5 expected-negative samples) and does not generalise to clinical prevalence. Sensitivity and specificity are the prevalence-independent metrics for inter-study comparison.*

**Organism-specific clinical sensitivity** ranged 0–100% (Table 1, body of `APPENDIX-VALIDATION-PANEL.md`). Notable patterns: (i) *E. coli* and *O. tsutsugamushi* achieved 100% sensitivity (with *Orientia* requiring two-tier rescue for several cases); (ii) three clinical samples (`00126_S6_L001`, `00369_S1_L001`, `25800370_S9_L001`) showed zero taxa, consistent with pre-sequencing failure rather than classification error; (iii) Streptococcus discordants reflect V1-V3 inability to resolve *S. pneumoniae* from *S. suis*; (iv) sample `00618_S7_L001` (Orientia-expected) carries a Rickettsia Tier-2 order-level rescue (`call = Probable`, breadth 0.2015) which is counted as concordant under the two-tier framework.

The V4 filter does not remove any organism from the validation panel that affects target detection. The 3 validation-panel *B. pseudomallei* samples (`09502813_S2_L001`, `09-0-02165`, `09700912_S3_L001`) trivially pass the species-level safeguard (species reads 13,744 / 20,061 / 65,016, all >> 500-read threshold). Without the positive-control spike-in bypass, the P-aeru_S5_L001 PC would have failed (since *Pseudomonas* is in Tier A); with the bypass, all 10 PCs are concordant — this is the bypass's specific contribution to the 72.1% sensitivity figure.

### 3.2 Study cohort: cohort size and pre-sequencing failures

The study cohort comprises **86 patient blood specimens** across 5 sequencing runs (1_and_2, 3, 4_and_5, 6_and_7, 8_and_9): 71/86 samples (82.6%) have ≥1 positive call; 15/86 (17.4%) have zero detections (pre-sequencing failures).

### 3.3 Study cohort: filter impact

The V4 filter was applied to all 217 genus-sample detections in the 71 samples with detections: 70 removed (32.3%), 147 retained. The largest contributors to removal were Tier A contaminants (61: predominantly *Cutibacterium*, *Staphylococcus*, *Brevundimonas*, *Acinetobacter*, *Corynebacterium*). Per-sample biomass distribution: 21 samples retain 100% of detected biomass; 2 retain 0% (everything removed as contaminants); 48 are mixed; 15 are empty (pre-seq failures).

Multiple NTCs in run 6_and_7 showed substantial reagent contamination (e.g., NTC2_ExDw_S13_L001: 78,691 *Leptospira*, 106,204 *Burkholderia* genus, 188,778 *Brevundimonas*). The pipeline reports `ntc_reads = 0` in study samples' `.calls.tsv` for these organisms in the same run, indicating the NCmax derivation excludes contaminated NTCs; the explicit rule should be documented in the SOP.

### 3.4 Study cohort: Rickettsiales rescue (key finding)

**Eleven of 86 study samples (12.8%; 11/71 of samples with detections = 15.5%) carry Rickettsiales rescue evidence** (Table 2). One sample (`16901195_S5_L001`) is a Tier-1 genus-level *Orientia* call (14,816 mapped reads, breadth 0.3235, alignment NTC 1,288 — ~12× headroom) — the cleanest Rickettsiales detection in the cohort. Ten additional samples are Tier-2 order-level "Rickettsiales detected" calls with breadth values 0.21–0.24 and per-sample confidence varying with alignment-NTC headroom.

| Sample | Run | Genus | Mapped reads | Breadth | Tier | Confidence |
|---|---|---|---|---|---|---|
| `16901195_S5_L001` | 8_and_9 | Orientia | 14,816 | 0.3235 | **Tier 1** | HIGH |
| `22600400_S6_L001` | 8_and_9 | Orientia | 9,747 | 0.2181 | Tier 2 | HIGH |
| `23900356_S5_L001` | 4_and_5 | Rickettsia | 441,173 | 0.3046 | Tier 2 | MODERATE |
| `23200430_S6_L001` | 4_and_5 | Rickettsia | 94,142 | 0.2189 | Tier 2 | LOW |
| `23200519_S12_L001` | 6_and_7 | Orientia | 2,492 | 0.2175 | Tier 2 | MODERATE |
| `09801652_S5_L001` | 6_and_7 | Orientia | 3,727 | 0.2282 | Tier 2 | MODERATE |
| `25800718_S8_L001` | 6_and_7 | Orientia | 4,143 | 0.2248 | Tier 2 | LOW |
| `10100409_S9_L001` | 8_and_9 | Orientia | 2,385 | 0.2221 | Tier 2 | LOW |
| `16601093_S7_L001` | 8_and_9 | Orientia | 2,284 | 0.2201 | Tier 2 | LOW |
| `23200736_S8_L001` | 8_and_9 | Orientia | 1,981 | 0.2141 | Tier 2 | LOW |
| `23900752_S9_L001` | 8_and_9 | Orientia | 1,684 | 0.2428 | Tier 2 | LOW |

Confidence reflects sample-vs-alignment-NTC headroom. All 11 samples warrant Rickettsiales-specific qPCR confirmation; doxycycline empirical coverage in this endemic context is standard practice and is supported by the alignment evidence.

### 3.5 Study cohort: species-level *B. pseudomallei* (corrected)

Four *Burkholderia* genus detections were observed in study samples. Species-level kreport parsing shows that the dominant Burkholderia species in all four are *B. cepacia* complex members (*B. contaminans*, *B. cenocepacia*, *B. multivorans*, *B. sola*, *B. cepacia*) — well-documented kit/water contaminants — while *B. pseudomallei* species reads are 0–8. In sample `23200430_S6_L001`, species-level *B. pseudomallei* reads are 8 against 12,789 reads of *B. cepacia* complex (run NTC carries 10 *B. pseudomallei* reads). **No study sample passes the species-level *B. pseudomallei* safeguard.** The 3 validation-panel *B. pseudomallei* samples carry 13,744 / 20,061 / 65,016 species reads, demonstrating that the assay detects melioidosis when present at meaningful abundance.

### 3.6 Study cohort: candidate fastidious-organism signals

- ***Mycoplasmopsis*** (cell-wall-deficient, fastidious organism class) was retained in **4 study samples** (mean abundance 39.78%, max 70.55%; read counts 537–6,568 across 5 samples). Mycoplasma-class organisms require sterol-supplemented media and 1–3 weeks of incubation; they would not be expected to grow on routine aerobic blood subculture within standard windows.
- ***Leptospira*** was detected in `09801652_S5_L001` (run 6_and_7) at 7,294 reads (22.78%). The same run's NTC2_ExDw_S13_L001 carries 78,691 *Leptospira* reads (10× the study sample). The pipeline reports `ntc_reads = 0` for this sample (excluding the contaminated NTC), but the NTC contamination is a real caveat. Treat as a candidate requiring orthogonal confirmation.
- ***Brucella*** was detected at sub-1% mean abundance (0.98% mean, max 1.96%) in 3 samples — near the noise floor. Candidate observations, not diagnoses.
- ***Streptococcus*** in 5 samples (mean 12.75%) is genus-level only; V1-V3 does not resolve *S. pneumoniae* / *S. suis* / fastidious vs. non-fastidious species.

### 3.7 Inter-run reproducibility

Inter-run reproducibility is 100% across 9 sequencing runs for the validation panel positive and negative controls. `[PLACEHOLDER - clinical team: per-run PC and NTC pass tables from deployment QC log]` for the study runs.

---

## 4. Discussion (~800 words)

### 4.1 The cohort: what the analysis shows

Three findings define the cohort's bacterial signal landscape:

**Rickettsiales involvement in ~13% of cases is the most prevalent identifiable etiology.** The 11 samples with Minimap2 rescue evidence — 1 Tier-1 *Orientia* and 10 Tier-2 order-level "Rickettsiales detected" — align directly with the expected endemic epidemiology of northeastern Thailand. Rickettsiales are obligate intracellular pathogens that cannot be cultured on routine blood agar, so their detection by 16S in a positive-bottle / no-subculture cohort is biologically coherent: the bottle signal could plausibly come from Rickettsiales metabolic activity in the broth, but no organism would grow on agar subculture. The Tier-1 detection (16901195_S5_L001, breadth 0.3235, ~12× NTC headroom) is the highest-confidence Rickettsiales finding. The 10 Tier-2 calls have variable confidence based on alignment-NTC headroom and warrant priority-stratified confirmation.

**No study sample carries *Burkholderia pseudomallei* above species-level background.** Species-level kreport parsing shows that the genus-level Burkholderia signal in the cohort is dominated by *B. cepacia* complex (kit contaminants); species-level *B. pseudomallei* reads are 0–8 per sample, below the run NTC. The 3 validation-panel *B. pseudomallei* samples carry 13,744–65,016 species reads, demonstrating that the assay detects melioidosis when present at meaningful abundance; this specific cohort simply does not contain it.

***Mycoplasmopsis* is the most prominent fastidious-organism-class signal in the cohort** (4 samples, mean abundance 39.78%, max 70.55%). Mycoplasma-class organisms are cell-wall-deficient by definition and require specialized media and extended incubation — exactly the organism profile that explains a positive-bottle / no-subculture phenotype. Cohort-level prevalence is ~5%, putting this in the hypothesis-generating category for follow-up testing.

### 4.2 Why filter design matters

The V4 filter's three core design choices each have outsized impact:

- **Aggressive Tier A** removes 61 of 70 detections in the cohort (predominantly *Cutibacterium*, *Staphylococcus*, *Brevundimonas*, *Acinetobacter*, *Corynebacterium*); these are universally documented kit/skin/water contaminants in low-biomass studies (1–5).
- **The *Burkholderia* species-level safeguard** is essential. Without species-rank parsing of the kreport, genus-level *Burkholderia* signals dominated by *B. cepacia* complex contamination could be misreported as melioidosis; with the safeguard, no study sample passes the threshold.
- **The positive-control spike-in bypass** restores PC sensitivity. Without the bypass, the P-aeru_S5_L001 PC would fail (because *Pseudomonas* is in Tier A); with the bypass, PC sensitivity is 10/10 and overall sensitivity (clinical + PC) is 31/43 = 72.1%.

### 4.3 What the data do not support

The analysis does not support cohort-level claims of "no Rickettsiales" (11 samples have rescue evidence), "strong support for fastidious-organism hypothesis" (the candidate signals affect ~10% of the cohort, not most of it), or "absence of anaerobes" (*Porphyromonas* and *Desulfovibrio* are retained; V1-V3 primer biases against many anaerobes preclude absolute negation). It also does not support causal claims linking 16S detections in patient blood to the bottle signals — the samples are not aliquots of the original bottles, and 16S does not distinguish viable from non-viable cells. The corrected reading frames the work as hypothesis-generating.

### 4.4 Clinical implications and deployment recommendations

For the 11 Rickettsiales-rescue samples, doxycycline empirical coverage in an endemic context is supported and species-specific qPCR is recommended. The 4 *Mycoplasmopsis* samples warrant Mycoplasma-specific PCR and sterol-supplemented culture; the *Leptospira* candidate warrants paired serology and species-specific qPCR (with the NTC contamination caveat); the 3 *Brucella* candidates warrant serology and Brucella-specific qPCR. The 15 zero-detection samples represent pre-sequencing failures requiring separate workup or repeat specimen collection.

Operational deployment should incorporate: (i) explicit documentation of the per-run NCmax derivation rule (especially in the presence of contaminated NTCs), (ii) a pre-sequencing QC step (Qubit + library qPCR) to distinguish biological negatives from extraction/library failures, (iii) continued species-level safeguards for clinically critical genera, and (iv) documentation of the positive-control spike-in bypass in the deployment SOP.

---

## 5. Limitations

1. **Sample source mismatch.** 16S samples are patient blood, not bottle aliquots — 16S detections are consistent with, but not causally linked to, the bottle signals.
2. **No viability assessment.** 16S detects DNA from viable, non-viable, and VBNC cells.
3. **V1-V3 primer biases.** Under-detection of certain Gram-positive anaerobes; limited Rickettsiales genus-level discrimination (mitigated by the two-tier framework but not eliminated).
4. **Genus-level resolution for most organisms.** Only *Burkholderia* receives species-level parsing; *S. pneumoniae* / *S. suis* and *Brucella* / *Leptospira* species-level discrimination is not resolved.
5. **NTC contamination in some study runs.** Pipeline NCmax derivation handles this implicitly; explicit documentation is needed.
6. **Cohort sample size (n=86) is a pilot.** Single-sample and few-sample findings are candidates for follow-up, not cohort-level prevalence claims.
7. **No clinical-outcome correlation.** Patient demographics, treatment, outcomes are not currently linked to 16S findings.

---

## 6. Conclusion

A 16S V1-V3 amplicon workflow with two-tier Rickettsiales rescue and species-aware decontamination achieves 72.1% sensitivity (31/43, 95% CI 57.3%–83.3%), 100% specificity (5/5 NTCs), and 100% positive predictive value (31/31) against a 48-sample reference panel. Applied to 86 AFI cases with positive blood culture / failed subculture, the workflow identifies Rickettsiales rescue evidence in 11 samples (~13% of cohort) — consistent with endemic epidemiology — and candidate fastidious-organism signals in an additional ~10%. The workflow is suitable as a complementary diagnostic for these cases; candidate detections require orthogonal confirmation before clinical reporting.

---

## References

1. Salter SJ, et al. Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. *BMC Biol.* 2014;12:87. PMID: 25387460.
2. Glassing A, Dowd SE, Galandiuk S, Davis B, Chiodini RJ. Inherent bacterial DNA contamination of extraction and sequencing reagents may affect interpretation of microbiota in low bacterial biomass samples. *Gut Pathog.* 2016;8:24. PMID: 27239228.
3. Lauder AP, Roche AM, Sherrill-Mix S, et al. Comparison of placenta samples with contamination controls does not provide evidence for a distinct placenta microbiota. *Microbiome.* 2016;4(1):29. PMID: 27338728.
4. de Goffau MC, Lager S, Salter SJ, et al. Recognizing the reagent microbiome. *Nat Microbiol.* 2018;3(8):851-853. PMID: 30046175.
5. Tan CCS, Ko KKK, Chen H, et al. No evidence for a common blood microbiome based on a population study of 9,770 healthy humans. *Nat Microbiol.* 2023;8(5):973-985. PMID: 36997797.

`[PLACEHOLDER - clinical team]`: software citations (Centrifuge, Minimap2, NCBI Scrubber), AFI epidemiology references for endemic Thailand, Rickettsiales clinical references.

---

## Appendix files (companion documents)

| File | Content |
|---|---|
| `APPENDIX-VALIDATION-PANEL.md` | Per-detection table for 48-sample validation panel (236 rows × 20 columns) + per-sample concordance summary. |
| `APPENDIX-STUDY-SAMPLES.md` | Per-detection table for 86-sample study cohort (336 rows × 22 columns) + per-sample biomass distribution summary. |
| `APPENDICES.xlsx` | Same content in Excel with 4 sheets. |
| `DECONTAMINATION-FILTER-REPORT-V4.txt` | Raw V4 filter output with summary statistics + Burkholderia species evidence per sample. |
| `afi_decontamination_filter_v4.py` | Executable V4 filter. |
| `generate_appendices.py` | Appendix generator (regenerates all tables from raw `.calls.tsv` and Centrifuge kreport files). |
| `MANUSCRIPT-FINAL-DRAFT.md` | Long-form (~13K-word) version of this manuscript for longer-venue submission. |
| `MANUSCRIPT-METHODS-DECONTAMINATION.md` | Extended methods detail. |
| `MANUSCRIPT-RESULTS-SECTION.md` | Extended results detail. |
| `MANUSCRIPT-DISCUSSION-SECTION.md` | Extended discussion detail. |

---

## Placeholders for clinical / wet-lab team

| Section | Content needed |
|---|---|
| Title page | Authors + affiliations; corresponding author; CoI; funding |
| §2.1 | IRB approval, ethics body, enrollment dates, geography, inclusion / exclusion criteria |
| §2.2 | Reference-laboratory diagnostic methods per pathogen |
| §2.2 wet-lab | Blood volume / container; DNA extraction kit; NCBI Scrubber version; 16S primers; library prep kit; sequencer + run config |
| §2.3 | Centrifuge version + reference DB; Minimap2 version + Rickettsiales reference panel; per-run NCmax derivation rule |
| §3.7 | Per-run PC and NTC pass tables for the deployment QC log |
| §4.4 | Confirmatory test results (if any) for the 11 Rickettsiales / 4 Mycoplasmopsis / 1 Leptospira / 3 Brucella samples; patient clinical outcomes |
| References | Software citations and AFI epidemiology references |

---

**Draft prepared 2026-05-12. Numerical results are derived from `.calls.tsv` and Centrifuge kreport files in `/Users/peerahemarajata/Downloads/AFI_P_Final/`.**
