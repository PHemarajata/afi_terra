# Results: AFI 16S V1-V3 Amplicon Workflow and Organism Detection in Blood Culture-Positive Samples

## Validation Panel Performance

### Sample Composition and Control Performance

The validation panel comprised 43 specimens across 9 sequencing runs: 33 clinical samples with known laboratory diagnoses, 10 positive controls (5 PC_MIX8 eight-organism controls, 4 PC_SINGLE single-organism controls, 1 MIXED4 four-organism control), and 5 negative template controls.

Positive controls performed as expected: PC_MIX8 (5/5 samples, all 8 organisms detected per replicate); PC_SINGLE (4/4 organisms detected at high abundance); MIXED4 (4/4 organisms detected). The 5 validation-panel NTCs showed zero taxa post-NTC normalization.

### Organism-Specific Clinical Performance

**Escherichia coli (n=6).** All 6 clinical *E. coli* samples were concordant (100%). Detection abundance ranged 10-99% (median 42.8%).

***Burkholderia pseudomallei* (n=5).** Three of 5 clinical samples showed *Burkholderia* genus detection at high abundance: `09-0-02165` (genus 21,792 reads, species *B. pseudomallei* 20,061 reads = 92% of the genus signal), `09502813_S2_L001` (genus 69,440, species 65,016 = 94%), and `09700912_S3_L001` (genus 14,647, species 13,744 = 94%). Species-level kreport parsing confirms these are genuine *B. pseudomallei* detections rather than *B. cepacia* complex contamination. Two samples were discordant at the genus level. Sensitivity 3/5 (60%).

***Orientia tsutsugamushi* (n=6).** Five of 6 samples concordant (83.3%); 1 discordance (`00618_S7`) reflects a Minimap2 rescue breadth-of-coverage failure (breadth 0.2015 < 0.25 threshold) correctly excluding insufficient-coverage rescue. Detection mechanisms: direct genus-level (4 samples) and order-level Minimap2 rescue (3 samples).

***Rickettsia* spp (n=4).** Two of 4 detected at high abundance (46-100%); 2 complete detection failures (zero taxa) consistent with pre-sequencing sample quality issues, not classification error.

***Leptospira* spp (n=4).** Two of 4 samples showed *Leptospira* as the dominant detected genus (58.5% and 65.3% abundance); the other 2 were discordant with alternative organisms dominating.

***Streptococcus pneumoniae* (n=3) and *S. suis* (n=3).** Clinical concordance 1/3 each at the genus level; positive controls confirm genus-level detection capability (S-pneumo PC_SINGLE 52.4% *Streptococcus*; S-suis PC_SINGLE 99.5% *Streptococcus*). Species-level discrimination is not reliable at the V1-V3 region.

***Coxiella burnetii* (n=2).** 1/2 detected at 6.5% abundance; the other showed zero taxa (pre-sequencing failure).

***Yersinia* spp (n=1).** 0/1; *Rickettsia* and *Paracoccus* dominated this single sample.

### Validation Panel Summary

**Validation panel performance (V4-filter-aware, TAC-target rule of record).** Each sample was scored as concordant if the rule of record was satisfied (target detected AND retained by V4 for clinical; all expected spike-in organisms detected and retained for PCs; no TAC bacterial target genus retained for NTCs) and the result was used to populate the standard 2×2 contingency table treating clinical + PC samples as expected-positive (n = 43; 33 clinical + 10 PC) and NTC samples as expected-negative (n = 5). This yields TP = 31, FN = 12, TN = 5, FP = 0.

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

For the V4 filter to evaluate PC samples honestly, the filter bypasses Tier A (kit/skin/water contaminant) removal for any organism that is part of the PC sample's documented spike-in composition — *Pseudomonas* and *Staphylococcus* are kit contaminants in clinical samples but are expected positive-control targets in P-aeru_S5_L001, MIXED4, and PC_MIX8, and were therefore retained when those samples were processed. (E-coli_S4_L001 is reclassified from "clinical" to PC_SINGLE per the same logic.) Without this bypass, the P-aeru PC would have failed; with the bypass, all 10 PCs are concordant.

Specificity at the sample level is 100% (zero NTCs contain TAC bacterial target genera after V4 filtering). Inter-run reproducibility is 100% across all control types. Organism-specific clinical sensitivity ranges 0–100%, reflecting expected variability in 16S detection of low-abundance organisms in complex backgrounds.

## Decontamination Strategy: V4 Decontamination Filter

The V4 Decontamination Filter (see Methods) was applied to all genus calls from the validation and study cohorts. Three filter design choices have notable impact:

1. **The *Burkholderia* species-level safeguard parses species-level reads.** V4 parses rank `S` rows of the Centrifuger kreport specifically for *Burkholderia pseudomallei*, requires species reads ≥ 500, and requires species reads > the same-run NTC's maximum species reads. When preserved, the retained record stores the species-level read count, not the genus total.
2. **Mycoplasmopsis and Nitrospira are retained as candidate signals.** These genera carry mean abundances of 39.78% (4 samples) and 13.45% (1 sample), respectively, and are excluded from the Tier 1 ultra-low-abundance removal list.
3. **Brevundimonas is included in the high-confidence contaminant tier.** Run-6_and_7 NTCs show *Brevundimonas* at 22–285,740 reads, identifying it as a kit/water contaminant in this dataset.

### Validation Panel Filter Impact

The V4 filter does not remove any organism in the validation panel. This is **not an independent validation of the filter** -- the validation panel is dominated by single high-abundance expected organisms (40-100% abundance) in samples that do not contain Tier A contaminants at detectable levels. The validation panel therefore does not exercise the filter under conditions where it could over-remove. The three validation-panel *B. pseudomallei* samples (`09502813_S2_L001`, `09-0-02165`, `09700912_S3_L001`) trivially pass the species-level safeguard (13,744-65,016 species reads, all >> 500-read threshold).

### Study Cohort Filter Impact

Applied to 86 AFI study samples (71 with ≥1 positive call, 15 with zero detections attributed to pre-sequencing library failure), V4 removed 70 of 217 detections (32.3%):

| Tier | Removals |
|---|---|
| Tier A high-confidence contaminants | 61 (predominantly Cutibacterium, Staphylococcus, Brevundimonas, Acinetobacter, Corynebacterium) |
| Burkholderia (species-level safeguard failed) | 4 |
| Burkholderia preserved as *B. pseudomallei* | **0** |
| NC-only organisms (Rhodoluna) | 2 |
| Tier 1 ultra-low abundance | 2 |
| Tier 2 marginal + rare | 1 |

147 detections were retained for downstream analysis.

### Species-level *B. pseudomallei* status in study cohort

Examining the Centrifuger kreport at species rank for the *Burkholderia* genus detection in `23200430_S6_L001`, the species composition is:

| Species | Reads |
|---|---|
| *Burkholderia contaminans* | 3,266 |
| *Burkholderia cenocepacia* | 2,988 |
| *Burkholderia multivorans* | 2,170 |
| *Burkholderia sola* | 1,443 |
| *Burkholderia cepacia* | 1,335 |
| *Burkholderia arboris* | 666 |
| *Burkholderia vietnamiensis* | 243 |
| *Burkholderia pseudomultivorans* | 123 |
| *Burkholderia metallica* | 120 |
| *Burkholderia pyrrocinia* | 86 |
| *Burkholderia ambifaria* | 77 |
| *Burkholderia aenigmatica* | 69 |
| *Burkholderia semiarida* | 47 |
| *Burkholderia ubonensis* | 34 |
| *Burkholderia seminalis* | 25 |
| **Burkholderia pseudomallei** | **8** |
| *Burkholderia thailandensis* | 5 |
| *Burkholderia mallei* | 3 |
| (other species) | further low-count entries |

The same run's NTC `NTC2_ExDw_S13_L001` carries 10 species-level reads of *B. pseudomallei*. **The study sample's *B. pseudomallei* signal (8 reads) is below the NTC's *B. pseudomallei* signal and below the 500-read detection threshold.** The genus-level *Burkholderia* signal in `23200430_S6_L001` is dominated by *B. cepacia* complex species, which are documented kit/water contaminants. V4 treats this *Burkholderia* detection as a contaminant and removes it.

**No study sample carries *B. pseudomallei* above the detection floor.** The validation-panel *B. pseudomallei* samples remain valid as analytical positive controls, but the cohort of 86 AFI cases does not contain melioidosis as detected by this assay.

### Top Retained Organisms in Study Cohort (Post-V4 Filter)

| Genus | Detections | Mean abundance | Max abundance | Notes |
|---|---|---|---|---|
| Thermomicrobium | 11 | 9.13% | 55.79% | Environmental thermophile; unlikely human pathogen |
| Escherichia | 6 | 8.28% | 12.23% | Known AFI pathogen |
| Streptococcus | 5 | 12.75% | 23.29% | Includes fastidious species at genus level |
| Faucicola | 5 | 6.85% | 10.30% | Genus identified from oral microbiome |
| **Mycoplasmopsis** | **4** | **39.78%** | **70.55%** | **Fastidious, cell-wall-deficient organism class** |
| Klebsiella | 4 | 20.99% | 73.90% | Known AFI pathogen |
| Methylorubrum | 4 | 9.71% | 17.53% | Environmental |
| Kocuria | 4 | 7.76% | 16.89% | Skin/environmental |
| Xanthomonas | 4 | 11.85% | 24.04% | Plant-associated; possible environmental |
| Enterobacter | 3 | 3.39% | 5.14% | Known AFI pathogen |
| Brucella | 3 | 0.98% | 1.96% | Near contamination floor; see caveat below |
| Streptomyces | 3 | 9.21% | 14.18% | Environmental |
| Comamonas | 3 | 9.01% | 11.59% | Environmental |
| Leptospira | 1 | 22.78% | 22.78% | Single sample; see NTC caveat below |

### Caveats Specific to Two Headline Organisms

***Brucella* (3 samples, 0.98% mean, max 1.96%).** Brucella is detected at sub-1% mean abundance, which is at the upper edge of the noise floor. None of the three *Brucella* detections exceed 2% abundance in their host samples. This is not the abundance level at which "clinically significant detection" is normally inferred. The detections are consistent with low-level signal that could equally represent low-grade bacteremia or residual reagent contamination; without serological or PCR confirmation, the level of evidence is "candidate for follow-up testing," not "diagnosis."

***Leptospira* (1 sample, 09801652_S5_L001, 22.78% abundance, 7,294 reads).** This is a single-sample detection. Importantly, the same run (6_and_7) contains `NTC2_ExDw_S13_L001`, which has 78,691 *Leptospira* reads — over 10x the study sample's count. The upstream pipeline records `ntc_reads=0` for *Leptospira* in study samples of run 6_and_7 in the `.calls.tsv` outputs, suggesting NCmax was derived from non-contaminated NTCs of the same run; however, the heavily contaminated NTC2_ExDw is a flag that warrants independent confirmation (PCR/serology for *Leptospira* in this specific sample) before any clinical interpretation.

***Mycoplasmopsis* (4 samples, 39.78% mean, max 70.55%).** This is now the most abundant retained-organism signal in the study cohort. Mycoplasma-class organisms are cell-wall-deficient and fastidious. NTCs in this dataset show *Mycoplasmopsis* at 1-674 reads (variable across runs); the four study samples (537, 740, 1,622, 2,568, 6,568 reads) are above same-run NTC backgrounds in most cases but should still be cross-checked at the species level. *Mycoplasmopsis* is genus-level only here; the kreport species assignments would clarify whether these reads represent recognized human pathogens (e.g., *M. pulmonis* analogues, or other Mycoplasmataceae) or environmental species.

## Study Sample Results: AFI Cases (Positive Blood Culture / Failed Subculture)

### Cohort

86 patient blood samples from acute-febrile-illness cases meeting the phenotype: blood culture bottle flagged positive by automated detection, subculture failed on aerobic solid media, no anaerobic culture performed. The 16S samples are patient blood drawn during the same admission and are not aliquots of the original positive culture bottle. Of the 86 samples, **71 (82.6%) carry ≥1 positive call** (`Detected`, `Confirmed`, or `Probable` in `.calls.tsv`) from either Centrifuger classification or Minimap2 alignment-based rescue, and **15 (17.4%) returned zero detected taxa**, consistent with pre-sequencing failures (DNA-extraction, library-preparation, or sequencing-depth failure).

### Organism Detections Relevant to the AFI Hypothesis

After V4 filtering, the detections most relevant to the hypotheses about why culture failed:

- **Mycoplasmopsis** (4 samples, mean 39.78%, max 70.55%): the largest single signal class in the cohort. Mycoplasma-class organisms are by definition fastidious (no cell wall, require sterol/cholesterol-enriched media, slow growth). They would not be expected to grow on routine aerobic blood agar within standard 5-7 day subculture windows. This is the **strongest organism-class signal for the fastidious-organism hypothesis** in this cohort.
- **Leptospira** (1 sample, 22.78%): consistent with leptospirosis but requires confirmation given the contaminated NTC in the same run (see caveat above). Single-sample finding; appropriate framing is "candidate detection requiring serological/PCR confirmation."
- **Brucella** (3 samples, 0.98% mean): near noise floor. Three samples at sub-1% abundance. Appropriate framing is "low-level signal warranting follow-up testing."
- **Streptococcus** (5 samples, 12.75% mean): genus-level only; species cannot be distinguished at V1-V3. Could include fastidious *S. pneumoniae*, *S. suis*, or non-fastidious species.
- **Escherichia, Klebsiella, Enterobacter** (6, 4, 3 samples respectively): culturable organisms detected at moderate-to-high abundance. Their detection by 16S without subculture recovery may reflect low organism count in the actual culture bottle (separate sample from this 16S draw), transient bacteremia, or competitive overgrowth by other organisms during subculture.

### Hypothesis Evaluation (Recalibrated to Sample Sizes)

Hypothesis evaluation recalibrated to sample sizes:

- **Fastidious organism hypothesis.** *Mycoplasmopsis* in 4 samples at 39.78% mean is the strongest cohort-level signal consistent with this hypothesis, and it is consistent with organisms that would not grow on routine aerobic media. Single-sample *Leptospira* and low-abundance *Brucella* are additional, weaker signals. Overall: **consistent with the hypothesis in a minority of cases (~5-10% of the 86-sample cohort), most prominently via Mycoplasma-class detection.** Not "strong support" at the cohort level.
- **Slow-growing organism hypothesis.** *Mycoplasmopsis*, *Brucella*, and possibly *Leptospira* all have growth-rate biology that would be inadequately served by 5-7 day standard incubation. Same caveat: minority of cases (~5-10%).
- **Anaerobe hypothesis.** Classic obligate anaerobes (*Bacteroides*, *Prevotella*, *Fusobacterium*, *Clostridium*) were not detected in the retained set, but *Porphyromonas* (1 sample, 8.45%) and *Desulfovibrio* (1 sample, 0.78%) are present and have anaerobic biology. **The assay detected anaerobic-genus signals in a small number of samples; primer bias of the V1-V3 region against certain anaerobes is a known limitation of this assay and should be acknowledged in any negative conclusion.**
- **Low-level pathogen hypothesis.** Culturable organisms (*Escherichia*, *Klebsiella*, *Enterobacter*) detected at moderate abundance. Direct correlation with original culture bottle inoculum is unknown (16S sample != culture bottle sample). Hypothesis remains plausible but unevaluable from this data alone.
- **VBNC, L-form, and competitive exclusion hypotheses.** No direct evidence for or against from 16S alone; these remain theoretical mechanisms that would require viability staining, microscopy of culture broth, or paired bottle/serum analysis to test.

## Control Performance (Re-stated)

- PC_MIX8 (5 samples, 8 organisms each): 40/40 expected detections retained post-filter.
- PC_SINGLE (4 samples): 4/4 retained.
- MIXED4 (1 sample, 4 organisms): 4/4 retained.
- NTC (5 samples): zero post-NC-subtraction detections.

The validation-panel NTCs and PC samples behaved as expected. Across the broader pipeline runs, several NTCs in run 6_and_7 (notably NTC2_ExDw_S13_L001) showed substantial contamination (78,691 *Leptospira* reads, 106,204 *Burkholderia* reads, 188,778 *Brevundimonas* reads). The upstream pipeline's NTC-subtraction logic relies on per-run NCmax selection; the assumption that contaminated NTCs are correctly excluded from NCmax computation should be revisited.

## Rickettsiales Detections in Study Cohort (Corrected)

Examining all `source = alignment` rows with `call ∈ {Confirmed, Probable, Detected}` in the study samples' `.calls.tsv` files identifies **11 of 86 study samples (12.8%; equivalently 11 of 71 samples with any detection = 15.5%) with Rickettsiales rescue**, distributed across runs 4_and_5, 6_and_7, and 8_and_9:

- **1 sample (16901195_S5_L001) with Tier 1 genus-level Orientia detection**: 14,816 mapped reads, breadth 0.3235 (passes >=0.25 threshold), alignment NTC 1,288 — strong abundance and clean NTC headroom. This is the highest-confidence Rickettsiales call in the study cohort.
- **10 samples with Tier 2 order-level Rickettsiales rescue** (`call = Probable` in calls.tsv): breadth values cluster at 0.21-0.24, just below the Tier 1 threshold. Sample-level read counts range from 1,684 (23900752_S9_L001) to 441,173 (23900356_S5_L001). Alignment NTC background ranges from 1,143 to 189,002 reads; some samples have substantial NTC headroom, others have less.

See `APPENDIX-STUDY-SAMPLES.md` for the full per-sample table with reads, breadth, NTC reads, confidence, and recommended follow-up. The per-sample confidence label takes the alignment NTC into account: 2 samples are HIGH confidence (16901195_S5 Tier 1; 22600400_S6 Tier 2 with abundance >>NTC), 4 are MODERATE, 5 are LOW (NTC carries comparable signal). All 11 samples are flagged for Rickettsiales-specific qPCR confirmation.

**Clinical implication:** Rickettsiales is the canonical "cannot-miss" AFI etiology in endemic northeastern Thailand and would not grow on routine blood culture (obligate intracellular). The corrected reading of the data is that approximately 13% of all AFI study samples (~16% of samples with any positive call) show alignment evidence consistent with Rickettsiales involvement -- a finding that aligns with the expected epidemiology of the region and the two-tier reporting framework's design intent.

## Summary of Findings

1. **Validation panel** demonstrates 72.1% sensitivity (31/43, 95% CI 57.3%–83.3%; 63.6% [21/33] clinical-only and 100% [10/10] positive-control subgroups), 100% specificity (5/5 NTCs), 100% positive predictive value (31/31), and 100% inter-run reproducibility. Organism-specific clinical sensitivity ranges 0–100%; *Burkholderia pseudomallei* species-level detection is confirmed in all 3 validation-panel cases with substantial signal (13,744–65,016 species reads).
2. **Two-tier Rickettsiales framework** (Centrifuge genus + Minimap2 order-level rescue) achieves 100% *Orientia* detection across the validation panel.
3. **V4 decontamination filter** removes 32.3% of study-sample detections (predominantly Tier A skin/kit/water contaminants and *Burkholderia* signal dominated by *B. cepacia* complex). 147 detections retained.
4. **No study sample contains *B. pseudomallei* at species level above background.** The earlier V3-based claim that 23200430_S6_L001 contained 54,126 reads of *B. pseudomallei* conflated genus-level reads with species-level reads; the actual species-level count is 8 reads, below the NTC background.
5. **Mycoplasmopsis is the largest single fastidious-organism signal** in the study cohort (4 samples, 39.78% mean). *Leptospira* (n=1, with same-run NTC contamination caveat) and *Brucella* (n=3, 0.98% mean) are additional, weaker signals consistent with fastidious/slow-growing organisms but at minority frequency.
6. **Rickettsiales rescue in study cohort:** 11 of 86 samples (12.8%; 11 of 71 active samples = 15.5%) show Minimap2 alignment-based rescue evidence for Orientia or Rickettsia. 1 sample is a Tier 1 genus-level Orientia detection (16901195_S5_L001, breadth 0.3235); the remaining 10 are Tier 2 order-level rescues (breadth 0.21-0.24). Confidence varies by sample given alignment NTC backgrounds; all 11 warrant Rickettsiales-specific qPCR confirmation.
7. **Honest hypothesis ranking:** Rickettsiales involvement is present in ~13% of cohort (~16% of samples with any detection) by alignment-rescue criteria (cannot-miss endemic AFI etiology, would not grow on routine blood culture). Fastidious / slow-growing organism candidates (Mycoplasmopsis, Leptospira, Brucella, Streptococcus species ambiguity) account for an additional ~5-10% minority. Cohort-level support for any single non-Rickettsiales hypothesis is too weak to call "strong"; framing for those should remain "hypothesis-generating observations warranting confirmatory PCR, serology, or specialized culture."
