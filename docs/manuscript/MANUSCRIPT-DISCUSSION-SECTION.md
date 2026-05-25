# Discussion: Interpreting 16S Detections in Blood Culture-Positive Samples with Failed Subculture Recovery

## Central Observation

This pilot analysis applied 16S rRNA gene sequencing to **86 acute-febrile-illness (AFI) cases** meeting a specific phenotype: positive automated blood culture bottle signal with subsequent failed subculture recovery on aerobic solid media (no anaerobic culture performed). Of these, 71 (82.6%) carry ≥1 positive call from either Centrifuger classification or Minimap2 alignment-based rescue, and 15 (17.4%) returned zero detected taxa attributed to pre-sequencing library failure. The 16S samples are patient blood drawn during the same admission; they are **not** aliquots of the original positive culture bottles. After applying the V4 decontamination filter (which removes well-documented kit/skin/water contaminants and gates *Burkholderia* genus retention on species-level *B. pseudomallei* evidence), 147 of 217 genus-sample detections were retained.

The signals most relevant to the AFI hypothesis (why did culture fail?) come from a small subset of samples:

- ***Mycoplasmopsis*** in 4 samples (mean 39.78%, max 70.55%) — the largest single fastidious-organism-class signal in the cohort
- ***Leptospira*** in 1 sample at 22.78%, with an important caveat about a contaminated NTC in the same run
- ***Brucella*** in 3 samples at 0.98% mean abundance, near the noise floor
- ***Streptococcus*** in 5 samples (genus-level only; V1-V3 cannot distinguish fastidious from non-fastidious species)

These detections are **consistent with**, but do not prove, the involvement of fastidious or slow-growing organisms in this cohort. They are best framed as candidate observations for confirmatory PCR, serology, or specialized culture, not as established diagnoses.

## What the Data Support

### *Mycoplasmopsis* as a Hypothesis-Generating Finding

The strongest cohort-level fastidious-organism signal in the V4-filtered results is *Mycoplasmopsis*: **4 of 86 samples (4.7%; 4 of 71 active samples = 5.6%)**, with mean abundance 39.78% and a maximum of 70.55%. Mycoplasma-class organisms are by definition fastidious: they lack a cell wall (cannot be Gram-stained reliably, do not retain peptidoglycan-targeted antibiotics like beta-lactams), require sterol-enriched media (e.g., Mycoplasma broth/agar with serum supplementation), and are slow-growing (visible colonies typically appear in 1-3 weeks, beyond the standard 5-7 day blood culture window). They would not be expected to grow on routine aerobic blood agar within the subculture window used by most hospital protocols.

This detection is plausibly relevant to the AFI phenotype: an organism present in blood, generating metabolic signal that could trip an automated bottle detector, but unable to grow on standard subculture. We emphasize that:

- The species-level identity of these *Mycoplasmopsis* detections has not been resolved in this analysis (genus-level only).
- *Mycoplasmopsis* read counts in same-run NTCs were variable (1-674 reads); each study-sample detection should be cross-checked against its specific run's NTC before clinical interpretation.
- The 4/86 (4.7%) detection rate is too small to support cohort-level causal claims; it is a hypothesis-generating observation.

### *Leptospira* as a Single-Case Candidate with NTC Caveat

One study sample (`09801652_S5_L001`, run 6_and_7) carries 7,294 *Leptospira* reads, corresponding to 22.78% of detected reads. This is consistent in abundance with a true clinical Leptospira detection. However, the same run contains a heavily contaminated NTC (`NTC2_ExDw_S13_L001`) with 78,691 *Leptospira* reads -- 10x the study sample's count. The pipeline's calls.tsv records `ntc_reads=0` for *Leptospira* in study samples of this run, suggesting NCmax was derived from a different (uncontaminated) NTC; but the existence of the contaminated NTC2_ExDw is a flag that should be addressed before any clinical claim.

The two validation-panel *Leptospira* samples (`00277_S4` and `09202659_S12`) with expected *Leptospira* detections are also in this overall dataset and were detected at 58.5% and 65.3% abundance — those serve as analytical positive controls. The study-sample candidate is suggestive but warrants independent serological/PCR confirmation before being reported.

### *Brucella* at Sub-1% Abundance

Three study samples contain *Brucella* genus signal at mean 0.98% abundance (max 1.96%). At this abundance level, the signal is at the upper edge of what could plausibly be reagent contamination. It is also at the lower edge of what could be true low-level bacteremia. Without serological confirmation or species-level resolution, this is a candidate observation, not a diagnosis, and does not constitute "strong support" for fastidious-organism involvement at the cohort level.

### *Streptococcus* as Ambiguous at V1-V3 Resolution

Five samples carry *Streptococcus* at mean 12.75%. The V1-V3 hypervariable region of 16S does not reliably distinguish *S. pneumoniae* and *S. suis* (well-documented in our validation panel: PC_SINGLE controls for both organisms hit the same *Streptococcus* genus call). Whether these *Streptococcus* detections include fastidious species relevant to the AFI hypothesis is not resolvable from this assay alone.

## What the Data Do Not Support

### Species-level *B. pseudomallei* status in the cohort

At the species level (Centrifuger kreport rank `S`), the *Burkholderia* genus signals in this cohort are dominated by *B. cepacia* complex species (e.g., 12,789 *B. cepacia* complex reads in `23200430_S6_L001` against 8 *B. pseudomallei* species reads). The same run's NTC contains 10 *B. pseudomallei* species reads, exceeding the study sample's count. The V4 filter, with species-level parsing and NTC comparison, removes this sample's *Burkholderia* detection as a *B. cepacia complex* contaminant. Across all 86 study samples, no detection passes the species-level *B. pseudomallei* safeguard. **The cohort does not contain melioidosis as detected by this assay.**

The 3 validation-panel *B. pseudomallei* samples (`09502813_S2`, `09-0-02165`, `09700912_S3`) carry 13,744–65,016 species-level *B. pseudomallei* reads each, demonstrating that the pipeline does detect *B. pseudomallei* when it is present at clinically meaningful abundance. The absence of such detections in the AFI cohort is itself a finding — it argues against melioidosis as a major contributor to the positive-culture/failed-subculture phenotype in this specific cohort.

### Anaerobes Cannot Be Excluded from the Differential

Two considerations make absolute exclusion of anaerobic contribution untenable:

1. **V1-V3 primer bias.** The 27F primer set is known to underdetect certain Gram-positive and Gram-negative anaerobic taxa. "Absent in this assay" is not "biologically absent."
2. **Anaerobic-genus signals are present in the retained set.** *Porphyromonas* (1 sample, 8.45%) and *Desulfovibrio* (1 sample, 0.78%) are anaerobic genera that were retained. Their presence, even at low frequency, undermines an absolute "no anaerobes" framing.

The assay's V1-V3 primer set is not optimal for anaerobe enumeration, and the available data neither rules anaerobes in nor out as contributors to the AFI phenotype in this cohort.

### Causal Claims Are Not Supported by the Design

Causal claims linking 16S detections to the bottle signal are not supported by the study design:

1. The 16S samples are **patient blood from the same admission**, not aliquots of the original blood culture bottles. We do not know whether the organisms detected by 16S in the patient's blood are the same organisms that triggered the bottle alarm.
2. 16S detects DNA from viable, non-viable, and viable-but-non-culturable cells alike. A high-abundance 16S detection does not establish that culturable, infection-causing organisms were present in the bottle.
3. We have no paired analysis of the original culture bottles (no Gram stain, no specialized-media subculture attempt, no 16S of the bottle broth itself).

The honest framing is: 16S of patient blood **identifies candidate organisms** that are detectable in the same patient at the time of admission. Whether these organisms caused the positive culture and whether they would be recoverable with appropriate culture conditions are questions for follow-up.

## Reframed Hypothesis Evaluation

The seven hypotheses we originally proposed (low-level pathogen, fastidious, anaerobes, slow-growing, VBNC, L-forms, polymicrobial competition) cannot be evaluated as binary cohort-level claims from this dataset. A more honest summary:

| Hypothesis | What this dataset supports | What it does not |
|---|---|---|
| **Fastidious organisms** | *Mycoplasmopsis* in 4 samples (4.7% of cohort; 5.6% of active samples) at 39.78% mean is suggestive of fastidious-organism involvement; *Leptospira* (n=1) and *Brucella* (n=3, sub-1%) are weaker candidate signals | Cohort-level claim that fastidious organisms explain >10% of cases; species-level confirmation; viability of detected organisms |
| **Slow-growing organisms** | Same organisms (Mycoplasma-class, Leptospira, Brucella) are also slow-growing; mechanism is coherent | Direct demonstration that bottle signal was generated by slow-growing organism specifically |
| **Low-level pathogen / inadequate inoculum** | Culturable organisms (E. coli, Klebsiella, Enterobacter) are detected at moderate abundance; abundance in blood != abundance in bottle | Direct correlation requires paired bottle/blood analysis |
| **Anaerobes** | Cannot be assessed: V1-V3 has primer bias against some anaerobes; a few anaerobic-genus signals (Porphyromonas, Desulfovibrio) are present | "Absent from cohort" claim is not supportable |
| **VBNC, L-forms** | Mechanism plausible (Mycoplasma-class are naturally cell-wall-deficient); no direct evidence in 16S data | Requires viability assays, microscopy of bottle broth |
| **Polymicrobial competition** | Mean 2.6 detections per sample post-filter; mechanism plausible | Requires paired bottle/blood analysis |

The fairest summary is that **fastidious / slow-growing / cell-wall-deficient organism involvement is plausible in a minority of cases (~5-10%) of this cohort, with Mycoplasma-class detection being the most prominent organism class.** Cohort-level "strong support" for any single hypothesis is not warranted by the data.

## Comparison with Expected AFI Pathogens

In northeastern Thailand, the major causes of acute febrile illness include rickettsial diseases (Orientia, Rickettsia), leptospirosis, melioidosis (*B. pseudomallei*), brucellosis, dengue and other arboviruses, enteric bacterial pathogens, and tuberculosis. Our validation panel demonstrates that the pipeline detects these organisms when they are present at adequate abundance:

- *Burkholderia pseudomallei*: detected at species level (13,744-65,016 species reads) in all 3 validation-panel cases with expected melioidosis
- *Orientia tsutsugamushi*: 100% detection across 7 expected cases (genus + order-level rescue combined)
- *Leptospira*: detected at 58-65% abundance in validation-panel positives

The study cohort itself yielded:
- No detectable *B. pseudomallei* at species level
- **Rickettsiales detected by Minimap2 alignment-based rescue in 11 of 86 study samples (12.8%; 11 of 71 samples with any positive call = 15.5%)**
- One candidate *Leptospira* with same-run NTC caveat
- Three low-abundance *Brucella* candidates
- Four *Mycoplasmopsis* detections at substantial abundance

The full per-sample Rickettsiales rescue table is:

| Sample | Run | Genus | Reads (alignment) | Breadth | Rescue tier | Alignment NTC reads | Confidence |
|---|---|---|---|---|---|---|---|
| 16901195_S5_L001 | 8_and_9 | Orientia | 14,816 | 0.3235 | **Tier 1: genus-level (Confirmed)** | 1,288 | HIGH |
| 22600400_S6_L001 | 8_and_9 | Orientia | 9,747 | 0.2181 | Tier 2: order-level | 1,288 | HIGH (abundance) |
| 23900356_S5_L001 | 4_and_5 | Rickettsia | 441,173 | 0.3046 | Tier 2: order-level | 189,002 | MODERATE (NTC ~43% of sample) |
| 23200430_S6_L001 | 4_and_5 | Rickettsia | 94,142 | 0.2189 | Tier 2: order-level | 53,908 | LOW (NTC ~57% of sample) |
| 23200519_S12_L001 | 6_and_7 | Orientia | 2,492 | 0.2175 | Tier 2: order-level | 1,143 | MODERATE |
| 23200736_S8_L001 | 8_and_9 | Orientia | 1,981 | 0.2141 | Tier 2: order-level | 1,288 | LOW |
| 23900752_S9_L001 | 8_and_9 | Orientia | 1,684 | 0.2428 | Tier 2: order-level | 1,288 | LOW |
| 25800718_S8_L001 | 6_and_7 | Orientia | 4,143 | 0.2248 | Tier 2: order-level | 2,202 | LOW |
| 10100409_S9_L001 | 8_and_9 | Orientia | 2,385 | 0.2221 | Tier 2: order-level | 1,288 | LOW |
| 16601093_S7_L001 | 8_and_9 | Orientia | 2,284 | 0.2201 | Tier 2: order-level | 1,288 | LOW |
| 09801652_S5_L001 | 6_and_7 | Orientia | 3,727 | 0.2282 | Tier 2: order-level | 1,143 | MODERATE |

**Of these, 1 sample (16901195_S5_L001) is a Tier 1 genus-level Orientia detection with strong abundance and NTC headroom — the cleanest Rickettsiales call in the study cohort.** The remaining 10 are Tier 2 order-level Rickettsiales detections; the breadth-of-coverage values (0.21-0.24) cluster just below the documented 0.25 Tier-1 threshold, which is why they appear as `Probable` rather than `Confirmed` in calls.tsv. Most carry meaningful NTC alignment background (1,143-189,002 NTC reads vs. sample reads), so absolute confidence varies from HIGH (for the one Tier 1 sample) to LOW (for samples where the NTC carries >=50% of the sample's signal).

**Clinical implication.** Rickettsiales is the canonical "cannot-miss" AFI etiology in endemic northeastern Thailand. Approximately 13% of the AFI cohort (~16% of samples with any positive detection) shows Rickettsiales alignment evidence, with one genus-level Orientia detection at high confidence and ten additional order-level rescues that warrant follow-up by Rickettsiales-specific qPCR (doxycycline empirical coverage is already standard in this clinical context).

This pattern -- one high-confidence Orientia detection, ~10 order-level Rickettsiales rescues, no melioidosis at species level, candidate Leptospira / Brucella / Mycoplasma in additional samples -- is consistent with a cohort in which a meaningful subset of patients (~13% of all samples, ~16% of samples with any detection) have Rickettsiales involvement that culture cannot recover (Rickettsiales are intracellular obligate parasites that do not grow on routine blood culture), alongside a smaller subset (~5-10%) showing fastidious / slow-growing organism candidate signals. The cohort is biologically consistent with expected AFI etiology in this region.

## Critical Limitations

1. **Sample size.** 86 cases (71 with detections) is a pilot. Single-sample detections (n=1 for *Leptospira*) cannot establish prevalence. Three-sample detections at sub-1% abundance (*Brucella*) cannot distinguish low-level signal from residual contamination.
2. **Sample source mismatch.** The 16S samples are patient blood, not the original culture bottles. We cannot directly link 16S detections to the bottles that flagged positive.
3. **No viability assessment.** 16S cannot distinguish viable cells from dead DNA or VBNC cells.
4. **V1-V3 primer biases.** The 27F primer set has known under-detection of certain anaerobes and some Gram-positive taxa.
5. **NTC contamination in run 6_and_7.** NTC2_ExDw_S13_L001 carried substantial contamination (78,691 *Leptospira* reads, 106,204 *Burkholderia* reads, 188,778 *Brevundimonas* reads). The upstream pipeline's NCmax derivation appears to exclude this NTC, but the assumption should be revisited.
6. **No serological / PCR / specialized-culture confirmation.** Every candidate organism in this analysis would need orthogonal confirmation before clinical reporting.
7. **No clinical-outcome correlation.** We do not know whether patients with *Mycoplasmopsis*, *Leptospira*, or *Brucella* candidate detections had clinical presentations or treatment responses consistent with these organisms.
8. **Filter dependence on dataset-specific tier assembly.** The Tier 1 and Tier 2 lists were assembled from this dataset's observed distributions. Performance on a new specimen type or kit batch is unknown.
9. **B. pseudomallei detection limit in this assay.** The validation-panel *B. pseudomallei* cases have species-level reads in the 13,744-65,016 range. The species-level detection floor (where the assay reliably distinguishes *B. pseudomallei* from *B. cepacia* complex contamination) is not formally established; 500 reads is used here as a practical threshold but may be too generous or too strict.

## Recommendations for Follow-Up

### Immediate Confirmatory Testing

1. **Mycoplasmopsis-positive samples (n=4):** PCR for *Mycoplasma* / Mycoplasmataceae targets on original blood samples; species-level kreport interrogation; consider Mycoplasma-specific culture if material remains.
2. **Leptospira candidate (`09801652_S5_L001`):** serology (IgM/IgG paired acute/convalescent if available) and Leptospira-specific qPCR; same-run NTC follow-up to clarify why NTC2_ExDw_S13_L001 had 78,691 *Leptospira* reads.
3. **Brucella candidates (n=3):** serology and Brucella-specific qPCR.

### Pipeline / Methods Improvements

4. **Run-specific NCmax should explicitly handle contaminated NTCs.** The presence of NTC2_ExDw_S13_L001 with 78,691 *Leptospira* reads but `ntc_reads=0` in study-sample `.calls.tsv` entries indicates the NCmax derivation skips this NTC -- this should be made transparent.
5. **Species-level reporting for clinically critical genera.** *Burkholderia* now has a species-level safeguard; analogous treatment for *Brucella* (species discrimination of *B. melitensis* / *B. abortus* / *B. suis* / *B. canis*) and *Leptospira* (pathogenic vs. saprophytic species) would strengthen interpretation.
6. **Cross-run negative-control profiling.** A formal NC vs. NC comparison across runs (which we did informally here) should be part of routine QC, with thresholds for flagging contaminated NTCs.

### Study Design Improvements

7. **Paired bottle and blood 16S.** Future cohorts should include 16S of the original positive bottle broth, not just patient blood, so direct comparison is possible.
8. **Specialized-culture attempts.** Subculture aliquots onto Mycoplasma broth (sterol-enriched), Fletcher/EMJH media (Leptospira), Brucella agar (Brucella) would convert "candidate detection" into "confirmed isolation."
9. **Clinical metadata.** Treatment history, serology, outcomes, and exposure history are needed to interpret which candidate organisms are biologically plausible in each case.

## Conclusion

This is a pilot 16S analysis of **86 AFI cases** (71 with ≥1 positive call, 15 zero-detection pre-sequencing failures) with positive blood cultures but failed subculture recovery. After correction of contaminant filtering, the most prominent organism-class signal of clinical interest is Rickettsiales alignment-rescue evidence in **11 of 86 samples (12.8%; 11/71 active = 15.5%)**, followed by *Mycoplasmopsis* in 4 samples (mean 39.78%), with weaker candidate signals from *Leptospira* (n=1, NTC caveat) and *Brucella* (n=3, sub-1% abundance). No study sample contains *B. pseudomallei* at species level above background.

The detections are consistent with the working hypothesis that fastidious, slow-growing, or cell-wall-deficient organisms contribute to the positive-culture / failed-subculture phenotype in some AFI cases, but the magnitude (~10-15% of cohort) and the absence of confirmatory testing preclude cohort-level "strong support" claims. The honest framing is **hypothesis-generating**: candidate organisms identified for confirmatory PCR, serology, and specialized culture in a future targeted study.

The validation panel demonstrates the underlying 16S pipeline detects expected AFI pathogens (notably *B. pseudomallei*, *Orientia*, *Leptospira*, *Escherichia*) when they are present at adequate abundance, with 100% specificity. The pipeline itself is suitable for analytical use. The interpretive narrative around the 86-sample study cohort -- and what it does or does not say about Rickettsiales and fastidious-organism involvement in failed-culture AFI -- should be reframed as pilot data rather than clinical validation.
