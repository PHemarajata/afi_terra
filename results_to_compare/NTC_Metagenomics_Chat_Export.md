# NTC Strategy in 16S Metagenomics — Full Chat Export
**Date:** 2026-05-08  
**Topic:** Negative Controls, NTC Design, and Background Subtraction in 16S rRNA Amplicon Sequencing

---

## 💬 Conversation 1

### 🧑 User Question
> Please research the literature re: the use of negative controls in metagenomics studies (especially in 16S metagenomics) and see what best practices are. I am using an accumulative NTC model, where all NTCs from all runs were accumulated to best capture the contamination that might be transient in the environment and apply the composite NTC on all samples for background subtraction. Is this an ok practice?

---

### 🤖 Agent Response

**Search Summary:**
- **SciSpace Basic Search**: 100 results
- **SciSpace Full Text Search**: 100 results
- **Google Scholar**: 20 results
- **PubMed**: 20 results
- **Total after deduplication**: 83 unique papers ranked by relevance
- **Zotero / Mendeley**: [Connect your account](https://www.scispace.com?settings=connected-apps)

---

#### TL;DR
Accumulative NTC models are used in practice and can capture transient contamination, but run-specific controls are the gold standard. The composite NTC approach is most defensible as a supplementary strategy, not a sole background subtraction method.

---

#### Key Takeaways
- Include negative controls at **every extraction and PCR batch** — contamination varies by kit lot, operator, and time
- **Extraction blanks** and **PCR blanks** capture different contamination sources and should both be used
- The **decontam** tool (prevalence + frequency methods) is the most widely validated statistical approach for contaminant removal
- **Reagent water and elution buffers** are known contamination sources — documented in multiple studies
- Accumulative NTCs can help detect **transient contamination** but risk masking run-specific signals or over-removing taxa

---

#### 1. Best Practices for NTC Design and Processing

Negative controls (NTCs) should be included at **every stage of the workflow** — extraction, PCR, and library preparation. The literature distinguishes between:

| Control Type | What It Captures |
|---|---|
| **Extraction blank** | Kit reagent DNA, environmental DNA introduced during extraction |
| **PCR blank / No-template control (NTC)** | Amplification reagent contaminants, cross-contamination during PCR setup |
| **Library prep blank** | Contamination introduced during indexing/library preparation |

Salter et al. (2014) [9] demonstrated that reagent and laboratory contamination can critically impact microbiome analyses, particularly in **low-biomass samples**, where contaminant reads can dominate the signal. The 2025 Nature Microbiology guidelines [1] explicitly recommend processing negative controls through the **entire workflow** (extraction → PCR → sequencing) to capture all contamination stages.

---

#### 2. Background Subtraction and Decontamination Methods

Several approaches are used in the literature:

- **Prevalence-based removal**: Taxa present in NTCs above a threshold are flagged as contaminants
- **Frequency-based removal (decontam)**: Uses the inverse relationship between ASV frequency and total sample DNA concentration to identify contaminants [5][7]
- **SourceTracker**: Probabilistic method that estimates the proportion of reads originating from known contamination sources — highly effective (>98% removal) when source environments are well-defined [5]
- **microDecon**: Specifically designed for low-biomass samples; uses NTC profiles directly for subtraction
- **Read-count thresholds**: Simple but not recommended as a sole method — the literature warns against arbitrary cutoffs [1][5]

---

#### 3. Is the Accumulative NTC Model Valid?

The accumulative/composite NTC approach — pooling NTCs across all runs to build a composite contamination profile — is used in practice and has logical merit:

**Arguments in favor:**
- Captures **transient contamination** that may not appear in every run's NTC ✅
- Builds a more **comprehensive contamination reference** when individual run NTCs have low read counts ✅
- Useful when the same kit and environment are used throughout ✅

**Caveats from the literature:**
- Contamination profiles can vary **by kit lot, day, and operator** even with the same kit [1][2][3]
- Pooling risks **masking run-specific contaminants** (diluted out of composite) ⚠️
- May cause **over-removal** if a taxon appeared transiently in one run's NTC but is biologically real in other runs ⚠️
- Run-specific controls are the **gold standard** per current guidelines [1][2]

**Verdict:** The accumulative NTC model is a **reasonable supplementary approach**, especially for capturing transient contamination. However, it should ideally be used **alongside** run-specific controls, not as a replacement.

---

#### 4. Recommended Tools

| Tool | Method | Best For |
|---|---|---|
| **decontam** (R package) | Prevalence + Frequency | General-purpose; widely validated |
| **SourceTracker** | Probabilistic source tracking | When source environments are known |
| **microDecon** | NTC-based direct subtraction | Low-biomass samples |
| **Recentrifuge** | Comparative analysis | Metagenomics (shotgun) |

---

#### 5. Known Limitations and Caveats

- NTCs can yield **false negatives** — not all contaminants appear in every NTC
- Low-read NTCs may not be representative of the contamination profile
- **Over-subtraction** is a real risk, especially in low-biomass samples where contaminants and true signal overlap
- Statistical methods require **DNA concentration data** (Qubit/qPCR) for the frequency approach
- Contaminant taxa (e.g., *Bradyrhizobium*, *Burkholderia*, *Ralstonia*, *Propionibacterium*) are well-documented from kit reagents [3][9]

---

#### References

[1] Nature Microbiology (2025). *Guidelines for preventing and reporting contamination in low-biomass microbiome studies.* https://doi.org/10.1038/s41564-025-02035-2

[2] Walker et al. (2020). *Sequence-Based Characterization of Intratumoral Bacteria — A Guide to Best Practice.* Frontiers in Oncology. https://doi.org/10.3389/FONC.2020.00179

[3] Salter et al. (2014). *Reagent and laboratory contamination can critically impact sequence-based microbiome analyses.* BMC Biology. https://doi.org/10.1186/S12915-014-0087-Z

[4] Karstens et al. (2019). *Controlling for Contaminants in Low-Biomass 16S rRNA Gene Sequencing Experiments.* mSystems. https://doi.org/10.1128/MSYSTEMS.00290-19

[5] Davis et al. (2018). *Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data.* Microbiome. https://doi.org/10.1186/S40168-018-0605-2

[6] Dyrhovden et al. (2021). *Managing Contamination and Diverse Bacterial Loads in 16S rRNA Deep Sequencing of Clinical Samples.* mBio. https://doi.org/10.1128/MBIO.00598-21

[7] Zozaya-Valdés et al. (2021). *Detection of cell-free microbial DNA using a contaminant-controlled analysis framework.* Genome Biology. https://doi.org/10.1186/S13059-021-02401-3

[8] Lejal et al. (2020). *Taxon Appearance From Extraction and Amplification Steps Demonstrates the Value of Multiple Controls in Tick Microbiota Analysis.* Frontiers in Microbiology. https://doi.org/10.3389/FMICB.2020.01093

[9] Díaz et al. (2021). *Identification and Removal of Potential Contaminants in 16S rRNA Gene Sequence Data Sets from Low-Microbial-Biomass Samples.* mSphere. https://doi.org/10.1128/MSPHERE.00506-21

---
---

## 💬 Conversation 2

### 🧑 User Question
> The kits for extraction and PCR were the same throughout. All work was performed in the same environment. What's the level of risk we should be willing to take? For 7 out of 9 runs, elution buffer alone was used as NTC for PCR. There were two out of 9 runs that had 4 extra NTCs:
> - **ExDwEB** = extract Dw (distilled water) + elute EB (elution buffer)
> - **ExDW** = extract Dw + elute Dw
> - **NExDw** = no extract (use Dw for 16S seq PCR)
> - **NExEB** = no extract (use EB for 16S seq PCR)
>
> We have been using run-specific NTCs for subtraction but for these two runs, we found some potential contamination in the water. We are at a dilemma:
> 1. Should we keep the by-run NTCs or accumulative since it's the same kit?
> 2. If we were to keep it by-run, should we only use one NTC (EB) as background or all concatenated (5) NTC types?

---

### 🤖 Agent Response

**Search Summary:**
- **SciSpace Basic Search (run 1)**: 100 results — pooling NTCs, batch effects, run-specific controls
- **SciSpace Full Text Search (run 1)**: 100 results
- **Google Scholar (run 1)**: 20 results
- **PubMed (run 1)**: 18 results
- **SciSpace Basic Search (run 2)**: 100 results — extraction blanks vs. PCR blanks, water contamination
- **SciSpace Full Text Search (run 2)**: 100 results
- **Google Scholar (run 2)**: 20 results
- **PubMed (run 2)**: 9 results
- **Total unique papers (combined, deduplicated)**: ~172 papers across both searches
- **Zotero / Mendeley**: [Connect your account](https://www.scispace.com?settings=connected-apps)

---

## 🧫 NTC Strategy for Your 9-Run 16S Study: Literature-Based Assessment

---

### 📋 TL;DR
Run-specific NTCs are the gold standard, but your asymmetric design (7 runs with EB-only vs. 2 runs with full NTC panels) requires a **tiered, hybrid strategy**: use the richest available NTC set per run, leverage the full 5-NTC panel as a "contamination reference map" across all runs, and apply statistical classifiers rather than hard cutoffs.

---

### 🔑 Key Takeaways

- **Run-specific NTCs are preferred** because contamination profiles can vary by day, operator, and kit lot — even when the same kit is used throughout [1][2]
- **Accumulative/pooled NTCs are not wrong** but carry the risk of masking or over-calling run-specific signals; they are a reasonable fallback when per-run controls are sparse [5]
- **Extraction blanks capture different contamination sources** than PCR blanks — excluding them loses sensitivity to kit-derived DNA [6][7][8]
- **Distilled water contamination is a well-documented, real phenomenon** — your finding in 2 runs is not unusual [3][6]
- **Statistical classifiers (e.g., decontam) outperform hard read-count cutoffs** and should guide your risk threshold decisions [5]

---

### 1. Run-Specific vs. Accumulative NTCs — What Should You Do?

#### What the literature says
The literature is **strongly in favor of run-specific (per-batch) controls**. Salter et al. [3] and the 2025 Nature Microbiology guidelines [1] both emphasize that contamination is **temporally and batch-dependent** — reagent lots, environment, and operator behavior all shift over time. Even when the same kit is used, contamination profiles are **not guaranteed to be stable** across runs.

Pooling/accumulating NTCs across runs risks two failure modes:
- **Under-removal**: a run-specific contaminant is diluted out of the composite profile and not flagged
- **Over-removal**: a taxon that appeared transiently in one run's NTC gets applied to all runs, removing legitimate signal

#### Assessment of your accumulative NTC model
> *"All NTCs from all runs were accumulated to best capture transient contamination"*

This is a **reasonable and defensible rationale**, especially given that:
- The same kits and environment were used throughout ✅
- You are trying to capture **transient** contamination that might not appear in every run's NTC ✅
- Your 7 EB-only runs have limited per-run NTC resolution ⚠️

**However**, the literature cautions that the composite NTC approach is better used as a **supplementary contamination reference map** rather than the sole basis for subtraction. The recommended approach is:

| Approach | Recommended Use |
|---|---|
| **Run-specific NTC** | Primary subtraction — apply per-run to catch run-level signal |
| **Accumulative/composite NTC** | Secondary check — flag taxa that appear across ≥2 runs' NTCs as globally suspect |
| **Statistical classifier (decontam)** | Arbitrates borderline cases using DNA concentration + prevalence data |

---

### 2. For the 2 Full-Panel Runs: Use EB-Only or All 5 NTC Types?

#### What the literature says
Studies that included multiple control types consistently show that **extraction blanks and PCR blanks capture non-overlapping contamination sources** [6][7][8]:

| NTC Type | What It Captures |
|---|---|
| **ExDwEB / ExDW** (extraction blanks) | Kit reagent DNA (columns, buffers, beads), environmental DNA introduced during extraction |
| **NExDw / NExEB** (PCR blanks) | Amplification reagent contaminants, cross-contamination during PCR setup |
| **EB only** (your 7-run NTC) | Only PCR-stage contamination from elution buffer |

Lejal et al. [7][8] showed that **extraction-step controls detected taxa not present in PCR blanks**, and vice versa — demonstrating that **using only one type misses real contamination**. Claassen-Weitz et al. [6] (low-biomass nasopharyngeal samples) showed extraction blanks often dominated the contamination signal.

#### Recommendation for your 2 full-panel runs
✅ **Use all 5 NTC types concatenated/pooled** for those 2 runs. The rationale:
1. You detected **water contamination** — this is captured by ExDwEB and ExDW, which would be **invisible** to EB-only subtraction
2. The literature supports using the richest available control set [1][6][7]
3. You can still **harmonize across runs** by tracking which ASVs flagged in the 5-NTC runs also appear in the EB-only runs — those are your highest-confidence contaminants to remove globally

#### What to do about the asymmetry (5 NTCs vs. 1 NTC)?

```
Step 1: For runs 1–7 (EB only)     → subtract using EB NTC (run-specific)
Step 2: For runs 8–9 (5 NTC types) → subtract using all 5 NTCs pooled (run-specific)
Step 3: Cross-run harmonization     → flag ASVs that appear in BOTH EB-only NTCs 
                                      AND in any of the 5 NTC types across runs 8–9
                                      → these are your "high-confidence global contaminants"
Step 4: Apply decontam              → use prevalence method (controls vs. samples) 
                                      across all 9 runs to statistically validate
```

---

### 3. Distilled Water Contamination — Is This a Real Risk?

**Yes — this is well-documented.** Salter et al. [3] and Davis et al. [5] both note that molecular-grade and distilled water are known contamination sources. Your ExDwEB and ExDW controls are doing exactly what they should — **they revealed contamination that the EB-only NTC would have missed**. This is a strong argument for:
- Retaining the full 5-NTC subtraction for runs 8–9
- Treating taxa found in your water-based NTCs as suspect across **all 9 runs**, not just runs 8–9

---

### 4. What Level of Risk is Acceptable? Statistical Thresholds

The literature **warns against arbitrary read-count cutoffs** [1][5]. The recommended framework:

#### Tiered confidence for contaminant flagging

| Criterion | Weight |
|---|---|
| Higher **prevalence in NTCs** than in samples | Strong evidence |
| **Inverse correlation** with sample DNA concentration (decontam frequency method) | Strong evidence |
| Appears in **extraction blank** (kit-origin) | Moderate-strong |
| Appears in **PCR blank only** | Moderate |
| **Biologically implausible** in your sample type | Supporting evidence |
| Present in only **1 of 9 runs' NTCs** | Weak — flag but do not remove without corroboration |

#### Recommended tool: `decontam` [5]
- **Prevalence method**: marks ASVs more prevalent in controls than samples → use threshold `p = 0.1` (liberal, catches more) to `p = 0.5` (conservative) depending on your biomass level
- **Frequency method**: requires DNA concentration data per sample (qPCR/Qubit); models inverse relationship between ASV frequency and total DNA → more powerful for low-biomass samples
- Benchmarking studies show decontam removes **74–91% of known contaminants** [5]

#### Risk tolerance guidance
- **Low-biomass samples** (e.g., tumor tissue, plasma, low-input clinical): use **liberal threshold** (p ≤ 0.5 in decontam prevalence), accept higher false-positive removal to protect signal integrity
- **High-biomass samples** (e.g., gut stool): use **conservative threshold** (p ≤ 0.1), prioritize retaining true signal

---

### 5. Summary Recommendation for Your Study

| Decision | Recommended Action | Rationale |
|---|---|---|
| Run-specific vs. accumulative | **Hybrid**: run-specific primary, accumulative as secondary reference | Literature favors per-run; accumulative catches transient signals |
| Runs 1–7 NTC | Use **EB NTC per run** for subtraction | Only available control; document limitation |
| Runs 8–9 NTC | Use **all 5 NTCs concatenated** per run | Maximizes contamination detection; water contamination confirmed |
| Cross-run harmonization | Flag ASVs in ≥2 runs' NTCs as **globally suspect** | Convergent evidence across runs = higher confidence |
| Statistical method | Apply **decontam** (prevalence + frequency if DNA conc. available) | Principled, validated, widely accepted |
| Water contamination finding | Treat water-detected taxa as suspect in **all 9 runs** | Same water source likely used across all runs |
| Reporting | **Fully disclose** the asymmetric NTC design and both strategies tested | Required for reproducibility; reviewers will ask |

---

### References

[1] Nature Microbiology (2025). *Guidelines for preventing and reporting contamination in low-biomass microbiome studies.* https://doi.org/10.1038/s41564-025-02035-2

[2] Aggarwal et al. (2023). *Optimization of high-throughput 16S rRNA gene amplicon sequencing.* https://doi.org/10.6084/m9.figshare.24188649

[3] Salter et al. (2014). *Reagent and laboratory contamination can critically impact sequence-based microbiome analyses.* BMC Biology. https://doi.org/10.1186/S12915-014-0087-Z

[4] Zozaya-Valdés et al. (2021). *Detection of cell-free microbial DNA using a contaminant-controlled analysis framework.* Genome Biology. https://doi.org/10.1186/S13059-021-02401-3

[5] Davis et al. (2018). *Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data.* Microbiome. https://doi.org/10.1186/S40168-018-0605-2

[6] Claassen-Weitz et al. (2020). *Optimizing 16S rRNA gene profile analysis from low biomass nasopharyngeal and induced sputum specimens.* BMC Microbiology. https://doi.org/10.1186/S12866-020-01795-7

[7] Lejal et al. (2019). *Taxon appearance from extraction and amplification steps demonstrates the value of multiple controls in tick microbiome analysis.* bioRxiv. https://doi.org/10.1101/714030

[8] Lejal et al. (2020). *Taxon Appearance From Extraction and Amplification Steps Demonstrates the Value of Multiple Controls in Tick Microbiota Analysis.* Frontiers in Microbiology. https://doi.org/10.3389/FMICB.2020.01093

[9] Díaz et al. (2021). *Identification and Removal of Potential Contaminants in 16S rRNA Gene Sequence Data Sets from Low-Microbial-Biomass Samples.* mSphere. https://doi.org/10.1128/MSPHERE.00506-21

---

*Export generated by SciSpace Research Agent — 2026-05-08*
