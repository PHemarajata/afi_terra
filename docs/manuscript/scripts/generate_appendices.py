#!/usr/bin/env python3
"""
Generate two manuscript appendices as markdown:

  APPENDIX-VALIDATION-PANEL.md
    - one row per (sample, detected genus, source) for the 43-sample
      validation panel
    - includes expected organism, concordance call, failure-mode tag,
      Minimap2 rescue evidence even when rescue failed

  APPENDIX-STUDY-SAMPLES.md
    - one row per (sample, detected genus, source) for the 56 AFI study
      samples
    - includes pipeline NCmax, same-run NTC cross-check, V4 filter tier
      decision, Burkholderia species-level evidence, confidence, and
      hypothesis-class tag

Reads:
  - .calls.tsv files (one per sample, per run, per shard)
  - .centrifuger.kreport.tsv files (Burkholderia species-level)
"""

import glob
import os
from collections import defaultdict

import pandas as pd

ROOT = '/Users/peerahemarajata/Downloads/AFI_P_Final'
OUT_DIR = '/Users/peerahemarajata/.claude/projects/-Users-peerahemarajata-afi-terra/session-data'

RUN_DIRS = ['validation_panel', '1_and_2', '3', '4_and_5', '6_and_7', '8_and_9']

# Expected organisms in the validation panel.
# Built from VALIDATION-PANEL-REPORT.md content.
EXPECTED_VALIDATION = {
    # Escherichia coli
    '09302710_S2_L001': 'Escherichia coli',
    '09602983_S7_L001': 'Escherichia coli',
    '21-0-01573': 'Escherichia coli',
    '23-0-00493': 'Escherichia coli',
    '23500805_S7_L001': 'Escherichia coli',
    'E-coli_S4_L001': 'Escherichia coli (PC_SINGLE)',

    # Orientia tsutsugamushi
    '00389_S5_L001': 'Orientia tsutsugamushi',
    '11800801_S7_L001': 'Orientia tsutsugamushi',
    '22900253_S4_L001': 'Orientia tsutsugamushi',
    '24500367_S5_L001': 'Orientia tsutsugamushi',
    '25900911_S4_L001': 'Orientia tsutsugamushi',
    '00618_S7_L001': 'Orientia tsutsugamushi',

    # Rickettsia
    '10900410_S10_L001': 'Rickettsia spp',
    '16401070_S11_L001': 'Rickettsia spp',
    '00126_S6_L001': 'Rickettsia spp',
    '00369_S1_L001': 'Rickettsia spp',

    # Leptospira
    '00277_S4_L001': 'Leptospira spp',
    '09202659_S12_L001': 'Leptospira spp',
    '00428_S8_L001': 'Leptospira spp',
    '02757_S3_L001': 'Leptospira spp',

    # Burkholderia pseudomallei
    '09-0-02165': 'Burkholderia pseudomallei',
    '09502813_S2_L001': 'Burkholderia pseudomallei',
    '09700912_S3_L001': 'Burkholderia pseudomallei',
    '09102966_S1_L001': 'Burkholderia pseudomallei',
    '01937_S2_L001': 'Burkholderia pseudomallei',

    # Streptococcus pneumoniae
    '09302890_S5_L001': 'Streptococcus pneumoniae',
    '09902043_S1_L001': 'Streptococcus pneumoniae',
    '23500382_S6_L001': 'Streptococcus pneumoniae',

    # Streptococcus suis
    '09402516_S8_L001': 'Streptococcus suis',
    '09300546_S4_L001': 'Streptococcus suis',
    '22900088_S3_L001': 'Streptococcus suis',

    # Coxiella burnetii
    '25300179_S6_L001': 'Coxiella burnetii',
    '25800370_S9_L001': 'Coxiella burnetii',

    # Yersinia
    '11200465_S2_L001': 'Yersinia spp',

    # Controls
    'PC-20251016_S7_L001': 'PC_MIX8 (8 organisms)',
    'PC_S10': 'PC_MIX8 (8 organisms)',
    'PC_S12': 'PC_MIX8 (8 organisms)',
    'PC_S13_L001': 'PC_MIX8 (8 organisms)',
    'PC_S8_L001': 'PC_MIX8 (8 organisms)',
    'P-aeru_S5_L001': 'Pseudomonas aeruginosa (PC_SINGLE)',
    'S-pneumo_S2_L001': 'Streptococcus pneumoniae (PC_SINGLE)',
    'S-suis_S3_L001': 'Streptococcus suis (PC_SINGLE)',
    'Mixed_S6_L001': 'MIXED4 (E.coli, P.aer, S.pneumo, S.suis)',
    'NC-20251016_S1_L001': 'NTC (no organisms)',
    'NC_S1': 'NTC (no organisms)',
    'NC_S14': 'NTC (no organisms)',
    'NC_S9': 'NTC (no organisms)',
    'NTC_S11': 'NTC (no organisms)',
}

# V4 filter contaminant tiers (kept consistent with afi_decontamination_filter_v4.py)
TIER_A = {
    'Pseudomonas', 'Ralstonia', 'Bradyrhizobium', 'Sphingomonas',
    'Stenotrophomonas', 'Methylobacterium', 'Acinetobacter',
    'Cutibacterium', 'Staphylococcus', 'Corynebacterium', 'Brevundimonas',
}
NC_ONLY = {
    'Cereibacter', 'Thioclava', 'Bdellovibrio', 'Saltatorellus',
    'Pseudogemmobacter', 'Minisyncoccus', 'Rhodoluna', 'Microbacterium',
    'Arcanobacterium',
}
TIER_1 = {
    'Shigella', 'Metapseudomonas', 'Stutzerimonas', 'Capsulimonas',
    'Chamaesiphon', 'Chloroflexus', 'Flavihumibacter', 'Hymenobacter',
    'Limnoglobus', 'Methylovirgula', 'Microvirga',
    'Pelagovum', 'Pseudonocardia', 'Rufibacter', 'Salmonella', 'Spirosoma',
}
TIER_2 = {
    'Klebsiella', 'Asticcacaulis', 'Methylobacterium', 'Xanthomonas',
    'Caulobacter', 'Fimbriimonas', 'Gemmatirosa', 'Roseateles',
    'Actinomycetospora', 'Algoriphagus', 'Chloroherpeton', 'Delftia',
    'Rhodococcus', 'Coxiella', 'Novosphingobium', 'Sphingopyxis',
    'Brasilonema', 'Candidatus_Amoebophilus', 'Chroococcidiopsis',
    'Dermatobacter', 'Erythrobacter', 'Leptolyngbya', 'Nostoc',
    'Pseudoluteimonas', 'Roseomonas', 'Rubrobacter', 'Variovorax',
}
B_PSEUDOMALLEI_MIN_SPECIES_READS = 500

# Hypothesis-class tagging for the AFI interpretation
HYPOTHESIS_CLASS = {
    # Fastidious / slow-growing / cell-wall-deficient organisms of interest
    'Mycoplasmopsis': 'fastidious / cell-wall-deficient',
    'Leptospira': 'fastidious / slow-growing',
    'Brucella': 'fastidious / slow-growing',
    # Known culturable AFI pathogens
    'Escherichia': 'culturable AFI pathogen',
    'Klebsiella': 'culturable AFI pathogen',
    'Enterobacter': 'culturable AFI pathogen',
    # Genus-level ambiguous
    'Streptococcus': 'fastidious-or-not (V1-V3 species-ambiguous)',
    'Burkholderia': 'see species-level evidence (B. pseudomallei vs cepacia complex)',
    # Anaerobic-ish
    'Porphyromonas': 'anaerobic',
    'Desulfovibrio': 'anaerobic',
    # Rickettsiales
    'Orientia': 'Rickettsiales (cannot-miss)',
    'Rickettsia': 'Rickettsiales (cannot-miss)',
    # Generally noise / environmental in this context
    'Thermomicrobium': 'environmental (likely noise)',
    'Brevundimonas': 'kit / water contaminant',
    'Faucicola': 'oral flora / contamination',
}


def get_run_for_sample(sample):
    for run in RUN_DIRS:
        hits = (
            glob.glob(f"{ROOT}/{run}/call-P2_Interpret/*/cacheCopy/{sample}.calls.tsv")
            + glob.glob(f"{ROOT}/{run}/call-P2_Interpret/*/{sample}.calls.tsv")
        )
        if hits:
            return run
    return None


def _find_calls_tsvs(run):
    """validation_panel stores calls.tsv directly under shard-N/; other runs
    use shard-N/cacheCopy/. Handle both."""
    patterns = [
        f"{ROOT}/{run}/call-P2_Interpret/*/cacheCopy/*.calls.tsv",
        f"{ROOT}/{run}/call-P2_Interpret/*/*.calls.tsv",
    ]
    found = []
    for p in patterns:
        found.extend(glob.glob(p))
    # De-duplicate (some runs may match both patterns if they have files at both levels)
    return sorted(set(found))


def load_calls(run):
    """Return a list of dicts: one per row of every calls.tsv in the run."""
    out = []
    files = _find_calls_tsvs(run)
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t')
        except Exception:
            continue
        if len(df) == 0:
            continue
        for _, row in df.iterrows():
            out.append({
                'run': run,
                'sample': str(row.get('sample', '')).strip(),
                'genus': str(row.get('genus', '')).strip(),
                'source': str(row.get('source', '')).strip(),
                'reads': int(row['reads']) if pd.notna(row.get('reads')) else 0,
                'breadth': float(row['breadth']) if pd.notna(row.get('breadth')) and str(row.get('breadth','')) != '' else None,
                'ntc_reads': int(row['ntc_reads']) if pd.notna(row.get('ntc_reads')) else 0,
                'call': str(row.get('call', '')).strip(),
                'cfr_reads': int(row['cfr_reads']) if pd.notna(row.get('cfr_reads')) else 0,
                'rescued': str(row.get('rescued', '')).strip().lower() == 'true',
                'align_confirmed': str(row.get('align_confirmed', '')).strip().lower() == 'true',
            })
    return out


POSITIVE_CALLS = ('Detected', 'Confirmed', 'Probable')

# Genus-level mapping of CDC AFI TaqMan Array Card (TAC) BACTERIAL targets.
# 16S V1-V3 detects bacteria only, so viral / protozoal TAC targets are not
# evaluable by this assay and are excluded here.
TAC_BACTERIAL_GENERA = {
    'Bartonella',   # TAC_BART
    'Brucella',     # TAC_BRUC
    'Rickettsia',   # TAC_RICK (also covered by Orientia)
    'Orientia',     # TAC_ORTS
    'Yersinia',     # TAC_YERS
    'Coxiella',     # TAC_COBU
    'Streptococcus',  # TAC_STSU, TAC_STPN
    'Salmonella',   # TAC_SATY, TAC_SALS
    'Escherichia',  # TAC_ESCH
    'Burkholderia',  # TAC_BUPS (species-level safeguard applied separately)
}

# Expected organisms per positive-control sample. The V4 contaminant filter
# is bypassed for these organisms in their respective PC samples (a spike-in
# organism is not a contaminant in its own control).
#
# PC_MIX8 is the ZymoBIOMICS Microbial Community Standard (8 bacterial organisms):
#   B. subtilis, E. faecalis, E. coli, L. fermentum (now Limosilactobacillus),
#   L. monocytogenes, P. aeruginosa, S. enterica, S. aureus.
PC_MIX8_GENERA = {
    'Bacillus', 'Enterococcus', 'Escherichia', 'Limosilactobacillus',
    'Listeria', 'Pseudomonas', 'Salmonella', 'Staphylococcus',
}
PC_EXPECTED = {
    # PC_SINGLE (one organism per sample)
    'E-coli_S4_L001':   {'Escherichia'},
    'P-aeru_S5_L001':   {'Pseudomonas'},
    'S-pneumo_S2_L001': {'Streptococcus'},
    'S-suis_S3_L001':   {'Streptococcus'},
    # PC_MIX8 replicates (ZymoBIOMICS Microbial Community Standard)
    'PC-20251016_S7_L001': PC_MIX8_GENERA,
    'PC_S8_L001':          PC_MIX8_GENERA,
    'PC_S10':              PC_MIX8_GENERA,
    'PC_S12':              PC_MIX8_GENERA,
    'PC_S13_L001':         PC_MIX8_GENERA,
    # MIXED4: E. coli + P. aeruginosa + S. pneumoniae + S. suis = 3 genera at 16S resolution
    'Mixed_S6_L001':       {'Escherichia', 'Pseudomonas', 'Streptococcus'},
}


def total_reads_by_sample(rows):
    """Sum the positively-called reads per sample.

    Strategy: take the union of (sample, genus) pairs from centrifuge and
    alignment positive calls. For each pair, use the max of centrifuge and
    alignment reads to avoid double-counting overlapping evidence. This
    correctly handles samples where Rickettsiales is detected only by
    Minimap2 alignment (e.g. 00618_S7_L001).
    """
    per_sample_genus = defaultdict(lambda: defaultdict(int))
    for r in rows:
        if r['call'] in POSITIVE_CALLS:
            key = (r['sample'], r['genus'])
            if r['reads'] > per_sample_genus[r['sample']][r['genus']]:
                per_sample_genus[r['sample']][r['genus']] = r['reads']
    totals = defaultdict(int)
    for sample, genera in per_sample_genus.items():
        totals[sample] = sum(genera.values())
    return totals


def parse_kreport_species_reads(kreport_file, species_name):
    if not os.path.exists(kreport_file):
        return 0
    total = 0
    try:
        with open(kreport_file, 'r') as f:
            for raw in f:
                cols = raw.rstrip('\n').split('\t')
                if len(cols) < 6 or cols[3].strip() != 'S':
                    continue
                if cols[5].strip() == species_name:
                    try:
                        total += int(cols[1])
                    except ValueError:
                        continue
    except Exception:
        return 0
    return total


def find_kreport(run, sample):
    hits = glob.glob(f"{ROOT}/{run}/call-P1_Centrifuger/*/cacheCopy/{sample}.centrifuger.kreport.tsv")
    return hits[0] if hits else None


def find_run_ntcs(run):
    """Return list of (ntc_name, calls.tsv path, kreport path) per run."""
    out = []
    calls = _find_calls_tsvs(run)
    for c in calls:
        name = os.path.basename(c).replace('.calls.tsv', '')
        if name.startswith('NTC') or name.startswith('NC_') or name.startswith('NC-'):
            kreport = find_kreport(run, name)
            out.append((name, c, kreport))
    return out


def same_run_ntc_max(run, genus):
    """Return the maximum reads of `genus` across the run's NTC calls.tsv."""
    out = 0
    for ntc_name, calls_path, _ in find_run_ntcs(run):
        try:
            df = pd.read_csv(calls_path, sep='\t')
        except Exception:
            continue
        for _, row in df.iterrows():
            if str(row.get('genus','')).strip() == genus and str(row.get('source','')).strip() != 'alignment':
                rds = int(row['reads']) if pd.notna(row.get('reads')) else 0
                if rds > out:
                    out = rds
    return out


def same_run_ntc_burkpseudomallei_max(run):
    out = 0
    for _, _, kreport in find_run_ntcs(run):
        if not kreport:
            continue
        rds = parse_kreport_species_reads(kreport, 'Burkholderia_pseudomallei')
        if rds > out:
            out = rds
    return out


def decide_filter(genus, reads, pct, source, sample, run, total_reads, total_genus_detections_in_cohort, expected_pc_set=None):
    """Replicate V4 filter decision for a single (sample, genus, source) row.

    expected_pc_set: optional set of genera that are expected to be present
    in this sample as positive-control spike-ins. A spike-in genus is NOT a
    contaminant in its own control, so it bypasses Tier A removal (and the
    Burkholderia species safeguard) when listed here.
    """
    expected_pc_set = expected_pc_set or set()
    if genus in expected_pc_set:
        return 'KEEP_PC_EXPECTED', 'Expected positive-control organism (filter bypassed)'
    if genus in TIER_A:
        return 'REMOVE_TIER_A', f'Tier A high-confidence contaminant'
    if genus == 'Burkholderia':
        kreport = find_kreport(run, sample) if run else None
        species_reads = parse_kreport_species_reads(kreport, 'Burkholderia_pseudomallei') if kreport else 0
        ntc_max = same_run_ntc_burkpseudomallei_max(run) if run else 0
        if species_reads >= B_PSEUDOMALLEI_MIN_SPECIES_READS and species_reads > ntc_max:
            return 'KEEP_BURK_PSEUDOMALLEI', f'B. pseudomallei species: {species_reads} reads (>= {B_PSEUDOMALLEI_MIN_SPECIES_READS}; > NTC {ntc_max})'
        else:
            return 'REMOVE_BURK', f'B. pseudomallei species reads {species_reads} below threshold or below NTC max {ntc_max}'
    if genus in NC_ONLY:
        return 'REMOVE_NC_ONLY', 'NC-only organism'
    if genus in TIER_1:
        return 'REMOVE_TIER_1', 'Ultra-low abundance (median <0.5% in dataset)'
    if genus in TIER_2 and pct < 1.0:
        if total_genus_detections_in_cohort <= 3:
            return 'REMOVE_TIER_2', f'Marginal abundance (<1%) in rare detection ({total_genus_detections_in_cohort} samples)'
    return 'KEEP', 'Passes filter'


def confidence_label(pct, sample_reads, same_run_ntc_max_val, genus, is_pc_or_ntc):
    """Confidence assignment based on abundance and NTC headroom.

    The strongest signal-killer is when the run's NTC carries MORE reads of
    the same organism than the sample itself does. This is treated as a
    severe red flag regardless of abundance.
    """
    if is_pc_or_ntc:
        return 'CONTROL'
    # Worst case: NTC carries more reads than sample (e.g. Leptospira in 09801652_S5_L001
    # where sample = 7,294 but same-run NTC = 78,691).
    if same_run_ntc_max_val > sample_reads and sample_reads > 0:
        return 'CANDIDATE (run NTC carries more reads than this sample)'
    # NTC has meaningful headroom-eating background
    if same_run_ntc_max_val > 0 and sample_reads > 0:
        ratio = same_run_ntc_max_val / sample_reads
        if ratio >= 0.5:
            return 'LOW (run NTC carries comparable signal)'
        if ratio >= 0.1:
            return 'MODERATE (NTC headroom modest)'
    if pct >= 30:
        return 'HIGH'
    if pct >= 5:
        return 'MODERATE'
    if pct >= 1:
        return 'LOW'
    return 'CANDIDATE (near noise floor)'


def hypothesis_tag(genus):
    return HYPOTHESIS_CLASS.get(genus, '-')


def main():
    # --- Load everything ----------------------------------------------------
    all_rows = []
    for run in RUN_DIRS:
        all_rows.extend(load_calls(run))

    # Cohort-level genus-detection counts (used for Tier 2 rare check)
    detected = [r for r in all_rows if r['call'] in POSITIVE_CALLS and r['source'] != 'alignment']
    sample_is_study = lambda r: not (r['sample'].startswith(('PC', 'NTC', 'NC_', 'NC-')) or r['run'] == 'validation_panel')
    study_detected = [r for r in detected if sample_is_study(r)]
    study_genus_counts = defaultdict(int)
    for r in study_detected:
        study_genus_counts[r['genus']] += 1

    # per-sample totals from centrifuge rows (matches V4 filter logic)
    sample_totals = total_reads_by_sample(all_rows)

    # =======================================================================
    # APPENDIX 1: VALIDATION PANEL
    # =======================================================================
    val_rows_out = []
    val_calls = [r for r in all_rows if r['run'] == 'validation_panel']
    val_samples = sorted({r['sample'] for r in val_calls})

    # build per-sample list of detections (all sources)
    by_sample = defaultdict(list)
    for r in val_calls:
        if r['call'] in POSITIVE_CALLS:
            by_sample[r['sample']].append(r)
        elif r['call'] == 'Negative' and r['reads'] >= 100:
            # Preserve failed-rescue alignment rows for transparency
            by_sample[r['sample']].append(r)

    for sample in val_samples:
        expected = EXPECTED_VALIDATION.get(sample, '(unknown - check manifest)')
        total_reads = sample_totals.get(sample, 0)
        sample_dets = by_sample.get(sample, [])
        is_pc = (
            sample.startswith(('PC', 'Mixed'))
            or '(PC_SINGLE)' in expected
            or 'MIXED4' in expected
        )
        is_ntc = sample.startswith(('NC_', 'NC-', 'NTC'))

        if not sample_dets:
            # No detections at all - flag as pre-sequencing failure or NTC
            val_rows_out.append({
                'sample': sample,
                'run': 'validation_panel',
                'category': 'NTC' if is_ntc else ('control' if is_pc else 'clinical'),
                'expected': expected,
                'total_reads': total_reads,
                'genus': '-',
                'source': '-',
                'reads': '-',
                'pct': '-',
                'ntc_reads': '-',
                'rescue_info': '-',
                'rescue_tier': '-',
                'filter_outcome': 'n/a',
                'concordance': 'Concordant' if is_ntc else ('Discordant' if expected != '(unknown - check manifest)' else 'n/a'),
                'failure_mode': '-' if is_ntc else 'pre-sequencing or expected-empty',
                'notes': 'No taxa detected' if not is_ntc else 'NTC expected empty - PASS',
            })
            continue

        for det in sample_dets:
            pct = (100.0 * det['reads'] / total_reads) if total_reads else 0.0
            # Rescue info: only for alignment source or when rescued / centrifuge for Rickettsiales
            rescue_info = '-'
            rescue_tier = '-'
            if det['source'] == 'alignment':
                rescue_info = f'mapped={det["reads"]}, breadth={det["breadth"]:.4f}, ntc={det["ntc_reads"]}'
                if det['call'] == 'Confirmed':
                    rescue_tier = 'Tier 1: genus-level (align_confirmed)'
                elif det['call'] == 'Probable':
                    rescue_tier = 'Tier 2: order-level (Rickettsiales detected)'
                elif det['call'] == 'Detected':
                    rescue_tier = 'genus-level (rescued)' if det['rescued'] else 'genus-level'
                elif det['call'] == 'Negative':
                    if det['breadth'] is not None and det['breadth'] < 0.25:
                        rescue_tier = 'FAIL: breadth < 0.25'
                    elif det['reads'] < 100:
                        rescue_tier = 'FAIL: reads < 100'
                    else:
                        rescue_tier = 'FAIL'
            # Filter outcome
            same_run_ntc = same_run_ntc_max('validation_panel', det['genus']) if det['call'] in POSITIVE_CALLS else 0
            if det['call'] not in POSITIVE_CALLS:
                # Failed-rescue alignment rows (call='Negative') are not subject
                # to V4 filtering. They were rejected by the pipeline upstream.
                filter_decision = 'n/a (pipeline rejected pre-filter)'
                filter_reason = 'Alignment-based rescue did not pass pipeline thresholds'
            else:
                filter_decision, filter_reason = decide_filter(
                    det['genus'], det['reads'], pct, det['source'], sample, 'validation_panel',
                    total_reads, study_genus_counts.get(det['genus'], 0),
                    expected_pc_set=PC_EXPECTED.get(sample, set()),
                )
            # Concordance
            if det['call'] not in POSITIVE_CALLS:
                concordance = 'n/a'
            elif is_ntc:
                concordance = 'Discordant (NTC should be empty)'
            elif is_pc:
                # PC samples: row-level concordance against the full expected
                # spike-in set (PC_EXPECTED dict), not just a single organism.
                pc_expected = PC_EXPECTED.get(sample, set())
                if det['genus'] in pc_expected:
                    concordance = 'PC TARGET (expected spike-in)'
                else:
                    concordance = 'PC co-detection (non-expected organism)'
            else:
                exp_genus = expected.split()[0]
                rickettsiales = {'Orientia', 'Rickettsia'}
                if det['genus'].lower().startswith(exp_genus.lower()):
                    concordance = 'TARGET match'
                elif exp_genus in rickettsiales and det['genus'] in rickettsiales:
                    # Within Rickettsiales, cross-genus rescue counts as
                    # clinically concordant (same doxycycline-treatable order).
                    if det['call'] == 'Probable':
                        concordance = 'TARGET match (order-level Rickettsiales)'
                    else:
                        concordance = 'TARGET match (cross-genus Rickettsiales)'
                else:
                    concordance = 'Co-detection (non-target organism)'

            # Failure mode (only meaningful where expected target was NOT recovered)
            failure_mode = '-'
            if concordance == 'Co-detection (non-target organism)' and not is_pc and not is_ntc:
                # This row represents a co-detection, not the target. Failure
                # mode tags should be reserved for sample-level interpretation
                # in the prose; per-row this is just "not the target".
                failure_mode = '-'
            if det['source'] == 'alignment' and det['call'] == 'Negative' and not is_pc and not is_ntc:
                failure_mode = 'Minimap2 rescue threshold failure'

            note = ''
            if 'Streptococcus' in det['genus'] and ('pneumoniae' in expected or 'suis' in expected):
                note = 'V1-V3 cannot distinguish S. pneumoniae from S. suis'
            if det['source'] == 'alignment' and det['call'] == 'Negative':
                note = 'Rickettsiales Minimap2 rescue attempted but failed thresholds'

            val_rows_out.append({
                'sample': sample,
                'run': 'validation_panel',
                'category': 'NTC' if is_ntc else ('control' if is_pc else 'clinical'),
                'expected': expected,
                'total_reads': total_reads,
                'genus': det['genus'],
                'source': det['source'],
                'reads': det['reads'],
                'pct': f'{pct:.2f}%',
                'ntc_reads': det['ntc_reads'],
                'rescue_info': rescue_info,
                'rescue_tier': rescue_tier,
                'filter_outcome': filter_decision,
                'concordance': concordance,
                'failure_mode': failure_mode,
                'notes': note,
            })

    # Compute per-sample biomass distribution (kept vs removed by V4 filter)
    # This sums reads across all rows that contributed to the sample, broken
    # down by filter outcome. Light double-counting can occur when a sample
    # has both centrifuge and alignment rows for the same genus; the % kept
    # ratio remains interpretable.
    def biomass_split(rows):
        kept = defaultdict(int)
        removed = defaultdict(int)
        for r in rows:
            if r.get('reads') in ('-', None) or r.get('reads') == '':
                continue
            try:
                rd = int(r['reads'])
            except (ValueError, TypeError):
                continue
            out = r.get('filter_outcome', '')
            sample = r['sample']
            if out == 'n/a':
                continue
            if out.startswith('KEEP'):
                kept[sample] += rd
            elif out.startswith('REMOVE'):
                removed[sample] += rd
        return kept, removed

    val_kept, val_removed = biomass_split(val_rows_out)
    study_kept_reads = defaultdict(int)
    study_removed_reads = defaultdict(int)
    # study split is computed after study_rows_out exists; placeholder for now

    def annotate_biomass(rows, kept_map, removed_map):
        for r in rows:
            sample = r['sample']
            k = kept_map.get(sample, 0)
            rem = removed_map.get(sample, 0)
            total = k + rem
            r['biomass_kept'] = k if total else '-'
            r['biomass_removed'] = rem if total else '-'
            r['pct_biomass_kept'] = f'{(100.0 * k / total):.1f}%' if total else '-'

    annotate_biomass(val_rows_out, val_kept, val_removed)

    # ----- Final concordance (TAC-aware, V4-filter-aware) ------------------
    def compute_final_concordance(rows_by_sample):
        """Return a dict {sample: final_concordance_string} for validation samples.

        Rules:
        - NTC: concordant if no TAC bacterial target genus is detected AND retained
          by the V4 filter; discordant otherwise (list the offending genera).
        - PC controls: report which TAC bacterial genera are retained vs not.
          Marked as 'Control' for manual review against the full PC organism mix.
        - Clinical sample (expected single organism):
            * Concordant if a TARGET match row exists AND is retained by V4
              (this includes Rickettsiales order-level / cross-genus rescue and
              the B. pseudomallei species-level safeguard).
            * Discordant (filter-removed) if target detected but removed by V4.
            * Discordant (target not detected) if no TARGET match row exists.
            * Discordant (pre-sequencing) if sample has zero total reads.
        """
        out = {}
        for sample, rows in rows_by_sample.items():
            if not rows:
                continue
            category = rows[0].get('category')
            expected = rows[0].get('expected', '')
            is_ntc = category == 'NTC'
            is_pc = category == 'control'

            # Retained TAC bacterial detections
            tac_retained = [
                r for r in rows
                if r.get('genus') in TAC_BACTERIAL_GENERA
                and str(r.get('filter_outcome', '')).startswith('KEEP')
            ]

            if is_ntc:
                if tac_retained:
                    organisms = ', '.join(sorted({r['genus'] for r in tac_retained}))
                    out[sample] = f'Discordant (NTC contains TAC target genus: {organisms})'
                else:
                    out[sample] = 'Concordant (NTC clear of TAC target genera)'
                continue

            if is_pc:
                # Positive controls: concordant iff all expected spike-in
                # organisms are detected AND retained by V4. PCs count
                # towards validation panel accuracy.
                expected_set = PC_EXPECTED.get(sample, set())
                if not expected_set:
                    out[sample] = 'Control (no expected organisms defined)'
                    continue
                retained_genera = {
                    r['genus'] for r in rows
                    if str(r.get('filter_outcome', '')).startswith('KEEP')
                }
                detected_genera = {r['genus'] for r in rows if r.get('genus') not in (None, '-', '')}
                missing_retained = expected_set - retained_genera
                missing_detected = expected_set - detected_genera
                if not missing_retained:
                    out[sample] = (
                        f'Concordant (PC: all {len(expected_set)}/{len(expected_set)} '
                        f'expected organisms detected and retained)'
                    )
                elif not missing_detected:
                    # detected but filter removed it
                    out[sample] = (
                        f'Discordant (PC: detected but removed by V4 filter — {", ".join(sorted(missing_retained))})'
                    )
                else:
                    out[sample] = (
                        f'Discordant (PC: missing expected organism(s): {", ".join(sorted(missing_detected))})'
                    )
                continue

            # Clinical samples
            # Detect target match rows (per-row Concordance already classifies)
            target_rows = [r for r in rows if str(r.get('concordance', '')).startswith('TARGET match')]
            target_retained = [r for r in target_rows if str(r.get('filter_outcome', '')).startswith('KEEP')]
            target_removed = [r for r in target_rows if str(r.get('filter_outcome', '')).startswith('REMOVE')]

            total_reads = rows[0].get('total_reads', 0)
            if target_retained:
                out[sample] = 'Concordant (target detected and retained by V4 filter)'
            elif target_removed:
                # Target organism was matched but V4 filter removed it
                removed_genera = ', '.join(sorted({r['genus'] for r in target_removed}))
                out[sample] = f'Discordant (target matched but removed by V4: {removed_genera})'
            elif total_reads == 0:
                out[sample] = 'Discordant (no detections - pre-sequencing failure suspected)'
            else:
                out[sample] = 'Discordant (target organism not detected)'
        return out

    val_by_sample = defaultdict(list)
    for r in val_rows_out:
        val_by_sample[r['sample']].append(r)
    val_final_concordance = compute_final_concordance(val_by_sample)
    for r in val_rows_out:
        r['final_concordance'] = val_final_concordance.get(r['sample'], '-')

    # Write the validation appendix
    val_md = ['# Appendix: Validation Panel Detection Table',
              '',
              '> **Format:** one row per (sample × detected genus × source). Rows where the row corresponds to a Minimap2 alignment-based rescue are marked `source = alignment`. Negative alignment rows are included where they document a rescue attempt that failed (reads >= 100 but breadth < 0.25, etc.) -- these are deliberately preserved to show why rescue did not trigger.',
              '> ',
              '> **Filter outcome** is the V4 decontamination filter\'s decision applied to this row. The validation panel rows uniformly show `KEEP` because validation samples do not contain Tier A contaminants at detectable levels; this is recorded for transparency, not as filter validation.',
              '',
              '| Sample | Run | Category | Expected | Total reads | Biomass kept (V4) | Biomass removed (V4) | % biomass kept | Detected genus | Source | Reads | % of sample | Pipeline NTC reads | Minimap2 rescue info | Rescue tier | V4 filter | Per-row concordance | Failure mode | Final concordance (TAC + V4) | Notes |',
              '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for row in val_rows_out:
        val_md.append('| ' + ' | '.join(str(row[c]) for c in [
            'sample','run','category','expected','total_reads','biomass_kept','biomass_removed','pct_biomass_kept',
            'genus','source','reads','pct','ntc_reads',
            'rescue_info','rescue_tier','filter_outcome','concordance','failure_mode','final_concordance','notes',
        ]) + ' |')
    val_md.append('')
    val_md.append('---')
    val_md.append('')
    val_md.append('## Concordance Vocabulary')
    val_md.append('')
    val_md.append('- **TARGET match:** detected genus matches the expected organism at genus level.')
    val_md.append('- **TARGET match (order-level Rickettsiales):** call=Probable in calls.tsv; the rescued evidence supports Rickettsiales presence at order level but V1-V3 sequence variation does not allow genus discrimination. Counts as clinically concordant in an endemic context.')
    val_md.append('- **TARGET match (cross-genus Rickettsiales):** Orientia-expected sample where the rescue produced a Rickettsia call (or vice versa). Same doxycycline-treatable order, counted as concordant.')
    val_md.append('- **Co-detection (non-target organism):** organism present in the sample that is NOT the expected target. These rows are informational, not failures per se. Whether a sample is concordant overall depends on whether at least one TARGET match row exists for that sample.')
    val_md.append('- **Discordant (no target detected):** the sample has detections but none of them match the expected organism. Per-row concordance is "Co-detection"; sample-level concordance is Discordant.')
    val_md.append('- **n/a (no detections at all):** sample showed zero taxa; consistent with DNA extraction / library prep / sequencing-depth failure rather than classifier error.')
    val_md.append('')
    val_md.append('## Failure Mode Taxonomy (sample-level, for samples with no TARGET match)')
    val_md.append('')
    val_md.append('- **pre-sequencing or expected-empty:** sample showed zero taxa. Examples in this panel: `00126_S6_L001`, `00369_S1_L001` (both expected Rickettsia), `25800370_S9_L001` (expected Coxiella).')
    val_md.append('- **V1-V3 species ambiguity:** the V1-V3 hypervariable region cannot reliably distinguish *S. pneumoniae* from *S. suis*; both are reported as genus *Streptococcus*. Affects validation samples expecting these species (genus-level detection is the correct call; species discrimination is a 16S region limitation).')
    val_md.append('- **Minimap2 rescue threshold failure:** alignment evidence is present (>=100 mapped reads) but below the breadth-of-coverage threshold (>=0.25). Documents that rescue logic was attempted and correctly declined.')
    val_md.append('- **Abundance-driven:** target organism is present but at low abundance; other organisms dominate the 16S signal. Reflects 16S genus-detection biology, not pipeline error.')
    val_md.append('')

    # Sample-level concordance summary
    val_md.append('## Per-Sample Concordance Summary')
    val_md.append('')
    val_md.append('> The **Sample-level concordance** column is a coarse internal label. The **Final concordance (TAC + V4)** column is the rule-of-record: clinical sample is Concordant iff its expected target was detected AND retained by the V4 filter; NTC is Concordant iff no TAC bacterial target genus is retained; PC sample is flagged as Control with the TAC genera it retained (manual review against the full PC mix).')
    val_md.append('')
    val_md.append('| Sample | Category | Expected | Sample-level concordance | Final concordance (TAC + V4) | Biomass kept | % biomass kept | Notes |')
    val_md.append('|---|---|---|---|---|---|---|---|')
    sample_target_match = defaultdict(bool)
    for row in val_rows_out:
        if str(row.get('concordance', '')).startswith('TARGET match'):
            sample_target_match[row['sample']] = True
    # Build summary rows from EXPECTED_VALIDATION
    seen = set()
    for row in val_rows_out:
        sample = row['sample']
        if sample in seen:
            continue
        seen.add(sample)
        expected = row['expected']
        is_pc = row['category'] == 'control'
        is_ntc = row['category'] == 'NTC'
        if is_ntc:
            level = 'Concordant (NTC empty - PASS)' if row.get('total_reads', 0) == 0 else 'CONTAMINATED NTC'
            note_summary = 'NTC expected to be empty'
        elif is_pc:
            # Use final concordance verdict (computed earlier) for PC samples,
            # which checks all expected spike-in organisms are detected + retained.
            level = val_final_concordance.get(sample, 'Control')
            note_summary = 'PC sample; concordance evaluated against full expected spike-in set'
        elif sample_target_match.get(sample):
            level = 'Concordant'
            note_summary = '-'
        elif row.get('total_reads', 0) == 0:
            level = 'Discordant (no detections - pre-sequencing failure suspected)'
            note_summary = 'No taxa detected'
        else:
            level = 'Discordant (target not in detections)'
            note_summary = 'Other organisms dominate; review failure mode'
        val_md.append(f'| {sample} | {row["category"]} | {expected} | {level} | {val_final_concordance.get(sample, "-")} | {row.get("biomass_kept", "-")} | {row.get("pct_biomass_kept", "-")} | {note_summary} |')
    val_md.append('')

    with open(os.path.join(OUT_DIR, 'APPENDIX-VALIDATION-PANEL.md'), 'w') as f:
        f.write('\n'.join(val_md))

    # =======================================================================
    # APPENDIX 2: STUDY SAMPLES
    # =======================================================================
    study_runs = [r for r in RUN_DIRS if r != 'validation_panel']
    study_rows_out = []

    for run in study_runs:
        run_calls = [r for r in all_rows if r['run'] == run]
        run_samples = sorted({r['sample'] for r in run_calls
                              if not r['sample'].startswith(('PC', 'NTC', 'NC_', 'NC-'))})
        for sample in run_samples:
            total_reads = sample_totals.get(sample, 0)
            sample_dets = [r for r in run_calls if r['sample'] == sample
                           and (r['call'] in POSITIVE_CALLS or (r['call'] == 'Negative' and r['source'] == 'alignment' and r['reads'] >= 100))]
            if not sample_dets:
                study_rows_out.append({
                    'sample': sample,
                    'run': run,
                    'total_reads': total_reads,
                    'genus': '-',
                    'source': '-',
                    'reads': '-',
                    'pct': '-',
                    'pipeline_ntc': '-',
                    'same_run_ntc_max': '-',
                    'rescue_info': '-',
                    'rescue_tier': '-',
                    'filter_tier': 'n/a',
                    'filter_outcome': 'n/a',
                    'final_reported': '(none)',
                    'burk_species': '-',
                    'confidence': 'n/a',
                    'hypothesis': '-',
                    'followup': '-',
                    'notes': 'No taxa detected; pre-sequencing failure suspected',
                })
                continue
            for det in sample_dets:
                pct = (100.0 * det['reads'] / total_reads) if total_reads else 0.0
                # Same-run NTC cross-check
                srnm = same_run_ntc_max(run, det['genus']) if det['call'] in POSITIVE_CALLS else 0
                # Rescue info
                rescue_info = '-'
                rescue_tier = '-'
                if det['source'] == 'alignment':
                    rescue_info = f'mapped={det["reads"]}, breadth={det["breadth"]:.4f}, ntc={det["ntc_reads"]}'
                    if det['call'] == 'Confirmed':
                        rescue_tier = 'Tier 1: genus-level'
                    elif det['call'] == 'Probable':
                        rescue_tier = 'Tier 2: order-level (Rickettsiales detected)'
                    elif det['call'] == 'Detected':
                        rescue_tier = 'genus-level (rescued)' if det['rescued'] else 'genus-level'
                    elif det['call'] == 'Negative':
                        if det['breadth'] is not None and det['breadth'] < 0.25:
                            rescue_tier = 'FAIL: breadth < 0.25'
                        elif det['reads'] < 100:
                            rescue_tier = 'FAIL: reads < 100'
                        else:
                            rescue_tier = 'FAIL'

                # Filter decision
                if det['call'] not in POSITIVE_CALLS:
                    filter_decision = 'n/a (pipeline rejected pre-filter)'
                    filter_reason = 'Alignment-based rescue did not pass pipeline thresholds'
                else:
                    filter_decision, filter_reason = decide_filter(
                        det['genus'], det['reads'], pct, det['source'], sample, run,
                        total_reads, study_genus_counts.get(det['genus'], 0)
                    )
                # Filter tier label
                if filter_decision == 'KEEP':
                    filter_tier = 'KEEP (no tier match)'
                else:
                    filter_tier = filter_decision.replace('REMOVE_', '').replace('KEEP_', '')

                # Burkholderia species evidence
                burk_species = '-'
                if det['genus'] == 'Burkholderia':
                    kreport = find_kreport(run, sample)
                    species_reads = parse_kreport_species_reads(kreport, 'Burkholderia_pseudomallei') if kreport else 0
                    ntc_max_sp = same_run_ntc_burkpseudomallei_max(run)
                    burk_species = f'B.pseudomallei={species_reads} reads | NTC species max={ntc_max_sp} | threshold={B_PSEUDOMALLEI_MIN_SPECIES_READS}'

                # Final reported taxa
                if filter_decision.startswith('KEEP'):
                    if det['genus'] == 'Burkholderia':
                        # Will only happen if species safeguard passed
                        final = 'Burkholderia pseudomallei (species-level confirmed)'
                    else:
                        final = f'{det["genus"]} (genus-level)'
                else:
                    final = '(removed by filter)'

                # Confidence
                # For alignment-source rows, the right "NTC max" is the alignment
                # ntc_reads from calls.tsv (per-run, per-organism Minimap2 NTC),
                # not the centrifuge-based same_run_ntc_max.
                ntc_for_conf = det['ntc_reads'] if det['source'] == 'alignment' else srnm
                conf = confidence_label(pct, det['reads'], ntc_for_conf, det['genus'],
                                        is_pc_or_ntc=False)

                # Hypothesis tag
                hyp = hypothesis_tag(det['genus'])

                # Follow-up
                follow = '-'
                if det['genus'] in ('Leptospira',) and filter_decision.startswith('KEEP'):
                    follow = 'qPCR + paired serology; cross-check same-run NTC'
                elif det['genus'] == 'Brucella' and filter_decision.startswith('KEEP'):
                    follow = 'Serology (IgM/IgG) + Brucella-specific qPCR'
                elif det['genus'] == 'Mycoplasmopsis' and filter_decision.startswith('KEEP'):
                    follow = 'Species-level kreport check; Mycoplasma-specific PCR'
                elif det['genus'] == 'Burkholderia' and filter_decision == 'KEEP_BURK_PSEUDOMALLEI':
                    follow = 'B. pseudomallei-specific qPCR; clinical correlation (melioidosis)'
                elif det['genus'] in ('Orientia', 'Rickettsia') and filter_decision.startswith('KEEP'):
                    follow = 'Doxycycline coverage; species-specific qPCR for confirmation'

                # Notes
                note = ''
                if det['genus'] == 'Leptospira' and run == '6_and_7':
                    note = 'Same-run NTC2_ExDw_S13_L001 carries 78,691 Leptospira reads -- pipeline reports ntc_reads=0 here; investigate NCmax derivation'
                if det['genus'] == 'Burkholderia' and burk_species != '-':
                    note = 'See species-level evidence column; genus signal does not equal B. pseudomallei'
                if 'Streptococcus' in det['genus']:
                    note = 'V1-V3 cannot resolve species; possible fastidious species present'
                if det['genus'] == 'Mycoplasmopsis':
                    note = 'Cell-wall-deficient organism class; would NOT grow on routine aerobic subculture'

                study_rows_out.append({
                    'sample': sample,
                    'run': run,
                    'total_reads': total_reads,
                    'genus': det['genus'],
                    'source': det['source'],
                    'reads': det['reads'],
                    'pct': f'{pct:.2f}%',
                    'pipeline_ntc': det['ntc_reads'],
                    'same_run_ntc_max': srnm,
                    'rescue_info': rescue_info,
                    'rescue_tier': rescue_tier,
                    'filter_tier': filter_tier,
                    'filter_outcome': filter_decision,
                    'final_reported': final,
                    'burk_species': burk_species,
                    'confidence': conf,
                    'hypothesis': hyp,
                    'followup': follow,
                    'notes': note,
                })

    # Annotate biomass on study rows now that they exist
    study_kept, study_removed = biomass_split(study_rows_out)
    annotate_biomass(study_rows_out, study_kept, study_removed)

    # Write the study-samples appendix
    n_study_samples = len({r['sample'] for r in study_rows_out})
    n_empty = sum(1 for r in study_rows_out if r.get('genus') == '-' and r.get('reads') == '-')
    study_md = ['# Appendix: AFI Study Samples Detection Table',
                '',
                f'> **Cohort size correction:** earlier reports cited "56 AFI study samples." The actual count of non-control samples with calls.tsv outputs across runs 1_and_2 / 3 / 4_and_5 / 6_and_7 / 8_and_9 is {n_study_samples}. Of these, {n_empty} samples have zero detected taxa (pre-sequencing or library-prep failures suspected). The "56" figure in prior reports appears to have counted only samples with >=1 centrifuge-Detected genus call.',
                '> ',
                f'> **Format:** one row per (sample × detected genus × source) for the {n_study_samples} AFI study samples. `source = alignment` rows document Minimap2 evidence (including failed rescues, preserved for transparency).',
                '> ',
                '> **Pipeline NTC reads** is the `ntc_reads` field from the pipeline\'s `.calls.tsv` -- this is the NCmax value the upstream pipeline used for its detection call.',
                '> ',
                '> **Same-run NTC max (cross-check)** is computed independently here from the run\'s NTC `.calls.tsv` files. Mismatches between this and the pipeline NTC reads (e.g. *Leptospira* in run 6_and_7 where same-run NTC2_ExDw_S13_L001 carries 78,691 reads but pipeline reports `ntc_reads=0`) indicate the upstream pipeline\'s NCmax derivation should be audited.',
                '> ',
                '> **Burkholderia species evidence** appears only on Burkholderia rows and is the species-rank read count for *Burkholderia pseudomallei* from the Centrifuge kreport, alongside the run\'s NTC species-level count and the V4 threshold (500 species reads).',
                '> ',
                '> **Confidence** is a coarse label: HIGH (high abundance with NTC headroom), MODERATE (mid-abundance with reasonable NTC headroom), LOW (low abundance OR same-run NTC carries non-trivial signal), CANDIDATE (near noise floor).',
                '',
                '| Sample | Run | Total reads | Biomass kept (V4) | Biomass removed (V4) | % biomass kept | Detected genus | Source | Reads | % of sample | Pipeline NTC reads | Same-run NTC max (cross-check) | Minimap2 rescue info | Rescue tier | V4 filter tier matched | V4 filter outcome | Final reported taxa | Burkholderia species evidence | Confidence | Hypothesis class | Recommended follow-up | Notes |',
                '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for row in study_rows_out:
        study_md.append('| ' + ' | '.join(str(row[c]) for c in [
            'sample','run','total_reads','biomass_kept','biomass_removed','pct_biomass_kept',
            'genus','source','reads','pct','pipeline_ntc','same_run_ntc_max',
            'rescue_info','rescue_tier','filter_tier','filter_outcome','final_reported','burk_species',
            'confidence','hypothesis','followup','notes',
        ]) + ' |')
    study_md.append('')
    study_md.append('---')
    study_md.append('')
    # Biomass distribution summary (compute inline from study_rows_out)
    sample_pct = {}
    for r in study_rows_out:
        sample_pct.setdefault(r['sample'], r.get('pct_biomass_kept'))
    full_kept = sum(1 for v in sample_pct.values() if v == '100.0%')
    zero_kept = sum(1 for v in sample_pct.values() if v == '0.0%')
    mixed = sum(1 for v in sample_pct.values() if v not in ('100.0%', '0.0%', '-'))
    empty = sum(1 for v in sample_pct.values() if v == '-')
    study_md.append('## Per-Sample Biomass Distribution (V4 Filter Impact)')
    study_md.append('')
    study_md.append(f'- **{full_kept} samples** retain 100% of detected biomass (no contaminants present at detection level)')
    study_md.append(f'- **{zero_kept} samples** retain 0% (all detected reads removed as Tier A / B / 1 / 2 contaminants)')
    study_md.append(f'- **{mixed} samples** have a mix of kept and removed signal (the most informative cases)')
    study_md.append(f'- **{empty} samples** had zero detections to begin with (pre-sequencing failures)')
    study_md.append(f'- See "Biomass kept (V4)" / "Biomass removed (V4)" / "% biomass kept" columns on each row for per-sample values.')
    study_md.append('')

    study_md.append('## Filter Tier Glossary')
    study_md.append('')
    study_md.append('- **TIER_A:** removed as high-confidence kit/skin/water contaminant (Pseudomonas, Ralstonia, Bradyrhizobium, Sphingomonas, Stenotrophomonas, Methylobacterium, Acinetobacter, Cutibacterium, Staphylococcus, Corynebacterium, Brevundimonas)')
    study_md.append('- **BURK_PSEUDOMALLEI / BURK (preserved/removed):** Burkholderia genus detection evaluated against species-level *B. pseudomallei* reads from the Centrifuge kreport. Preservation requires species reads >= 500 AND > same-run NTC species max.')
    study_md.append('- **NC_ONLY:** removed because the genus only appears in NTC samples in this dataset.')
    study_md.append('- **TIER_1:** removed as ultra-low-abundance environmental noise (median <0.5%). Note: *Mycoplasmopsis* and *Nitrospira* were previously in this tier in V3 but have been MOVED OUT in V4 after raw-data inspection showed substantial abundance.')
    study_md.append('- **TIER_2:** removed only if abundance <1% AND detected in <=3 samples in the cohort.')
    study_md.append('- **KEEP:** passes filter; appears in final reported taxa.')
    study_md.append('')
    study_md.append('## Confidence Glossary')
    study_md.append('')
    study_md.append('- **HIGH:** post-filter abundance >=30% AND same-run NTC carries minimal background for this organism.')
    study_md.append('- **MODERATE:** post-filter abundance 5-30%, NTC headroom acceptable.')
    study_md.append('- **LOW:** post-filter abundance 1-5%, OR same-run NTC carries non-trivial signal for this organism.')
    study_md.append('- **CANDIDATE (near noise floor):** post-filter abundance <1%; treat as a candidate observation, not a diagnosis.')
    study_md.append('')

    with open(os.path.join(OUT_DIR, 'APPENDIX-STUDY-SAMPLES.md'), 'w') as f:
        f.write('\n'.join(study_md))

    print(f"Wrote APPENDIX-VALIDATION-PANEL.md ({len(val_rows_out)} rows)")
    print(f"Wrote APPENDIX-STUDY-SAMPLES.md ({len(study_rows_out)} rows)")

    # =======================================================================
    # EXCEL EXPORT
    # =======================================================================
    val_columns = [
        ('sample', 'Sample'),
        ('run', 'Run'),
        ('category', 'Category'),
        ('expected', 'Expected organism'),
        ('total_reads', 'Total reads in sample'),
        ('biomass_kept', 'Biomass kept (V4 filter)'),
        ('biomass_removed', 'Biomass removed (V4 filter)'),
        ('pct_biomass_kept', '% biomass kept'),
        ('genus', 'Detected genus'),
        ('source', 'Source'),
        ('reads', 'Reads'),
        ('pct', '% of sample'),
        ('ntc_reads', 'Pipeline NTC reads'),
        ('rescue_info', 'Minimap2 rescue info'),
        ('rescue_tier', 'Rescue tier'),
        ('filter_outcome', 'V4 filter outcome'),
        ('concordance', 'Per-row concordance'),
        ('failure_mode', 'Failure mode'),
        ('final_concordance', 'Final concordance (TAC + V4)'),
        ('notes', 'Notes'),
    ]
    study_columns = [
        ('sample', 'Sample'),
        ('run', 'Run'),
        ('total_reads', 'Total reads in sample'),
        ('biomass_kept', 'Biomass kept (V4 filter)'),
        ('biomass_removed', 'Biomass removed (V4 filter)'),
        ('pct_biomass_kept', '% biomass kept'),
        ('genus', 'Detected genus'),
        ('source', 'Source'),
        ('reads', 'Reads'),
        ('pct', '% of sample'),
        ('pipeline_ntc', 'Pipeline NTC reads (calls.tsv)'),
        ('same_run_ntc_max', 'Same-run NTC max (cross-check)'),
        ('rescue_info', 'Minimap2 rescue info'),
        ('rescue_tier', 'Rescue tier'),
        ('filter_tier', 'V4 filter tier matched'),
        ('filter_outcome', 'V4 filter outcome'),
        ('final_reported', 'Final reported taxa'),
        ('burk_species', 'Burkholderia species evidence'),
        ('confidence', 'Confidence'),
        ('hypothesis', 'Hypothesis class'),
        ('followup', 'Recommended follow-up'),
        ('notes', 'Notes'),
    ]
    val_df = pd.DataFrame([{disp: r.get(key, '') for key, disp in val_columns} for r in val_rows_out])
    study_df = pd.DataFrame([{disp: r.get(key, '') for key, disp in study_columns} for r in study_rows_out])

    # Per-sample validation concordance summary (reuse the earlier logic)
    val_summary_rows = []
    seen2 = set()
    for row in val_rows_out:
        sample = row['sample']
        if sample in seen2:
            continue
        seen2.add(sample)
        is_pc = row['category'] == 'control'
        is_ntc = row['category'] == 'NTC'
        any_target_match = any(
            r['sample'] == sample and str(r.get('concordance', '')).startswith('TARGET match')
            for r in val_rows_out
        )
        if is_ntc:
            level = 'Concordant (NTC empty - PASS)' if row.get('total_reads', 0) == 0 else 'CONTAMINATED NTC'
        elif is_pc:
            level = 'Control - review against full PC list'
        elif any_target_match:
            level = 'Concordant'
        elif row.get('total_reads', 0) == 0:
            level = 'Discordant (no detections - pre-sequencing failure suspected)'
        else:
            level = 'Discordant (target not in detections)'
        val_summary_rows.append({
            'Sample': sample,
            'Category': row['category'],
            'Expected': row['expected'],
            'Sample-level concordance': level,
            'Final concordance (TAC + V4)': val_final_concordance.get(sample, '-'),
            'Biomass kept (V4 filter)': row.get('biomass_kept', '-'),
            'Biomass removed (V4 filter)': row.get('biomass_removed', '-'),
            '% biomass kept': row.get('pct_biomass_kept', '-'),
        })
    val_summary_df = pd.DataFrame(val_summary_rows)

    # Per-sample study summary
    study_summary_rows = []
    seen3 = set()
    for row in study_rows_out:
        sample = row['sample']
        if sample in seen3:
            continue
        seen3.add(sample)
        # Aggregate by listing kept genera and any Rickettsiales rescue
        kept_genera = sorted({
            r['genus'] for r in study_rows_out
            if r['sample'] == sample and str(r.get('filter_outcome', '')).startswith('KEEP')
            and r.get('genus') not in (None, '-', '')
        })
        rick_evidence = [
            f"{r['genus']} ({r.get('rescue_tier', '?')}, {r.get('reads', '?')} reads, breadth {r.get('rescue_info','').split('breadth=')[-1].split(',')[0] if 'breadth=' in r.get('rescue_info','') else '?'})"
            for r in study_rows_out
            if r['sample'] == sample and r.get('genus') in ('Orientia', 'Rickettsia')
            and str(r.get('filter_outcome', '')).startswith('KEEP')
        ]
        # Highest-confidence retained organism in this sample
        retained_rows = [r for r in study_rows_out
                         if r['sample'] == sample and str(r.get('filter_outcome', '')).startswith('KEEP')]
        if retained_rows:
            top = max(retained_rows, key=lambda r: float(r.get('pct', '0').rstrip('%')) if isinstance(r.get('pct'), str) and r['pct'].endswith('%') else 0)
            top_summary = f"{top.get('genus','')} ({top.get('pct','-')}, {top.get('confidence','-')})"
        else:
            top_summary = '(none retained)'
        study_summary_rows.append({
            'Sample': sample,
            'Run': row['run'],
            'Total reads in sample': row.get('total_reads', '-'),
            'Biomass kept (V4 filter)': row.get('biomass_kept', '-'),
            'Biomass removed (V4 filter)': row.get('biomass_removed', '-'),
            '% biomass kept': row.get('pct_biomass_kept', '-'),
            'Retained genera (count)': len(kept_genera),
            'Retained genera (list)': ', '.join(kept_genera) if kept_genera else '(none)',
            'Top retained organism': top_summary,
            'Rickettsiales rescue evidence': ' | '.join(rick_evidence) if rick_evidence else '(none)',
        })
    study_summary_df = pd.DataFrame(study_summary_rows)

    xlsx_path = os.path.join(OUT_DIR, 'APPENDICES.xlsx')
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
        val_df.to_excel(writer, sheet_name='Validation - detections', index=False)
        val_summary_df.to_excel(writer, sheet_name='Validation - sample summary', index=False)
        study_df.to_excel(writer, sheet_name='Study - detections', index=False)
        study_summary_df.to_excel(writer, sheet_name='Study - sample summary', index=False)
        # Light formatting: freeze top row on each sheet
        for ws in writer.book.worksheets:
            ws.freeze_panes = 'A2'
            # Auto-set column widths from header text length (capped)
            for col_idx, col_cells in enumerate(ws.iter_cols(min_row=1, max_row=1), start=1):
                from openpyxl.utils import get_column_letter
                header = str(col_cells[0].value or '')
                ws.column_dimensions[get_column_letter(col_idx)].width = min(max(len(header) + 2, 12), 60)
    print(f"Wrote APPENDICES.xlsx (4 sheets: validation detections + summary, study detections + summary)")


if __name__ == '__main__':
    main()
