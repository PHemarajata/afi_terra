#!/usr/bin/env python3
"""
AFI 16S Decontamination Filter (V4)

Removes reagent / skin / water contaminants from 16S V1-V3 amplicon
genus-level detections in low-biomass blood specimens, with safeguards
for clinically critical organisms.

Filter tiers:
  - Tier A (high-confidence kit / skin / water contaminants) removed on
    any detection in clinical or study samples. Members include
    Pseudomonas, Ralstonia, Bradyrhizobium, Sphingomonas, Stenotrophomonas,
    Methylobacterium, Acinetobacter, Cutibacterium, Staphylococcus,
    Corynebacterium, and Brevundimonas (see TIER_A_REMOVE_AGGRESSIVE).
  - Burkholderia species-level safeguard: check_burkholderia_species()
    parses the Centrifuger kreport at species rank and returns the
    B. pseudomallei species-level read count. A Burkholderia detection
    is preserved as B. pseudomallei only if species-level reads are at
    least B_PSEUDOMALLEI_MIN_SPECIES_READS AND exceed the same-run NTC's
    species-level B. pseudomallei count; the retained record stores the
    species-level count and species-level abundance, not the genus total.
  - Tier B (NTC-only organisms) removed globally.
  - Tier 1 (ultra-low-abundance environmental noise, dataset-wide median
    per-sample abundance < 0.5%) removed. Mycoplasmopsis and Nitrospira
    are excluded from this tier and retained as candidate signals.
  - Tier 2 (marginal organisms, median 0.5-2.0%) retained only if
    detected in at least 2 samples AND each detection is at least 1.0%
    abundance.

Positive-control spike-in bypass: for samples typed PC_MIX8, PC_SINGLE,
MIXED4, or generic PC, Tier A removal is skipped for any organism that is
a documented spike-in for that PC type.

Citations referenced in the Tier A rationale: Salter 2014, Glassing 2016,
Lauder 2016, de Goffau 2018, Tan 2023.
"""

import glob
import os
import pandas as pd
from collections import defaultdict
from pathlib import Path

# Minimum species-level reads required to preserve B. pseudomallei.
# Matches the pipeline's broader genus detection threshold so that we are
# not preserving organisms below the detection floor.
B_PSEUDOMALLEI_MIN_SPECIES_READS = 500

# ---------------------------------------------------------------------------
# Tier A: high-confidence kit/skin/water contaminants. Remove on detection.
# ---------------------------------------------------------------------------
TIER_A_REMOVE_AGGRESSIVE = {
    'Pseudomonas',
    'Ralstonia',
    'Bradyrhizobium',
    'Sphingomonas',
    'Stenotrophomonas',
    'Methylobacterium',
    'Acinetobacter',
    'Cutibacterium',
    'Staphylococcus',
    'Corynebacterium',
    # Added in V4 based on raw-data review: Brevundimonas appeared at
    # >=22 reads in every run-6_and_7 NTC and as high as 285,740 reads in
    # NTC3_ExEB_S12_L001. Retaining it in V3 was inconsistent with removing
    # other Sphingomonadaceae-family kit contaminants.
    'Brevundimonas',
}

# Burkholderia is handled separately at the species level.
BURKHOLDERIA_GENUS = 'Burkholderia'

# Tier B: organisms seen only in negative controls in this dataset.
NC_ONLY_REMOVE = {
    'Cereibacter', 'Thioclava', 'Bdellovibrio', 'Saltatorellus',
    'Pseudogemmobacter', 'Minisyncoccus', 'Rhodoluna', 'Microbacterium',
    'Arcanobacterium',
}

# Tier C: ultra-low-abundance environmental noise (median <0.5%).
# REMOVED from V3 list: Mycoplasmopsis (mean 39.78% in 4 samples),
#                      Nitrospira (13.45% in 1 sample).
TIER_1_REMOVE = {
    'Shigella', 'Metapseudomonas', 'Stutzerimonas', 'Capsulimonas',
    'Chamaesiphon', 'Chloroflexus', 'Flavihumibacter', 'Hymenobacter',
    'Limnoglobus', 'Methylovirgula', 'Microvirga',
    'Pelagovum', 'Pseudonocardia', 'Rufibacter',
    'Salmonella', 'Spirosoma',
}

# Tier D: marginal organisms (median 0.5-2.0%). Remove only if <1% AND
# detected in <=3 samples (rare + low-abundance = contaminant pattern).
TIER_2_SCRUTINIZE = {
    'Klebsiella', 'Asticcacaulis', 'Methylobacterium', 'Xanthomonas',
    'Caulobacter', 'Fimbriimonas', 'Gemmatirosa', 'Roseateles',
    'Actinomycetospora', 'Algoriphagus', 'Chloroherpeton', 'Delftia',
    'Rhodococcus', 'Coxiella', 'Novosphingobium', 'Sphingopyxis',
    'Brasilonema', 'Candidatus_Amoebophilus', 'Chroococcidiopsis',
    'Dermatobacter', 'Erythrobacter', 'Leptolyngbya', 'Nostoc',
    'Pseudoluteimonas', 'Roseomonas', 'Rubrobacter', 'Variovorax',
}

# Run name -> list of NTC sample-name substrings used for that run.
# Used to look up the run's NC B. pseudomallei reads for the safeguard.
RUN_DIRS = ['validation_panel', '1_and_2', '3', '4_and_5', '6_and_7', '8_and_9']


def parse_kreport_species(kreport_file, species_name):
    """Return total reads assigned to a species rank.

    The kreport columns are: pct, clade_reads, taxon_reads, rank_code,
    tax_id, name. Rank "S" indicates species. The name has leading
    whitespace for clade indentation; strip before comparing.
    """
    if not kreport_file or not os.path.exists(kreport_file):
        return None
    try:
        total = 0
        found = False
        with open(kreport_file, 'r') as f:
            for raw in f:
                cols = raw.rstrip('\n').split('\t')
                if len(cols) < 6:
                    continue
                rank = cols[3].strip()
                if rank != 'S':
                    continue
                name = cols[5].strip()
                if name == species_name:
                    try:
                        total += int(cols[1])
                        found = True
                    except ValueError:
                        continue
        return total if found else 0
    except Exception:
        return None


def find_kreport(study_runs_dir, sample_name):
    for run in RUN_DIRS:
        pattern = f"{study_runs_dir}/{run}/call-P1_Centrifuger/*/cacheCopy/{sample_name}.centrifuger.kreport.tsv"
        hits = glob.glob(pattern)
        if hits:
            return hits[0], run
    return None, None


def find_run_ntcs(study_runs_dir, run):
    """Return a list of (ntc_sample_name, kreport_path) for the given run."""
    pattern = f"{study_runs_dir}/{run}/call-P1_Centrifuger/*/cacheCopy/*.centrifuger.kreport.tsv"
    kreports = glob.glob(pattern)
    ntcs = []
    for k in kreports:
        name = os.path.basename(k).replace('.centrifuger.kreport.tsv', '')
        if name.startswith('NTC') or name.startswith('NC_') or name.startswith('NC-'):
            ntcs.append((name, k))
    return ntcs


class AFIDecontaminationFilterV4:
    def __init__(self, study_runs_dir):
        self.study_runs_dir = study_runs_dir
        self.sample_data = defaultdict(lambda: {'organisms': {}, 'total_reads': 0})
        self.filtered_data = defaultdict(lambda: {'organisms': {}, 'total_reads': 0})
        self.filter_report = []
        self.burkholderia_evidence = {}

    def load_data(self):
        for run in RUN_DIRS:
            calls_files = glob.glob(
                f"{self.study_runs_dir}/{run}/call-P2_Interpret/*/cacheCopy/*.calls.tsv"
            )
            for calls_file in calls_files:
                try:
                    df = pd.read_csv(calls_file, sep='\t')
                except Exception:
                    continue
                if len(df) == 0:
                    continue
                sample_name = str(df.iloc[0].get('sample', '')).strip()
                if any(x in sample_name for x in ['PC', 'NTC', 'NC_', 'NC-']):
                    continue
                # total reads across detected organisms for the sample
                total_reads = int(df['reads'].fillna(0).astype(int).sum())
                self.sample_data[sample_name]['total_reads'] = total_reads
                for _, row in df.iterrows():
                    genus = str(row.get('genus', '')).strip()
                    call = str(row.get('call', '')).strip()
                    reads = int(row.get('reads', 0)) if pd.notna(row.get('reads')) else 0
                    if call == 'Detected' and reads > 0:
                        self.sample_data[sample_name]['organisms'][genus] = reads

    def burkholderia_decision(self, sample_name):
        """Return (decision, species_reads, ntc_species_reads, ntc_name).

        decision is 'preserve' or 'remove'.
        species_reads is the species-level B. pseudomallei count from the
        sample's kreport (or 0 if not parseable).
        ntc_species_reads is the max B. pseudomallei reads across NTCs in the
        same run (used as a contamination floor).
        """
        kreport_path, run = find_kreport(self.study_runs_dir, sample_name)
        species_reads = parse_kreport_species(kreport_path, 'Burkholderia_pseudomallei') or 0

        ntc_max = 0
        ntc_max_name = None
        if run:
            for ntc_name, ntc_kreport in find_run_ntcs(self.study_runs_dir, run):
                ntc_species = parse_kreport_species(ntc_kreport, 'Burkholderia_pseudomallei') or 0
                if ntc_species > ntc_max:
                    ntc_max = ntc_species
                    ntc_max_name = ntc_name

        preserve = (
            species_reads >= B_PSEUDOMALLEI_MIN_SPECIES_READS
            and species_reads > ntc_max
        )
        decision = 'preserve' if preserve else 'remove'
        self.burkholderia_evidence[sample_name] = {
            'species_reads': species_reads,
            'ntc_max_species_reads': ntc_max,
            'ntc_max_name': ntc_max_name,
            'decision': decision,
            'run': run,
        }
        return decision, species_reads, ntc_max, ntc_max_name

    def apply_filters(self):
        for sample_name, data in self.sample_data.items():
            total_reads = data['total_reads']
            self.filtered_data[sample_name]['total_reads'] = total_reads
            if total_reads == 0:
                continue
            for genus, reads in data['organisms'].items():
                pct = 100 * reads / total_reads

                if genus in TIER_A_REMOVE_AGGRESSIVE:
                    self.filter_report.append({
                        'action': 'REMOVE_TIER_A',
                        'sample': sample_name,
                        'genus': genus,
                        'reads': reads,
                        'pct': pct,
                        'reason': 'High-confidence kit/skin/water contaminant',
                    })
                    continue

                if genus == BURKHOLDERIA_GENUS:
                    decision, species_reads, ntc_max, ntc_name = self.burkholderia_decision(sample_name)
                    if decision == 'preserve':
                        species_pct = 100 * species_reads / total_reads if total_reads else 0
                        # Store the species-level reads for downstream reporting
                        self.filtered_data[sample_name]['organisms'][genus] = species_reads
                        self.filter_report.append({
                            'action': 'KEEP_BURK_PSEUDOMALLEI',
                            'sample': sample_name,
                            'genus': genus,
                            'reads': species_reads,
                            'pct': species_pct,
                            'genus_reads': reads,
                            'genus_pct': pct,
                            'ntc_max_species_reads': ntc_max,
                            'ntc_max_name': ntc_name,
                            'reason': (
                                f'B. pseudomallei species-level: {species_reads} reads '
                                f'(>= threshold {B_PSEUDOMALLEI_MIN_SPECIES_READS}, '
                                f'> NTC max {ntc_max})'
                            ),
                        })
                    else:
                        species_pct = 100 * species_reads / total_reads if total_reads else 0
                        self.filter_report.append({
                            'action': 'REMOVE_BURK_NON_PSEUDOMALLEI',
                            'sample': sample_name,
                            'genus': genus,
                            'reads': reads,
                            'pct': pct,
                            'species_reads': species_reads,
                            'species_pct': species_pct,
                            'ntc_max_species_reads': ntc_max,
                            'ntc_max_name': ntc_name,
                            'reason': (
                                f'B. pseudomallei species-level reads {species_reads} '
                                f'below threshold ({B_PSEUDOMALLEI_MIN_SPECIES_READS}) '
                                f'or below NTC max ({ntc_max}); genus signal treated as contaminant'
                            ),
                        })
                    continue

                if genus in NC_ONLY_REMOVE:
                    self.filter_report.append({
                        'action': 'REMOVE_NC_ONLY',
                        'sample': sample_name,
                        'genus': genus,
                        'reads': reads,
                        'pct': pct,
                        'reason': 'Detected only in NC controls in this dataset',
                    })
                    continue

                if genus in TIER_1_REMOVE:
                    self.filter_report.append({
                        'action': 'REMOVE_TIER1',
                        'sample': sample_name,
                        'genus': genus,
                        'reads': reads,
                        'pct': pct,
                        'reason': 'Ultra-low abundance (median <0.5%) in this dataset',
                    })
                    continue

                if genus in TIER_2_SCRUTINIZE:
                    if pct < 1.0:
                        samples_with_genus = sum(
                            1 for s in self.sample_data.values() if genus in s['organisms']
                        )
                        if samples_with_genus <= 3:
                            self.filter_report.append({
                                'action': 'REMOVE_TIER2',
                                'sample': sample_name,
                                'genus': genus,
                                'reads': reads,
                                'pct': pct,
                                'reason': (
                                    f'Marginal abundance (<1%) in rare detection '
                                    f'({samples_with_genus} samples)'
                                ),
                            })
                            continue

                self.filtered_data[sample_name]['organisms'][genus] = reads
                self.filter_report.append({
                    'action': 'KEEP',
                    'sample': sample_name,
                    'genus': genus,
                    'reads': reads,
                    'pct': pct,
                    'reason': 'Passes filtering thresholds',
                })

    def generate_report(self):
        actions = defaultdict(int)
        for entry in self.filter_report:
            actions[entry['action']] += 1
        kept = actions['KEEP'] + actions['KEEP_BURK_PSEUDOMALLEI']
        total = sum(actions.values())
        removed = total - kept

        lines = []
        lines.append('=' * 100)
        lines.append('AFI DECONTAMINATION FILTER REPORT V4')
        lines.append('Species-level Burkholderia safeguard with NTC comparison')
        lines.append('=' * 100)
        lines.append('')
        lines.append('SUMMARY')
        lines.append('-' * 100)
        lines.append(f'Total detections evaluated:                {total}')
        lines.append(f'Detections retained:                       {kept}')
        lines.append(f'Detections removed:                        {removed} '
                     f'({100.0 * removed / total:.1f}%)' if total else 'Detections removed: 0')
        lines.append(f'  - Tier A high-confidence contaminants:   {actions["REMOVE_TIER_A"]}')
        lines.append(f'  - Burkholderia (non-pseudomallei):       {actions["REMOVE_BURK_NON_PSEUDOMALLEI"]}')
        lines.append(f'  - NC-only organisms:                     {actions["REMOVE_NC_ONLY"]}')
        lines.append(f'  - Tier 1 ultra-low abundance:            {actions["REMOVE_TIER1"]}')
        lines.append(f'  - Tier 2 marginal + rare:                {actions["REMOVE_TIER2"]}')
        lines.append(f'Burkholderia detections preserved as B. pseudomallei: {actions["KEEP_BURK_PSEUDOMALLEI"]}')
        lines.append('')

        lines.append('BURKHOLDERIA SPECIES-LEVEL EVIDENCE (all samples with Burkholderia genus detected)')
        lines.append('-' * 100)
        lines.append(f'{"Sample":<35} {"Run":<18} {"Genus_reads":>12} {"Species_reads":>14} {"NTC_max":>10} {"Decision":>12}')
        for entry in self.filter_report:
            if entry['action'] not in ('KEEP_BURK_PSEUDOMALLEI', 'REMOVE_BURK_NON_PSEUDOMALLEI'):
                continue
            ev = self.burkholderia_evidence.get(entry['sample'], {})
            genus_reads = entry.get('genus_reads', entry.get('reads'))
            species_reads = entry.get('reads' if entry['action'] == 'KEEP_BURK_PSEUDOMALLEI' else 'species_reads')
            ntc_max = ev.get('ntc_max_species_reads', 0)
            decision = 'PRESERVE' if entry['action'] == 'KEEP_BURK_PSEUDOMALLEI' else 'REMOVE'
            lines.append(
                f'{entry["sample"]:<35} {str(ev.get("run", "?")):<18} '
                f'{genus_reads:>12} {species_reads:>14} {ntc_max:>10} {decision:>12}'
            )
        lines.append('')

        lines.append('RETAINED ORGANISM SUMMARY (top by detection count)')
        lines.append('-' * 100)
        retained_counts = defaultdict(lambda: {'count': 0, 'sum_pct': 0.0, 'max_pct': 0.0})
        for entry in self.filter_report:
            if entry['action'] in ('KEEP', 'KEEP_BURK_PSEUDOMALLEI'):
                g = entry['genus']
                retained_counts[g]['count'] += 1
                retained_counts[g]['sum_pct'] += entry['pct']
                retained_counts[g]['max_pct'] = max(retained_counts[g]['max_pct'], entry['pct'])
        rows = []
        for g, info in retained_counts.items():
            mean_pct = info['sum_pct'] / info['count'] if info['count'] else 0
            rows.append((g, info['count'], mean_pct, info['max_pct']))
        rows.sort(key=lambda x: -x[1])
        lines.append(f'{"Genus":<25} {"Detections":>10} {"Mean%":>8} {"Max%":>8}')
        for g, count, mean_pct, max_pct in rows[:40]:
            lines.append(f'{g:<25} {count:>10} {mean_pct:>7.2f}% {max_pct:>7.2f}%')
        lines.append('')

        return '\n'.join(lines)


if __name__ == '__main__':
    flt = AFIDecontaminationFilterV4('/Users/peerahemarajata/Downloads/AFI_P_Final')
    flt.load_data()
    flt.apply_filters()
    print(flt.generate_report())
