#!/usr/bin/env python3
"""
AFI Terra Sheet Builder
-----------------------
Two-screen wizard that produces a Terra-compatible TSV for the AFI_16S_Batch
workflow.

Screen 1  — How many runs? → enter run names → browse per-run FASTQ folders
Screen 2  — Sample metadata table with run_id dropdown, sample_type combo,
            auto-fill of mode/expected_taxa, batch editing, Validate, Export TSV
"""

import csv
import difflib
import io
import os
import re
import sys
from datetime import date
from pathlib import Path

from PySide6.QtCore import Qt, QDate, QModelIndex, QAbstractTableModel, QSortFilterProxyModel
from PySide6.QtGui import QColor, QBrush, QFont, QPixmap, QIcon
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QStackedWidget,
    QVBoxLayout, QHBoxLayout, QGridLayout, QFormLayout,
    QLabel, QPushButton, QLineEdit, QSpinBox, QDateEdit,
    QTableWidget, QTableWidgetItem, QHeaderView,
    QFileDialog, QMessageBox, QComboBox, QScrollArea,
    QGroupBox, QFrame, QSizePolicy, QAbstractItemView,
    QStyledItemDelegate, QStyle,
)

# ---------------------------------------------------------------------------
# Asset path helper (works both in dev and PyInstaller --onefile bundles)
# ---------------------------------------------------------------------------

def _bundle_path(relative: str) -> Path:
    """Resolve a path relative to the script; works in dev and PyInstaller."""
    base = Path(getattr(sys, "_MEIPASS", Path(__file__).parent))
    return base / relative


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SAMPLE_TYPES = ["clinical", "NTC", "NC", "PC_MIX8", "PC_SINGLE", "MIXED4", "PC"]

PC_MIX8_EXPECTED = (
    "Bacillus;Listeria;Staphylococcus;Enterococcus;"
    "Limosilactobacillus;Salmonella;Escherichia;Pseudomonas"
)
MIXED4_EXPECTED = (
    "Escherichia coli;Pseudomonas aeruginosa;"
    "Streptococcus pneumoniae;Streptococcus suis"
)
PC_SINGLE_OPTIONS = [
    "Streptococcus pneumoniae",
    "Streptococcus suis",
    "Escherichia coli",
    "Pseudomonas aeruginosa",
]

AUTO_FILL = {
    "clinical":  ("routine",     ""),
    "NTC":       ("routine",     ""),
    "NC":        ("routine",     ""),
    "PC_MIX8":   ("validation",  PC_MIX8_EXPECTED),
    "PC_SINGLE": ("validation",  PC_SINGLE_OPTIONS[0]),
    "MIXED4":    ("validation",  MIXED4_EXPECTED),
    "PC":        ("routine",     ""),
}

# Row highlight colors (APHL palette)
COLOR_NTC = QColor("#e8f5e2")   # green tint  — NTC / NC
COLOR_PC  = QColor("#fff8e0")   # yellow tint — positive controls
COLOR_CLN = QColor("#ffffff")   # white       — clinical / routine

NTC_TYPES = {"NTC", "NC"}
PC_TYPES  = {"PC_MIX8", "PC_SINGLE", "MIXED4", "PC"}

# Columns in the metadata table
COL_SAMPLE_ID    = 0
COL_RUN_ID       = 1
COL_SAMPLE_TYPE  = 2
COL_MODE         = 3
COL_EXPECTED     = 4
NUM_COLS         = 5
COL_HEADERS      = ["sample_id", "run_id", "sample_type", "mode", "expected_taxa"]

R1_PATTERNS = re.compile(r"_R1[_.]|_R1_001\.")
R2_PATTERNS = re.compile(r"_R2[_.]|_R2_001\.")
FASTQ_SUFFIXES = {".fastq", ".fastq.gz", ".fq", ".fq.gz"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def is_fastq(p: Path) -> bool:
    return "".join(p.suffixes[-2:]) in {".fastq.gz", ".fq.gz"} or p.suffix in {".fastq", ".fq"}


def discover_pairs(folder: str) -> list[tuple[str, str, str]]:
    """Return (sample_id, r1_path, r2_path) for every R1/R2 pair in a folder."""
    p = Path(folder)
    r1_files = sorted(f for f in p.iterdir() if is_fastq(f) and R1_PATTERNS.search(f.name))
    return _pair_r1_list(r1_files, p)


def pair_fastqs_from_files(file_paths: list[str]) -> list[tuple[str, str, str]]:
    """Auto-pair a flat list of FASTQ paths by R1/R2 name convention."""
    paths = [Path(f) for f in file_paths]
    r1_files = sorted(f for f in paths if R1_PATTERNS.search(f.name))
    r2_lookup: dict[str, Path] = {f.name: f for f in paths if R2_PATTERNS.search(f.name)}
    pairs = []
    for r1 in r1_files:
        r2_name = R1_PATTERNS.sub(lambda m: m.group().replace("R1", "R2"), r1.name)
        r2 = r2_lookup.get(r2_name) or r2_lookup.get(r2_name + ".gz")
        if r2:
            stem = R1_PATTERNS.split(r1.name)[0]
            pairs.append((stem, str(r1), str(r2)))
    return pairs


def _pair_r1_list(r1_files: list[Path], folder: Path) -> list[tuple[str, str, str]]:
    pairs = []
    for r1 in r1_files:
        r2_name = R1_PATTERNS.sub(lambda m: m.group().replace("R1", "R2"), r1.name)
        r2 = folder / r2_name
        if not r2.exists():
            r2_gz = folder / (r2_name + ".gz")
            r2 = r2_gz if r2_gz.exists() else None
        if r2:
            stem = R1_PATTERNS.split(r1.name)[0]
            pairs.append((stem, str(r1), str(r2)))
    return pairs


def fuzzy_match_sample(
    sample_id: str,
    candidates: list[Path],
    cutoff: float = 0.6,
) -> tuple[Path | None, Path | None]:
    r1_cands = [f for f in candidates if R1_PATTERNS.search(f.name)]
    if not r1_cands:
        return None, None
    sid_lower = sample_id.lower()

    def stem_of(f: Path) -> str:
        return R1_PATTERNS.split(f.name)[0].lower()

    for r1 in r1_cands:
        st = stem_of(r1)
        if st == sid_lower or st.startswith(sid_lower + "_") or st.startswith(sid_lower + "-"):
            r2 = _find_r2(r1, candidates)
            if r2:
                return r1, r2
    for r1 in r1_cands:
        if sid_lower in stem_of(r1):
            r2 = _find_r2(r1, candidates)
            if r2:
                return r1, r2
    stems = [stem_of(r1) for r1 in r1_cands]
    matches = difflib.get_close_matches(sid_lower, stems, n=1, cutoff=cutoff)
    if matches:
        best = r1_cands[stems.index(matches[0])]
        r2 = _find_r2(best, candidates)
        if r2:
            return best, r2
    return None, None


def _find_r2(r1: Path, candidates: list[Path]) -> Path | None:
    r2_name = R1_PATTERNS.sub(lambda m: m.group().replace("R1", "R2"), r1.name)
    r2_lookup = {f.name: f for f in candidates if R2_PATTERNS.search(f.name)}
    return r2_lookup.get(r2_name) or r2_lookup.get(r2_name + ".gz")


def parse_mapping_tsv(path: str) -> tuple[list[str], list[dict]]:
    ALIASES = {
        "run_id":    {"run_id", "run", "run_name"},
        "sample_id": {"sample_id", "sample", "sample_name"},
        "r1":        {"r1_fastq", "r1", "fastq_r1", "read1", "r1_path"},
        "r2":        {"r2_fastq", "r2", "fastq_r2", "read2", "r2_path"},
    }
    with open(path, newline="", encoding="utf-8") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.DictReader(fh, dialect=dialect)
        raw_cols = reader.fieldnames or []
        col_map: dict[str, str] = {}
        for raw in raw_cols:
            for canon, aliases in ALIASES.items():
                if raw.strip().lower() in aliases:
                    col_map[raw] = canon
                    break
        rows = []
        for raw_row in reader:
            row: dict[str, str] = {"run_id": "", "sample_id": "", "r1": "", "r2": ""}
            for raw_col, canon in col_map.items():
                row[canon] = (raw_row.get(raw_col) or "").strip()
            rows.append(row)
    return raw_cols, rows


def unique_sample_id(base: str, existing: set[str]) -> str:
    if base not in existing:
        return base
    i = 2
    while f"{base}_{i}" in existing:
        i += 1
    return f"{base}_{i}"


def parse_terra_tsv(path: str) -> dict:
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        headers = reader.fieldnames or []
        raw_rows = list(reader)

    if not headers:
        raise ValueError("TSV has no column headers.")

    entity_col = headers[0]
    m = re.match(r"^entity:(.+?)_id$", entity_col)
    table_name = m.group(1) if m else ""

    rows: list[dict] = []
    run_names_seen: list[str] = []
    errors: list[str] = []

    for i, raw in enumerate(raw_rows, start=1):
        sid  = (raw.get(entity_col)           or "").strip()
        rid  = (raw.get("run_id")              or "").strip()
        st   = (raw.get("sample_type")         or "").strip()
        mode = (raw.get("mode")                or "").strip()
        exp  = (raw.get("expected_taxa")        or "").strip()
        r1   = (raw.get("r1_fastq")            or "").strip()
        r2   = (raw.get("r2_fastq")            or "").strip()
        cmt  = (raw.get("analysis_comments")   or "").strip()

        if rid and rid not in run_names_seen:
            run_names_seen.append(rid)

        label = f"Row {i}" + (f" ({sid})" if sid else "")

        if st not in SAMPLE_TYPES:
            errors.append(f"{label}: unrecognized sample_type '{st}'. Expected one of: {', '.join(SAMPLE_TYPES)}")
        if mode not in ("routine", "validation", ""):
            errors.append(f"{label}: unrecognized mode '{mode}'. Expected 'routine' or 'validation'.")
        if mode == "validation" and not exp:
            errors.append(f"{label}: mode=validation but expected_taxa is empty.")
        if mode == "routine" and exp:
            errors.append(f"{label}: mode=routine but expected_taxa is '{exp}'. This value will be cleared on export.")

        rows.append({
            "sample_id":         sid,
            "run_id":            rid,
            "sample_type":       st,
            "mode":              mode,
            "expected_taxa":     exp,
            "r1":                r1,
            "r2":                r2,
            "analysis_comments": cmt,
        })

    sids = [r["sample_id"] for r in rows]
    dupes = sorted({s for s in sids if sids.count(s) > 1 and s})
    if dupes:
        errors.append(f"Duplicate sample_ids: {', '.join(dupes)}")

    analysis_date = ""
    initials = ""
    if rows:
        cmt = rows[0]["analysis_comments"]
        dm = re.search(r"(\d{4}-\d{2}-\d{2})", cmt)
        if dm:
            analysis_date = dm.group(1)
            after = cmt[dm.end():].lstrip("_")
            initials = after.split("_")[0] if after else ""

    return {
        "table_name":    table_name,
        "analysis_date": analysis_date,
        "initials":      initials,
        "run_names":     run_names_seen,
        "rows":          rows,
        "errors":        errors,
    }


# ---------------------------------------------------------------------------
# APHL Header bar widget
# ---------------------------------------------------------------------------

class APHLHeader(QWidget):
    """Branded teal header bar with APHL logo and app title."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedHeight(52)
        self.setStyleSheet("background-color: #006E79;")

        layout = QHBoxLayout(self)
        layout.setContentsMargins(20, 0, 20, 0)
        layout.setSpacing(14)

        # Logo — try several candidate paths for dev + all bundle types
        logo_pix = None
        for candidate in [
            _bundle_path("assets/aphl-logo-white.png"),
            Path(__file__).parent / "assets" / "aphl-logo-white.png",
            Path(sys.executable).parent / "assets" / "aphl-logo-white.png",
        ]:
            if candidate.exists():
                logo_pix = QPixmap(str(candidate))
                break

        if logo_pix and not logo_pix.isNull():
            logo_lbl = QLabel()
            logo_lbl.setPixmap(
                logo_pix.scaledToHeight(28, Qt.TransformationMode.SmoothTransformation)
            )
            logo_lbl.setStyleSheet("background: transparent;")
            layout.addWidget(logo_lbl)
        else:
            fallback = QLabel("APHL")
            fallback.setStyleSheet(
                "color: white; font-weight: 800; font-size: 18px; "
                "letter-spacing: -1px; background: transparent;"
            )
            layout.addWidget(fallback)

        # Vertical separator
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.VLine)
        sep.setFixedHeight(30)
        sep.setStyleSheet("color: rgba(255,255,255,0.3); background: rgba(255,255,255,0.3);")
        layout.addWidget(sep)

        # Title block
        title_block = QVBoxLayout()
        title_block.setSpacing(1)

        title = QLabel("AFI Terra Sheet Builder")
        title.setStyleSheet(
            "color: #ffffff; font-weight: 700; font-size: 14px; "
            "letter-spacing: -0.3px; background: transparent;"
        )
        title_block.addWidget(title)

        subtitle = QLabel("AFI 16S Batch Pipeline · APHL")
        subtitle.setStyleSheet(
            "color: rgba(255,255,255,0.65); font-size: 11px; background: transparent;"
        )
        title_block.addWidget(subtitle)

        layout.addLayout(title_block)
        layout.addStretch()

        # Version chip
        ver = QLabel("v2.0")
        ver.setStyleSheet(
            "color: rgba(255,255,255,0.75); font-size: 11px; "
            "background: rgba(255,255,255,0.15); padding: 2px 8px;"
        )
        layout.addWidget(ver)


# ---------------------------------------------------------------------------
# Screen 1
# ---------------------------------------------------------------------------

class Screen1(QWidget):
    """Three-part: run count → run names → per-run folder selection."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._run_name_edits: list[QLineEdit] = []
        self._run_folder_sections: list[dict] = []
        self._run_names: list[str] = []

        self._main_layout = QVBoxLayout(self)
        self._main_layout.setSpacing(12)
        self._main_layout.setContentsMargins(16, 16, 16, 16)

        # ── Part A: how many runs ──
        self._part_a = QGroupBox("Step 1 — How many sequencing runs are in this batch?")
        a_layout = QHBoxLayout(self._part_a)
        a_layout.addWidget(QLabel("Number of runs:"))
        self._run_count = QSpinBox()
        self._run_count.setMinimum(1)
        self._run_count.setMaximum(50)
        self._run_count.setValue(1)
        a_layout.addWidget(self._run_count)
        self._btn_continue_a = QPushButton("Continue →")
        self._btn_continue_a.clicked.connect(self._show_part_b)
        a_layout.addWidget(self._btn_continue_a)
        a_layout.addStretch()
        self._main_layout.addWidget(self._part_a)

        # ── Part B: run name entry ──
        self._part_b = QGroupBox("Step 2 — Enter a name for each run")
        self._part_b.hide()
        self._b_form = QFormLayout(self._part_b)
        self._btn_next_b = QPushButton("Next →")
        self._btn_next_b.clicked.connect(self._show_part_c)
        self._main_layout.addWidget(self._part_b)

        # ── Part C: per-run folder selection ──
        self._part_c_outer = QGroupBox("Step 3 — Select FASTQ files for each run")
        self._part_c_outer.hide()
        c_outer_layout = QVBoxLayout(self._part_c_outer)

        # Global TSV import row
        tsv_row = QHBoxLayout()
        lbl_tsv = QLabel("Import all runs at once:")
        btn_tsv = QPushButton("Import from TSV…")
        btn_tsv.setProperty("secondary", True)
        btn_tsv.setToolTip(
            "TSV must have a 'run_id' column and a 'sample_id' column.\n"
            "Optionally include 'r1_fastq' / 'r2_fastq' columns for direct paths,\n"
            "or choose a FASTQ folder for fuzzy filename matching."
        )
        btn_tsv.clicked.connect(self._import_all_from_tsv)
        tsv_row.addWidget(lbl_tsv)
        tsv_row.addWidget(btn_tsv)
        tsv_row.addStretch()
        c_outer_layout.addLayout(tsv_row)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setFrameShadow(QFrame.Shadow.Sunken)
        c_outer_layout.addWidget(sep)

        self._part_c_scroll = QScrollArea()
        self._part_c_scroll.setWidgetResizable(True)
        self._part_c_inner = QWidget()
        self._part_c_layout = QVBoxLayout(self._part_c_inner)
        self._part_c_scroll.setWidget(self._part_c_inner)
        c_outer_layout.addWidget(self._part_c_scroll)

        self._btn_next_c = QPushButton("Next →  (proceed to sample metadata)")
        self._btn_next_c.clicked.connect(self._go_to_screen2)
        c_outer_layout.addWidget(self._btn_next_c)
        self._main_layout.addWidget(self._part_c_outer)

        self._main_layout.addStretch()

    # ── Part A → B ──────────────────────────────────────────────────────────

    def _show_part_b(self):
        n = self._run_count.value()
        for i in reversed(range(self._b_form.rowCount())):
            self._b_form.removeRow(i)
        self._run_name_edits.clear()
        for i in range(n):
            edit = QLineEdit()
            edit.setPlaceholderText(f"e.g. Run{i+1}")
            edit.setText(f"Run{i+1}")
            self._run_name_edits.append(edit)
            self._b_form.addRow(f"Run {i+1} name:", edit)
        self._b_form.addRow("", self._btn_next_b)
        self._part_b.show()

    # ── Part B → C ──────────────────────────────────────────────────────────

    def _show_part_c(self):
        names = [e.text().strip() for e in self._run_name_edits]
        if any(not n for n in names):
            QMessageBox.warning(self, "Validation", "All run names must be non-empty.")
            return
        if len(set(names)) != len(names):
            QMessageBox.warning(self, "Validation", "Run names must be unique.")
            return
        if any(re.search(r"\s", n) for n in names):
            QMessageBox.warning(self, "Validation", "Run names must not contain whitespace.")
            return
        self._run_names = names
        while self._part_c_layout.count():
            item = self._part_c_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self._run_folder_sections.clear()
        for name in names:
            section = self._build_run_section(name)
            self._part_c_layout.addWidget(section["widget"])
        self._part_c_layout.addStretch()
        self._part_c_outer.show()

    def _build_run_section(self, run_name: str) -> dict:
        box = QGroupBox(f"Run: {run_name}")
        box_layout = QVBoxLayout(box)

        table = QTableWidget(0, 4)
        table.setHorizontalHeaderLabels(["Sample name", "R1 file", "R2 file", ""])
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        table.setAlternatingRowColors(True)
        table.verticalHeader().setDefaultSectionSize(28)
        table.setMinimumHeight(120)

        btn_layout = QHBoxLayout()
        btn_browse = QPushButton("Browse folder…")
        btn_browse.setProperty("secondary", True)
        btn_add = QPushButton("Add files…")
        btn_add.setProperty("secondary", True)
        btn_layout.addWidget(btn_browse)
        btn_layout.addWidget(btn_add)
        btn_layout.addStretch()
        box_layout.addLayout(btn_layout)
        box_layout.addWidget(table)

        section = {"widget": box, "run_name": run_name, "table": table}
        self._run_folder_sections.append(section)

        btn_browse.clicked.connect(lambda checked=False, s=section: self._browse_folder(s))
        btn_add.clicked.connect(lambda checked=False, s=section: self._add_files_multi(s))

        return section

    def _browse_folder(self, section: dict):
        folder = QFileDialog.getExistingDirectory(self, "Select FASTQ folder")
        if not folder:
            return
        pairs = discover_pairs(folder)
        if not pairs:
            QMessageBox.information(
                self, "No pairs found",
                "No R1/R2 FASTQ pairs found in the selected folder.\n"
                "Expected filenames containing _R1_ and _R2_."
            )
            return
        for stem, r1, r2 in pairs:
            self._add_row(section, stem, r1, r2)

    def _add_files_multi(self, section: dict):
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select FASTQ files (R1 and/or R2)", "",
            "FASTQ files (*.fastq *.fastq.gz *.fq *.fq.gz)",
        )
        if not files:
            return
        pairs = pair_fastqs_from_files(files)
        if not pairs:
            if len(files) == 2:
                f0, f1 = Path(files[0]), Path(files[1])
                if R1_PATTERNS.search(f0.name):
                    r1, r2 = files[0], files[1]
                elif R1_PATTERNS.search(f1.name):
                    r1, r2 = files[1], files[0]
                else:
                    r1, r2 = files[0], files[1]
                stem = R1_PATTERNS.split(Path(r1).name)[0] or Path(r1).name
                self._add_row(section, stem, r1, r2)
            else:
                QMessageBox.warning(
                    self, "No pairs found",
                    "Could not detect R1/R2 pairs from the selected files.\n"
                    "Expected filenames containing _R1_ or _R2_."
                )
            return
        for stem, r1, r2 in pairs:
            self._add_row(section, stem, r1, r2)

    def _add_row(self, section: dict, sample_id: str, r1: str, r2: str):
        table: QTableWidget = section["table"]
        row = table.rowCount()
        table.insertRow(row)
        table.setItem(row, 0, QTableWidgetItem(sample_id))
        table.setItem(row, 1, QTableWidgetItem(os.path.basename(r1)))
        table.setItem(row, 2, QTableWidgetItem(os.path.basename(r2)))
        btn_rm = QPushButton("✕")
        btn_rm.setProperty("iconOnly", True)
        btn_rm.setFixedWidth(32)

        def _make_remover(t: QTableWidget, b: QPushButton):
            def remove():
                for r in range(t.rowCount()):
                    if t.cellWidget(r, 3) is b:
                        t.removeRow(r)
                        break
            return remove

        btn_rm.clicked.connect(_make_remover(table, btn_rm))
        table.setCellWidget(row, 3, btn_rm)

    # ── TSV global import ────────────────────────────────────────────────────

    def _import_all_from_tsv(self):
        tsv_path, _ = QFileDialog.getOpenFileName(
            self, "Select mapping TSV", "", "TSV / CSV files (*.tsv *.csv *.txt)"
        )
        if not tsv_path:
            return
        try:
            _, rows = parse_mapping_tsv(tsv_path)
        except Exception as exc:
            QMessageBox.critical(self, "TSV parse error", str(exc))
            return
        if not rows:
            QMessageBox.warning(self, "Empty TSV", "No data rows found in the file.")
            return

        has_paths = any(r["r1"] for r in rows)
        folder_candidates: list[Path] = []
        if not has_paths:
            folder = QFileDialog.getExistingDirectory(
                self, "Select FASTQ folder for fuzzy matching"
            )
            if not folder:
                return
            folder_candidates = [f for f in Path(folder).iterdir() if is_fastq(f)]
            if not folder_candidates:
                QMessageBox.warning(self, "No FASTQ files",
                                    "No FASTQ files found in the selected folder.")
                return

        section_map = {s["run_name"]: s for s in self._run_folder_sections}
        unmatched: list[str] = []
        unknown_runs: set[str] = set()
        matched = 0

        for row in rows:
            run_id = row["run_id"]
            sample_id = row["sample_id"]
            if not run_id or not sample_id:
                continue
            section = section_map.get(run_id)
            if section is None:
                unknown_runs.add(run_id)
                continue
            if has_paths and row["r1"] and row["r2"]:
                self._add_row(section, sample_id, row["r1"], row["r2"])
                matched += 1
            elif has_paths and row["r1"]:
                r1 = Path(row["r1"])
                r2_name = R1_PATTERNS.sub(lambda m: m.group().replace("R1", "R2"), r1.name)
                r2 = r1.parent / r2_name
                if r2.exists():
                    self._add_row(section, sample_id, str(r1), str(r2))
                    matched += 1
                else:
                    unmatched.append(f"{run_id}/{sample_id} (R2 not found)")
            else:
                r1_path, r2_path = fuzzy_match_sample(sample_id, folder_candidates)
                if r1_path and r2_path:
                    self._add_row(section, sample_id, str(r1_path), str(r2_path))
                    matched += 1
                else:
                    unmatched.append(f"{run_id}/{sample_id}")

        lines = [f"Imported {matched} sample(s)."]
        if unknown_runs:
            lines.append("\nRun IDs not found in this batch:\n  " + "\n  ".join(sorted(unknown_runs)))
        if unmatched:
            lines.append("\nNo FASTQ match found for:\n  " + "\n  ".join(unmatched))
        if unknown_runs or unmatched:
            QMessageBox.warning(self, "Import complete with warnings", "\n".join(lines))
        else:
            QMessageBox.information(self, "Import complete", "\n".join(lines))

    # ── Part C → Screen 2 ───────────────────────────────────────────────────

    def _go_to_screen2(self):
        all_rows = []
        for section in self._run_folder_sections:
            table: QTableWidget = section["table"]
            run_name = section["run_name"]
            for r in range(table.rowCount()):
                def cell(col: int) -> str:
                    item = table.item(r, col)
                    return item.text().strip() if item else ""
                sample_id = cell(0)
                r1 = cell(1)
                r2 = cell(2)
                if sample_id and r1 and r2:
                    all_rows.append({
                        "sample_id": sample_id,
                        "run_id":    run_name,
                        "r1":        r1,
                        "r2":        r2,
                    })
        if not all_rows:
            QMessageBox.warning(self, "No samples",
                                "Please add at least one sample before continuing.")
            return
        self.window().go_to_screen2(self._run_names, all_rows)


# ---------------------------------------------------------------------------
# ComboBox delegate
# ---------------------------------------------------------------------------

class ComboDelegate(QStyledItemDelegate):
    def __init__(self, choices: list[str], parent=None):
        super().__init__(parent)
        self._choices = choices

    def createEditor(self, parent, option, index):
        combo = QComboBox(parent)
        combo.addItems(self._choices)
        return combo

    def setEditorData(self, editor, index):
        val = index.data(Qt.ItemDataRole.EditRole) or ""
        idx = editor.findText(val)
        if idx >= 0:
            editor.setCurrentIndex(idx)

    def setModelData(self, editor, model, index):
        model.setData(index, editor.currentText(), Qt.ItemDataRole.EditRole)

    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)


# ---------------------------------------------------------------------------
# Screen 2
# ---------------------------------------------------------------------------

class Screen2(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._run_names: list[str] = []
        self._rows: list[dict] = []
        self._validated = False

        layout = QVBoxLayout(self)
        layout.setSpacing(10)
        layout.setContentsMargins(16, 16, 16, 16)

        # ── Header fields ──
        header_box = QGroupBox("Run metadata")
        form = QFormLayout(header_box)

        self._date_edit = QDateEdit()
        self._date_edit.setCalendarPopup(True)
        self._date_edit.setDate(QDate.currentDate())
        form.addRow("Analysis date:", self._date_edit)

        self._table_name_edit = QLineEdit()
        self._table_name_edit.setPlaceholderText("e.g. afi_run_1  (lowercase, no spaces)")
        self._table_name_edit.textChanged.connect(self._validate_table_name_live)
        form.addRow("Data table name:", self._table_name_edit)

        self._initials_edit = QLineEdit()
        self._initials_edit.setPlaceholderText("e.g. JS")
        form.addRow("Operator initials:", self._initials_edit)

        layout.addWidget(header_box)

        # ── Sample metadata table ──
        self._table = QTableWidget(0, NUM_COLS)
        self._table.setHorizontalHeaderLabels(COL_HEADERS)
        self._table.horizontalHeader().setSectionResizeMode(COL_SAMPLE_ID, QHeaderView.ResizeMode.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(COL_RUN_ID, QHeaderView.ResizeMode.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(COL_SAMPLE_TYPE, QHeaderView.ResizeMode.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(COL_MODE, QHeaderView.ResizeMode.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(COL_EXPECTED, QHeaderView.ResizeMode.Stretch)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self._table.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked |
            QAbstractItemView.EditTrigger.EditKeyPressed
        )
        self._table.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self._table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self._table.setAlternatingRowColors(True)
        self._table.verticalHeader().setDefaultSectionSize(28)
        self._table.itemChanged.connect(self._on_item_changed)
        layout.addWidget(self._table)

        hint = QLabel("Tip: double-click a cell (or select + press F2) to edit it.")
        hint.setObjectName("hint")
        layout.addWidget(hint)

        # ── Buttons ──
        btn_row = QHBoxLayout()
        self._btn_back = QPushButton("← Back")
        self._btn_back.setProperty("secondary", True)
        self._btn_back.clicked.connect(self._go_back)

        self._btn_load = QPushButton("📂 Load existing TSV…")
        self._btn_load.setProperty("secondary", True)
        self._btn_load.setToolTip("Load a previously-exported Terra TSV for editing.")
        self._btn_load.clicked.connect(self._load_existing_tsv)

        self._btn_validate = QPushButton("✓ Validate")
        self._btn_validate.clicked.connect(self._validate)

        self._btn_export = QPushButton("⬇ Export TSV")
        self._btn_export.setEnabled(False)
        self._btn_export.clicked.connect(self._export)

        btn_row.addWidget(self._btn_back)
        btn_row.addWidget(self._btn_load)
        btn_row.addStretch()
        btn_row.addWidget(self._btn_validate)
        btn_row.addWidget(self._btn_export)
        layout.addLayout(btn_row)

        self._blocking_itemChanged = False

    # ── Populate ─────────────────────────────────────────────────────────────

    def populate(self, run_names: list[str], rows: list[dict]):
        self._run_names = run_names
        self._rows = rows
        self._validated = False
        self._btn_export.setEnabled(False)

        run_delegate  = ComboDelegate(run_names, self._table)
        type_delegate = ComboDelegate(SAMPLE_TYPES, self._table)
        self._table.setItemDelegateForColumn(COL_RUN_ID, run_delegate)
        self._table.setItemDelegateForColumn(COL_SAMPLE_TYPE, type_delegate)

        seen: set[str] = set()
        self._blocking_itemChanged = True
        self._table.setRowCount(0)
        for row_data in rows:
            sid = unique_sample_id(row_data["sample_id"], seen)
            seen.add(sid)
            r = self._table.rowCount()
            self._table.insertRow(r)
            self._set_cell(r, COL_SAMPLE_ID, sid,
                           file_data={"r1": row_data["r1"], "r2": row_data["r2"]})
            self._set_cell(r, COL_RUN_ID, row_data["run_id"])
            self._set_cell(r, COL_SAMPLE_TYPE, "clinical")
            self._set_cell(r, COL_MODE, AUTO_FILL["clinical"][0])
            self._set_cell(r, COL_EXPECTED, AUTO_FILL["clinical"][1])
            self._apply_row_color(r, "clinical")
        self._blocking_itemChanged = False

    def _apply_row_color(self, row: int, sample_type: str):
        """Color-code rows by sample type using APHL palette."""
        if sample_type in NTC_TYPES:
            color = COLOR_NTC
        elif sample_type in PC_TYPES:
            color = COLOR_PC
        else:
            color = COLOR_CLN
        for col in range(NUM_COLS):
            item = self._table.item(row, col)
            if item:
                item.setBackground(QBrush(color))

    def _set_cell(self, row: int, col: int, value: str, file_data: dict | None = None):
        item = QTableWidgetItem(value)
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)
        if col == COL_SAMPLE_ID and file_data is not None:
            item.setData(Qt.ItemDataRole.UserRole, file_data)
        self._table.setItem(row, col, item)

    # ── Live validation ───────────────────────────────────────────────────────

    def _validate_table_name_live(self, text: str):
        valid = bool(re.match(r"^[a-z][a-z0-9_]{0,31}$", text))
        self._table_name_edit.setProperty("invalid", not valid and bool(text))
        self._table_name_edit.style().unpolish(self._table_name_edit)
        self._table_name_edit.style().polish(self._table_name_edit)

    def _on_item_changed(self, item: QTableWidgetItem):
        if self._blocking_itemChanged:
            return
        col = item.column()
        row = item.row()
        if col == COL_SAMPLE_TYPE:
            st = item.text().strip()
            if st in AUTO_FILL:
                mode, expected = AUTO_FILL[st]
                self._blocking_itemChanged = True
                self._set_cell(row, COL_MODE, mode)
                self._set_cell(row, COL_EXPECTED, expected)
                self._blocking_itemChanged = False
            self._apply_row_color(row, st)
        self._validated = False
        self._btn_export.setEnabled(False)

    # ── Validate ──────────────────────────────────────────────────────────────

    def _validate(self):
        errors: list[str] = []

        table_name = self._table_name_edit.text().strip()
        if not re.match(r"^[a-z][a-z0-9_]{0,31}$", table_name):
            errors.append(
                "Data table name must start with a lowercase letter, "
                "contain only a-z, 0-9, and _, and be at most 32 characters."
            )
        if not self._initials_edit.text().strip():
            errors.append("Operator initials must not be empty.")

        n = self._table.rowCount()
        sample_ids, run_ids, sample_types, modes, expected_taxa = [], [], [], [], []
        for r in range(n):
            sample_ids.append(self._cell_text(r, COL_SAMPLE_ID))
            run_ids.append(self._cell_text(r, COL_RUN_ID))
            sample_types.append(self._cell_text(r, COL_SAMPLE_TYPE))
            modes.append(self._cell_text(r, COL_MODE))
            expected_taxa.append(self._cell_text(r, COL_EXPECTED))

        seen: set[str] = set()
        dupes: set[str] = set()
        for sid in sample_ids:
            if sid in seen:
                dupes.add(sid)
            seen.add(sid)
        if dupes:
            errors.append(f"Duplicate sample_ids: {', '.join(sorted(dupes))}")

        from collections import defaultdict
        run_ntc: dict[str, bool] = defaultdict(bool)
        run_pc:  dict[str, bool] = defaultdict(bool)
        for rid, st in zip(run_ids, sample_types):
            if st in NTC_TYPES:
                run_ntc[rid] = True
            if st in PC_TYPES:
                run_pc[rid] = True
        for rname in self._run_names:
            if not run_ntc.get(rname):
                errors.append(f"Run '{rname}' has no NTC or NC sample.")
            if not run_pc.get(rname):
                errors.append(f"Run '{rname}' has no positive control (PC_MIX8/PC_SINGLE/MIXED4/PC).")

        for r in range(n):
            if modes[r] == "validation" and not expected_taxa[r].strip():
                errors.append(
                    f"Row {r+1} ({sample_ids[r]}): mode=validation but expected_taxa is empty."
                )

        if errors:
            QMessageBox.warning(self, "Validation failed",
                                "\n\n".join(f"• {e}" for e in errors))
            return

        self._validated = True
        self._btn_export.setEnabled(True)
        QMessageBox.information(self, "Validation passed", "All checks passed. Ready to export.")

    def _cell_text(self, row: int, col: int) -> str:
        item = self._table.item(row, col)
        return item.text().strip() if item else ""

    # ── Export ────────────────────────────────────────────────────────────────

    def _export(self):
        if not self._validated:
            QMessageBox.warning(self, "Export", "Please validate first.")
            return

        table_name  = self._table_name_edit.text().strip()
        analysis_dt = self._date_edit.date().toString("yyyy-MM-dd")
        initials    = self._initials_edit.text().strip()
        comment     = f"{table_name}_{analysis_dt}_{initials}"

        save_path, _ = QFileDialog.getSaveFileName(
            self, "Save TSV", f"{table_name}.tsv", "TSV files (*.tsv)"
        )
        if not save_path:
            return

        n = self._table.rowCount()
        fieldnames = [
            f"entity:{table_name}_id",
            "run_id", "sample_type", "mode", "expected_taxa",
            "r1_fastq", "r2_fastq", "analysis_comments",
        ]
        with open(save_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for r in range(n):
                sid  = self._cell_text(r, COL_SAMPLE_ID)
                rid  = self._cell_text(r, COL_RUN_ID)
                st   = self._cell_text(r, COL_SAMPLE_TYPE)
                mode = self._cell_text(r, COL_MODE)
                exp  = self._cell_text(r, COL_EXPECTED)
                if mode == "routine":
                    exp = ""
                file_data = {}
                id_item = self._table.item(r, COL_SAMPLE_ID)
                if id_item:
                    file_data = id_item.data(Qt.ItemDataRole.UserRole) or {}
                r1 = file_data.get("r1", "")
                r2 = file_data.get("r2", "")
                writer.writerow({
                    f"entity:{table_name}_id": sid,
                    "run_id":            rid,
                    "sample_type":       st,
                    "mode":              mode,
                    "expected_taxa":     exp,
                    "r1_fastq":          r1,
                    "r2_fastq":          r2,
                    "analysis_comments": comment,
                })

        QMessageBox.information(self, "Exported", f"Saved {n} rows to:\n{save_path}")

    # ── Load existing TSV ─────────────────────────────────────────────────────

    def _load_existing_tsv(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load existing Terra TSV", "",
            "TSV files (*.tsv *.txt);;All files (*)"
        )
        if not path:
            return
        try:
            data = parse_terra_tsv(path)
        except Exception as exc:
            QMessageBox.critical(self, "Load error", f"Could not parse TSV:\n{exc}")
            return
        if not data["rows"]:
            QMessageBox.warning(self, "Empty TSV", "No data rows found in the file.")
            return
        if data["errors"]:
            msg = (
                "The following issues were found in the TSV:\n\n"
                + "\n".join(f"  • {e}" for e in data["errors"])
                + "\n\nYou can still load and fix these in the editor.\nProceed?"
            )
            reply = QMessageBox.warning(
                self, "TSV issues found", msg,
                QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Ok,
            )
            if reply != QMessageBox.StandardButton.Ok:
                return

        if data["table_name"]:
            self._table_name_edit.setText(data["table_name"])
        if data["analysis_date"]:
            qd = QDate.fromString(data["analysis_date"], "yyyy-MM-dd")
            if qd.isValid():
                self._date_edit.setDate(qd)
        if data["initials"]:
            self._initials_edit.setText(data["initials"])

        run_names = data["run_names"]
        self._run_names = run_names
        self._table.setItemDelegateForColumn(COL_RUN_ID, ComboDelegate(run_names, self._table))
        self._table.setItemDelegateForColumn(COL_SAMPLE_TYPE, ComboDelegate(SAMPLE_TYPES, self._table))

        self._blocking_itemChanged = True
        self._table.setRowCount(0)
        self._rows = []
        for r_data in data["rows"]:
            self._rows.append({
                "sample_id": r_data["sample_id"],
                "run_id":    r_data["run_id"],
                "r1":        r_data["r1"],
                "r2":        r_data["r2"],
            })
            r = self._table.rowCount()
            self._table.insertRow(r)
            self._set_cell(r, COL_SAMPLE_ID, r_data["sample_id"],
                           file_data={"r1": r_data["r1"], "r2": r_data["r2"]})
            self._set_cell(r, COL_RUN_ID,     r_data["run_id"])
            self._set_cell(r, COL_SAMPLE_TYPE, r_data["sample_type"])
            self._set_cell(r, COL_MODE,        r_data["mode"])
            self._set_cell(r, COL_EXPECTED,    r_data["expected_taxa"])
            self._apply_row_color(r, r_data["sample_type"])
        self._blocking_itemChanged = False

        self._validated = False
        self._btn_export.setEnabled(False)

        n = len(data["rows"])
        suffix = (f"\n\nPlease fix {len(data['errors'])} issue(s) before exporting."
                  if data["errors"] else "")
        QMessageBox.information(self, "Loaded",
                                f"Loaded {n} row(s) from {Path(path).name}.{suffix}")

    # ── Back ──────────────────────────────────────────────────────────────────

    def _go_back(self):
        self.window().stack.setCurrentIndex(0)


# ---------------------------------------------------------------------------
# Main Window
# ---------------------------------------------------------------------------

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("AFI Terra Sheet Builder")
        self.resize(1100, 720)
        self.setMinimumSize(800, 560)

        # ── Window icon ──
        icon_path = _bundle_path("assets/aphl-icon.png")
        if icon_path.exists():
            self.setWindowIcon(QIcon(str(icon_path)))

        # ── Root widget: header + stacked screens ──
        root = QWidget()
        root_layout = QVBoxLayout(root)
        root_layout.setContentsMargins(0, 0, 0, 0)
        root_layout.setSpacing(0)

        root_layout.addWidget(APHLHeader())

        self.stack = QStackedWidget()
        root_layout.addWidget(self.stack)
        self.setCentralWidget(root)

        self._screen1 = Screen1(self)
        self._screen2 = Screen2(self)
        self.stack.addWidget(self._screen1)
        self.stack.addWidget(self._screen2)
        self.stack.setCurrentIndex(0)

    def go_to_screen2(self, run_names: list[str], rows: list[dict]):
        self._screen2.populate(run_names, rows)
        self.stack.setCurrentIndex(1)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Embedded fallback QSS — applied if aphl_style.qss is not found in bundle.
# Keeps branding intact even if the asset path resolution fails.
# ---------------------------------------------------------------------------

_FALLBACK_QSS = """
QMainWindow, QDialog { background-color: #E2E9EC; }
QWidget { font-family: Arial, sans-serif; font-size: 13px; color: #404040; }
QGroupBox { background-color: #ffffff; border: 1px solid #d0d9dc; margin-top: 18px; padding-top: 6px; }
QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left; left: 0px; top: 0px; background-color: #006E79; color: #ffffff; font-weight: bold; font-size: 12px; padding: 5px 14px; }
QLabel { color: #404040; background: transparent; }
QPushButton { background-color: #006E79; color: #ffffff; border: none; padding: 6px 18px; font-weight: bold; font-size: 13px; min-height: 28px; min-width: 80px; }
QPushButton:hover { background-color: #005057; }
QPushButton:pressed { background-color: #003d42; }
QPushButton:disabled { background-color: #A3CCCC; color: #ffffff; }
QPushButton[secondary="true"] { background-color: #ffffff; color: #006E79; border: 1.5px solid #006E79; }
QPushButton[secondary="true"]:hover { background-color: #f0f7f8; }
QPushButton[danger="true"] { background-color: #B42E34; color: #ffffff; border: none; }
QPushButton[danger="true"]:hover { background-color: #862226; }
QPushButton[iconOnly="true"] { background-color: transparent; color: #B42E34; border: none; padding: 2px 6px; min-width: 28px; min-height: 20px; font-size: 14px; font-weight: bold; }
QPushButton[iconOnly="true"]:hover { background-color: #f3d1d2; }
QLineEdit, QSpinBox, QDateEdit, QComboBox { border: 1px solid #c4d0d4; background-color: #ffffff; color: #404040; padding: 4px 8px; min-height: 26px; }
QLineEdit:focus, QSpinBox:focus, QDateEdit:focus, QComboBox:focus { border-color: #00A0AF; }
QLineEdit[invalid="true"] { border-color: #B42E34; }
QTableWidget { background-color: #ffffff; alternate-background-color: #f4f7f8; gridline-color: #E2E9EC; border: 1px solid #d0d9dc; selection-background-color: #e5f4f5; selection-color: #004048; }
QTableWidget::item { padding: 4px 8px; border: none; }
QHeaderView { background-color: #006E79; }
QHeaderView::section { background-color: #006E79; color: #ffffff; font-weight: bold; font-size: 11px; padding: 5px 10px; border: none; border-right: 1px solid #005057; text-transform: uppercase; }
QScrollBar:vertical { background-color: #E2E9EC; width: 10px; margin: 0; }
QScrollBar::handle:vertical { background-color: #A3CCCC; min-height: 20px; }
QScrollBar::handle:vertical:hover { background-color: #006E79; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
QScrollBar:horizontal { background-color: #E2E9EC; height: 10px; margin: 0; }
QScrollBar::handle:horizontal { background-color: #A3CCCC; min-width: 20px; }
QScrollBar::handle:horizontal:hover { background-color: #006E79; }
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }
QScrollArea { border: none; background-color: transparent; }
QComboBox QAbstractItemView { background-color: #ffffff; border: 1px solid #00A0AF; selection-background-color: #006E79; selection-color: #ffffff; }
QToolTip { background-color: #005057; color: #ffffff; border: 1px solid #006E79; padding: 4px 8px; }
"""


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("AFI Terra Sheet Builder")
    app.setOrganizationName("APHL")

    # ── Force Fusion style for consistent cross-platform QSS rendering ──
    # Without this, macOS/Windows native styles partially override QSS rules
    # (especially QGroupBox titles and input borders).
    app.setStyle("Fusion")

    # ── Apply APHL stylesheet to the application (cascades to all widgets) ──
    qss_path = _bundle_path("aphl_style.qss")
    if qss_path.exists():
        try:
            app.setStyleSheet(qss_path.read_text(encoding="utf-8"))
        except Exception:
            app.setStyleSheet(_FALLBACK_QSS)
    else:
        app.setStyleSheet(_FALLBACK_QSS)

    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
