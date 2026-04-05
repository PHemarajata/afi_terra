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
from PySide6.QtGui import QColor, QBrush, QFont
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
    # Build a lookup for R2 candidates
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
    """
    Find the best-matching R1/R2 pair for *sample_id* among *candidates*.

    Strategy (in order):
      1. Exact prefix match:  filename starts with sample_id + separator
      2. Substring match:     sample_id appears anywhere in the filename stem
      3. difflib fuzzy match: highest ratio among R1-candidate stems
    Returns (r1_path, r2_path) or (None, None) if nothing scores above *cutoff*.
    """
    r1_cands = [f for f in candidates if R1_PATTERNS.search(f.name)]
    if not r1_cands:
        return None, None

    sid_lower = sample_id.lower()

    def stem_of(f: Path) -> str:
        return R1_PATTERNS.split(f.name)[0].lower()

    # 1. Exact prefix
    for r1 in r1_cands:
        st = stem_of(r1)
        if st == sid_lower or st.startswith(sid_lower + "_") or st.startswith(sid_lower + "-"):
            r2 = _find_r2(r1, candidates)
            if r2:
                return r1, r2

    # 2. Substring
    for r1 in r1_cands:
        if sid_lower in stem_of(r1):
            r2 = _find_r2(r1, candidates)
            if r2:
                return r1, r2

    # 3. difflib
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
    """
    Parse a mapping TSV into (column_names, rows).

    Accepted column names (case-insensitive, flexible):
      run_id / run
      sample_id / sample_name / sample
      r1_fastq / r1 / fastq_r1 / read1
      r2_fastq / r2 / fastq_r2 / read2

    Returns normalised dicts with keys: run_id, sample_id, r1 (may be ""), r2 (may be "").
    """
    ALIASES = {
        "run_id":    {"run_id", "run", "run_name"},
        "sample_id": {"sample_id", "sample", "sample_name"},
        "r1":        {"r1_fastq", "r1", "fastq_r1", "read1", "r1_path"},
        "r2":        {"r2_fastq", "r2", "fastq_r2", "read2", "r2_path"},
    }

    with open(path, newline="", encoding="utf-8") as fh:
        # Sniff delimiter
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.DictReader(fh, dialect=dialect)
        raw_cols = reader.fieldnames or []

        col_map: dict[str, str] = {}   # raw header → canonical key
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


# ---------------------------------------------------------------------------
# Screen 1
# ---------------------------------------------------------------------------

class Screen1(QWidget):
    """Three-part: run count → run names → per-run folder selection."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._run_name_edits: list[QLineEdit] = []
        self._run_folder_sections: list[dict] = []   # {run_name, table, rows}
        self._run_names: list[str] = []

        self._main_layout = QVBoxLayout(self)
        self._main_layout.setSpacing(12)

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

        # ── Part B: run name entry (hidden until Part A done) ──
        self._part_b = QGroupBox("Step 2 — Enter a name for each run")
        self._part_b.hide()
        self._b_form = QFormLayout(self._part_b)
        self._btn_next_b = QPushButton("Next →")
        self._btn_next_b.clicked.connect(self._show_part_c)
        self._main_layout.addWidget(self._part_b)

        # ── Part C: per-run folder selection (hidden until Part B done) ──
        self._part_c_outer = QGroupBox("Step 3 — Select FASTQ files for each run")
        self._part_c_outer.hide()
        c_outer_layout = QVBoxLayout(self._part_c_outer)

        # Global TSV import row
        tsv_row = QHBoxLayout()
        lbl_tsv = QLabel("Import all runs at once:")
        btn_tsv = QPushButton("Import from TSV…")
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
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)
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

    # ── Part A → B ──

    def _show_part_b(self):
        n = self._run_count.value()
        # Rebuild Part B form
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

    # ── Part B → C ──

    def _show_part_c(self):
        names = [e.text().strip() for e in self._run_name_edits]
        # Validate
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
        # Rebuild Part C
        # Remove old sections
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

        # Table of discovered pairs
        table = QTableWidget(0, 4)
        table.setHorizontalHeaderLabels(["Sample name", "R1 file", "R2 file", "Remove"])
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        table.setMinimumHeight(120)

        btn_layout = QHBoxLayout()
        btn_browse = QPushButton("Browse folder…")
        btn_add    = QPushButton("Add files…")
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
        """Multi-select any number of FASTQ files; auto-pair by R1/R2 naming."""
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Select FASTQ files (R1 and/or R2)",
            "",
            "FASTQ files (*.fastq *.fastq.gz *.fq *.fq.gz)",
        )
        if not files:
            return
        pairs = pair_fastqs_from_files(files)
        if not pairs:
            # Fallback: if user selected files that don't match R1/R2 patterns,
            # let them pair one R1 + one R2 manually.
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
        table.setItem(row, 1, QTableWidgetItem(r1))
        table.setItem(row, 2, QTableWidgetItem(r2))
        # Remove button: look up own row index at click time so it survives
        # other rows being deleted above it.
        btn_rm = QPushButton("✕")
        btn_rm.setFixedWidth(30)
        def _make_remover(t: QTableWidget, b: QPushButton):
            def remove():
                for r in range(t.rowCount()):
                    if t.cellWidget(r, 3) is b:
                        t.removeRow(r)
                        break
            return remove
        btn_rm.clicked.connect(_make_remover(table, btn_rm))
        table.setCellWidget(row, 3, btn_rm)

    def _remove_row(self, table: QTableWidget, row: int):
        table.removeRow(row)

    # ── TSV global import ──

    def _import_all_from_tsv(self):
        """
        Import samples for ALL runs from a mapping TSV.

        Accepted columns (case-insensitive):
          run_id / run
          sample_id / sample_name / sample
          r1_fastq / r1 / fastq_r1 / read1     (optional)
          r2_fastq / r2 / fastq_r2 / read2     (optional)

        If r1/r2 columns are absent or empty the user is asked to choose
        a FASTQ folder and sample names are matched to filenames with fuzzy
        matching (exact prefix → substring → difflib).
        """
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

        # Determine whether direct paths or fuzzy matching needed
        has_paths = any(r["r1"] for r in rows)

        folder_candidates: list[Path] = []
        if not has_paths:
            folder = QFileDialog.getExistingDirectory(
                self, "Select FASTQ folder for fuzzy matching"
            )
            if not folder:
                return
            folder_candidates = [
                f for f in Path(folder).iterdir() if is_fastq(f)
            ]
            if not folder_candidates:
                QMessageBox.warning(
                    self, "No FASTQ files",
                    "No FASTQ files found in the selected folder."
                )
                return

        # Build a mapping run_name → section
        section_map = {s["run_name"]: s for s in self._run_folder_sections}

        unmatched: list[str] = []
        unknown_runs: set[str] = set()
        matched = 0

        for row in rows:
            run_id    = row["run_id"]
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
                # Only R1 supplied — try to derive R2
                r1 = Path(row["r1"])
                r2_name = R1_PATTERNS.sub(
                    lambda m: m.group().replace("R1", "R2"), r1.name
                )
                r2 = r1.parent / r2_name
                if r2.exists():
                    self._add_row(section, sample_id, str(r1), str(r2))
                    matched += 1
                else:
                    unmatched.append(f"{run_id}/{sample_id} (R2 not found)")
            else:
                # Fuzzy match against folder
                r1_path, r2_path = fuzzy_match_sample(sample_id, folder_candidates)
                if r1_path and r2_path:
                    self._add_row(section, sample_id, str(r1_path), str(r2_path))
                    matched += 1
                else:
                    unmatched.append(f"{run_id}/{sample_id}")

        # Report
        lines = [f"Imported {matched} sample(s)."]
        if unknown_runs:
            lines.append(
                f"\nRun IDs in TSV not found in this batch "
                f"(check spelling):\n  " + "\n  ".join(sorted(unknown_runs))
            )
        if unmatched:
            lines.append(
                f"\nNo FASTQ match found for:\n  " + "\n  ".join(unmatched)
            )
        if unknown_runs or unmatched:
            QMessageBox.warning(self, "Import complete with warnings", "\n".join(lines))
        else:
            QMessageBox.information(self, "Import complete", "\n".join(lines))

    # ── Part C → Screen 2 ──

    def _go_to_screen2(self):
        # Collect all rows from all section tables
        all_rows = []
        for section in self._run_folder_sections:
            table: QTableWidget = section["table"]
            run_name = section["run_name"]
            for r in range(table.rowCount()):
                def cell(col: int) -> str:
                    item = table.item(r, col)
                    return item.text().strip() if item else ""
                sample_id = cell(0)
                r1        = cell(1)
                r2        = cell(2)
                if sample_id and r1 and r2:
                    all_rows.append({
                        "sample_id": sample_id,
                        "run_id":    run_name,
                        "r1":        r1,
                        "r2":        r2,
                    })
        if not all_rows:
            QMessageBox.warning(
                self, "No samples",
                "Please add at least one sample before continuing."
            )
            return
        # Find parent MainWindow and navigate
        main = self.window()
        main.go_to_screen2(self._run_names, all_rows)


# ---------------------------------------------------------------------------
# ComboBox delegate for sample_type / run_id columns
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
        val = index.data(Qt.EditRole) or ""
        idx = editor.findText(val)
        if idx >= 0:
            editor.setCurrentIndex(idx)

    def setModelData(self, editor, model, index):
        model.setData(index, editor.currentText(), Qt.EditRole)

    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)


# ---------------------------------------------------------------------------
# Screen 2
# ---------------------------------------------------------------------------

class Screen2(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._run_names: list[str] = []
        self._rows: list[dict] = []   # {"sample_id", "run_id", "r1", "r2"}
        self._validated = False

        layout = QVBoxLayout(self)
        layout.setSpacing(10)

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
        self._table_name_label = QLabel("Data table name:")
        form.addRow(self._table_name_label, self._table_name_edit)

        self._initials_edit = QLineEdit()
        self._initials_edit.setPlaceholderText("e.g. JS")
        form.addRow("Operator initials:", self._initials_edit)

        layout.addWidget(header_box)

        # ── Sample metadata table ──
        self._table = QTableWidget(0, NUM_COLS)
        self._table.setHorizontalHeaderLabels(COL_HEADERS)
        self._table.horizontalHeader().setSectionResizeMode(COL_SAMPLE_ID, QHeaderView.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(COL_RUN_ID, QHeaderView.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(COL_SAMPLE_TYPE, QHeaderView.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(COL_MODE, QHeaderView.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(COL_EXPECTED, QHeaderView.Stretch)
        self._table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self._table.setSelectionBehavior(QAbstractItemView.SelectItems)
        self._table.itemChanged.connect(self._on_item_changed)
        layout.addWidget(self._table)

        # ── Buttons ──
        btn_row = QHBoxLayout()
        self._btn_back = QPushButton("← Back")
        self._btn_back.clicked.connect(self._go_back)
        self._btn_validate = QPushButton("✓ Validate")
        self._btn_validate.clicked.connect(self._validate)
        self._btn_export = QPushButton("⬇ Export TSV")
        self._btn_export.setEnabled(False)
        self._btn_export.clicked.connect(self._export)
        btn_row.addWidget(self._btn_back)
        btn_row.addStretch()
        btn_row.addWidget(self._btn_validate)
        btn_row.addWidget(self._btn_export)
        layout.addLayout(btn_row)

        self._blocking_itemChanged = False

    # ── Populate from Screen 1 data ──

    def populate(self, run_names: list[str], rows: list[dict]):
        self._run_names = run_names
        self._rows = rows
        self._validated = False
        self._btn_export.setEnabled(False)

        # Install delegates
        run_delegate = ComboDelegate(run_names, self._table)
        type_delegate = ComboDelegate(SAMPLE_TYPES, self._table)
        self._table.setItemDelegateForColumn(COL_RUN_ID, run_delegate)
        self._table.setItemDelegateForColumn(COL_SAMPLE_TYPE, type_delegate)

        # Deduplicate sample_ids
        seen: set[str] = set()
        self._blocking_itemChanged = True
        self._table.setRowCount(0)
        for row_data in rows:
            sid = unique_sample_id(row_data["sample_id"], seen)
            seen.add(sid)
            r = self._table.rowCount()
            self._table.insertRow(r)
            self._set_cell(r, COL_SAMPLE_ID, sid)
            self._set_cell(r, COL_RUN_ID, row_data["run_id"])
            self._set_cell(r, COL_SAMPLE_TYPE, "clinical")
            self._set_cell(r, COL_MODE, AUTO_FILL["clinical"][0])
            self._set_cell(r, COL_EXPECTED, AUTO_FILL["clinical"][1])
        self._blocking_itemChanged = False

    def _set_cell(self, row: int, col: int, value: str):
        item = QTableWidgetItem(value)
        if col in (COL_SAMPLE_ID, COL_EXPECTED):
            item.setFlags(item.flags() | Qt.ItemIsEditable)
        elif col in (COL_RUN_ID, COL_SAMPLE_TYPE, COL_MODE):
            item.setFlags(item.flags() | Qt.ItemIsEditable)
        self._table.setItem(row, col, item)

    # ── Live table name validation ──

    def _validate_table_name_live(self, text: str):
        valid = bool(re.match(r"^[a-z][a-z0-9_]{0,31}$", text))
        self._table_name_edit.setStyleSheet(
            "" if valid or not text else "border: 2px solid red;"
        )

    # ── Handle cell changes: auto-fill mode/expected on sample_type change ──

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
        self._validated = False
        self._btn_export.setEnabled(False)

    # ── Validate ──

    def _validate(self):
        errors: list[str] = []

        # Header checks
        table_name = self._table_name_edit.text().strip()
        if not re.match(r"^[a-z][a-z0-9_]{0,31}$", table_name):
            errors.append(
                "Data table name must start with a lowercase letter, contain only "
                "a-z, 0-9, and _, and be at most 32 characters."
            )
        if not self._initials_edit.text().strip():
            errors.append("Operator initials must not be empty.")

        # Collect table data
        n = self._table.rowCount()
        sample_ids: list[str] = []
        run_ids: list[str] = []
        sample_types: list[str] = []
        modes: list[str] = []
        expected_taxa: list[str] = []

        for r in range(n):
            sample_ids.append(self._cell_text(r, COL_SAMPLE_ID))
            run_ids.append(self._cell_text(r, COL_RUN_ID))
            sample_types.append(self._cell_text(r, COL_SAMPLE_TYPE))
            modes.append(self._cell_text(r, COL_MODE))
            expected_taxa.append(self._cell_text(r, COL_EXPECTED))

        # Duplicate sample_ids
        seen: set[str] = set()
        dupes: set[str] = set()
        for sid in sample_ids:
            if sid in seen:
                dupes.add(sid)
            seen.add(sid)
        if dupes:
            errors.append(f"Duplicate sample_ids: {', '.join(sorted(dupes))}")

        # Per-run checks: needs ≥1 NTC/NC and ≥1 positive control
        from collections import defaultdict
        run_ntc: dict[str, bool] = defaultdict(bool)
        run_pc:  dict[str, bool] = defaultdict(bool)
        PC_TYPES = {"PC_MIX8", "PC_SINGLE", "MIXED4", "PC"}
        NTC_TYPES = {"NTC", "NC"}
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

        # Validation-mode rows must have non-empty expected_taxa
        for r in range(n):
            if modes[r] == "validation" and not expected_taxa[r].strip():
                errors.append(
                    f"Row {r+1} ({sample_ids[r]}): mode=validation but expected_taxa is empty."
                )

        if errors:
            QMessageBox.warning(
                self, "Validation failed",
                "\n\n".join(f"• {e}" for e in errors)
            )
            return

        self._validated = True
        self._btn_export.setEnabled(True)
        QMessageBox.information(self, "Validation passed", "All checks passed.")

    def _cell_text(self, row: int, col: int) -> str:
        item = self._table.item(row, col)
        return item.text().strip() if item else ""

    # ── Export ──

    def _export(self):
        if not self._validated:
            QMessageBox.warning(self, "Export", "Please validate first.")
            return

        table_name  = self._table_name_edit.text().strip()
        analysis_dt = self._date_edit.date().toString("yyyy-MM-dd")
        initials    = self._initials_edit.text().strip()
        comment     = f"{table_name}_{analysis_dt}_{initials}"

        # Build a lookup from sample_id → r1/r2
        r1_map: dict[str, str] = {}
        r2_map: dict[str, str] = {}
        for row_data in self._rows:
            sid = row_data["sample_id"]
            r1_map[sid] = row_data["r1"]
            r2_map[sid] = row_data["r2"]

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
                # Look up r1/r2 by original sample_id (before dedup rename)
                r1 = r1_map.get(sid, "")
                r2 = r2_map.get(sid, "")
                writer.writerow({
                    f"entity:{table_name}_id": sid,
                    "run_id":           rid,
                    "sample_type":      st,
                    "mode":             mode,
                    "expected_taxa":    exp,
                    "r1_fastq":         r1,
                    "r2_fastq":         r2,
                    "analysis_comments": comment,
                })

        QMessageBox.information(
            self, "Exported",
            f"Saved {n} rows to:\n{save_path}"
        )

    # ── Back ──

    def _go_back(self):
        self.window().stack.setCurrentIndex(0)


# ---------------------------------------------------------------------------
# Main Window
# ---------------------------------------------------------------------------

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("AFI Terra Sheet Builder")
        self.resize(1000, 680)

        self.stack = QStackedWidget()
        self.setCentralWidget(self.stack)

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

def main():
    app = QApplication(sys.argv)
    app.setApplicationName("AFI Terra Sheet Builder")
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
