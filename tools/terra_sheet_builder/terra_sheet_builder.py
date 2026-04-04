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

def discover_pairs(folder: str) -> list[tuple[str, str, str]]:
    """Return list of (sample_id, r1_path, r2_path) from a folder."""
    p = Path(folder)
    r1_files = sorted(
        f for f in p.iterdir()
        if f.suffix in FASTQ_SUFFIXES or "".join(f.suffixes) in FASTQ_SUFFIXES
        if R1_PATTERNS.search(f.name)
    )
    pairs = []
    for r1 in r1_files:
        r2_name = R1_PATTERNS.sub(lambda m: m.group().replace("R1", "R2"), r1.name)
        r2 = p / r2_name
        if not r2.exists():
            # Try .gz variant
            r2_gz = p / (r2_name + ".gz")
            if r2_gz.exists():
                r2 = r2_gz
            else:
                continue
        stem = R1_PATTERNS.split(r1.name)[0]
        pairs.append((stem, str(r1), str(r2)))
    return pairs


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
        self._part_c_scroll = QScrollArea()
        self._part_c_scroll.setWidgetResizable(True)
        self._part_c_inner = QWidget()
        self._part_c_layout = QVBoxLayout(self._part_c_inner)
        self._part_c_scroll.setWidget(self._part_c_inner)
        c_outer_layout = QVBoxLayout(self._part_c_outer)
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
        btn_add = QPushButton("Add files manually…")
        btn_layout.addWidget(btn_browse)
        btn_layout.addWidget(btn_add)
        btn_layout.addStretch()
        box_layout.addLayout(btn_layout)
        box_layout.addWidget(table)

        section = {"widget": box, "run_name": run_name, "table": table, "rows": []}
        self._run_folder_sections.append(section)

        btn_browse.clicked.connect(lambda checked=False, s=section: self._browse_folder(s))
        btn_add.clicked.connect(lambda checked=False, s=section: self._add_files_manually(s))

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

    def _add_files_manually(self, section: dict):
        r1, _ = QFileDialog.getOpenFileName(
            self, "Select R1 FASTQ", "", "FASTQ files (*.fastq *.fastq.gz *.fq *.fq.gz)"
        )
        if not r1:
            return
        r2, _ = QFileDialog.getOpenFileName(
            self, "Select R2 FASTQ", "", "FASTQ files (*.fastq *.fastq.gz *.fq *.fq.gz)"
        )
        if not r2:
            return
        stem = R1_PATTERNS.split(Path(r1).name)[0]
        self._add_row(section, stem, r1, r2)

    def _add_row(self, section: dict, sample_id: str, r1: str, r2: str):
        table: QTableWidget = section["table"]
        row = table.rowCount()
        table.insertRow(row)
        table.setItem(row, 0, QTableWidgetItem(sample_id))
        table.setItem(row, 1, QTableWidgetItem(r1))
        table.setItem(row, 2, QTableWidgetItem(r2))
        btn_rm = QPushButton("✕")
        btn_rm.setFixedWidth(30)
        btn_rm.clicked.connect(lambda checked=False, r=row, t=table: self._remove_row(t, r))
        table.setCellWidget(row, 3, btn_rm)
        section["rows"].append({"sample_id": sample_id, "r1": r1, "r2": r2})

    def _remove_row(self, table: QTableWidget, row: int):
        table.removeRow(row)

    # ── Part C → Screen 2 ──

    def _go_to_screen2(self):
        # Collect all rows from all sections
        all_rows = []
        for section in self._run_folder_sections:
            table: QTableWidget = section["table"]
            run_name = section["run_name"]
            for r in range(table.rowCount()):
                sample_id = (table.item(r, 0) or QTableWidgetItem("")).text().strip()
                r1 = (table.item(r, 1) or QTableWidgetItem("")).text().strip()
                r2 = (table.item(r, 2) or QTableWidgetItem("")).text().strip()
                if sample_id and r1 and r2:
                    all_rows.append({
                        "sample_id": sample_id,
                        "run_id": run_name,
                        "r1": r1,
                        "r2": r2,
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
