#!/usr/bin/env python3
# =============================================================================
# Title:        Extract Gene Signatures from van Galen 2019 Supplement (MMC4)
# Project:      AML_Cellecta
# Author:       Ayeh Sadr (updated)
# Created:      2026-05-23
# Description:  Parses the supplementary table Excel file mmc4(1).xlsx
#               directly using python's standard library (zipfile/xml) to
#               ensure zero external dependencies (no pandas or openpyxl needed).
#               Extracts high-confidence cell state signature genes and
#               updates config/signatures.yml.
# =============================================================================

import os
import sys
import zipfile
import xml.etree.ElementTree as ET

# Number of top genes to extract per signature
NUM_GENES = 40

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
HPC_SETUP_DIR = os.path.dirname(SCRIPT_DIR)
GLOBAL_ROOT = os.path.dirname(HPC_SETUP_DIR)

XLSX_PATH = os.path.join(GLOBAL_ROOT, "Data", "mmc4(1).xlsx")
SIGS_PATH = os.path.join(HPC_SETUP_DIR, "config", "signatures.yml")

def get_shared_strings(z):
    try:
        with z.open('xl/sharedStrings.xml') as f:
            tree = ET.parse(f)
            root = tree.getroot()
            ns = {'n': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
            return [t.text if t.text is not None else "" for t in root.findall('.//n:t', ns)]
    except KeyError:
        return []

def parse_sheet(z, sheet_rel_path, shared_strings):
    with z.open(sheet_rel_path) as f:
        tree = ET.parse(f)
        root = tree.getroot()
        ns = {'n': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        
        rows = []
        for row_elem in root.findall('.//n:row', ns):
            row_idx = int(row_elem.attrib['r'])
            row_cells = {}
            for cell_elem in row_elem.findall('n:c', ns):
                ref = cell_elem.attrib['r']
                col_letter = ''.join([c for c in ref if c.isalpha()])
                val_elem = cell_elem.find('n:v', ns)
                val = ""
                if val_elem is not None:
                    val = val_elem.text
                    t_attr = cell_elem.attrib.get('t')
                    if t_attr == 's':  # shared string
                        val = shared_strings[int(val)]
                row_cells[col_letter] = val
            rows.append((row_idx, row_cells))
        return rows

def clean_float(val):
    if not val:
        return 0.0
    try:
        return float(val)
    except ValueError:
        return 0.0

def main():
    print(f"=== van Galen Signature Extractor ===")
    print(f"Excel file: {XLSX_PATH}")
    print(f"YAML file:  {SIGS_PATH}")

    if not os.path.exists(XLSX_PATH):
        print(f"ERROR: Excel file not found at {XLSX_PATH}", file=sys.stderr)
        print("Please copy/ensure 'mmc4(1).xlsx' is placed in your Data/ folder.", file=sys.stderr)
        sys.exit(1)

    print("Parsing Excel workbook structure...")
    with zipfile.ZipFile(XLSX_PATH) as z:
        shared_strings = get_shared_strings(z)
        print(f"Loaded {len(shared_strings)} shared strings.")
        
        # 1. Parse Sheet 1 (Table S4A) - Stem, GMP, Promono (Myeloid)
        print("Extracting HSC, GMP, Promono markers from Table S4A...")
        rows_a = parse_sheet(z, "xl/worksheets/sheet1.xml", shared_strings)
        
        # Verify columns
        # Col H = HSC/Prog correlation, Col I = GMP correlation, Col J = Myeloid correlation
        hsc_candidates = []
        gmp_candidates = []
        promono_candidates = []
        
        for r_idx, r_cells in rows_a[2:]:  # Skip row 1 descriptive and row 2 header
            gene = r_cells.get('A', '').strip()
            if not gene or gene == "Gene":
                continue
            
            corr_h = clean_float(r_cells.get('H', '0'))
            corr_i = clean_float(r_cells.get('I', '0'))
            corr_j = clean_float(r_cells.get('J', '0'))
            
            hsc_candidates.append((gene, corr_h))
            gmp_candidates.append((gene, corr_i))
            promono_candidates.append((gene, corr_j))
            
        # 2. Parse Sheet 2 (Table S4B) - Monocyte (Mono)
        print("Extracting Monocyte markers from Table S4B...")
        rows_b = parse_sheet(z, "xl/worksheets/sheet2.xml", shared_strings)
        
        mono_candidates = []
        for r_idx, r_cells in rows_b[2:]:  # Skip headers
            gene = r_cells.get('A', '').strip()
            if not gene or gene == "Gene":
                continue
            # Col B = Normal BM average expression
            expr_b = clean_float(r_cells.get('B', '0'))
            mono_candidates.append((gene, expr_b))

    # Sort descending and slice top N
    hsc_candidates.sort(key=lambda x: x[1], reverse=True)
    gmp_candidates.sort(key=lambda x: x[1], reverse=True)
    promono_candidates.sort(key=lambda x: x[1], reverse=True)
    mono_candidates.sort(key=lambda x: x[1], reverse=True)

    hsc_genes = [x[0] for x in hsc_candidates[:NUM_GENES]]
    gmp_genes = [x[0] for x in gmp_candidates[:NUM_GENES]]
    promono_genes = [x[0] for x in promono_candidates[:NUM_GENES]]
    mono_genes = [x[0] for x in mono_candidates[:NUM_GENES]]

    print(f"\nExtracted Signatures (top {NUM_GENES} genes):")
    print(f"  vgalen_hsc_mmc4:     {', '.join(hsc_genes[:5])} ...")
    print(f"  vgalen_gmp_mmc4:     {', '.join(gmp_genes[:5])} ...")
    print(f"  vgalen_promono_mmc4: {', '.join(promono_genes[:5])} ...")
    print(f"  vgalen_mono_mmc4:    {', '.join(mono_genes[:5])} ...")

    # 3. Read signatures.yml and update/append these entries
    if not os.path.exists(SIGS_PATH):
        print(f"ERROR: {SIGS_PATH} not found. Please verify project structure.", file=sys.stderr)
        sys.exit(1)

    with open(SIGS_PATH, "r") as f:
        sigs_lines = f.readlines()

    # We will strip out any previous instances of our custom mmc4 keys to avoid duplicates
    # YAML parser block identification (simple prefix check)
    custom_keys = ["vgalen_hsc_mmc4:", "vgalen_gmp_mmc4:", "vgalen_promono_mmc4:", "vgalen_mono_mmc4:"]
    cleaned_lines = []
    skip_mode = False
    
    for line in sigs_lines:
        line_strip = line.strip()
        # If we see one of our custom keys at the root level, enter skip mode
        if any(line.startswith(key) for key in custom_keys):
            skip_mode = True
            continue
        # If we see another root level key or a divider starting at column 0 while skipping, exit skip mode
        if skip_mode and line.strip() and not line.startswith(" ") and not line.startswith("-") and not line.startswith("#"):
            skip_mode = False
            
        if not skip_mode:
            cleaned_lines.append(line)

    # Clean trailing empty lines
    while cleaned_lines and cleaned_lines[-1].strip() == "":
        cleaned_lines.pop()

    # Append our new signatures
    new_sigs_block = f"""

# --- van Galen Cell 2019 supplementary Table S4 (mmc4) derived signatures ---

vgalen_hsc_mmc4:
  description: HSC-like malignant blast markers (top {NUM_GENES} by correlation in Table S4A)
  citation: van Galen P et al., Cell 176(6):1265-1281, 2019
  genes: [{", ".join(hsc_genes)}]

vgalen_gmp_mmc4:
  description: GMP-like malignant blast markers (top {NUM_GENES} by correlation in Table S4A)
  citation: van Galen P et al., Cell 176(6):1265-1281, 2019
  genes: [{", ".join(gmp_genes)}]

vgalen_promono_mmc4:
  description: Promonocyte-like malignant blast markers (top {NUM_GENES} by correlation in Table S4A)
  citation: van Galen P et al., Cell 176(6):1265-1281, 2019
  genes: [{", ".join(promono_genes)}]

vgalen_mono_mmc4:
  description: Monocyte-like malignant blast markers (top {NUM_GENES} by Normal BM expression in Table S4B)
  citation: van Galen P et al., Cell 176(6):1265-1281, 2019
  genes: [{", ".join(mono_genes)}]
"""
    
    with open(SIGS_PATH, "w") as f:
        f.writelines(cleaned_lines)
        f.write(new_sigs_block)

    print("\nconfig/signatures.yml updated successfully!")

if __name__ == "__main__":
    main()
