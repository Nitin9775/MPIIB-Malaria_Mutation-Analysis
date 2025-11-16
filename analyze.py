#!/usr/bin/env python3
"""
analyze.py
Author: Nitin Sharma 
Affiliation: BTech Biotechnology, NIT Allahabad
Date: 2025-11-16

Purpose:
  Simple script to scan a VCF for the PfK13 C580Y mutation and write a short report.
  (Designed for the MPIIB application mini-project.)

How to run:
  python3 analyze.py round_1.vcf Pfalciparum_3D7.fasta

Notes:
  - This file is intentionally kept lightweight and uses only Python standard library.
  - Contact: https://github.com/Nitin9775
"""
import sys
import os

# -------------------- Configuration --------------------
TARGET_CONTIG = 'NC_004331.3'               
CODON580_POS = (1725260, 1725261, 1725262)  
REPORT_FILE = "C580Y_report.txt"
# -------------------------------------------------------

CODON_TO_AA = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def revcomp(seq: str) -> str:
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]

def read_fasta_get_contig_seq(fasta_path: str, contig_name: str):
    """
    Return sequence string for contig_name.
    Looks for a header line that contains contig_name.
    """
    seq_lines = []
    found = False
    with open(fasta_path, 'r') as fh:
        header = None
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                hdr = line[1:].split()[0]
                if contig_name == hdr or contig_name in line[1:]:
                    found = True
                    header = line
                    seq_lines = []
                else:
                    if found:
                        break
                    seq_lines = []
                continue
            if found:
                seq_lines.append(line.strip())
    if not found:
        raise ValueError(f"Contig '{contig_name}' not found in FASTA '{fasta_path}'.\n"
                         "Available headers (first words) may be different; open the FASTA and check.")
    seq = ''.join(seq_lines).upper()
    return seq

def parse_vcf(vcf_path: str):
    variants = {}
    with open(vcf_path, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3].upper()
            alt_field = parts[4].upper()
            alts = [a for a in alt_field.split(',') if a != '.']
            variants[pos] = {
                'ref': ref,
                'alts': alts,
                'raw': line.rstrip('\n')
            }
    return variants

def translate_codon(codon: str) -> str:
    codon = codon.upper().replace('U','T')
    return CODON_TO_AA.get(codon, '?')

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 analyze.py round_1.vcf Pfalciparum_3D7.fasta")
        sys.exit(2)

    vcf_path = sys.argv[1]
    fasta_path = sys.argv[2]

    if not os.path.exists(vcf_path):
        print(f"VCF not found: {vcf_path}")
        sys.exit(1)
    if not os.path.exists(fasta_path):
        print(f"FASTA not found: {fasta_path}")
        sys.exit(1)

    #  parse VCF
    variants = parse_vcf(vcf_path)

    # read reference contig sequence
    try:
        ref_seq = read_fasta_get_contig_seq(fasta_path, TARGET_CONTIG)
    except Exception as e:
        print("Error reading FASTA:", e)
        sys.exit(1)

    # extract reference codon (1-based coordinates)
    p1, p2, p3 = CODON580_POS
    # convert to 0-based indices
    idx1 = p1 - 1
    idx2 = p2 - 1
    idx3 = p3 - 1
    seq_len = len(ref_seq)
    if idx3 >= seq_len or idx1 < 0:
        print(f"Requested codon positions {CODON580_POS} are out of range for contig (length {seq_len})")
        sys.exit(1)
    ref_codon = ref_seq[idx1] + ref_seq[idx2] + ref_seq[idx3]
    ref_aa = translate_codon(ref_codon)

    #  collect variants overlapping the codon positions
    overlapping = {pos: variants[pos] for pos in CODON580_POS if pos in variants}

    # Build possible mutated codons from SNP alts (only handle simple single-base substitutions here)
    mutated_codons = []
    mutation_records = []  # record strings for report
    # start with ref codon
    mutated_codons.append(('REF', ref_codon, ref_aa, 'reference'))

    
    pos_to_alts = {}
    for pos in CODON580_POS:
        if pos in variants:
            var = variants[pos]
            ref_base = var['ref']
            # Only consider SNP alts (len==1)
            snp_alts = [a for a in var['alts'] if len(a) == 1 and len(ref_base) == 1]
            if snp_alts:
                pos_to_alts[pos] = (ref_base, snp_alts)
                mutation_records.append(f"Pos {pos}: {ref_base} -> {','.join(snp_alts)}")

    # generate combinations (cartesian product) of alts across positions
    import itertools
    positions = list(CODON580_POS)
    # Build lists of possible bases at each position
    base_options = []
    for pos in positions:
        if pos in pos_to_alts:
            refb, alts = pos_to_alts[pos]
            base_options.append([refb] + alts)  
        else:
            # no variant -> only ref base
            base_options.append([ref_seq[pos - 1]])

    # build all possible codons from base_options
    for combo in itertools.product(*base_options):
        cod = ''.join(combo)
        aa = translate_codon(cod)
        # Skip the pure reference (already added)
        if cod == ref_codon:
            continue
        mutated_codons.append(('ALT', cod, aa, ' / '.join(mutation_records) if mutation_records else 'variant'))

    #  Check specifically for C580Y signature:
    # We expect reference AA == 'C' (Cysteine). The C580Y mutation is 'C' -> 'Y'
    found_C580Y = False
    evidence = []
    for kind, cod, aa, note in mutated_codons:
        if ref_aa == 'C' and aa == 'Y':
            # Identify which nucleotide change produced it: compare to ref_codon
            diffs = []
            for i,(r,b) in enumerate(zip(ref_codon, cod), start=1):
                if r != b:
                    diffs.append((i, CODON580_POS[i-1], r, b))
            evidence.append({
                'codon': cod,
                'aa': aa,
                'diffs': diffs,
                'note': note
            })
            found_C580Y = True

    # Write report
    with open(REPORT_FILE, 'w') as out:
        out.write(f"Analysis of VCF: {vcf_path}\n")
        out.write(f"Reference FASTA: {fasta_path}\n\n")
        out.write(f"Target contig: {TARGET_CONTIG}\n")
        out.write(f"Codon 580 genomic positions: {CODON580_POS} (1-based)\n")
        out.write(f"Reference codon (genome): {ref_codon} -> {ref_aa}\n\n")

        if overlapping:
            out.write("Variants overlapping codon positions:\n")
            for pos, var in overlapping.items():
                out.write(f"  {pos}: REF={var['ref']} ALT={','.join(var['alts'])}\n")
        else:
            out.write("No variants found overlapping the codon positions.\n")

        out.write("\nPossible mutated codons (from VCF SNP alts):\n")
        for kind,cod,aa,note in mutated_codons:
            out.write(f"  {kind}: {cod} -> {aa}   ({note})\n")

        out.write("\nResult:\n")
        if found_C580Y:
            out.write("SUCCESS: C580Y (Cysteine -> Tyrosine) IS PRESENT.\n")
            out.write("Evidence:\n")
            for ev in evidence:
                out.write(f"  Mutated codon: {ev['codon']} -> {ev['aa']}\n")
                for idx, gpos, rbase, abase in ev['diffs']:
                    out.write(f"    genomic pos {gpos} (codon base {idx}): {rbase} -> {abase}\n")
        else:
            out.write("C580Y NOT FOUND in this VCF.\n")
            if overlapping:
                out.write("Variants present do not change codon to Tyrosine.\n")
            else:
                out.write("No relevant SNPs were present at the middle base (expected G->A) of the codon.\n")

    #  Print summary to console
    print(f"Reference codon at {CODON580_POS}: {ref_codon} -> {ref_aa}")
    if overlapping:
        print("Overlapping variants:")
        for pos,var in overlapping.items():
            print(f"  {pos}: REF={var['ref']} ALT={','.join(var['alts'])}")
    else:
        print("No overlapping variants found for codon positions.")

    if found_C580Y:
        print("SUCCESS: C580Y IS PRESENT. See", REPORT_FILE)
    else:
        print("C580Y NOT FOUND. See", REPORT_FILE)

if __name__ == '__main__':
    main()
