#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ContaminOFF: Fast Metagenomic Contamination Filter
Developed for human sequencing data (WGS), specifically targeting saliva and high-noise samples.
"""

import os
import sys
import argparse
import subprocess
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description="ContaminOFF: Fast Metagenomic Contamination Filter")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-r1', '--read1', help="Input FASTQ R1 (.gz) [Single Mode]")
    group.add_argument('-i', '--input-dir', help="Directory containing FASTQ pairs (.gz) [Batch Mode]")
    
    parser.add_argument('-r2', '--read2', help="Input FASTQ R2 (.gz) [Required if using -r1]")
    parser.add_argument('-o', '--out-dir', required=True, help="Output directory for clean FASTQs and reports")
    parser.add_argument('-d', '--db', required=True, help="Path to the Kraken2 database")
    parser.add_argument('-t', '--threads', type=int, default=8, help="Number of CPU threads to use")
    parser.add_argument('-p', '--prefix', default=None, help="Output prefix (Auto-detected if not provided)")
    
    parser.add_argument('-id', '--taxid', default="9606", help="Target TaxID to keep (default: 9606 - Homo sapiens)")
    parser.add_argument('-l', '--tax-levels', nargs='+', default=["G"], help="Taxonomic levels for QC bar plots (e.g., P C O F G)")
    parser.add_argument('-n', '--top-taxa', type=int, default=10, help="Number of top contaminants to show in QC plots")
    parser.add_argument('--keep-kraken', action='store_true', help="Do not delete the temporary .kraken output file")
    
    return parser.parse_args()

def run_kraken2(r1, r2, db, threads, out_kraken, out_report):
    print(f"[*] Running Kraken2 classification...")
    cmd = [
        "kraken2", "--db", db, "--threads", str(threads),
        "--paired", "--gzip-compressed",
        "--output", out_kraken, "--report", out_report,
        r1, r2
    ]
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        print(f"\n[FATAL ERROR] Kraken2 failed to load the database.")
        print(f"Please verify if the path '{db}' is correct and contains the .k2d files.")
        sys.exit(1)

def extract_target_reads(kraken_out, target_taxid):
    print(f"[*] Parsing Kraken2 output for TaxID {target_taxid}...")
    valid_ids = set()
    with open(kraken_out, 'r') as f:
        for line in f:
            parts = line.split('\t')
            if parts[2] == target_taxid:
                valid_ids.add(parts[1].encode('ascii'))
    print(f"    -> Found {len(valid_ids)} target read pairs.")
    return valid_ids

def filter_fastq_pigz(r1_in, r2_in, r1_out, r2_out, valid_ids, threads):
    print(f"[*] Filtering FASTQ files using high-speed pigz I/O...")
    
    p1_in = subprocess.Popen(['pigz', '-dc', r1_in], stdout=subprocess.PIPE)
    p2_in = subprocess.Popen(['pigz', '-dc', r2_in], stdout=subprocess.PIPE)
    
    t_out = str(max(1, threads // 2))
    p1_out = subprocess.Popen(['pigz', '-p', t_out, '-c'], stdin=subprocess.PIPE, stdout=open(r1_out, 'wb'))
    p2_out = subprocess.Popen(['pigz', '-p', t_out, '-c'], stdin=subprocess.PIPE, stdout=open(r2_out, 'wb'))
    
    while True:
        id1 = p1_in.stdout.readline()
        if not id1: break
        seq1 = p1_in.stdout.readline()
        plus1 = p1_in.stdout.readline()
        qual1 = p1_in.stdout.readline()
        
        id2 = p2_in.stdout.readline()
        seq2 = p2_in.stdout.readline()
        plus2 = p2_in.stdout.readline()
        qual2 = p2_in.stdout.readline()
        
        read_id = id1.split(b' ')[0].split(b'/')[0][1:]
        
        if read_id in valid_ids:
            p1_out.stdin.write(id1 + seq1 + plus1 + qual1)
            p2_out.stdin.write(id2 + seq2 + plus2 + qual2)
            
    p1_out.stdin.close()
    p2_out.stdin.close()
    p1_out.wait()
    p2_out.wait()

def generate_qc_report(report_file, prefix, out_dir, target_taxid, top_n, tax_levels):
    print(f"[*] Generating QC Reports and Plots for {prefix}...")
    
    col_names = ['Pct', 'Reads_rooted', 'Reads_direct', 'Rank', 'TaxID', 'Name']
    df = pd.read_csv(report_file, sep='\t', names=col_names)
    df['Name'] = df['Name'].str.strip()
    
    try: pct_unclass = df.loc[df['TaxID'] == 0, 'Pct'].values[0]
    except IndexError: pct_unclass = 0.0
    
    try:
        pct_target = df.loc[df['TaxID'] == int(target_taxid), 'Pct'].values[0]
        reads_target = df.loc[df['TaxID'] == int(target_taxid), 'Reads_rooted'].values[0]
    except IndexError: 
        pct_target, reads_target = 0.0, 0
        
    try: pct_bact = df.loc[df['TaxID'] == 2, 'Pct'].values[0]
    except IndexError: pct_bact = 0.0
    
    summary_path = os.path.join(out_dir, f"{prefix}_summary.tsv")
    with open(summary_path, 'w') as f:
        f.write("Sample\tTarget_Reads\tTarget_Pct\tBacteria_Pct\tUnclassified_Pct\n")
        f.write(f"{prefix}\t{reads_target}\t{pct_target}\t{pct_bact}\t{pct_unclass}\n")
        
    sns.set_theme(style="whitegrid", context="paper")
    
    # 1. Generate Overall Composition Pie Chart
    fig, ax = plt.subplots(figsize=(8, 8))
    labels = ['Target (Human)', 'Bacteria', 'Unclassified']
    sizes = [pct_target, pct_bact, pct_unclass]
    colors = ['#1F77B4', '#D62728', '#A0A0A0']
    
    ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors, 
           wedgeprops={'edgecolor': 'black'}, normalize=True)
    ax.set_title('Overall Sample Composition', fontweight='bold', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}_QC_Composition.png"), dpi=300)
    plt.close()

    # Define full human lineage to exclude from contamination plots
    human_lineage = ['Eukaryota', 'Chordata', 'Mammalia', 'Primates', 'Hominidae', 'Homo', 'Homo sapiens']

    # 2. Generate Separate Bar Plots for each requested taxonomic level
    for level in tax_levels:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        df_contam = df[(df['Rank'] == level) & 
                       (df['TaxID'] != int(target_taxid)) & 
                       (~df['Name'].isin(human_lineage))].copy()
                       
        df_contam = df_contam.sort_values(by='Pct', ascending=False).head(top_n)
        
        if not df_contam.empty:
            sns.barplot(data=df_contam, y='Name', x='Pct', ax=ax, hue='Name', palette='Reds_r', legend=False, edgecolor='black')
            ax.set_title(f'Top {top_n} Contaminants (Level: {level})', fontweight='bold', fontsize=14)
            ax.set_xlabel('Percentage of Total Reads (%)', fontsize=12)
            ax.set_ylabel('')
        else:
            ax.text(0.5, 0.5, 'No significant contaminants found.', ha='center', va='center', fontsize=12)
            ax.set_title(f'Top {top_n} Contaminants (Level: {level})', fontweight='bold', fontsize=14)
            
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{prefix}_QC_Contaminants_Level_{level}.png"), dpi=300)
        plt.close()

def compile_final_reports(processed_prefixes, out_dir, target_taxid):
    """Compiles individual summaries and taxonomy profiles into final cohort tables."""
    print(f"\n{'='*50}\n[*] Compiling Final Reports...\n{'='*50}")
    
    summary_files = [os.path.join(out_dir, f"{p}_summary.tsv") for p in processed_prefixes]
    valid_summaries = [f for f in summary_files if os.path.exists(f)]
    
    if valid_summaries:
        df_list = [pd.read_csv(f, sep='\t') for f in valid_summaries]
        df_merged = pd.concat(df_list, ignore_index=True)
        
        out_name = "samples_summary.tsv" if len(processed_prefixes) > 1 else f"{processed_prefixes[0]}_summary_final.tsv"
        cohort_summary_path = os.path.join(out_dir, out_name)
        
        df_merged.to_csv(cohort_summary_path, sep='\t', index=False)
        for f in valid_summaries: os.remove(f)
        print(f"[+] Saved merged summary: {cohort_summary_path}")

    taxa_dict = {}
    tax_hierarchy = ['D', 'P', 'C', 'O', 'F', 'G']
    human_lineage = ['Eukaryota', 'Chordata', 'Mammalia', 'Primates', 'Hominidae', 'Homo', 'Homo sapiens']
    
    for prefix in processed_prefixes:
        rep = os.path.join(out_dir, f"{prefix}.report")
        if not os.path.exists(rep): continue
        
        lineage = {r: "" for r in tax_hierarchy}
        
        with open(rep, 'r') as f:
            for line in f:
                parts = line.strip('\n').split('\t')
                if len(parts) < 6: continue
                
                reads_rooted = int(parts[1])
                rank = parts[3]
                taxid = int(parts[4])
                name = parts[5].strip()
                
                if taxid == 0 or taxid == int(target_taxid) or name in human_lineage:
                    continue
                    
                if rank in tax_hierarchy:
                    lineage[rank] = name
                    idx = tax_hierarchy.index(rank)
                    for lower_rank in tax_hierarchy[idx+1:]:
                        lineage[lower_rank] = ""
                        
                if rank == 'G':
                    tax_path = " | ".join([lineage[r] for r in tax_hierarchy if lineage[r]])
                    
                    if tax_path not in taxa_dict:
                        taxa_dict[tax_path] = {'Taxonomy': tax_path}
                    
                    taxa_dict[tax_path][prefix] = reads_rooted

    if taxa_dict:
        df_contam = pd.DataFrame.from_dict(taxa_dict, orient='index')
        
        samples = [col for col in df_contam.columns if col != 'Taxonomy']
        df_contam[samples] = df_contam[samples].fillna(0).astype(int)
        
        df_contam['Total'] = df_contam[samples].sum(axis=1)
        df_contam = df_contam.sort_values(by='Total', ascending=False).drop(columns=['Total'])
        
        out_name_contam = "contaminants_table.tsv" if len(processed_prefixes) > 1 else f"{processed_prefixes[0]}_contaminants_table.tsv"
        out_contam = os.path.join(out_dir, out_name_contam)
        
        cols = ['Taxonomy'] + samples
        df_contam = df_contam[cols]
        df_contam.to_csv(out_contam, sep='\t', index=False)
        print(f"[+] Saved contamination table: {out_contam}")

def process_sample(r1, r2, prefix, args):
    print(f"\n{'='*50}\n[ContaminOFF] Processing: {prefix}\n{'='*50}")
    
    out_kraken = os.path.join(args.out_dir, f"{prefix}.kraken")
    out_report = os.path.join(args.out_dir, f"{prefix}.report")
    out_r1 = os.path.join(args.out_dir, f"{prefix}_clean_R1.fq.gz")
    out_r2 = os.path.join(args.out_dir, f"{prefix}_clean_R2.fq.gz")
    
    run_kraken2(r1, r2, args.db, args.threads, out_kraken, out_report)
    valid_ids = extract_target_reads(out_kraken, args.taxid)
    filter_fastq_pigz(r1, r2, out_r1, out_r2, valid_ids, args.threads)
    generate_qc_report(out_report, prefix, args.out_dir, args.taxid, args.top_taxa, args.tax_levels)
    
    if not args.keep_kraken:
        os.remove(out_kraken)
        print(f"[*] Cleaned up temporary .kraken file.")
        
    print(f"[+] Finished {prefix} successfully!")

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    
    processed_prefixes = []
    
    if args.read1:
        if not args.read2:
            print("[ERROR] -r2/--read2 is required in single-sample mode.")
            sys.exit(1)
            
        if args.prefix is None:
            prefix = os.path.basename(args.read1).split('_R1')[0].split('_1')[0]
            prefix = prefix.replace('.fastq.gz', '').replace('.fq.gz', '')
        else:
            prefix = args.prefix
            
        process_sample(args.read1, args.read2, prefix, args)
        processed_prefixes.append(prefix)
        
    elif args.input_dir:
        print(f"[*] Initializing Batch Mode in directory: {args.input_dir}")
        search_pattern = os.path.join(args.input_dir, "*_R1*.fq.gz")
        r1_files = glob.glob(search_pattern)
        if not r1_files:
            search_pattern = os.path.join(args.input_dir, "*_1*.fastq.gz")
            r1_files = glob.glob(search_pattern)
            
        for r1 in sorted(r1_files):
            r2 = r1.replace('_R1', '_R2').replace('_1.fastq.gz', '_2.fastq.gz')
            if os.path.exists(r2):
                prefix = os.path.basename(r1).split('_R1')[0].split('_1')[0]
                process_sample(r1, r2, prefix, args)
                processed_prefixes.append(prefix)
            else:
                print(f"[WARNING] Could not find matching R2 for {r1}. Skipping.")
                
    if processed_prefixes:
        compile_final_reports(processed_prefixes, args.out_dir, args.taxid)

if __name__ == "__main__":
    main()
