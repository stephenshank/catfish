import os
import csv
import json

import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

from catfish import *


GENES = read_lines('tables/genes.txt')
TREES = read_lines('tables/trees.txt')
TRAITS = ['B2P', 'M2F', 'S2E']

HYPHY_ROOT = os.environ.get('HYPHY_ROOT')
DEV_ROOT = os.path.join(HYPHY_ROOT, 'hyphy', 'hyphy-dev')
HA_ROOT = os.path.join(HYPHY_ROOT, 'hyphy-analyses')

rule genes:
  input:
    "tables/non-dupes.tsv"
  output:
    "tables.genes.txt"
  shell:
    "tail -n + 2 {input} | cut -f 1 > {output}"

rule bealign:
  input:
    unaligned='data/{gene}/unaligned.fasta',
    reference='data/{gene}/reference.fasta',
  output:
    bam='data/{gene}/codons.bam',
    fasta='data/{gene}/codons.fasta'
  shell:
    '''
      NCPU=1 bealign -r {input.reference} {input.unaligned} {output.bam}
      bam2msa {output.bam} {output.fasta}
    '''

rule absrel_root_to_tip:
  input:
    alignment=rules.bealign.output.fasta,
    tree='data/{gene}/root_to_tip/{tree}/tree.nwk'
  params:
    stdout='data/{gene}/root_to_tip/{tree}/stdout.txt',
    stderr='data/{gene}/root_to_tip/{tree}/stderr.txt'
  output:
    'data/{gene}/root_to_tip/{tree}/absrel.json'
  shell:
    'mpirun -np 8 HYPHYMPI absrel --alignment {input.alignment} --tree {input.tree} --branches Foreground --output {output}> {params.stdout} 2> {params.stderr}'

rule absrel:
  input:
    alignment=rules.bealign.output.fasta,
    tree='data/{gene}/aBSREL/{trait}/tree.nwk'
  output:
    json='data/{gene}/aBSREL/{trait}/absrel.json'
  params:
    stdout='data/{gene}/aBSREL/{trait}/stdout.txt',
    stderr='data/{gene}/aBSREL/{trait}/stderr.txt'
  shell:
    'mpirun -np 8 HYPHYMPI absrel --alignment {input.alignment} --tree {input.tree} --branches Foreground --output {output.json} > {params.stdout} 2> {params.stderr}'

rule busted_e:
  input:
    alignment=rules.bealign.output.fasta,
    tree='data/{gene}/BUSTED/{trait}/tree.nwk'
  output:
    json='data/{gene}/BUSTED/{trait}/busted_e.json'
  params:
    stdout='data/{gene}/aBSREL/{trait}/stdout.txt',
    stderr='data/{gene}/aBSREL/{trait}/stderr.txt'
  shell:
    '%s/hyphy busted CPU=8 --alignment {input.alignment} --tree {input.tree} --branches Foreground --output {output.json} --error-sink Yes --starting-points 5 > {params.stdout} 2> {params.stderr}' % DEV_ROOT

rule busted_ph:
  input:
    alignment=rules.bealign.output.fasta,
    tree='data/{gene}/BUSTED-PH/{trait}/tree.nwk'
  output:
    json='data/{gene}/BUSTED-PH/{trait}/busted_e.json'
  params:
    stdout='data/{gene}/aBSREL/{trait}/stdout.txt',
    stderr='data/{gene}/aBSREL/{trait}/stderr.txt'
  shell:
    'hyphy %s/BUSTED-PH/BUSTED-PH.bf --alignment {input.alignment} --tree {input.tree} --srv No --branches Primates > {params.stdout} 2> {params.stderr}' % HA_ROOT

rule relax:
  input:
    alignment=rules.bealign.output.fasta,
    tree='data/{gene}/RELAX/{trait}/tree.nwk'
  output:
    json='data/{gene}/RELAX/{trait}/relax.json'
  params:
    stdout='data/{gene}/aBSREL/{trait}/stdout.txt',
    stderr='data/{gene}/aBSREL/{trait}/stderr.txt'
  shell:
    'hyphy relax CPU=8 --alignment {input.alignment} --tree {input.tree} --test TEST --reference REFERENCE --output {output.json} > {params.stdout} 2> {params.stderr}'

rule bh_extraction:
  input:
    expand(
      'data/{gene}/{tree}/absrel.json',
      gene=GENES,
      tree=TREES
    )
  output:
    'tables/bh.tsv'
  run:
    table = []
    for absrel_filepath in input:
      with open(absrel_filepath) as json_file:
        absrel = json.load(json_file)
        _, gene, tree, _ = absrel_filepath.split('/')
        branch_attributes = absrel['branch attributes']['0']
        for branch, attributes in branch_attributes.items():
          rds = attributes['Rate Distributions']
          table.append({
            'branch': branch,
            'tree': tree,
            'gene': gene,
            'pvalue': attributes['Uncorrected P-value'],
            'rate_1': rds[0][0],
            'frac_1': rds[0][1],
            'rate_2': None if len(rds) < 2 else rds[1][0],
            'frac_2': None if len(rds) < 2 else rds[1][1],
            'rate_3': None if len(rds) < 3 else rds[2][0],
            'frac_3': None if len(rds) < 3 else rds[2][1],
            'full_adaptive_model_bl': attributes['Full adaptive model']
          })
    df = pd.DataFrame(table)
    rejected, corrected = fdrcorrection(df.pvalue)
    df['rejected'] = rejected
    df['qvalue'] = corrected
    df.to_csv(output[0], sep='\t', index=False)

rule tip_data_rows_for_gene:
  input:
    bh=rules.bh_extraction.output[0],
    func='tables/Functional_categories.csv',
    absrels=expand(
      'data/{{gene}}/{tree}/absrel.json',
      tree=TREES
    )
  output:
    'data/{gene}/rows.csv'
  run:
    functional = None
    with open(input.func) as csv_file:
      func_reader = csv.DictReader(csv_file)
      for row in func_reader:
        if row['File name'] == wildcards.gene:
          functional = row
      if not functional:
        functional = {
          'Functional Category': 'UNKNOWN',
          'Functional Category Analyses': 'UNKNOWN'
        }

    full_bh_hash = build_bh_hash(input.bh, wildcards.gene)
        
    rows = []
    for absrel_filepath in input.absrels:
      _, _, tree, _ = absrel_filepath.split('/')
      absrel = read_json(absrel_filepath)
      all_mean_pss = calculate_mean_pss(absrel, full_bh_hash[tree])
      for tip, mean_pss in all_mean_pss.items():
        rows.append({
          'gene': wildcards.gene,
          'tip': tip,
          'functional_category': functional['Functional Category'],
          'functional_category_analyses': functional['Functional Category Analyses'],
          'mean_pss': mean_pss,
          'tree': tree
        })
    pd.DataFrame(rows).to_csv(output[0], index=False)


rule make_full_csv:
  input:
    expand(
      'data/{gene}/rows.csv',
      gene=GENES
    )
  output:
    'tables/full_absrel.csv'
  run:
    all_tables = [
      pd.read_csv(single_gene_rows)
      for single_gene_rows in input
    ]
    pd.concat(all_tables).to_csv(output[0], index=False)

rule make_tip_csv:
  input:
    rules.make_full_csv.output[0]
  output:
    tip='tables/tip_absrel.csv',
    func='tables/func_absrel.csv'
  run:
    df = pd.read_csv(input[0])
    df['mean_psg'] = df['mean_pss'] > 0

    mean_pss = df.groupby('tip')['mean_pss'].mean()
    mean_psg = df.groupby('tip')['mean_psg'].sum() / 10
    tip_csv = pd.concat([mean_pss, mean_psg], axis=1)
    tip_csv.to_csv(output.tip)

    mean_pss = df.groupby('functional_category')['mean_pss'].mean()
    mean_psg = df.groupby('functional_category')['mean_psg'].sum() / 10
    func_csv = pd.concat([mean_pss, mean_psg], axis=1)
    func_csv.to_csv(output.func)

rule all_root_to_tip:
  input:
    expand(
      'data/{gene}/root_to_tip/{tree}/absrel.json',
      gene=GENES,
      tree=TREES,
    )

rule all_hyphy:
  input:
    expand(
      'data/{gene}/BUSTED/{trait}/busted_e.json',
      gene=GENES,
      trait=TRAITS,
    ),
    expand(
      'data/{gene}/aBSREL/{trait}/absrel.json',
      gene=GENES,
      trait=TRAITS,
    ),
    expand(
      'data/{gene}/RELAX/{trait}/relax.json',
      gene=GENES,
      trait=TRAITS,
    )
