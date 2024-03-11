import os
import csv
import json
import sys

import pandas as pd
from Bio import SeqIO
from statsmodels.stats.multitest import fdrcorrection

from catfish import *


GENES = read_lines('tables/genes.txt')
TREES = read_lines('tables/trees.txt')
TRAITS = ['B2P', 'M2F', 'S2E']

HYPHY_ROOT = os.environ.get('HYPHY_ROOT')
DEV_ROOT = os.path.join(HYPHY_ROOT, 'hyphy', 'hyphy-dev')
HA_ROOT = os.path.join(HYPHY_ROOT, 'hyphy-analyses')

functional_categories = [
  'UNCLASSIFIED',
  'GO:0009987',
  'GO:0065007',
  'GO:0032502',
  'GO:0032501',
  'GO:0008152',
  'GO:0051179',
  'GO:0050896',
  'GO:0042592',
  'GO:0022414',
  'transmembrane transport',
  'GO:0002376',
  'GO:0098772',
  'GO:0005488',
  'GO:0040011',
  'GO:0140110',
  'GO:0005198',
  'GO:0003824',
  'GO:0005215',
  'GO:0044419',
  'GO:0048511',
  'GO:0040007',
  'GO:0098754',
  'GO:0000003',
  'signaling',
  'GO:0043473',
  'locomotion'
]


rule bh_extraction:
  input:
    expand(
      'data/{gene}/root_to_tip/{tree}/absrel.json',
      gene=GENES,
      tree=TREES
    )
  output:
    'tables/bh.tsv'
  run:
    table = []
    for absrel_filepath in input:
      with open(absrel_filepath) as json_file:
        try:
          absrel = json.load(json_file)
        except:
          with open('error.txt', 'w') as f:
            f.write('could not load %s' % absrel_filepath)
          sys.exit(1)
        _, gene, _, tree, _ = absrel_filepath.split('/')
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
    func='tables/function.tsv',
    absrels=expand(
      'data/{{gene}}/root_to_tip/{tree}/absrel.json',
      tree=TREES
    )
  output:
    'data/{gene}/rows.csv'
  run:
    functional = None
    with open(input.func) as csv_file:
      func_reader = csv.DictReader(csv_file, delimiter='\t')
      for row in func_reader:
        if row['Gene_name'] == wildcards.gene:
          functional = {
            fc: row[fc]
            for fc in functional_categories
          }

    full_bh_hash = build_bh_hash(input.bh, wildcards.gene)
        
    rows = []
    for absrel_filepath in input.absrels:
      _, _, _, tree, _ = absrel_filepath.split('/')
      absrel = read_json(absrel_filepath)
      all_mean_pss = calculate_mean_pss(absrel, full_bh_hash[tree])
      for tip, mean_pss in all_mean_pss.items():
        row = {
          'gene': wildcards.gene,
          'tip': tip,
          'mean_pss': mean_pss,
          'tree': tree
        }
        rows.append({**row, **functional})
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
    pss='tables/pss-functional.csv',
    psg='tables/psg-functional.csv'
  run:
    df = pd.read_csv(input[0])
    df['mean_psg'] = df['mean_pss'] > 0

    all_pss_columns = []
    all_psg_columns = []
    for fc in functional_categories:
      subset = df.loc[df[fc]]
      mean_pss = subset.groupby('tip')['mean_pss'].mean().rename(fc)
      mean_psg = (subset.groupby('tip')['mean_psg'].sum() / 10).rename(fc)
      all_pss_columns.append(mean_pss)
      all_psg_columns.append(mean_psg)
    full_pss = pd.concat(all_pss_columns, axis=1)
    full_psg = pd.concat(all_psg_columns, axis=1)
    full_pss.to_csv(output.pss)
    full_psg.to_csv(output.psg)
