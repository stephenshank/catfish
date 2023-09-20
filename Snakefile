import os
import csv
import json

import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

from catfish import calculate_mean_pss


def read_lines(input_filepath):
  with open(input_filepath) as input_file:
    lines = [line.strip() for line in input_file.readlines()]
  return lines


def read_json(input_filepath):
  with open(input_filepath) as input_file:
    result = json.load(input_file)
  return result


GENES = read_lines('tables/genes.txt')
TREES = read_lines('tables/trees.txt')

rule absrel:
  input:
    alignment='data/{gene}/codons.fasta',
    tree='data/{gene}/{tree}/tree.nwk'
  output:
    'data/{gene}/{tree}/absrel.json'
  shell:
    'mpirun -np 16 HYPHYMPI absrel --alignment {input.alignment} --tree {input.tree} --branches Foreground --output {output}'

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
          table.append({
            'branch': branch,
            'tree': tree,
            'gene': gene,
            'pvalue': attributes['Uncorrected P-value'],
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
    'data/{gene}/rows.tsv'
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

    full_bh_hash = {}
    with open(input.bh) as tsv_file:
      bh_reader = csv.DictReader(tsv_file, delimiter='\t')
      for row in bh_reader:
        if row['gene'] == wildcards.gene:
          if row['tree'] in full_bh_hash:
            full_bh_hash[row['tree']][row['branch']] = row['rejected']
          else:
            full_bh_hash[row['tree']] = {
              row['branch']: row['rejected']
            }
        
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
      'data/{gene}/row.tsv',
      gene=GENES
    )
  output:
    'tables/full_absrel.csv'
  run:
    # TODO: redo with new data model!
    pass

rule make_tip_csv:
  input:
    rules.make_full_csv.output[0]
  output:
    tip='data/tip_absrel.csv',
    func='data/func_absrel.csv'
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
