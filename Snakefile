import os
import csv

import pandas as pd

from catfish import extract


rule extract_fasta_names:
  input:
    'data/fasta'
  output:
    'data/fasta_names.txt'
  shell:
    'ls {input} | sed "s/\.fasta//g" > {output}'

rule absrel:
  input:
    alignment='data/Stephen/aBSREL_Input/{name}.fasta',
    tree='data/Stephen/aBSREL_Input/{name}.nwk'
  output:
    'data/absrel/{name}.ABSREL.json' 
  shell:
    'mpirun -np 16 HYPHYMPI absrel --alignment {input.alignment} --tree {input.tree} --branches Foreground --output {output}'

rule extract:
  input:
    rules.absrel.output[0]
  output:
    'data/absrel/{name}.extract.json' 
  run:
    extract(input[0], output[0])

rule make_full_csv:
  input:
    expand(
      'data/absrel/{name}.extract.json',
      name=[
        f.split('.')[0]
        for f in os.listdir('./data/Stephen/aBSREL_Input/')
        if f[-6:] == '.fasta'
      ])
  output:
    'data/full_absrel.csv'
  run:
    csv_file = open(output[0], 'w')
    field_names = [
      'id',
      'tip',
      'functional_category',
      'functional_category_analyses',
      'mean_pss'
    ]
    csv_writer = csv.DictWriter(csv_file, field_names)
    csv_writer.writeheader()
    for extract_filename in input:
      extract_id = extract_filename.split('/')[-1].split('.')[0]
      with open(extract_filename) as extract_file:
        extracted = json.load(extract_file)
      for key, value in extracted.items():
        row = {
          'id': extract_id,
          'tip': key,
        }
        row.update(value)
        csv_writer.writerow(row)
    csv_file.close()

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
