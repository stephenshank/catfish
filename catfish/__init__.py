import json
import csv

import ete3


def read_lines(input_filepath):
  with open(input_filepath) as input_file:
    lines = [line.strip() for line in input_file.readlines()]
  return lines


def read_json(input_filepath):
  with open(input_filepath) as input_file:
    result = json.load(input_file)
  return result



def build_bh_hash(bh_tsv_filepath, gene):
    full_bh_hash = {}
    with open(bh_tsv_filepath) as tsv_file:
      bh_reader = csv.DictReader(tsv_file, delimiter='\t')
      for row in bh_reader:
        value = row['rejected'] == 'True'
        if row['gene'] == gene:
          if row['tree'] in full_bh_hash:
            full_bh_hash[row['tree']][row['branch']] = value
          else:
            full_bh_hash[row['tree']] = {
              row['branch']: value
            }
    return full_bh_hash


def calculate_mean_pss(absrel, bh_hash):
    try:
        tree = ete3.Tree(absrel['input']['trees']['0'] + ';', format=8)
    except:
        print('ERROR PARSING:')
        print(absrel['input']['trees']['0'])
    tip_hash = {}
    for leaf in tree.get_leaves():
        node = leaf
        pss_sum = 0.0
        pss_count = 0.0
        while not node.is_root():
            pss_count += 1.0
            if bh_hash[node.name]:
                branch_attributes = absrel['branch attributes']['0'][node.name]
                rate_distributions = branch_attributes['Rate Distributions']
                for rd in rate_distributions:
                    rate = float(rd[0])
                    if rate > 1.0:
                        pss_sum += float(rd[1])
            node = node.up
        if pss_count == 0.0:
            pss_average = 0.0
        else:
            pss_average = pss_sum / pss_count
        tip_hash[leaf.name] = pss_average
    return tip_hash


if __name__ == '__main__':
    gene = '117480at7898'
    tree = 'allastral'
    full_bh_hash = build_bh_hash('tables/bh.tsv', gene)
    absrel = read_json('data/%s/%s/absrel.json' % (gene, tree))
    all_mean_pss = calculate_mean_pss(absrel, full_bh_hash[tree])
