import json
import csv

import ete3


def extract(input_absrel, output_json):
    name = input_absrel.split('/')[-1].split('.')[0]
    attributes = None
    with open('data/Functional_categories.csv') as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            if row['File name'] == name:
                attributes = row
        if not attributes:
            attributes = {
                'Functional Category': 'UNKNOWN',
                'Functional Category Analyses': 'UNKNOWN'
            }
        
    with open(input_absrel) as json_file:
        absrel = json.load(json_file)
    tree = ete3.Tree(absrel['input']['trees']['0'] + ';', format=8)    
    full_hash = {}
    for leaf in tree.get_leaves():
        node = leaf
        padding = ''
        pss_sum = 0
        pss_count = 0
        while not node.is_root():
            branch_attributes = absrel['branch attributes']['0'][node.name]
            p_value = branch_attributes['Corrected P-value']
            if p_value < .05:
                rate_distributions = branch_attributes['Rate Distributions']
                for rd in rate_distributions:
                    rate = float(rd[0])
                    if rate > 1:
                        pss_sum += float(rd[1])
                        pss_count += 1
            node = node.up
        if pss_count == 0:
            pss_average = 0
        else:
            pss_average = pss_sum / pss_count
        full_hash[leaf.name] = {
            'mean_pss': pss_average,
            'functional_category': attributes['Functional Category'],
            'functional_category_analyses': attributes['Functional Category Analyses'],
        }
        with open(output_json, 'w') as json_file:
            json.dump(full_hash, json_file, indent=2)


