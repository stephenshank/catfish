import json
import csv

import ete3


def calculate_mean_pss(absrel, bh_hash):
    tree = ete3.Tree(absrel['input']['trees']['0'] + ';', format=8)    
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
