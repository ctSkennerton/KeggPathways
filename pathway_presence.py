#!/usr/bin/env python
import argparse
import pickle
import pandas as pd
from plotnine import *
import os


if __name__ == '__main__':
    from kegg_module import KeggModule, KeggReaction, Catalyst, KeggEnzyme
    parser = argparse.ArgumentParser()
    parser.add_argument('-2', '--2', dest='two_column', action='store_true', default=False, help='input is a two column file, the first column is the organism name and the second column is the gene identifier')
    parser.add_argument("gene_list", help="A file containing gene identifiers, one per line")
    parser.add_argument("pathway", nargs='+', help="A pickled KeggModule object")
    args = parser.parse_args()
    gene_list = {}
    with open(args.gene_list) as fp:
        for n, line in enumerate(fp):
            if args.two_column:
                try:
                    organism, gene = line.strip().split('\t')
                except ValueError as e:
                    raise ValueError("problem on line {}\n{}".format(n, e))
                try:
                    gene_list[organism].add(gene)
                except KeyError:
                    gene_list[organism] = set([gene])
            else:
                gene = line.strip()
                try:
                    gene_list[''].add(gene)
                except KeyError:
                    gene_list[''] = set([gene])

    for pathway in args.pathway:
        data = {'organism': [], 'reaction': [], 'completeness': []}
        reaction_steps = []
        with open(pathway, 'rb') as pickled:
            kegg_module = pickle.load(pickled)
            for reaction in kegg_module.reactions:
                reaction_steps.append(str(reaction))
                for organism, genes in gene_list.items():
                    data['organism'].append(organism)
                    data['reaction'].append(str(reaction))
                    data['completeness'].append(reaction.completeness(genes))
                    #print(organism, reaction, reaction.completeness(genes), sep="\t")
        data = pd.DataFrame(data)
        p = (ggplot(data, aes('reaction', 'organism', fill='completeness')) + geom_tile() + scale_x_discrete(limits=reaction_steps) + theme_minimal() + theme(axis_text_x=element_text(rotation=45, hjust=1)))
        p.save("{}.pdf".format(os.path.basename(pathway)))
