#!/usr/bin/env python
import argparse
import pickle
from Bio.KEGG.REST import kegg_get


if __name__ == '__main__':
    from kegg_module import KeggModule, KeggReaction, Catalyst, KeggEnzyme, parse
    parser = argparse.ArgumentParser()
    parser.add_argument("pathway", help="An accession to a kegg module")
    parser.add_argument("outfile", help="output file for the pickled KeggModule object")
    args = parser.parse_args()

    for module in parse(kegg_get(args.pathway)):
        with open(args.outfile, 'wb') as pickled:
            pickle.dump(file=pickled, obj=module, protocol=-1)
