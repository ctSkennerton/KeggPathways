#!/usr/bin/env python
import argparse
import pickle


if __name__ == '__main__':
    from kegg_module import KeggModule, KeggReaction, Catalyst, KeggEnzyme
    import IPython
    parser = argparse.ArgumentParser()
    parser.add_argument("pathway", help="A pickled KeggModule object")
    args = parser.parse_args()

    with open(args.pathway, 'rb') as pickled:
        kegg_module = pickle.load(pickled)
        IPython.embed()
