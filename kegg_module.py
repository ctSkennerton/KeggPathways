# Copyright 2017 by Connor T. Skennerton. All rights reserved.
'''Parser for a KEGG module entry

'''

from __future__ import print_function
import re
import pickle

from Bio.KEGG import _write_kegg, _wrap_kegg
from Bio.KEGG.REST import kegg_get
from Bio.KEGG.Compound import Record as KeggCompound, parse as kegg_compound_parse
from Bio.Pathway import Reaction
from urllib.error import HTTPError


def kegg_compound_read(handle):
    iterator = kegg_compound_parse(handle)
    try:
        first = next(iterator)
    except StopIteration:
        first = None
    if first is None:
        raise ValueError("No records found in handle")
    try:
        second = next(iterator)
    except StopIteration:
        second = None
    if second is not None:
        raise ValueError("More than one record found in handle")
    return first

class KeggParseError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class KeggModule(object):
    '''Object representation of a KEGG Module record
    '''
    def __init__(self):
        self.entry = ""
        self.name = []
        self.definition = []
        self.orthologs = []
        self.classname = []
        self.pathway = []
        self.reactions = []
        self.compounds = {}
        self.definition = ""

    def present(self, gene_list):
        '''Check that all of the reactions can proceed

        The will test to make sure that all of the reactions
        have at least one enzyme catalyst that has a gene
        listed in gene_list
        '''
        for r in self.reactions:
            if not r.present(gene_list):
                return False
        return True

    def completeness(self, gene_list):
        '''return the mean of the completeness of all reactions
        '''
        total = 0.0
        for r in self.reactions:
            total += r.completeness(gene_list)

        return total / len(self.reactions)

    def _parse_kegg_module_definition(self):
        '''convert a kegg module into a directed graph

        The M number entry is defined by a logical expression of K numbers
        (and other M numbers), allowing automatic evaluation of whether
        the gene set is complete, i.e., the module is present, in a given
        genome. A space or a plus sign represents an AND operation, and a
        comma sign represents an OR operation in this expression. A plus
        sign is used for a molecular complex and a minus sign designates an
        optional item in the complex.
        '''
        bracket = 0
        current_KO = ''
        current_enzyme = []
        step_definition = ''
        steps = []
        step_enzymes = []
        current_enzyme = []
        for c in self.definition:
            if c == '(':
                bracket += 1
            elif c == ' ':
                current_enzyme.append(current_KO)
                current_KO = ''
                if bracket == 0:
                    # we have reached the end of this step

                    #step_definition = ''
                    step_enzymes.append(current_enzyme)
                    current_enzyme = []
                    steps.append(step_enzymes)
                    step_enzymes = []
                    continue
                else:
                    # we are still internal to this step
                    pass
            elif c == ')':
                bracket -= 1
            elif c == '+':
                # we have a complex (multiple subunits)
                current_enzyme.append(current_KO)
                current_KO = ''

            elif c == '-':
                # we have an optional subunit comming up
                current_enzyme.append(current_KO)
                current_KO = ''
            elif c == ',':
                # this is an alternation for this step
                # this should be the end of this enzyme
                current_enzyme.append(current_KO)
                step_enzymes.append(current_enzyme)
                current_enzyme = []
                current_KO = ''
            elif c in ['K', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                # this char is part of a KO number
                current_KO += c

            #step_definition += c
        #print(steps)

class KeggOrtholog(object):
    '''Representation of a kegg ortholog
    '''
    def __init__(self, accession):
        self.entry = accession
        self.definition = []
        self.name = []
        self.enzyme = []
        self.reaction = []

    def __hash__(self):
        return hash(self.entry)

    def __eq__(self, other):
        return self.entry == other.entry

    def __str__(self):
        return str(self.entry)

class KeggEnzyme(object):
    '''Representation of a kegg enzyme entry

    This is a container class for holding multiple KEGG orthologs
    that each act as a subunit in the enzyme.
    '''

    def __init__(self, accession=None):
        '''Initialize the object

        Arguments:
            - accession: the EC number for the enzyme. eg. 1.1.1.1
        '''
        self.entry = accession
        self.orthologs = set()

    def __hash__(self):
        return hash(self.entry) + hash(''.join(self.orthologs))

    def __eq__(self, other):
        return self.entry == other.entry and self.orthologs == other.orthologs

    def __repr__(self):
        return '<KeggEnzyme {} composed of {} gene products>'.format(self.entry, len(self.orthologs))

    def __str__(self):
        info = []
        for i in self.orthologs:
            info.append(str(i))
        return self.entry + '\n\t' + '\n\t'.join(info)

    def add(self, ortholog):
        #if not isinstance(KeggOrtholog, ortholog):
        #    raise RuntimeError("You must provide a KeggOrtholog object to KeggEnzyme")
        self.orthologs.add(ortholog)

    def remove(self, ortholog):
        self.orthologs.remove(ortholog)

    def present(self, gene_list):
        '''Check if the gene_list contains all of the orthologs for this enzyme
        '''
        if len(self.orthologs) == 0:
            return False
        return self.orthologs.issubset(gene_list)

    def completeness(self, gene_list):
        orthos = len(self.orthologs)
        if orthos == 0:
            return 0.0
        p = self.orthologs - gene_list
        diff = len(self.orthologs - gene_list)
        return (orthos - diff) / orthos


class Catalyst(object):
    '''this is a container class for genes involved in a reaction

    Often a reaction can be performed by multiple different enzymes
    that may have multiple subunits. When we decide whether an
    organism can perform a certain reaction we want to know whether
    all of the subunits are present in at least one of the enzymes
    '''

    def __init__(self):
        self.enzymes = set()

    def add(self, enzyme):
        '''add a KeggEnzyme object to the catalyst
        '''
        if isinstance(enzyme, list):
            for e in enzyme:
                self.enzymes.add(e)
        else:
            self.enzymes.add(enzyme)

    def __repr__(self):
        orthologs = 0
        enzymes = 0
        for i in self.enzymes:
            orthologs += len(i.orthologs)
            enzymes += 1
        return '<Catalyst containing {} different enzymes and {} gene products>'.format(enzymes, orthologs)

    def __str__(self):
        info = []
        for i in self.enzymes:
            info.append(str(i))
        return '\n'.join(info)

    def present(self, gene_list):
        '''Given a list of gene identifiers this method will check to see
        if any of the catalysts given have all of their required
        subunits present in the list
        '''
        for i in self.enzymes:
            if i.present(gene_list):
                return True
        return False

    def completeness(self, gene_list):
        max_comp = 0.0
        for i in self.enzymes:
            m = i.completeness(gene_list)
            if m > max_comp:
                max_comp = m
        return max_comp


class KeggReaction(Reaction):
    '''Representation of a kegg reaction
    '''
    def __init__(self,reactants=None, catalysts=(),
                             reversible=0, data=None):
        super(KeggReaction, self).__init__(reactants=reactants, catalysts=catalysts,
                                 reversible=reversible, data=data)

    def present(self, gene_list):
        '''Check if this reaction can be performed using the specified genes.

        Given a list of gene identifiers this method will check to see
        if any of the catalysts given have all of their required
        subunits present in the list
        '''
        for c in self.catalysts:
            if c.present(gene_list):
                return True
        return False

    def completeness(self, gene_list):
        '''Return the maximal percentage of present genes.

        Will look through all of the enzymes and return a
        number between 0.0 and 1.0 representing the fraction
        of that enzyme's orthologs that are present in gene_list
        '''
        max_gene = 0.0
        for c in self.catalysts:
            m = c.completeness(gene_list)
            if m > max_gene:
                max_gene = m
        return max_gene


def parse_enzyme_file(handle, enzyme):
    record = enzyme
    state = ''
    def parse_orthology(line, record):
        accession, _ = line.split(None, 1)
        record.add(accession)
    for line in handle:
        if line[:3] == '///':
            return
        elif line[:5] == 'GENES':
            state = 'GENES'
        elif line[:9] == 'ORTHOLOGY':
            state = 'ORTHOLOGY'
            parse_orthology(line[9:], record)
        elif line[0] == ' ':
            if state == 'ORTHOLOGY':
                parse_orthology(line, record)


def parse_reaction_file(handle):
    state = ''
    ret = []
    def parse_enzyme(line):
        enzymes = line.strip().split()
        for e in enzymes:
            enzyme = KeggEnzyme(e)
            try:
                parse_enzyme_file(kegg_get(e), enzyme)
            except HTTPError as error:
                print("cannot get information for enzyme: {}\nsource line:\n{}\nenzyme list:\n{}".format(e, line, str(enzymes)))
            else:
                ret.append(enzyme)

    for line in handle:
        if line[:3] == '///':
            return ret
        elif line[:6] == 'ENZYME':
            state = 'ENZYME'
            parse_enzyme(line[6:])

def parse(handle):
    def parse_orthology(string, record):
        return
        match = re.search(r'([\dK,]+)\s+([\w,:/() -]+)\s+\[EC:(.*)\] \[RN:(.*)\]', string)
        if match:
            #rxns = match.group(3).split()
            #ecs = match.group(2).split()
            if ',' in match.group(1):
                # there are multiple orthologs given on this line
                kegg_ids = match.group(1).split(',')
                for i in kegg_ids:
                    kegg_ortholog = i #KeggOrtholog(i)
                    record.orthologs.append(kegg_ortholog)
            else:
                #record.orthologs.append(KeggOrtholog(match.group(1)))
                record.orthologs.append(match.group(1))
        else:
            raise KeggParseError("The line:\n{}\ncould not be parsed properly".format(line))

    def parse_pathway(string, record):
        _, map_num, name = string.strip().split(None, 2)
        pathway = ('PATH', map_num, name)
        record.pathway.append(pathway)

    def parse_reaction(string, record):
        # a map of reaction to set of orthologs
        rxn_enzymes = {}
        enzymes = {}

        rxns, eq = string.strip().split(None, 1)
        substrates, products = eq.split(' -> ')
        substrates = substrates.split(' + ')
        products = products.split(' + ')
        reactants = {}
        compounds = set()
        for substrate in substrates:
            #if substrate not in record.compounds:
            #    try:
            #        #sub = kegg_compound_read(kegg_get(substrate))
            #    except IndexError as e:
            #        print(substrate)
            #        raise e
            #    record.compounds[substrate] = sub
            reactants[substrate] = -1
        for product in products:
            #if product not in record.compounds:
            #    prod = kegg_compound_read(kegg_get(product))
            #    record.compounds[product] = prod
            reactants[product] = 1

        catalyst = Catalyst()
        for rxn in rxns.split(','):
            enzymes = []
            for r in rxn.split('+'):
                enzymes.extend(parse_reaction_file(kegg_get(rxn)))
            catalyst.add(enzymes)

        #print(catalyst)
        record.reactions.append( KeggReaction(reactants, catalysts=(catalyst,), reversible=1, data = rxns))

    def parse_compound(string, record):
        accession, name = string.strip().split(None, 2)
        comp = KeggCompound(accession)
        comp.name.append(name)
        self.compounds.append(comp)

    state = ''
    record = KeggModule()
    for line in handle:
        if line[:3] == '///':
            yield record
            record = KeggModule()
        elif line[:5] == 'ENTRY':
            try:
                pathway, accession, entry_class = line.strip().split(None, 2)
            except ValueError as e:
                print(line)
                raise e
            record.entry = accession
        elif line[:4] == 'NAME':
            record.name.append(line[4:].strip())
        elif line[:10] == 'DEFINITION':
            record.definition = line[10:].strip()
        elif line[:9] == 'ORTHOLOGY':
            state = 'ORTHOLOGY'
            parse_orthology(line[9:], record)
        elif line[:5] == 'CLASS':
            record.classname.append(line[5:].strip())
        elif line[:7] == 'PATHWAY':
            state = 'PATHWAY'
            parse_pathway(line[7:], record)
        elif line[:8] == 'REACTION':
            state = 'REACTION'
            parse_reaction(line[8:], record)
        elif line[:8] == 'COMPOUND':
            state = 'COMPOUND'
            # this should already be taken care of in the reactions
            #parse_compound(line[8:], record)

        elif line[0] == ' ':
            if state == 'ORTHOLOGY':
                parse_orthology(line, record)
            elif state == 'PATHWAY':
                parse_pathway(line, record)
            elif state == 'REACTION':
                parse_reaction(line, record)

