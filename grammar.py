import pyparsing as pp

LPAR,RPAR = map(pp.Suppress, "()")
ortholog = pp.Combine('K' + pp.Word(pp.nums, exact=5))

ortholog_group = pp.Forward()
ortholog_group <<= pp.Group(LPAR + pp.delimitedList(pp.Group(ortholog_group*(1,) & ortholog*(0,))) + RPAR) | pp.delimitedList(ortholog)
#ortholog_group <<= pp.Group(LPAR + pp.OneOrMore(ortholog_group | pp.delimitedList(ortholog)) + RPAR)
expr = pp.OneOrMore(ortholog_group)

tests = """\
        ((K00134,K00150) K00927,K11389) (K00234,K00235)
        """
expr.runTests(tests)

#steps = OneOrMore(alternating_ortholog | multi_ortholog | ortholog)

#steps.runTests('''
#        K12345 (K12345,K12345,K12345+K12345) K54321+K56789-K11111
#        ((K00134,K00150) K00927,K11389)
#(K00234+K00235+K00236+K00237,K00239+K00240+K00241-(K00242,K18859,K18860),K00244+K00245+K00246-K00247)
#''')
