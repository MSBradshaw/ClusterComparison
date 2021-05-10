import unittest
from BOCC import BOCC

SMALL_TEST_COMS = 'Data/three_communities.txt'
BIG_TEST_COMS = 'Data/example_communities.txt'
EXAMPLE_MEMBERS = ['SEMA4A', 'ABCA4', 'MERTK', 'HP:0000608']
EXAMPLE_TYPES = ['gene', 'gene', 'gene', 'hpo']

class ClusterTests(unittest.TestCase):
    def test_blank(self):
        c = BOCC()
        self.assertEqual(c.name, None, 'Name should be None when first initialized')
        self.assertEqual(len(c.members), 0, 'Members should be initialized as an empty list')
        self.assertEqual(len(c.members), len(c.types), 'length of members and types should always be equal')

    def test_add_members(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        self.assertEqual(len(c.members), len(EXAMPLE_MEMBERS), 'Number of members should be 4')
        self.assertEqual(len(c.types), len(EXAMPLE_TYPES), 'Number of types should be 4')
        # add more to the members
        c.add_members(['gene1','HP:000007'], ['gene','hpo'])
        self.assertEqual(len(c.members), 6, 'Number of members should be 6')
        self.assertEqual(len(c.types), 6, 'Number of types should be 6')
        # add more but without types listed
        c.add_members(['gene2', 'HP:000008'])
        self.assertEqual(len(c.members), 8, 'Number of members should be 8')
        self.assertEqual(len(c.types), 8, 'Number of types should be 8')
        del c
        for line in open(SMALL_TEST_COMS, 'r'):
            c = BOCC()
            row = line.strip().split('\t')
            c.add_members(row[1:])
            self.assertEqual(len(c.members), len(row) - 1, 'Number of members in BOCC obj doesn\'t match expectation')

    def test_get_genes(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        self.assertEqual(len(c.get_genes()), 3, 'Number of genes should be 3')
        self.assertEqual(len(c.genes), 3, 'Number of genes should be 3')
        self.assertEqual(c.genes, EXAMPLE_MEMBERS[:-1], 'Given and recieved genes do not match')

    def test_reset(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        c.get_genes()
        self.assertEqual(len(c.genes), 3, 'Number of genes should be 3')
        c.reset()
        self.assertEqual(c.genes, None, 'Genes should be set to None')

    def test_go_enrichment(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        df = c.go_enrichment()
        self.assertNotEqual(df.shape[0], 0, 'Results should not be empty')
        self.assertEqual(df.iloc[0, 6], 'GO:0006649', 'Results to not match expected')
        self.assertEqual(df.iloc[0, 5], 0.0002912762741932548, 'Results to not match expected')


if __name__ == '__main__':
    unittest.main()

