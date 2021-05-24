import unittest
from BOCC import BOCC
from BOCC import load_clusters, plot_basic_com_stats, summarize_clusters
import pickle


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
        c.add_members(['gene1', 'HP:000007'], ['gene', 'hpo'])
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

    def test_load_file(self):
        coms = load_clusters(SMALL_TEST_COMS)
        self.assertEqual(len(coms), 3, 'Incorrect number of communities loaded')

    def test_plotting(self):
        # count lines in BIG_TEST_COMS
        num_coms = sum(1 for x in open(BIG_TEST_COMS, 'r') if x[0] != '#')
        coms = load_clusters(BIG_TEST_COMS)
        # res = plot_basic_com_stats(coms, output='del.png', logx=True)
        res = plot_basic_com_stats(coms)
        self.assertEqual(num_coms, len(res['ratio']), 'Number of communities and number of ratios do not match')
        self.assertEqual(num_coms, len(res['size']), 'Number of communities and number of sizes do not match')

    def test_get_ts(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        c.get_genes()
        stuff = c.get_gene_tissue_specificities()
        self.assertEqual(len(stuff[0]), len(c.get_genes()),
                         'Length of genes and tissue specific information should match')
        self.assertEqual(len(stuff[1]), len(c.get_genes()),
                         'Length of genes and cell-type group information should match')
        self.assertEqual(len(stuff[2]), len(c.get_genes()),
                         'Length of genes and cell-type specific information should match')

    def test_disease_associated_genes(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        c.get_genes()
        stuff = c.get_diseases_associated_with_genes()
        self.assertEqual(len(stuff), len(c.get_genes()),
                         'Length of genes and tissue specific information should match')

    def test_get_disease_counts(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        c.get_genes()
        stuff = c.get_disease_counts()
        self.assertEqual(len(stuff), 6, 'Number of returned associations does not match expectation')
        self.assertEqual(stuff['Retinitis pigmentosa'], 3, 'All gene in the community (3) should be associated with RP')

    def test_summarize_disease_associations(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        stuff = c.summarize_disease_associations()
        self.assertEqual(stuff[1], len(c.get_genes()),
                         'The max association count should not be larger than the number of genes in the com')

    def test_get_cell_type_counts(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        stuff = c.get_cell_type_counts()
        self.assertEqual(stuff['Cone photoreceptor cells'], 1, 'There should be 1 Cone cell occurrence')

    def test_summarize_cell_type_specificity(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        stuff = c.summarize_cell_type_specificity()
        self.assertEqual(stuff[1], 1, 'The max association count should not be 1')

    def test_get_summary_stats(self):
        c = BOCC()
        c.add_members(EXAMPLE_MEMBERS, EXAMPLE_TYPES)
        df = c.get_summary_stats(0.001)
        self.assertEqual(df.iloc[0, :]['cluster_size'], 4, 'Cluster size does not match expectation')
        self.assertEqual(df.iloc[0, :]['gene_ratio'], 0.75, 'gene ratio not as expected')
        self.assertEqual(df.iloc[0, :]['HPO_ratio'], 0.25, 'HPO ratio not as expected')
        self.assertEqual(df.iloc[0, :]['num_sig_go_enrichment_terms'], 8, 'wrong number of expected significant terms')
        self.assertEqual(df.iloc[0, :]['go_sig_threshold'], 0.001, 'wrong threshold for significance')
        self.assertEqual(df.iloc[0, :]['max_norm_cell_type_specificity'], 1/3, 'wrong max norm cell type specificity')
        self.assertEqual(df.iloc[0, :]['max_norm_disease_specificity'], 3/3, 'wrong max norm cell type specificity')

    def test_summarize_clusters(self):
        coms = load_clusters(SMALL_TEST_COMS)
        df = summarize_clusters(coms)
        self.assertEqual(df.shape[0], len(coms), 'Incorrect number of communities loaded')
        self.assertEqual(df.shape[1], 12, 'Incorrect number of features')

if __name__ == '__main__':
    unittest.main()
