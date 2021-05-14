import requests
import json
import pandas as pd
import typing
import matplotlib.pyplot as plt


class BOCC:
    """
    Biological Ontology Cluster Comparison (BOCC) (Pronounced like Bach the musician)
    """

    def __init__(self):
        self.members = []
        self.types = []
        self.name = None
        self.genes = None
        self.go_results = None

    def reset(self) -> None:
        """
        Resets certain object values that are dependent on the content of the self.members list
        Anytime self.members is changed these values should be reset
        :return: None
        """
        self.genes = None
        self.go_results = None

    def add_members(self, mems: typing.List[str], types=None) -> None:
        """

        :param mems: List of cluster members to be appended
        :param types: List of node types associated with each of the items mems
        :return: None
        """
        self.reset()
        # if no types are given, create a list of equal length with mems to add
        if types is None:
            types = [None] * len(mems)
        self.members += mems
        self.types += types

    def get_genes(self) -> typing.List[str]:
        """
        Sets self.genes equal to a list of all the items in self.members that are listed as 'gene' in self.types
        If self.types contains only None then anything without the 'HP:' prefix is considered a gene
        :return: a list the genes in the self.members
        """
        if self.genes is not None:
            return self.genes
        # if there are no types listed assuming that anything without the HP: (denoting an HPO term) prefix is a gene
        elif self.types is None or sum(0 if x is None else 1 for x in self.types) == 0:
            self.genes = [x for x in self.members if 'HP:' not in x]
        else:
            self.genes = [self.members[i] for i in range(len(self.types)) if self.types[i] == 'gene']
        return self.genes

    def go_enrichment(self) -> pd.DataFrame:
        """
        Uses send the content of self.genes to the panther API for over representation analysis
        :return: pandas DataFrame of the panther API results, not all returned rows are significant
        """
        if self.go_results is not None:
            print('Returning Cached results')
            return self.go_results
        self.get_genes()

        g_string = ','.join(self.genes)
        url = 'http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=' + g_string + \
              '&organism=9606&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR'

        resp = requests.get(url)
        resp_obj = json.loads(resp.content)

        results = {'number_in_list': [],
                   'fold_enrichment': [],
                   'fdr': [],
                   'expected': [],
                   'number_in_reference': [],
                   'pValue': [],
                   'id': [],
                   'label': [],
                   'plus_minus': []}
        try:
            for i in range(len(resp_obj['results']['result'])):
                for key in resp_obj['results']['result'][i].keys():
                    if key != 'term':
                        results[key].append(resp_obj['results']['result'][i][key])
                    else:
                        try:
                            results['id'].append(resp_obj['results']['result'][i]['term']['id'])
                        except KeyError:
                            results['id'].append('None')
                        try:
                            results['label'].append(resp_obj['results']['result'][i]['term']['label'])
                        except KeyError:
                            results['label'].append('None')

        except KeyError:
            # there are no results, return empty dataframe
            self.go_results = pd.DataFrame()
            return pd.DataFrame()

        self.go_results = pd.DataFrame(results)
        return self.go_results


def load_clusters(input_file: str) -> typing.List[BOCC]:
    """
    Load a community file into a list of BOCC objects
    :param input_file: path to file with each community listed on one tab separated line with the first item being the
    community name/ID
    :return: list of BOCC objects
    """
    clusters = []
    for line in open(input_file, 'r'):
        # check if line is commented out
        if line[0] == '#':
            continue
        row = line.strip().split('\t')
        clusters.append(BOCC())
        clusters[-1].add_members(row[1:])
        clusters[-1].name = row[0]
    return clusters


def plot_basic_com_stats(coms: typing.List[BOCC], output: str = None, logx: bool = False) -> typing.Dict[str, typing.List]:
    """
    Take a list of communities and plot the size vs the HPO:Gene ratio of each community
    :param coms: list of BOCC objects to plot
    :param output: where to save figure, if None fig is not created or saved
    :param logx: if the x scale should be log
    :return: dictionary with keys 'ratio' and 'size' all values are lists
    """
    results = {'ratio': [], 'size': []}
    for com in coms:
        gl = len(com.get_genes())
        tl = len(com.members)
        # using pseudo counts to avoid a dividing by zero issue
        hpo_gene_ratio = ((tl + 1) - (gl + 1)) / float((gl + 1))
        results['ratio'].append(hpo_gene_ratio)
        results['size'].append(tl)
    if output is not None:
        plt.scatter(results['size'],results['ratio'])
        plt.xlabel('Size')
        if logx:
            plt.xscale('log')
        plt.ylabel('Ratio (HPO:Gene)')
        plt.savefig(output)
    return results
