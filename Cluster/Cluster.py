import requests
import json
import pandas as pd


class Cluster:
    def __init__(self):
        self.members = []
        self.types = []
        self.name = None
        self.genes = None
        self.go_results = None

    def get_genes(self):
        if self.genes is not None:
            return self.genes
        elif self.types is None:
            self.genes = [ x for x in self.members if 'HP:' not in x]
        else:
            self.genes = [ self.members[i] for i in range(len(self.types)) if self.types[i] == 'gene']
        return self.genes

    def go_enrichment(self):
        if self.go_results is not None:
            print('Returning Cached results')
            return self.go_results
        self.get_genes()
        """
        :param genes: list of genes
        :return: pandas DataFrame of the panther API results, not all returned rows are significant
        """
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
