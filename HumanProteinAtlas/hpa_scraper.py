import urllib.request
import pandas as pd
from bs4 import BeautifulSoup
import sys

start = int(sys.argv[1])
end = int(sys.argv[2])

norm_tissues = pd.read_csv('proteins.txt', sep='\t')
url_parts = list(norm_tissues['Gene-Gene'])
# url_parts = list(set(norm_tissues['Gene'] + '-' + norm_tissues['Gene name']))
base_url = 'https://www.proteinatlas.org/'

pre_gathered_info = pd.read_csv('tissue_specificity.tsv',sep='\t',header=None)
pg_parts = set(pre_gathered_info.iloc[:,-1])

genes = []
infos = {'tissue_specificity':[],
         'single_cell_tissue_specificity':[],
         'blood_specificity':[],
         'cancer_prognosis':[],
         'predicted_location':[]}
skip_count = 0
for i,part in enumerate(url_parts[start:end]):
    if part in pg_parts:
        skip_count += 1
        continue
    if skip_count % 1000 == 0:
        print('Skip count: ' + str(skip_count))
    try:
        url = base_url + part
        uf = urllib.request.urlopen(url)
        html = uf.read()
        soup = BeautifulSoup(html, 'lxml')
        # the second table, 4th row, first value, first div's text (how to find the tissue specificity table row)
        # the second table, 4th row, first value, first div's text (how to find the tissue specificity table row)
        ts = soup.find_all('table')[1].find_all('tr')[3].find_all('td')[0].find_all('div')[0].text
        # the second table, 5th row, first value, first div's text (how to find the single cell tissue specificity row)
        scts = soup.find_all('table')[1].find_all('tr')[4].find_all('td')[0].find_all('div')[0].text
        # the second table, 6th row, first value, first div's text (how to find the blood specificity)
        bs = soup.find_all('table')[1].find_all('tr')[5].find_all('td')[0].find_all('div')[0].text
        # the second table, 7th row, first value, first div's text (how to find the cancer prognosis)
        cp = soup.find_all('table')[1].find_all('tr')[6].find_all('td')[0].find_all('div')[0].text
        # the second table, 8th row, first value, first div's text (how to find the predicted location)
        pl = soup.find_all('table')[1].find_all('tr')[7].find_all('td')[0].find_all('div')[0].text

        infos['tissue_specificity'].append(ts)
        infos['single_cell_tissue_specificity'].append(scts)
        infos['blood_specificity'].append(bs)
        infos['cancer_prognosis'].append(cp)
        infos['predicted_location'].append(pl)
        genes.append(part)
    except:
        print(part)
    if i % 1000 == 0:
        print('Saving at page: ' + str(i))
        df = pd.DataFrame(infos)
        df['Gene'] = genes
        # append to file
        df.to_csv('tissue_specificity.tsv', sep='\t', index=False, mode='a', header=False)
        genes = []
        infos = {'tissue_specificity': [],
                 'single_cell_tissue_specificity': [],
                 'blood_specificity': [],
                 'cancer_prognosis': [],
                 'predicted_location': []}



# with open('gene_names.txt','w') as file:
#     file.write(','.join(list(set(norm_tissues['Gene name']))))


quit()

genes.append(part.split('-')[1])
url = base_url + part
uf = urllib.request.urlopen(url)
html = uf.read()
soup = BeautifulSoup(html, 'lxml')
# the second table, 4th row, first value, first div's text (how to find the tissue specificity table row)
ts = soup.find_all('table')[1].find_all('tr')[3].find_all('td')[0].find_all('div')[0].text
# the second table, 5th row, first value, first div's text (how to find the single cell tissue specificity table row)
scts = soup.find_all('table')[1].find_all('tr')[4].find_all('td')[0].find_all('div')[0].text
# the second table, 6th row, first value, first div's text (how to find the blood specificity)
bs = soup.find_all('table')[1].find_all('tr')[5].find_all('td')[0].find_all('div')[0].text
# the second table, 7th row, first value, first div's text (how to find the cancer prognosis)
cp = soup.find_all('table')[1].find_all('tr')[6].find_all('td')[0].find_all('div')[0].text
# the second table, 8th row, first value, first div's text (how to find the predicted location)
pl = soup.find_all('table')[1].find_all('tr')[7].find_all('td')[0].find_all('div')[0].text

infos['tissue_specificity'].append(ts)
infos['single_cell_tissue_specificity'].append(scts)
infos['blood_specificity'].append(bs)
infos['cancer_prognosis'].append(cp)
infos['predicted_location'].append(pl)


