import json
import sys
import pandas as pd
import requests

api_url = "https://www.pgscatalog.org/rest/"
traits = pd.read_csv("traits.csv")

out={
    'trait_name':[],
    'trait_id':[],
    'pgs_id':[],
    'nvar':[],
    'pgs_ftp':[],
    'r2_max': []
  #  'rsid_check':[]
  #  'cohort':[],
  #  'r2':[]
}


for x in traits['id']:
    print(x)
    response = requests.get(f"https://www.pgscatalog.org/rest/score/search?trait_id={x}")
    response = response.json()

    #'ancestry_distribution': {'dev': {'dist': {'EUR': 100}, 'count': 408031}, 'eval': {'dist': {'EUR': 100}, 'count': 1}}
    ## check for 'EUR':100
    # get all unique publications in response
    #if 'results' not in response.keys():
    #    print(response)
    for pgs in response['results']:
        category = list(pgs['ancestry_distribution'].keys())[0]
        if 'EUR' in pgs['ancestry_distribution'][category]['dist']:
            if pgs['ancestry_distribution'][category]['dist']['EUR'] == 100 or pgs['samples_variants'][0]['cohorts'][0]['name_short'] == 'UKB':
               #https://www.pgscatalog.org/rest/performance/search?pgs_id=PGS000001
                performance = requests.get(f"https://www.pgscatalog.org/rest/performance/search?pgs_id={pgs['id']}")
                performance = performance.json()

                check = False
                if 'results' in performance.keys():
                    r2s = []
                    for perf in performance['results']:
                        r2s.append(perf['performance_metrics']['othermetrics'][0]['estimate'])
                        if any(r2 >= 0.3 for r2 in r2s):
                            check = True
                if check:
#                    if 'rsID' in pgs_df.columns:
                #    rsid = True
#                    rsid_list.append(pgs_df['rsID'].tolist())
                    out['trait_id'].append(x)
                    out['trait_name'].append(pgs['trait_reported'])
                    out['pgs_id'].append(pgs['id'])
                    out['nvar'].append(pgs['variants_number'])
                    out['pgs_ftp'].append(pgs['ftp_harmonized_scoring_files']['GRCh38']['positions'])
                    out['r2_max'].append(max(r2s))

                    #                out['cohort'].append(perf['sampleset']['samples'][0]['cohorts'][0]['name_short'])
#                out['r2'].append(perf['performance_metrics']['othermetrics'][0]['estimate'])

df = pd.DataFrame(out)
df.to_csv('pgs_results.csv', index=False)

filtered_df = pd.DataFrame()

for group_name, group_df in df.groupby('trait_id'):
    
    if len(group_df[group_df['r2_max'] > 0.4]) > 0:
        filtered_df = pd.concat([filtered_df, group_df[group_df['r2_max'] > 0.4]])
    else:
        filtered_df = pd.concat([filtered_df, pd.DataFrame([group_df.loc[group_df['r2_max'].idxmax()]])])

filtered_df.to_csv('pgs_results_filtered.csv', index=False)


#                    ftp = pgs['ftp_harmonized_scoring_files']['GRCh38']['positions']
#                    pgs_df = pd.read_csv(ftp, compression='gzip', sep = '\t', comment = '#', low_memory=False)

rsid_list = []
for i,row in filtered_df.iterrows():
    pgs_df = pd.read_csv(row['pgs_ftp'], compression='gzip', sep = '\t', comment = '#', low_memory = False)
    if 'rsID' in pgs_df.columns:
        rsid_list.append(pgs_df['rsID'].tolist())

rsid_list = [id for sublist in rsid_list for id in sublist]
print(len(rsid_list))
rsid_set = list(set(rsid_list))
print(len(rsid_set))

with open("union_rsIDs.txt", "w") as file:
  for item in rsid_set:
    file.write(str(item) + "\n")
### download scores
#
#for scorefile in df['pgs_ftp']:
#    print(f'downloading {scorefile}')
#    pgs_df = pd.read_csv(scorefile, compression='gzip', sep = '\t', comment = '#')
#    print(pgs_df)
#    break
