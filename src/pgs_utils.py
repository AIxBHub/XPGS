import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry


def setup_session():
    """Set up requests session with retry strategy."""
    retries = Retry(total=20, backoff_factor=0.1)
    session = requests.Session()
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session


def query_pgs_catalog(trait_id, session):
    """Query PGS catalog for models associated with a trait ID."""
    url = f"https://www.pgscatalog.org/rest/score/search?trait_id={trait_id}"
    response = session.get(url)
    return response.json()


def filter_by_ancestry(pgs_models):
    """
    Filter PGS models for European ancestry (100% EUR) or UKB cohorts.
    """
    filtered = []
    for pgs in pgs_models.get('results', []): ## details of each model are denoted by "results" in json list
        category = list(pgs['ancestry_distribution'].keys())[0] ## need the name of first element in ancestry dist
        ancestry_dist = pgs['ancestry_distribution'][category]['dist'] ## eg. {'dist': {'EUR':100}}
        
        ## check for 'EUR':100
        is_eur = 'EUR' in ancestry_dist and ancestry_dist['EUR'] == 100
        
        ## or UKB in cohorts - only check samples_variants (original behavior)
        is_ukb = False
        if pgs['samples_variants'] and pgs['samples_variants'][0]['cohorts']:
            is_ukb = pgs['samples_variants'][0]['cohorts'][0]['name_short'] == 'UKB'

        ## keep only Eur and UKB denoted ancestry/cohorts 
        if is_eur or is_ukb:
            filtered.append(pgs)
    
    return filtered


def filter_by_performance(pgs_models, session, min_r2=0.3):
    """
    Takes a list of PGS model json data and filters by performance metric R2
    """
    filtered = []
    
    for pgs in pgs_models:
        ## Get performance metrics (R2) for each PGS model, extracts ID from model json
        url = f"https://www.pgscatalog.org/rest/performance/search?pgs_id={pgs['id']}"
        response = session.get(url)
        performance_data = response.json()

        r2_values = []
        if 'results' in performance_data.keys():
            for perf in performance_data['results']:
                r2_values.append(perf['performance_metrics']['othermetrics'][0]['estimate'])
        
        ## Some models will have multiple R2 for tests on different cohorts
        ## we want any model that passes filtering in at least one test
        ## and the we can store the max R2 for our records
        if r2_values and any(r2 >= min_r2 for r2 in r2_values):
            pgs['r2_max'] = max(r2_values)
            filtered.append(pgs)
    
    return filtered


def extract_variants_from_scoring_file(ftp_url):
    """
    Extract variant information from PGS scoring file.
    """
    dtype_dict = {
        'rsID': str,
        'chr_name': str,
        'chr_position': 'Int64',
        'hm_chr': str,
        'hm_pos': 'Int64',
        'hm_rsID': str,
    }
    
    try:
        pgs_df = pd.read_csv(ftp_url, 
                           compression='gzip', 
                           sep='\t', 
                           comment='#', 
                           low_memory=False, 
                           dtype=dtype_dict, 
                           na_values=['NA', 'na', ''])
        
        #TODO this drops models that are missing rsIDs - which there are at least one, 
        # would be better to solve to keep models that have only coordinates, and map rsID?
        if 'rsID' not in pgs_df.columns: 
            return pd.DataFrame()
        
        variant_df = pd.DataFrame()
        variant_df['rsid'] = pgs_df['rsID']
        ## nested get, if chr_name isn't there, get for 'hm_chr'
        ## some score files use different column names
        variant_df['chr'] = pgs_df.get('chr_name', pgs_df.get('hm_chr', pd.NA))
        variant_df['start_GRCh37'] = pgs_df.get('chr_position', pd.NA)
        variant_df['start_GRCh38'] = pgs_df.get('hm_pos', pd.NA)
        
        return variant_df
    
    except Exception as e:
        print(f"Error processing {ftp_url}: {e}")
        return pd.DataFrame()


def process_traits(trait_csv_path, output_dir="data"):
    """Main function to process all traits and generate outputs."""
    session = setup_session()
    traits = pd.read_csv(trait_csv_path)
    
    results = []
    all_variants = []
    
    for trait_id in traits['id']:
        """
        Loop all trait ids in the sample sheet and pull matching pgs models.
        """
        print(f"Processing: {trait_id}")

        ## gets all models for a specific trait id and returns parsed json 
        pgs_data = query_pgs_catalog(trait_id, session)

        ## keep only models trained with Eur ancest or UKB cohorts 
        ancestry_filtered = filter_by_ancestry(pgs_data)
        
        ## keep only high performing models based on R2
        performance_filtered = filter_by_performance(ancestry_filtered, session)
        
        ## Store results
        for pgs in performance_filtered:
            results.append({
                'trait_id': trait_id,
                'trait_name': pgs['trait_reported'],
                'pgs_id': pgs['id'],
                'nvar': pgs['variants_number'],
                'pgs_ftp': pgs['ftp_harmonized_scoring_files']['GRCh38']['positions'],
                'r2_max': pgs['r2_max']
            })
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(f'{output_dir}/pgs_results.csv', index=False)
    
    ## apply a second round of R2 filtering to get a smaller union set of variants for bct
    ## keeps anything over R2 0.4 or the best model for each BCT
    filtered_df = pd.DataFrame()
    for trait_id, group in results_df.groupby('trait_id'):
        high_r2 = group[group['r2_max'] > 0.4]
        if len(high_r2) > 0:
            filtered_df = pd.concat([filtered_df, high_r2])
        else:
            best = group.loc[group['r2_max'].idxmax()]
            filtered_df = pd.concat([filtered_df, pd.DataFrame([best])])
    
    filtered_df.to_csv(f'{output_dir}/pgs_results_filtered.csv', index=False)
    
    ## The next step is to actually download the scoring files for each of the filtered models
    print("Extracting variants from scoring files...")
    for _, row in filtered_df.iterrows():
        variant_df = extract_variants_from_scoring_file(row['pgs_ftp'])
        if not variant_df.empty:
            all_variants.append(variant_df)
    
    ## Create union of all variants
    if all_variants:
        ## Filter out empty DataFrames and those with all-NA entries to avoid FutureWarning -- but still getting warning
        valid_variants = [df for df in all_variants if not df.empty and 'rsid' in df.columns and df['rsid'].notna().any()]
        
        union_variants = pd.concat(valid_variants, ignore_index=True)
        union_variants = union_variants.drop_duplicates(subset='rsid').reset_index(drop=True)
        union_variants = union_variants.sort_values(['chr', 'start_GRCh38'])
        union_variants.to_csv(f'{output_dir}/union_rsID.csv', index=False)
        
        print(f"Total unique variants: {len(union_variants)}")
    
    return results_df, filtered_df



if __name__ == "__main__":
    import os

    outdir = 'test'
    traits_sheet = 'sheets/traits.csv'

    os.makedirs(outdir, exist_ok=True)
    
    results_df, filtered_df = process_traits(traits_sheet, outdir)