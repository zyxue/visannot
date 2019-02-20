import pandas as pd
from tqdm import tqdm

from google.cloud import firestore

# Project ID is determined by the GCLOUD_PROJECT environment variable
db = firestore.Client()


# path to csv file produced by gtf2csv. ref: https://gitlab.com/zyxue/gtf2csv-csvs
pkl = '/path/to/gtf2csv/data/Homo_sapiens.GRCh38.92.pkl'
df_gtf = pd.read_pickle(pkl)
df_gtf = df_gtf.sort_values(['seqname', 'start', 'end']).reset_index(drop=True)
df_gtf['seqname'] = df_gtf['seqname'].astype(str)
print(f'loaded {pkl}')


batch = db.batch()

batch = None
for key, row in tqdm(df_gtf.iterrows()):
    if key % 400 == 0:
        if batch is not None:
            batch.commit()
        batch = db.batch()
    doc = db.collection('gtf_entries').document(str(key))
    batch.set(doc, row.dropna().to_dict())
