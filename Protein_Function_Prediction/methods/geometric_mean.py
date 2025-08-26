import pandas as pd
import numpy as np
import pickle
from scipy.stats import gmean

def process_blast_for_direct_lookup(input_file, output_file, output_pfile):
    columns = ['qseqid', 'sseqid', 'qstart', 'qend', 'qlen', 'sstart', 'send',
               'slen', 'pident', 'evalue', 'bitscore']
    
    dtype = {
        'qseqid': str,
        'sseqid': str,
        'qstart': np.uint32,
        'qend': np.uint32,
        'qlen': np.uint32,
        'sstart': np.uint32,
        'send': np.uint32,
        'slen': np.uint32,
        'pident': np.float32,
        'evalue': str,
        'bitscore': np.float32
    }

    print("🧪 Reading file...")
    df = pd.read_csv(input_file, sep='\t', names=columns, dtype=dtype)

    print("⚙️ Calculating overlap...")
    df['overlap'] = ((df['qend'] - df['qstart']) + (df['send'] - df['sstart'])) / (df['qlen'] + df['slen']) * 100

    print("📇 Creating index...")
    df.set_index(['qseqid', 'sseqid'], inplace=True)
    reverse_index = df.copy()

    print("📊 Calculating geometric means...")
    geometric_dict = {}
    processed = set()

    for (q, s), row in df.iterrows():
        if (s, q) in reverse_index.index and (q, s) not in processed:
            rev_row = reverse_index.loc[(s, q)]
            values = [row['overlap'], rev_row['overlap'], row['pident'], rev_row['pident']]
            gm = gmean(values)
            geometric_dict[(q, s)] = gm
            geometric_dict[(s, q)] = gm
            processed.add((q, s))
            processed.add((s, q))

    print("💾 Saving outputs...")
    with open(output_pfile, 'wb') as f:
        pickle.dump(geometric_dict, f)

    gm_df = pd.DataFrame([(q, s, gm) for (q, s), gm in geometric_dict.items()],
                         columns=['qseqid', 'sseqid', 'geometric_mean'])
    gm_df.to_csv(output_file, index=False, sep='\t')

    print(f"✅ Done! {len(geometric_dict)} pairwise entries saved.")

# Uso del script:
process_blast_for_direct_lookup(
    input_file="results/diamond_uniref50_morethan1.txt",
    output_file="results/geometric_mean_pairs.tsv",
    output_pfile="results/geometric_mean_pairs.pkl"
)
