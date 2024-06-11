# %%
# Get files
import os
import csv
import pandas as pd
dir = '../icp_runs/pod_vals/An'
files = os.listdir(f'{dir}')
files
# %%

df_full = None
names = []
data = []
for f in files:
    df = pd.read_csv(f"{dir}/{f}", names=['name', 'pody'])
    df_full = df if df_full is None else df_full
    print(df_full.head())
    df_full.append(df)
    names = names + df['name'].values.tolist()
    data = data + df['pody'].values.tolist()
df_full.head(24)
# %%

# df_full.to_csv(f"{dir}/pods.csv")
# %%
data
out = dict(zip(names, data))
out
# df_out = pd.DataFrame(out)
# df_out.head()

with open(f"{dir}/pods.csv", 'w', newline='') as out_file:
    spamwriter = csv.writer(out_file, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for k, v in out.items():
        spamwriter.writerow([k, v])
# %%
