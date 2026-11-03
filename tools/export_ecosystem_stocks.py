"""Create the single ecosystem stock CSV from the prepared reference values."""
from pathlib import Path
import csv
import gzip
import hashlib
import json

folder=Path(__file__).resolve().parent.parent/'data/reference'
source=folder/'reference_depth_values.csv.gz'
with (folder/'resources.csv').open(encoding='utf-8-sig',newline='') as f:
    resources={row['dataset']:row for row in csv.DictReader(f)}
with gzip.open(source,'rt',encoding='utf-8-sig',newline='') as f:
    original=[row for row in csv.DictReader(f) if row['variable']=='stock_kgm2']
rows=[]
for row in original:
    resource=resources[row['dataset']]
    rows.append(dict(dataset=row['dataset'],source=resource['resource_name'],source_url=resource['source_url'],
        profile_uid=row['profile_uid'],ecosystem=row['ecosystem'],interval=row['interval'],
        stock_kg_m2=row['value'],unit=row['unit'],datum=row['datum'],top_cm=row['top_cm'],
        bottom_cm=row['bottom_cm'],sampled_top_cm=row['sampled_top_cm'],sampled_bottom_cm=row['sampled_bottom_cm'],
        complete_target_depth=row['complete_target_depth'],n_layers=row['n_layers'],
        general_region=row['general_region'],source_region=row['source_region'],study_group=row['study_group'],
        longitude=row['lon'],latitude=row['lat'],excluded_invalid_depth_rows=row['excluded_invalid_depth_rows']))
target=folder/'ecosystem_carbon_stocks.csv'
with target.open('w',encoding='utf-8',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
report=dict(source='data/reference/reference_depth_values.csv.gz',rows=len(rows),
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    output='data/reference/ecosystem_carbon_stocks.csv',output_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
    operation='Select stock values; preserve sources, ecosystems, depth groups and sampled depth flags')
(folder/'provenance/ecosystem_stock_export.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'Created {target.name}: {len(rows):,} stock records.')
