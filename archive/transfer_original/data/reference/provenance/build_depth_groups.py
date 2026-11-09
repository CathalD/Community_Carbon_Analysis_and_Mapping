"""Rebuild the four reference depth groups from the bundled prepared layers.

Run: python3 data/reference/provenance/build_depth_groups.py
No network or third-party Python libraries are needed.
"""
from pathlib import Path
from collections import defaultdict, Counter
import csv
import gzip
import hashlib
import json
import math
import statistics

ROOT = Path(__file__).resolve().parent.parent
TOL = 1e-8
VARIABLES = {'stock_kgm2': 'kg C/m2', 'soc_gkg': 'g C/kg', 'bd': 'g/cm3'}

def number(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None

def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        yield from csv.DictReader(handle)

def write(path, data, fields=None):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(data[0]))
        writer.writeheader()
        writer.writerows(data)

def scale(row, key, multiplier=1):
    value = number(row.get(key))
    return None if value is None else value * multiplier

specs = [
    ('NRCan peat', 'canpeat_layers_v1.csv.gz', 'upper_depth', 'lower_depth', 'OrgC', 10, 'BDOD', 'TOTC_Stock', 0.1),
    ('WoSIS', 'wosis_canada_layers_v2.csv.gz', 'upper_depth', 'lower_depth', 'OrgC', 1, 'BDOD', 'OrgC_stock_kgm2', 1),
    ('NPDB', 'npdb_canada_layers_v1_upload_fixed.csv.gz', 'top_depth_cm', 'bot_depth_cm', 'carbon_pct', 10, 'bulk_density', 'stock_Mgha', 0.1),
    ('Janousek', 'janousek_layers_v1.csv.gz', 'upper_depth', 'lower_depth', 'OrgC_pct', 10, 'BDOD', 'OrgC_Stock_kgm2', 1),
    ('Riley', 'riley_layers_v1.csv.gz', 'upper_depth_cm', 'lower_depth_cm', 'organic_carbon_pct', 10, None, None, 1),
]
profiles = defaultdict(list)
bad_profiles = {}
invalid_depth_rows = Counter()
for dataset, filename, top_key, bottom_key, soc_key, soc_scale, bd_key, stock_key, stock_scale in specs:
    for source in rows(ROOT / 'cores' / filename):
        uid = source['profile_uid']
        profiles[(dataset, uid)]
        top, bottom = number(source[top_key]), number(source[bottom_key])
        if top is None or bottom is None or bottom <= top:
            invalid_depth_rows[(dataset, uid)] += 1
            continue
        if dataset == 'Janousek':
            thickness = number(source.get('layer_thickness_cm'))
            if thickness is None or abs(thickness - (bottom-top)) > 0.01:
                bad_profiles[(dataset, uid)] = 'Source layer thickness disagrees with boundaries'
        profiles[(dataset, uid)].append(dict(top=top, bottom=bottom,
            soc_gkg=scale(source, soc_key, soc_scale),
            bd=scale(source, bd_key) if bd_key else None,
            stock_kgm2=scale(source, stock_key, stock_scale) if stock_key else None))

metadata = {r['profile_uid']: r for r in rows(ROOT / 'profile_metadata.csv.gz')}
values, eligibility = [], []
for (dataset, uid), layers in sorted(profiles.items()):
    if not layers:
        for interval in ['0-30','0-50','0-100','total']:
            for variable in VARIABLES:
                eligibility.append(dict(dataset=dataset,profile_uid=uid,interval=interval,variable=variable,
                    eligible=False,reason='No valid source layer boundaries',datum='surface',
                    sampled_top_cm=None,sampled_bottom_cm=None,complete_target_depth=False))
        continue
    # Retain the previous SOD surface convention: shift NPDB organic-inclusive
    # profiles to their recorded profile top. Other collections keep supplied depths.
    if dataset == 'NPDB':
        shift = min(layer['top'] for layer in layers)
        layers = [dict(layer, top=layer['top']-shift, bottom=layer['bottom']-shift) for layer in layers]
    layers.sort(key=lambda layer: (layer['top'], layer['bottom']))
    top = min(layer['top'] for layer in layers)
    bottom = max(layer['bottom'] for layer in layers)
    meta = metadata[uid]
    for interval, start, target in [('0-30',0,30),('0-50',0,50),('0-100',0,100),('total',top,bottom)]:
        end=min(target,bottom)
        complete_target=interval=='total' or bottom>=target-TOL
        reason = bad_profiles.get((dataset, uid), '')
        if interval != 'total' and not reason:
            reaches = bottom>=30-TOL if interval=='0-30' else bottom>({'0-50':30,'0-100':50}[interval])+TOL
            if not reaches:
                reason = 'Profile does not meet the group entry depth'
        selected = [(layer, max(start,layer['top']), min(end,layer['bottom']))
                    for layer in layers if min(end,layer['bottom']) > max(start,layer['top'])]
        if not reason:
            continuous = bool(selected) and abs(selected[0][1]-start)<=TOL and abs(selected[-1][2]-end)<=TOL
            continuous = continuous and all(abs(selected[i][1]-selected[i-1][2])<=TOL for i in range(1,len(selected)))
            if not continuous:
                reason = 'Gap, overlap, or missing coverage in reporting interval'
        for variable, unit in VARIABLES.items():
            why = reason
            if not why:
                measurements=[part[0][variable] for part in selected]
                if any(value is None for value in measurements):
                    why='Unverified Riley density; no stock calculation' if dataset=='Riley' and variable in ('bd','stock_kgm2') else 'Missing required layer measurement'
                elif any(value<0 or (variable=='bd' and value<=0) for value in measurements):
                    why='Nonphysical required layer measurement'
            eligibility.append(dict(dataset=dataset,profile_uid=uid,interval=interval,variable=variable,
                eligible=not bool(why),reason=why,datum='surface',sampled_top_cm=top,sampled_bottom_cm=bottom,
                complete_target_depth=complete_target))
            if why:
                continue
            if variable=='stock_kgm2':
                value=sum(layer[variable]*(hi-lo)/(layer['bottom']-layer['top']) for layer,lo,hi in selected)
            else:
                value=sum(layer[variable]*(hi-lo) for layer,lo,hi in selected)/(end-start)
            values.append(dict(dataset=dataset,profile_uid=uid,interval=interval,datum='surface',
                variable=variable,unit=unit,value=value,top_cm=start,bottom_cm=end,
                sampled_top_cm=top,sampled_bottom_cm=bottom,n_layers=len(selected),
                complete_target_depth=complete_target,
                excluded_invalid_depth_rows=invalid_depth_rows[(dataset,uid)],
                ecosystem=meta['ecosystem'],general_region=meta['general_region'],
                source_region=meta['source_region'],study_group=meta['study_group'],
                lon=meta['lon'],lat=meta['lat']))

assert len({(r['profile_uid'],r['interval'],r['variable'],r['datum']) for r in values})==len(values)
assert all(math.isfinite(r['value']) for r in values)
write(ROOT/'reference_depth_values.csv.gz', values)
write(ROOT/'profile_eligibility.csv.gz', eligibility)

def quantile(numbers, p):
    ordered=sorted(numbers)
    index=(len(ordered)-1)*p
    low=math.floor(index);high=math.ceil(index)
    return ordered[low]+(ordered[high]-ordered[low])*(index-low)

for dimension in ['dataset','study_group','ecosystem','general_region','source_region']:
    groups=defaultdict(list)
    for value in values:
        group=value[dimension] or 'Not recorded'
        groups[(value['dataset'],value['interval'],value['variable'],group)].append(value)
    summary=[]
    for (dataset,interval,variable,group), parts in sorted(groups.items()):
        x=[part['value'] for part in parts]
        coordinates={(round(float(part['lon']),5),round(float(part['lat']),5)) for part in parts if number(part['lon']) is not None and number(part['lat']) is not None}
        sd=statistics.stdev(x) if len(x)>1 else None
        mean=statistics.mean(x)
        summary.append(dict(dataset=dataset,variable=variable,interval=interval,datum='surface',
            unit=VARIABLES[variable],group=group,n_profile_records=len(x),
            recorded_coordinate_pairs=len(coordinates),mean=mean,median=statistics.median(x),
            sample_sd=sd,cv_pct=100*sd/mean if sd is not None and mean else None,
            q25=quantile(x,0.25),q75=quantile(x,0.75),min=min(x),max=max(x),
            n_complete_target=sum(part['complete_target_depth'] for part in parts),
            n_partial_target=sum(not part['complete_target_depth'] for part in parts),
            min_sampled_bottom_cm=min(part['sampled_bottom_cm'] for part in parts),
            max_sampled_bottom_cm=max(part['sampled_bottom_cm'] for part in parts)))
    write(ROOT/f'summary_by_{dimension}.csv', summary)

counts=Counter((r['dataset'],r['interval'],r['variable']) for r in values)
count_rows=[dict(dataset=dataset,interval=interval,variable=variable,eligible_profiles=counts[(dataset,interval,variable)])
    for dataset in sorted(r[0] for r in specs) for interval in ['0-30','0-50','0-100','total'] for variable in VARIABLES]
write(ROOT/'depth_group_counts.csv',count_rows)

settings=dict(depth_groups=['0-30','0-50','0-100','total'],datum='surface',
    eligible_depth_rule={'0-30':'Reaches or exceeds 30 cm; capped at 30',
        '0-50':'Extends beyond 30 cm; capped at 50 or sampled bottom',
        '0-100':'Extends beyond 50 cm; capped at 100 or sampled bottom',
        'total':'All sampled depths; valid measurements required'},
    coverage_rule='Continuous measured coverage; no gaps/overlaps beyond floating point tolerance',
    coverage_tolerance_cm=TOL,no_extrapolation=True,stock_method='Source layer stock allocated by depth overlap',
    total_definition='All continuous sampled material from recorded top to sampled bottom; not confirmed whole-soil-profile stock',
    source_profile_count=len(profiles),eligible_value_rows=len(values),invalid_source_profiles=len(bad_profiles),
    excluded_invalid_depth_rows=sum(invalid_depth_rows.values()),
    snapshot_date='2026-10-30')
(ROOT/'provenance/depth_group_settings.json').write_text(json.dumps(settings,indent=2)+'\n')
print(json.dumps(settings,indent=2))
for dataset in sorted(r[0] for r in specs):
    print(dataset,{interval:counts[(dataset,interval,'stock_kgm2')] for interval in ['0-30','0-50','0-100','total']})
