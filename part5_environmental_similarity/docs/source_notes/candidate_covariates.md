# Candidate covariates to explore

This is a hypothesis inventory, not a list of established local predictors or a requirement to include every variable. Relevance can differ with carbon pool, landscape, observation depth and spatial support. Current modelling comparisons remain organic-layer depth and organic-layer carbon stock; other pools below provide context for future work.

| Candidate family | Potential relevance to explore | Pool/setting where the question is useful | First decision |
|---|---|---|---|
| AlphaEarth similarity | Do similar multisensor surface environments identify transferable cores? | Organic depth/carbon; vegetation context across settings | In default draft; use full vectors |
| NDVI and woody vegetation | Does vegetation context distinguish accumulation environments or carbon-density profiles? | Organic soils; aboveground biomass would be a separate target | NDVI in draft; review local vegetation distinctions |
| NDMI, green–NIR NDWI | Do surface/canopy moisture and open-water contrast distinguish core environments? | Wetland organic depth/carbon | In default draft; inspect season and clouds |
| Dynamic World class probabilities | Do soft tree, shrub, grass, water and flooded-vegetation distinctions improve transfer? | Landscape transitions; neither peat type nor stock directly | Optional in draft; keep all nine for QA, choose a small subset for modelling |
| Sentinel-1 VV/VH | Does radar add information beyond optical features and embeddings? | Wetness/vegetation structure in organic-soil landscapes | Optional; inspect pass/orbit effects and available scenes |
| Channel/coast distance, connectivity | Is proximity alone useful, or does hydrological connection matter more? | Floodplains, coastal landscapes, peatlands | Verify local hydrography and regional model overlap first |
| Water occurrence/seasonality | Does repeated inundation distinguish environments missed by a single season? | Wetland/floodplain soils | Candidate; historic composites may not match field year |
| Wetland/landform classes and local knowledge | Do locally recognized bog/fen/swamp, forest/open and drainage distinctions organize the evidence? | Organic soils, landscape-specific | Prioritize locally defensible mapping; broad classes are not ground truth |
| Elevation and relative landscape position | Are broad gradients informative independently of distance? | Regional context across pools | Elevation optional; test overlap with geography |
| LiDAR microrelief, valid terrain derivatives | Can reliable small elevation differences resolve drainage or erosion environments? | Peatland microforms; mineral topsoil/erosion settings | Conditional on DEM accuracy and process-appropriate flow calculations |
| Seasonal optical/radar contrasts | Does early wetness versus summer growth distinguish sites better than a pooled season? | Organic soils; crop/management systems separately | Candidate after scene counts and field dates are reviewed |
| Fire, harvest, drainage and land-use history | Does disturbance history explain differences between similar present-day surfaces? | Aboveground and soil pools, with different response times | Need dated, locally reliable records; no flux inference from stock differences |
| Parent material, texture, mineralogy | Do mineral-soil properties help explain density or stabilization? | Primarily mineral soil carbon; possible indirect hydrological context | Future mineral phase, not added to the organic-only pilot |
| Climate, permafrost, thermal regime | Do broad environmental gradients or thaw conditions matter within this domain? | Regional organic/mineral stocks; not necessarily locally constant | Inspect actual variation/support; do not resample coarse layers and call them local detail |
| Productivity, canopy height, biomass | Do inputs/vegetation structure add explanatory context? | Aboveground biomass directly; soil carbon only indirectly | Separate targets; assess redundancy and measurement support |

## How to explore without overfitting

1. Inventory coverage, source dates, native support, uncertainty and ecological rationale.
2. Examine spatial variation and redundancy without using carbon outcomes to select a winner. A variable that scarcely varies locally cannot distinguish cores; a geographic ramp may simply duplicate distance.
3. Compare core environments with the AOI distribution. Identify conditions absent from the samples, not just which covariates correlate with the few observations.
4. Prespecify a few alternative similarity definitions. Hold out whole cores from all contributing components and show every error. Variable, scale and bandwidth choices still add model-selection flexibility even without regression coefficients.
5. Retain variables only if they offer useful, interpretable and reasonably stable additional information. Do not interpret exploratory correlations as causal or statistically established controls.

With eight depth contacts and four carbon profiles, this inventory should inform a small experiment and future sampling, not a many-variable fitted model.

## Scale choices

Distinguish native sensor resolution, output grid, neighborhood radius and process scale. The draft retains 30 m cells and can append circular optical means with 100 m or 300 m radii. Those are approximately 200 m or 600 m diameters, not 100 m/300 m pixels. The neighborhood is calculated before clipping to the export area and requires a minimum valid-area fraction. Compare scales as alternatives; do not automatically include every correlated scale in a similarity metric. Fine-resolution terrain is useful only if its vertical quality supports the intended process.

## Source interpretation

The online multiscale note is inspiration, not an adopted model specification. SCORPAN is an organizing framework, not evidence that particular variables dominate here. Remote-sensing carbon flux, aboveground biomass, mineral topsoil concentration and whole organic-layer stock studies address different targets.

Dynamic World provides nine class-probability bands. Seasonal means of those probabilities are not measured fractional cover, peat thickness or carbon uncertainty. It shares Sentinel-2 information with other predictors, so it is not independent evidence. [Official product documentation](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_DYNAMICWORLD_V1).

See `alphaearth_local_pilot.md` for the reviewed AlphaEarth studies and their limits.
