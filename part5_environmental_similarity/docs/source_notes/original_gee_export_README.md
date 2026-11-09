# Route L GEE covariate imports

Place the completed Earth Engine GeoTIFF exports for the 2024–2025 optical /
Dynamic World period and 2024 AlphaEarth layer in this directory.

Expected logical groups:

- `surface_context`
- `embeddings`
- `terrain`
- `QA`

An export group may consist of several numbered GeoTIFF tiles because the GEE
workflow uses `fileDimensions` for large outputs. Keep every tile and its
original filename. These files are external inputs to the next local checking
stage; they do not replace the Li rasters or existing Route L outputs.
