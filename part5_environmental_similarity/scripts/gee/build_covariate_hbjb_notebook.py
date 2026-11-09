"""Rebuild Covariate_HBJB.ipynb from reviewed source cells."""

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
NOTEBOOK_PATH = HERE / "Covariate_HBJB.ipynb"
POINT_SOURCE_PATH = HERE / "canada_external_core_sampling.py"
COMMUNITY_SOURCE_PATH = HERE / "aoi_community_core_sampling.py"


def clean_cell(cell):
    cleaned = dict(cell)
    cleaned["execution_count"] = None
    cleaned["outputs"] = []
    cleaned["metadata"] = {}
    return cleaned


def code_cell(source):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source,
    }


def markdown_cell(source):
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": source,
    }


notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))

# Retain the first ten reviewed code cells: imports, configuration, helpers,
# AOI stack and guarded AOI export definition. Markdown part headings added by
# this builder are deliberately ignored so repeat builds stay stable.
base_cells = [
    clean_cell(cell)
    for cell in notebook["cells"]
    if cell["cell_type"] == "code"
][:10]
base_cells[0]["source"] = [
    '"""Route L covariate workflow in three independent parts.\n',
    "\n",
    "The AOI raster stack is retained for visualization and reproducibility,\n",
    "but its completed raster exports are not resubmitted unless\n",
    "RUN_RASTER_EXPORTS is explicitly changed to True. Part 2 samples external\n",
    "Canadian cores with its wider lazy stack. Part 3 samples community cores\n",
    "from the AOI stack already constructed in Part 1.\n",
    '"""\n',
    "\n",
    "# In a fresh Colab runtime, run once if geemap is unavailable:\n",
    "# !pip install geemap --quiet\n",
    "\n",
    "import ee\n",
    "import geemap\n",
]

def source_cells(path):
    """Split a percent-format Python source into clean notebook code cells."""
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    cells = []
    current = []
    for line in lines:
        if line.startswith("# %%"):
            if current:
                cells.append(code_cell(current))
                current = []
            continue
        current.append(line)
    if current:
        cells.append(code_cell(current))
    return cells


point_cells = source_cells(POINT_SOURCE_PATH)
community_cells = source_cells(COMMUNITY_SOURCE_PATH)

part1_guide = markdown_cell(
    [
        "## Part 1 — AOI covariate stack and raster exports\n",
        "\n",
        "Build and inspect the AOI predictor and QA images. The raster export "
        "switch remains off because those files already exist in Cloud Storage.\n",
    ]
)

part2_guide = markdown_cell(
    [
        "## Part 2 — Canadian external-core point extraction\n",
        "\n",
        "This section never exports a Canada-wide raster. Earth Engine builds "
        "a lazy image expression and evaluates it only at the selected core "
        "points and their focal neighborhoods. Its export switch is currently "
        "off so the already-started external tasks are not duplicated. Set "
        "`RUN_POINT_EXPORTS = True` only if those tasks must be resubmitted.\n",
    ]
)

part3_guide = markdown_cell(
    [
        "## Part 3 — Community-core extraction from the AOI stack\n",
        "\n",
        "Run Part 1 before this section. The eight community cores are known "
        "to be inside the AOI, so this part samples the existing AOI-clipped "
        "`covariate_stack` and `qa_stack` directly. It does not run Part 2, "
        "rebuild a Canada-wide image, or submit a raster export. Upload "
        "`community_cores_for_gee.csv` when prompted.\n",
    ]
)

notebook["cells"] = (
    [part1_guide]
    + base_cells
    + [part2_guide]
    + point_cells
    + [part3_guide]
    + community_cells
)
notebook["nbformat"] = 4
notebook["nbformat_minor"] = 0
notebook["metadata"] = {
    "colab": {"provenance": []},
    "kernelspec": {"display_name": "Python 3", "name": "python3"},
    "language_info": {"name": "python"},
}
NOTEBOOK_PATH.write_text(
    json.dumps(notebook, indent=1, ensure_ascii=False) + "\n",
    encoding="utf-8",
)
print(f"Rebuilt {NOTEBOOK_PATH} with {len(notebook['cells'])} clean cells")
