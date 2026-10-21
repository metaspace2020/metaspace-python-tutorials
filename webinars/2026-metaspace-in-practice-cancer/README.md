# METASPACE in Practice: Cancer Research with Public Data

Notebooks for the webinar *METASPACE in Practice: Cancer Research with Public Data*. They show how public cancer datasets in the [METASPACE](https://metaspace2020.org) knowledgebase can be found, compared and analyzed with the METASPACE Python client and the METASPACE converter.

## Notebooks

| Notebook | What you do |
|---|---|
| [`01_find_public_cancer_datasets.ipynb`](01_find_public_cancer_datasets.ipynb) | Search the knowledgebase for cancer datasets with keywords and collect them in a table with dataset links, project IDs and publication DOIs. |
| [`02_cancer_project_rois_anndata.ipynb`](02_cancer_project_rois_anndata.ipynb) | Download six datasets of a cancer project with the converter, label pixels with ROIs from GeoJSON files, keep the ROI pixels and merge everything into one AnnData object. |

For a general introduction to the Python client, see [`python-client/`](../../python-client/metaspace_python_client.ipynb).

## Setup

Python 3.11 or newer is required.

```bash
python -m venv .venv && source .venv/bin/activate   # or use conda
pip install -r requirements.txt
jupyter lab
```

No METASPACE account or API key is needed: all datasets used here are public.

## Files

- `roi_utils.py`: `add_rois_from_geojson` adds an `ROI` column to the `obs` table of an AnnData object, using ROIs exported from METASPACE as GeoJSON (ion image pixel coordinates).
- `data/*.geojson`: ROIs of the six datasets in notebook 2 (tumor `T`, paratumor `PT` and normal `N` regions), exported from the METASPACE web app.
- `data/cancer_datasets.csv`: the table produced by notebook 1.
- `data/cancer_rois_merged.h5ad`: the merged AnnData object produced by notebook 2.
