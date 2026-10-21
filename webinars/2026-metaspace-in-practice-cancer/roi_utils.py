"""Helpers for labelling AnnData pixels with ROIs exported from METASPACE as GeoJSON."""

import json
import warnings
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import shapely
from anndata import AnnData
from shapely.geometry import shape
from shapely.validation import make_valid

X_COL = "ion_image_pixel_x"
Y_COL = "ion_image_pixel_y"
METASPACE_KEY = "metaspace"


def load_roi_geometries(geojson: Union[str, Path, dict], label_property: str = "name") -> dict:
    """
    Read a METASPACE ROI export and return ``{label: shapely geometry}``.

    Features that share a label are merged. Invalid geometries (for example self-intersecting
    polygons drawn by hand) are repaired. Non-polygon features are skipped with a warning.
    """
    if isinstance(geojson, (str, Path)):
        with open(geojson) as f:
            geojson = json.load(f)

    geometries: dict = {}
    for i, feature in enumerate(geojson.get("features", [])):
        geometry = shape(feature["geometry"])
        if geometry.geom_type not in ("Polygon", "MultiPolygon"):
            warnings.warn(f"Skipping feature {i}: unsupported geometry {geometry.geom_type}.", stacklevel=3)
            continue
        if not geometry.is_valid:
            geometry = make_valid(geometry)
        label = (feature.get("properties") or {}).get(label_property) or f"ROI_{i + 1}"
        geometries[label] = geometries[label].union(geometry) if label in geometries else geometry
    return geometries


def add_rois_from_geojson(
    adata: AnnData,
    geojson: Union[str, Path, dict],
    key: str = "ROI",
    label_property: str = "name",
    pixel_offset: float = 0.0,
    overlap: str = "error",
) -> AnnData:
    """
    Add a categorical ``adata.obs[key]`` column labelling each pixel with the ROI containing it.

    The GeoJSON must be in ion image pixel coordinates (as exported by METASPACE), where x is
    the column and y the row of ``adata.obs[["ion_image_pixel_x", "ion_image_pixel_y"]]``.
    Pixels outside every ROI are NaN. Pixels on a ROI border count as inside. All pixels stay
    in ``adata``; subset on the column afterwards, e.g. ``adata[adata.obs["ROI"] == "T_1"]``.

    Args:
        adata: An AnnData from ``metaspace_to_anndata``. Modified in place.
        geojson: Path to a GeoJSON file or an already loaded dict.
        key: Name of the ``obs`` column to create.
        label_property: Feature property holding the ROI name.
        pixel_offset: Added to the pixel indices before testing them against the polygons. Use
            0.5 if the file's coordinates refer to pixel corners rather than pixel indices.
            METASPACE exports pixel indices, so the default 0.0 is right for its files.
        overlap: What to do with pixels inside more than one ROI. "error" raises. "first" gives
            them to the first ROI in the file. Adjacent ROIs can share a border pixel because
            METASPACE counts border pixels as inside.

    Returns:
        The same AnnData, for chaining.

    Raises:
        ValueError: If ROIs overlap and ``overlap="error"``, or the file's image size does not
            match the AnnData.
    """
    if overlap not in ("error", "first"):
        raise ValueError("overlap must be 'error' or 'first'")
    if isinstance(geojson, (str, Path)):
        with open(geojson) as f:
            geojson = json.load(f)
    _check_matches_adata(adata, geojson.get("metadata", {}))
    geometries = load_roi_geometries(geojson, label_property=label_property)

    points = shapely.points(
        adata.obs[X_COL].to_numpy() + pixel_offset, adata.obs[Y_COL].to_numpy() + pixel_offset
    )
    inside = {}
    for label, geometry in geometries.items():
        shapely.prepare(geometry)
        inside[label] = shapely.covers(geometry, points)

    n_rois_per_pixel = np.sum(list(inside.values()), axis=0) if inside else np.zeros(adata.n_obs)
    n_shared = int((n_rois_per_pixel > 1).sum())
    if n_shared and overlap == "error":
        overlapping = [(a, b) for a in inside for b in inside if a < b and (inside[a] & inside[b]).any()]
        raise ValueError(
            f"ROIs overlap in {n_shared} pixels: {overlapping}. "
            "Use overlap='first' to give shared pixels to the first ROI in the file."
        )
    elif n_shared:
        warnings.warn(
            f"{n_shared} pixel(s) lie in more than one ROI; assigned to the first ROI.", stacklevel=2
        )

    labels = pd.Series(np.nan, index=adata.obs.index, dtype=object)
    # Assign in reverse so that the first ROI in the file wins for shared pixels.
    for label, mask in reversed(list(inside.items())):
        labels[mask] = label
    adata.obs[key] = pd.Categorical(labels, categories=list(geometries))
    return adata


def _check_matches_adata(adata: AnnData, metadata: dict) -> None:
    # The METASPACE export records which dataset and image size it was drawn on.
    uns = adata.uns.get(METASPACE_KEY, {})
    dataset_id, size = uns.get("id"), uns.get("image_size")
    if dataset_id and metadata.get("datasetId") and metadata["datasetId"] != dataset_id:
        warnings.warn(
            f"GeoJSON was exported for dataset {metadata['datasetId']}, AnnData is {dataset_id}.",
            stacklevel=3,
        )
    if size and metadata.get("imageWidth") and metadata.get("imageHeight"):
        if (metadata["imageWidth"], metadata["imageHeight"]) != (size["x"], size["y"]):
            raise ValueError(
                f"GeoJSON image size {metadata['imageWidth']}x{metadata['imageHeight']} does not "
                f"match the AnnData image size {size['x']}x{size['y']}."
            )
