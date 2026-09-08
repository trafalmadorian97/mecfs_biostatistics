from pathlib import Path

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_plot_task import (
    PLOT_PNG_FILENAME,
    PLOT_SVG_FILENAME,
    PolyfunExplainPlotTask,
    _wrap_callout_label,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.genetic_map_constants import (
    GMAP_CM_COL,
    GMAP_POS_COL,
    GMAP_RATE_COL,
)
from test_mecfs_bio.unit.build_system.task.polyfun_explain.test_polyfun_explain_contrast_task import (
    build_synthetic_explain_inputs,
)


def test_wrap_callout_label_short_unchanged_long_wrapped():
    fonts = (10.0, 9.5)
    short = "173855298:A:T (conserved ++, coding +)"
    text, _size = _wrap_callout_label(short, fonts)
    assert text == short

    long = "47731228:A:C (coding ++, ld related continuous +, open chromatin +)"
    text, _size = _wrap_callout_label(long, fonts)
    assert text == (
        "47731228:A:C\ncoding ++\nld related continuous +\nopen chromatin +"
    )


def test_plot_writes_png_and_svg(tmp_path: Path):
    inputs = build_synthetic_explain_inputs(tmp_path)

    gene_info = pl.DataFrame(
        {
            GENE_INFO_CHROM_COL: [1],
            GENE_INFO_START_COL: [5],
            GENE_INFO_END_COL: [65],
            GENE_INFO_STRAND_COL: ["+"],
            GENE_INFO_NAME_COL: ["GENE1"],
        }
    )
    gene_path = tmp_path / "genes.parquet"
    gene_info.write_parquet(gene_path)
    gene_task = FakeTask(
        SimpleFileMeta("genes", read_spec=DataFrameReadSpec(DataFrameParquetFormat()))
    )

    # hg19 genetic map covering the locus (POS 10-60); the recomb track reads
    # the rate column directly.
    genetic_map = pl.DataFrame(
        {
            "CHR": [1, 1, 1],
            GMAP_POS_COL: [10, 35, 60],
            GMAP_RATE_COL: [0.5, 2.0, 1.0],
            GMAP_CM_COL: [0.0, 0.5, 1.2],
        }
    )
    gmap_path = tmp_path / "genetic_map.parquet"
    genetic_map.write_parquet(gmap_path)
    gmap_task = FakeTask(
        SimpleFileMeta(
            "genetic_map", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
        )
    )

    plot_task = PolyfunExplainPlotTask.create(
        asset_id="plot",
        susie_uniform_task=inputs.uni_task,
        susie_polyfun_task=inputs.pf_task,
        contrast_task=inputs.contrast_task,
        ridge_weights_task=inputs.weights_task,
        annotation_parquet_task=inputs.annot_task,
        gene_info_task=gene_task,
        genetic_map_task=gmap_task,
        gene_info_pipe=IdentityPipe(),
    )

    def fetch(asset_id: AssetId) -> Asset:
        mapping = dict(inputs.fetch_map)
        mapping["genes"] = FileAsset(gene_path)
        mapping["genetic_map"] = FileAsset(gmap_path)
        return mapping[str(asset_id)]

    scratch = tmp_path / "plot_scratch"
    scratch.mkdir()
    result = plot_task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    png_path = result.path / PLOT_PNG_FILENAME
    svg_path = result.path / PLOT_SVG_FILENAME
    assert png_path.is_file()
    assert svg_path.is_file()
    assert png_path.stat().st_size > 0
    assert svg_path.stat().st_size > 0
