import json
from pathlib import Path, PurePath

import numpy as np
import polars as pl
import pytest
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    FAMILY_COL,
    WEIGHTS_PARQUET_FILENAME,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    DETAILED_DISPLAY_TABLE_FILENAME,
    DISP_ANNOT_PREFIX,
    DISP_CS_PF,
    DISP_CS_U,
    DISP_LIFT,
    DISP_PIP_PF,
    PER_FAMILY_CONTRAST_FILENAME,
    PER_VARIANT_ANNOTATION_TABLE_FILENAME,
    SELECTION_JSON_FILENAME,
    TOP_LINE_DISPLAY_TABLE_FILENAME,
    PolyfunExplainContrastTask,
    SecondaryPositionFromSnpid,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    CS_COLUMN,
    FILTERED_GWAS_FILENAME,
    FILTERED_LD_FILENAME,
    PIP_COLUMN,
    PIP_FILENAME,
    PRIOR_FILENAME,
    PRIOR_WEIGHT_COLUMN,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_SE_COL,
    GWASLAB_SNPID_COL,
)
from mecfs_bio.constants.polyfun_annotation_families import FAMILY_SHORT_LABELS

# Two real baseline-LF annotations from different families so family aggregation
# is exercised: Coding_UCSC_common -> coding, GERP.NS -> conserved.
_ANNOT_A = "Coding_UCSC_common"
_ANNOT_B = "GERP.NS"

_N_VARIANTS = 6
# Diffuse uniform-run PIP used by the closed-form test.
_DEFAULT_UNIFORM_PIP: tuple[float, ...] = (0.2, 0.2, 0.2, 0.2, 0.1, 0.1)
# Polyfun-run PIP is concentrated on variant 0 (focal), for both tests.
_POLYFUN_PIP: tuple[float, ...] = (0.8, 0.05, 0.05, 0.05, 0.03, 0.02)
_PRIOR_WEIGHTS: tuple[float, ...] = (8.0, 1.0, 1.0, 1.0, 1.0, 1.0)
# BETA/SE so the focal variant (row 0) is the Manhattan-panel lead; used by
# the plot task (Task 3), harmless extra columns for the contrast task.
_BETAS: tuple[float, ...] = (2.0, 0.1, 0.1, 0.1, 0.1, 0.1)
_SES: tuple[float, ...] = (0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
# The synthetic SNPID encodes a "secondary build" position offset from the hg19
# POS, mimicking gwaslab's SNPID keeping the pre-liftover coordinate.
_SECONDARY_POS_OFFSET = 1000


def _write_run_dir(
    directory: Path,
    variants: pl.DataFrame,  # CHR, POS, EA, NEA
    pip: np.ndarray,
    cs_members: dict[str, list[int]],  # cs label -> row indices
    prior_weights: np.ndarray | None,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    # pip.parquet: one PIP column, row order == variants
    pl.DataFrame({PIP_COLUMN: pip}).write_parquet(directory / PIP_FILENAME)
    gwas = variants.with_columns(
        pl.Series(name=GWASLAB_BETA_COL, values=_BETAS[: variants.height]),
        pl.Series(name=GWASLAB_SE_COL, values=_SES[: variants.height]),
        # SNPID as CHR:POS:NEA:EA with the position field offset from the hg19
        # POS, so a parsed secondary position is distinguishable from POS.
        (
            pl.col("CHR").cast(pl.String)
            + ":"
            + (pl.col("POS") + _SECONDARY_POS_OFFSET).cast(pl.String)
            + ":"
            + pl.col("NEA")
            + ":"
            + pl.col("EA")
        ).alias(GWASLAB_SNPID_COL),
    )
    gwas.write_parquet(directory / FILTERED_GWAS_FILENAME)
    # Identity LD matrix (only the plot task, Task 3, reads this).
    np.save(directory / FILTERED_LD_FILENAME, np.eye(variants.height))
    if prior_weights is not None:
        variants.with_columns(
            pl.Series(name=PRIOR_WEIGHT_COLUMN, values=prior_weights)
        ).write_parquet(directory / PRIOR_FILENAME)
    else:
        variants.with_columns(
            pl.Series(name=PRIOR_WEIGHT_COLUMN, values=np.ones(variants.height))
        ).write_parquet(directory / PRIOR_FILENAME)
    cs_rows = []
    for label, idxs in cs_members.items():
        for i in idxs:
            row = variants.row(i, named=True)
            # A failed uniform run still lists its credible-set variants even
            # when their PIP is 0 - the combined_cs writer never drops rows on
            # PIP value, so the fixture must not either.
            cs_rows.append({**row, CS_COLUMN: label, PIP_COLUMN: float(pip[i])})
    pl.DataFrame(cs_rows).write_parquet(directory / COMBINED_CS_FILENAME)


def _make_contrast_fixture(
    tmp_path: Path, uniform_pip: tuple[float, ...] = _DEFAULT_UNIFORM_PIP
) -> tuple[Path, Path, Path, Path]:
    """Build a two-run (uniform + polyfun), two-annotation locus fixture on disk.

    Returns (uni_dir, pf_dir, weights_dir, annot_path). Shared by every test in
    this module (and available for other tests in this package to build on) so
    the fixture shape stays in one place. weights_dir mirrors the real shape
    RidgeAnnotationWeightsTask.execute produces (a DirectoryAsset containing
    weights.parquet), not a bare file, so the production DirectoryAsset branch
    of _load_weights is what the tests actually exercise.
    """
    assert len(uniform_pip) == _N_VARIANTS
    variants = pl.DataFrame(
        {
            "CHR": [1] * _N_VARIANTS,
            "POS": [10, 20, 30, 40, 50, 60],
            "EA": ["A"] * _N_VARIANTS,
            "NEA": ["C"] * _N_VARIANTS,
        }
    )
    # Annotation values (by CHR/BP). BP == POS. Alleles (A1/A2) match the run
    # variants' EA/NEA so the allele-aware join lines up.
    a = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # coding: focal only
    b = np.array([2.0, 1.0, 1.0, 1.0, 1.0, 1.0])  # conserved: focal higher
    annot = pl.DataFrame(
        {
            "CHR": [1] * _N_VARIANTS,
            "BP": [10, 20, 30, 40, 50, 60],
            "SNP": [f"rs{i}" for i in range(_N_VARIANTS)],
            "A1": ["A"] * _N_VARIANTS,
            "A2": ["C"] * _N_VARIANTS,
            _ANNOT_A: a,
            _ANNOT_B: b,
        }
    )
    annot_path = tmp_path / "annot.parquet"
    annot.write_parquet(annot_path)

    weights = pl.DataFrame(
        {
            "annotation": [_ANNOT_A, _ANNOT_B],
            "gamma_raw": [3.0, 0.5],
            "gamma_standardized": [0.0, 0.0],
            "family": ["coding", "conserved"],
        }
    )
    weights_dir = tmp_path / "weights"
    weights_dir.mkdir()
    weights.write_parquet(weights_dir / WEIGHTS_PARQUET_FILENAME)

    pip_u = np.array(uniform_pip)
    pip_pf = np.array(_POLYFUN_PIP)
    uni_dir = tmp_path / "uniform"
    pf_dir = tmp_path / "polyfun"
    _write_run_dir(uni_dir, variants, pip_u, {"L1": [0, 1, 2, 3]}, prior_weights=None)
    _write_run_dir(
        pf_dir,
        variants,
        pip_pf,
        {"L1": [0]},
        prior_weights=np.array(_PRIOR_WEIGHTS),
    )
    return uni_dir, pf_dir, weights_dir, annot_path


def _build_contrast_task_and_fetch_map(
    uni_dir: Path,
    pf_dir: Path,
    weights_dir: Path,
    annot_path: Path,
    n_important_families: int = 2,
    secondary_position: SecondaryPositionFromSnpid | None = None,
) -> tuple[PolyfunExplainContrastTask, dict[str, Asset]]:
    """Build a PolyfunExplainContrastTask plus its {asset_id: Asset} fetch map
    over a fixture produced by _make_contrast_fixture. Shared by every test in
    this module and by build_synthetic_explain_inputs."""
    uni_task = FakeTask(ResultDirectoryMeta(id=AssetId("uni"), trait="t", project="p"))
    pf_task = FakeTask(ResultDirectoryMeta(id=AssetId("pf"), trait="t", project="p"))
    weights_task = FakeTask(
        ReferenceFileMeta(
            group="polyfun",
            sub_group="annotations",
            sub_folder=PurePath("ridge"),
            id=AssetId("weights"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
    )
    annot_task = FakeTask(
        SimpleFileMeta("annot", read_spec=DataFrameReadSpec(DataFrameParquetFormat()))
    )

    task = PolyfunExplainContrastTask.create(
        asset_id="contrast",
        susie_uniform_task=uni_task,
        susie_polyfun_task=pf_task,
        ridge_weights_task=weights_task,
        annotation_parquet_task=annot_task,
        n_important_families=n_important_families,
        secondary_position=secondary_position,
    )

    fetch_map: dict[str, Asset] = {
        "uni": DirectoryAsset(uni_dir),
        "pf": DirectoryAsset(pf_dir),
        # Mirrors RidgeAnnotationWeightsTask.execute's real return type.
        "weights": DirectoryAsset(weights_dir),
        "annot": FileAsset(annot_path),
    }
    return task, fetch_map


def _run_contrast_task(
    tmp_path: Path,
    uni_dir: Path,
    pf_dir: Path,
    weights_dir: Path,
    annot_path: Path,
    n_important_families: int = 2,
    secondary_position: SecondaryPositionFromSnpid | None = None,
) -> DirectoryAsset:
    """Build and execute a PolyfunExplainContrastTask over a fixture produced by
    _make_contrast_fixture. Shared by every test in this module."""
    task, fetch_map = _build_contrast_task_and_fetch_map(
        uni_dir,
        pf_dir,
        weights_dir,
        annot_path,
        n_important_families,
        secondary_position=secondary_position,
    )

    def fetch(asset_id: AssetId) -> Asset:
        return fetch_map[str(asset_id)]

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    return result


@frozen
class _ExplainInputs:
    """Shared synthetic fixture for the polyfun-explain contrast and plot
    tasks: the task objects (so a caller can wire them as deps into a
    downstream task) plus a fetch map for the contrast task's own dep set and
    the already-executed contrast task directory."""

    uni_task: Task
    pf_task: Task
    weights_task: Task
    annot_task: Task
    contrast_task: Task
    fetch_map: tuple  # tuple of (str asset_id, Asset)


def build_synthetic_explain_inputs(tmp_path: Path) -> _ExplainInputs:
    """Build the two-run, two-annotation synthetic locus fixture, run the
    contrast task once over it, and return the tasks + fetch map + executed
    contrast directory so Task 3's plot task can be wired against the exact
    same data the contrast task's tables were computed from."""
    uni_dir, pf_dir, weights_dir, annot_path = _make_contrast_fixture(tmp_path)
    contrast_task, fetch_map = _build_contrast_task_and_fetch_map(
        uni_dir, pf_dir, weights_dir, annot_path
    )

    def fetch(asset_id: AssetId) -> Asset:
        return fetch_map[str(asset_id)]

    scratch = tmp_path / "contrast_scratch"
    scratch.mkdir()
    contrast_dir = contrast_task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(contrast_dir, DirectoryAsset)

    full_fetch_map = {**fetch_map, "contrast": contrast_dir}
    return _ExplainInputs(
        uni_task=contrast_task.susie_uniform_task,
        pf_task=contrast_task.susie_polyfun_task,
        weights_task=contrast_task.ridge_weights_task,
        annot_task=contrast_task.annotation_parquet_task,
        contrast_task=contrast_task,
        fetch_map=tuple(full_fetch_map.items()),
    )


def test_contrast_closed_form(tmp_path: Path):
    inputs = build_synthetic_explain_inputs(tmp_path)
    result = dict(inputs.fetch_map)["contrast"]

    # Focal variant is the max-PIP-polyfun variant (POS 10).
    selection = json.loads((result.path / SELECTION_JSON_FILENAME).read_text())
    assert selection["focal_variant"]["pos"] == 10

    # abar_c (uniform PIP-weighted over ALL variants):
    #   sum(pip_u)=1.0; abar_A = 0.2*1/1.0 = 0.2; abar_B = (0.2*2 + 0.8*1)/1.0 = 1.2
    # Focal contrast: A: 3.0*(1-0.2)=2.4 ; B: 0.5*(2-1.2)=0.4. Coding wins.
    assert selection["important_families"][0] == "coding"

    per_family = pl.read_parquet(result.path / PER_FAMILY_CONTRAST_FILENAME)
    focal_coding = per_family.filter(
        (pl.col("POS") == 10) & (pl.col("family") == "coding")
    )["family_contrast"][0]
    assert abs(focal_coding - 2.4) < 1e-9

    # Top-line table: union of both CS (rows 0,1,2,3), sorted desc by pip_pf, and
    # carries no annotation-family columns.
    top_line = pl.read_parquet(result.path / TOP_LINE_DISPLAY_TABLE_FILENAME)
    assert top_line.height == 4
    assert top_line[DISP_PIP_PF].to_list() == sorted(
        top_line[DISP_PIP_PF].to_list(), reverse=True
    )
    assert top_line.schema["chr"] == pl.Int32
    assert top_line.schema["pos"] == pl.Int32
    assert not any(c.startswith(DISP_ANNOT_PREFIX) for c in top_line.columns)
    # Focal lift = m * pi = 6 * (8/13) ~= 3.692
    focal_lift = top_line.filter(pl.col("pos") == 10)[DISP_LIFT][0]
    assert abs(focal_lift - 6.0 * (8.0 / 13.0)) < 1e-6
    # cs columns present, focal in both runs' CS
    focal_row = top_line.filter(pl.col("pos") == 10)
    assert focal_row[DISP_CS_PF][0] == 1
    assert focal_row[DISP_CS_U][0] == 1

    # Detailed table: same rows, plus a prefixed contrast column per family. The
    # family columns carry the local contrast gamma_raw_c*(a_ic-abar_c), NOT the
    # raw scaled value, so the focal coding column equals the closed-form 2.4.
    detailed = pl.read_parquet(result.path / DETAILED_DISPLAY_TABLE_FILENAME)
    assert detailed.height == 4
    # All eleven families are emitted regardless of what this locus carries.
    annot_cols = [c for c in detailed.columns if c.startswith(DISP_ANNOT_PREFIX)]
    assert len(annot_cols) == len(FAMILY_SHORT_LABELS)
    coding_col = f"{DISP_ANNOT_PREFIX}coding"
    focal_detail = detailed.filter(pl.col("pos") == 10)
    assert abs(focal_detail[coding_col][0] - 2.4) < 1e-9
    # A family with no annotations at this locus is an all-null column.
    assert (
        detailed[f"{DISP_ANNOT_PREFIX}non_synonymous"].null_count() == detailed.height
    )


def test_contrast_uniform_all_zero_pip_uses_equal_weights(tmp_path: Path):
    # When the uniform run finds no signal (all PIPs 0), abar_c falls back to an
    # unweighted mean, so contrasts are still well-defined (not NaN).
    uni_dir, pf_dir, weights_dir, annot_path = _make_contrast_fixture(
        tmp_path, uniform_pip=(0.0,) * _N_VARIANTS
    )
    result = _run_contrast_task(tmp_path, uni_dir, pf_dir, weights_dir, annot_path)

    per_family = pl.read_parquet(result.path / PER_FAMILY_CONTRAST_FILENAME)
    assert per_family["family_contrast"].is_nan().sum() == 0
    # abar_A now = unweighted mean of a = (1,0,0,0,0,0) over ALL 6 locus
    # variants = 1/6, so focal A contrast = 3*(1 - 1/6) = 2.5.
    focal_coding = per_family.filter(
        (pl.col("POS") == 10) & (pl.col("family") == "coding")
    )["family_contrast"][0]
    assert abs(focal_coding - 2.5) < 1e-9


def test_contrast_raises_on_duplicate_annotation_allele_key(tmp_path: Path):
    # The annotation matrix is built unique on (CHR, BP, unordered-allele-key).
    # A duplicate at the same position AND same alleles would cross-multiply a
    # run's variant rows into doubled/misattributed contrast values, so the task
    # must fail fast. (A genuine multiallelic site has a distinct allele key and
    # is fine - see test_contrast_resolves_co_located_variants_by_allele.)
    uni_dir, pf_dir, weights_dir, annot_path = _make_contrast_fixture(tmp_path)
    annot = pl.read_parquet(annot_path)
    dup_row = annot.filter(pl.col("BP") == 10).with_columns(
        pl.lit("rs0_dup").alias("SNP")  # same A1/A2 -> same allele key
    )
    pl.concat([annot, dup_row], how="vertical").write_parquet(annot_path)

    with pytest.raises(ValueError):
        _run_contrast_task(tmp_path, uni_dir, pf_dir, weights_dir, annot_path)


def test_contrast_resolves_co_located_variants_by_allele(tmp_path: Path):
    # Two variants at the same (CHR, POS) with different alleles must each pick
    # up their OWN annotation row (allele-aware join), not cross-multiply or
    # mispair. Focal is variant 0 (A/C, coding a=1); a second variant (A/G) at
    # the same position carries a different coding value.
    n = 2
    variants = pl.DataFrame(
        {
            "CHR": [1, 1],
            "POS": [10, 10],
            "EA": ["A", "A"],
            "NEA": ["C", "G"],
        }
    )
    annot = pl.DataFrame(
        {
            "CHR": [1, 1],
            "BP": [10, 10],
            "SNP": ["rsAC", "rsAG"],
            "A1": ["A", "A"],
            "A2": ["C", "G"],
            _ANNOT_A: [1.0, 4.0],  # coding: A/C -> 1, A/G -> 4
            _ANNOT_B: [0.0, 0.0],
        }
    )
    annot_path = tmp_path / "annot.parquet"
    annot.write_parquet(annot_path)

    weights = pl.DataFrame(
        {
            "annotation": [_ANNOT_A, _ANNOT_B],
            "gamma_raw": [3.0, 0.5],
            "gamma_standardized": [0.0, 0.0],
            "family": ["coding", "conserved"],
        }
    )
    weights_dir = tmp_path / "weights"
    weights_dir.mkdir()
    weights.write_parquet(weights_dir / WEIGHTS_PARQUET_FILENAME)

    pip = np.array([0.6, 0.4])
    uni_dir = tmp_path / "uniform"
    pf_dir = tmp_path / "polyfun"
    _write_run_dir(uni_dir, variants, pip, {"L1": [0, 1]}, prior_weights=None)
    _write_run_dir(
        pf_dir, variants, pip, {"L1": [0, 1]}, prior_weights=np.array([2.0, 1.0])
    )

    result = _run_contrast_task(tmp_path, uni_dir, pf_dir, weights_dir, annot_path)

    # family_scaled coding = gamma_raw_coding * a. A/C variant -> 3*1 = 3;
    # A/G variant -> 3*4 = 12. Correct allele pairing gives distinct values.
    per_family = pl.read_parquet(result.path / PER_FAMILY_CONTRAST_FILENAME)
    ac = per_family.filter((pl.col("NEA") == "C") & (pl.col("family") == "coding"))
    ag = per_family.filter((pl.col("NEA") == "G") & (pl.col("family") == "coding"))
    # abar_coding (uniform PIP-weighted over both) = (0.6*1 + 0.4*4)/1.0 = 2.2
    # contrast A/C = 3*(1-2.2) = -3.6 ; A/G = 3*(4-2.2) = 5.4
    assert abs(ac["family_contrast"][0] - (-3.6)) < 1e-9
    assert abs(ag["family_contrast"][0] - 5.4) < 1e-9


def test_secondary_position_from_snpid_adds_pos_column(tmp_path: Path):
    # With a SecondaryPositionFromSnpid config, both display tables gain a
    # build-labelled position column parsed from the SNPID's position field,
    # placed immediately after pos, holding POS + the synthetic offset.
    uni_dir, pf_dir, weights_dir, annot_path = _make_contrast_fixture(tmp_path)
    result = _run_contrast_task(
        tmp_path,
        uni_dir,
        pf_dir,
        weights_dir,
        annot_path,
        secondary_position=SecondaryPositionFromSnpid(build_label="hg38"),
    )
    for filename in (TOP_LINE_DISPLAY_TABLE_FILENAME, DETAILED_DISPLAY_TABLE_FILENAME):
        table = pl.read_parquet(result.path / filename)
        assert table.columns.index("pos_hg38") == table.columns.index("pos") + 1
        for row in table.iter_rows(named=True):
            assert row["pos_hg38"] == row["pos"] + _SECONDARY_POS_OFFSET


def test_per_variant_annotation_table(tmp_path: Path):
    # The per-variant annotation table characterizes the top variants of each
    # polyfun credible set: rows are the detailed annotations, columns are the
    # selected variants (chr:pos:nea:ea, hg19), cells the raw annotation values.
    variants = pl.DataFrame(
        {
            "CHR": [1] * _N_VARIANTS,
            "POS": [10, 20, 30, 40, 50, 60],
            "EA": ["A"] * _N_VARIANTS,
            "NEA": ["C"] * _N_VARIANTS,
        }
    )
    a = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # coding: variant 0 only
    b = np.array([2.0, 1.0, 1.0, 1.0, 1.0, 1.0])  # conserved
    annot = pl.DataFrame(
        {
            "CHR": [1] * _N_VARIANTS,
            "BP": [10, 20, 30, 40, 50, 60],
            "SNP": [f"rs{i}" for i in range(_N_VARIANTS)],
            "A1": ["A"] * _N_VARIANTS,
            "A2": ["C"] * _N_VARIANTS,
            _ANNOT_A: a,
            _ANNOT_B: b,
        }
    )
    annot_path = tmp_path / "annot.parquet"
    annot.write_parquet(annot_path)

    weights = pl.DataFrame(
        {
            "annotation": [_ANNOT_A, _ANNOT_B],
            "gamma_raw": [3.0, 0.5],
            "gamma_standardized": [0.0, 0.0],
            "family": ["coding", "conserved"],
        }
    )
    weights_dir = tmp_path / "weights"
    weights_dir.mkdir()
    weights.write_parquet(weights_dir / WEIGHTS_PARQUET_FILENAME)

    # Polyfun PIP + two credible sets chosen to exercise every selection branch:
    #   L1 top=var0 (0.5): var1 (0.4) within the 0.2 gap and above the 0.1 floor
    #   -> kept; var2 (0.2) outside the gap -> dropped; var3 (0.05) below the floor
    #   -> dropped. L2 top=var4 (0.05) kept despite being below the floor (it is the
    #   set's max); var5 (0.03) dropped. Selected: var0, var1, var4.
    pip_pf = np.array([0.5, 0.4, 0.2, 0.05, 0.05, 0.03])
    pf_dir = tmp_path / "polyfun"
    _write_run_dir(
        pf_dir,
        variants,
        pip_pf,
        {"L1": [0, 1, 2, 3], "L2": [4, 5]},
        prior_weights=np.ones(_N_VARIANTS),
    )
    uni_dir = tmp_path / "uniform"
    _write_run_dir(
        uni_dir, variants, np.array(_DEFAULT_UNIFORM_PIP), {"L1": [0, 1]}, None
    )

    result = _run_contrast_task(tmp_path, uni_dir, pf_dir, weights_dir, annot_path)
    table = pl.read_parquet(result.path / PER_VARIANT_ANNOTATION_TABLE_FILENAME)

    # Leading columns are the family/annotation labels plus the per-annotation
    # gamma and alpha_bar context; the variant columns follow in (credible set,
    # descending PIP) order. Gap-/floor-excluded variants (POS 30, 40, 60) absent.
    assert table.columns[:4] == [FAMILY_COL, ANNOTATION_COL, "gamma", "alpha_bar"]
    assert table.columns[4:] == ["1:10:C:A", "1:20:C:A", "1:50:C:A"]
    # One row per detailed annotation.
    assert table.height == len(weights)
    # Raw annotation values (no gamma) per selected variant, plus gamma (the ridge
    # coefficient) and alpha_bar (uniform PIP-weighted mean over the locus). With
    # uniform PIP summing to 1: alpha_bar_A = 0.2*1 = 0.2; alpha_bar_B = 0.2*2 +
    # (0.2+0.2+0.2+0.1+0.1)*1 = 1.2.
    coding = table.filter(pl.col(ANNOTATION_COL) == _ANNOT_A)
    assert coding[FAMILY_COL][0] == "coding"
    assert coding["gamma"][0] == 3.0
    assert abs(coding["alpha_bar"][0] - 0.2) < 1e-9
    assert coding.select("1:10:C:A", "1:20:C:A", "1:50:C:A").row(0) == (1.0, 0.0, 0.0)
    conserved = table.filter(pl.col(ANNOTATION_COL) == _ANNOT_B)
    assert conserved["gamma"][0] == 0.5
    assert abs(conserved["alpha_bar"][0] - 1.2) < 1e-9
    assert conserved.select("1:10:C:A", "1:20:C:A", "1:50:C:A").row(0) == (
        2.0,
        1.0,
        1.0,
    )


def test_secondary_position_malformed_snpid_raises(tmp_path: Path):
    # A SNPID whose position field is not an integer must fail fast rather than
    # silently produce a null (mislabelled) secondary position.
    uni_dir, pf_dir, weights_dir, annot_path = _make_contrast_fixture(tmp_path)
    corrupted = pl.read_parquet(pf_dir / FILTERED_GWAS_FILENAME).with_columns(
        pl.lit("not_a_valid_snpid").alias(GWASLAB_SNPID_COL)
    )
    corrupted.write_parquet(pf_dir / FILTERED_GWAS_FILENAME)
    with pytest.raises(ValueError):
        _run_contrast_task(
            tmp_path,
            uni_dir,
            pf_dir,
            weights_dir,
            annot_path,
            secondary_position=SecondaryPositionFromSnpid(build_label="hg38"),
        )
