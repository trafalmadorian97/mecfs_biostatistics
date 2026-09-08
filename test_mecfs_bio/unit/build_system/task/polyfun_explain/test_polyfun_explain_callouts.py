from pathlib import Path

import polars as pl

from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    FAMILY_COL,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    _CS_NUMBER_COL,
    FAMILY_CONTRAST_COL,
    _callout_families,
    _format_callout_label,
    _load_cs_numbers,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)


def _key(pos: int) -> dict:
    return {
        GWASLAB_CHROM_COL: 1,
        GWASLAB_POS_COL: pos,
        GWASLAB_EFFECT_ALLELE_COL: "T",
        GWASLAB_NON_EFFECT_ALLELE_COL: "A",
    }


def test_callout_families_buckets_and_orders_by_z():
    focal = _key(123)
    per_family = pl.DataFrame(
        [
            {**focal, FAMILY_COL: "conserved", FAMILY_CONTRAST_COL: 3.0},  # z=3 -> ++
            {**focal, FAMILY_COL: "coding", FAMILY_CONTRAST_COL: 1.5},  # z=1.5 -> +
            {
                **focal,
                FAMILY_COL: "repressed",
                FAMILY_CONTRAST_COL: 0.5,
            },  # z=0.5 -> drop
            {
                **focal,
                FAMILY_COL: "histone_marks",
                FAMILY_CONTRAST_COL: -4.0,
            },  # neg -> drop
        ]
    )
    family_sd = {
        "conserved": 1.0,
        "coding": 1.0,
        "repressed": 1.0,
        "histone_marks": 1.0,
    }
    result = _callout_families(per_family, focal, family_sd, max_families=3)
    assert result == [("conserved", "++"), ("coding", "+")]


def test_callout_families_drops_negligible_runner_up():
    # When one family dominates a variant's PIP lift, a runner-up whose per-family
    # contrast is below _CALLOUT_FAMILY_RELATIVE_MIN of the top contrast is dropped
    # even if it clears the background-SD bar on its own.
    focal = _key(123)
    per_family = pl.DataFrame(
        [
            {**focal, FAMILY_COL: "conserved", FAMILY_CONTRAST_COL: 10.0},
            {**focal, FAMILY_COL: "coding", FAMILY_CONTRAST_COL: 0.3},
        ]
    )
    # coding's z = 0.3 / 0.1 = 3 would normally earn a "++", but 0.3 is below 5% of
    # the dominant 10.0 contrast, so only conserved survives.
    family_sd = {"conserved": 1.0, "coding": 0.1}
    assert _callout_families(per_family, focal, family_sd, max_families=3) == [
        ("conserved", "++")
    ]


def test_callout_families_skips_degenerate_sd():
    focal = _key(123)
    per_family = pl.DataFrame(
        [{**focal, FAMILY_COL: "conserved", FAMILY_CONTRAST_COL: 3.0}]
    )
    assert _callout_families(per_family, focal, {"conserved": 0.0}, 3) == []


def test_load_cs_numbers_empty_yields_canonical_key_types(tmp_path: Path):
    # A SUSIE run that finds no credible set writes an empty combined_cs whose
    # parquet key columns default to Float64. _load_cs_numbers must still return
    # canonical Int64 keys (+ a cs_number column) so union_keys/joins line up
    # with runs that did find a credible set.
    pl.DataFrame(
        schema={
            GWASLAB_CHROM_COL: pl.Float64,
            GWASLAB_POS_COL: pl.Float64,
            GWASLAB_EFFECT_ALLELE_COL: pl.String,
            GWASLAB_NON_EFFECT_ALLELE_COL: pl.String,
        }
    ).write_parquet(tmp_path / COMBINED_CS_FILENAME)
    out = _load_cs_numbers(tmp_path)
    assert out.schema[GWASLAB_CHROM_COL] == pl.Int64
    assert out.schema[GWASLAB_POS_COL] == pl.Int64
    assert _CS_NUMBER_COL in out.columns


def test_format_label_with_and_without_families():
    focal = _key(174128548)
    assert (
        _format_callout_label(focal, [("conserved", "++"), ("coding", "+")])
        == "174128548:A:T (conserved ++, coding +)"
    )
    assert _format_callout_label(focal, []) == "174128548:A:T"
