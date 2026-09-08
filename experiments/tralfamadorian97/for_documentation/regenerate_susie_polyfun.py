from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr15_54_925_638 import \
    POLYFUN_EXPLAIN_CHR15_54
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr17_50_237_377 import \
    POLYFUN_EXPLAIN_CHR17_50
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr1_174_128_548 import \
    POLYFUN_EXPLAIN_CHR1_174
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [

            # POLYFUN_EXPLAIN_CHR1_174.upset_all_polyfun,
            # POLYFUN_EXPLAIN_CHR1_174.upset_cs50_polyfun,
            # POLYFUN_EXPLAIN_CHR1_174.groups_by_label["l10"].plot_svg,
            # POLYFUN_EXPLAIN_CHR1_174.groups_by_label["l10"].detailed_table,
            # POLYFUN_EXPLAIN_CHR1_174.groups_by_label["l10"].per_variant_annotation_table,
            # POLYFUN_EXPLAIN_CHR15_54.upset_all_polyfun,
            # POLYFUN_EXPLAIN_CHR15_54.upset_cs50_polyfun,
            # POLYFUN_EXPLAIN_CHR15_54.groups_by_label["l10"].plot_svg,
            # POLYFUN_EXPLAIN_CHR15_54.groups_by_label["l10"].detailed_table,
            # POLYFUN_EXPLAIN_CHR15_54.groups_by_label["l10"].per_variant_annotation_table,

            POLYFUN_EXPLAIN_CHR17_50.upset_all_polyfun,
            POLYFUN_EXPLAIN_CHR17_50.upset_cs50_polyfun,
            POLYFUN_EXPLAIN_CHR17_50.groups_by_label["l10"].plot_svg,
            POLYFUN_EXPLAIN_CHR17_50.groups_by_label["l10"].detailed_table,
            POLYFUN_EXPLAIN_CHR17_50.groups_by_label["l10"].per_variant_annotation_table,

            POLYFUN_EXPLAIN_CHR17_50.upset_all_polyfun,
            POLYFUN_EXPLAIN_CHR17_50.upset_cs50_polyfun,
            POLYFUN_EXPLAIN_CHR17_50.groups_by_label["l10"].plot_svg,
            POLYFUN_EXPLAIN_CHR17_50.groups_by_label["l10"].detailed_table,
            POLYFUN_EXPLAIN_CHR17_50.groups_by_label["l10"].per_variant_annotation_table,

        ]
    )
    # push_figures()

if __name__ == '__main__':
    go()