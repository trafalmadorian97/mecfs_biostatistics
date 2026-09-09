from mecfs_bio.build_system.task.susie_stacked_plot_task import _credible_set_label


def test_credible_set_label_positional_numbers_by_plot_order():
    # Positional mode ignores the credible-set value and numbers 1..N in order.
    assert _credible_set_label(("L2",), 0, "positional") == "CS 1"
    assert _credible_set_label(("L3",), 1, "positional") == "CS 2"


def test_credible_set_label_susie_index_uses_the_l_index():
    # susie_index mode uses SUSIE's own L-index (parsed from the group key polars
    # hands back as a 1-tuple), so a locus whose surviving sets are L2/L3 reads
    # "CS 2"/"CS 3" -- matching the contrast task's detailed-table cs_pf/cs_u.
    assert _credible_set_label(("L2",), 0, "susie_index") == "CS 2"
    assert _credible_set_label(("L3",), 1, "susie_index") == "CS 3"
