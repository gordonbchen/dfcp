import matplotlib.pyplot as plt
from dfcp import build_dfcp, run_dfcp, get_dfcp_parser


if __name__ == "__main__":
    p = get_dfcp_parser()
    args = p.parse_args()

    build_dfcp()

    INIT_METHODS = (
        "block",
        # ("pbwt", 2), ("pbwt", 5), ("pbwt", 10), ("pbwt", 20), ("pbwt", 50),
        ("pbwt", 5),
        "viterbi"
    )
    METRICS = (
        "t_init",
        "dfcp_impute_acc",
        "clade_iou", "mean_excess_parsimony",
        "mean_iou",
        "mean_clusters",
    )
    data = {i: {m: [] for m in METRICS} for i in INIT_METHODS}

    MASKS = (0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99)

    for init_method in INIT_METHODS:
        for mask in MASKS:
            res = run_dfcp(
                args.seq_file,

                val=0.2, mask=mask,

                tree=args.tree, variant_pos_fname=args.variant_pos_fname,
                variant_start_pos=args.variant_start_pos,
                clade_beta=args.clade_beta,

                block_init=int(init_method == "block"),
                pbwt_init=int(init_method[0] == "pbwt"),
                pbwt_match_len=init_method[1] if init_method[0] == "pbwt" else None,
                init_only=1,
            )
            for metric in METRICS:
                data[init_method][metric].append(res[metric])

    fig, axs = plt.subplots(nrows=len(METRICS), ncols=1, figsize=(12, len(METRICS) * 5),
                            layout="constrained")
    fig.get_layout_engine().set(h_pad=0.3, hspace=0.1)
    for metric, ax in zip(METRICS, axs.flat):
        for init_method in INIT_METHODS:
            ax.plot(MASKS, data[init_method][metric], label=f"{init_method}")
        ax.legend(loc="best")
        ax.set_xlabel("mask frac")
        ax.set_ylabel(metric)
    fig.suptitle(f"DFCP init methods: pre-training metrics \n seq_file: {args.seq_file}")
    fig.savefig("init.png")

