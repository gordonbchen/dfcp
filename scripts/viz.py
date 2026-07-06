import sys
import json
import matplotlib.pyplot as plt
from argparse import ArgumentParser


p = ArgumentParser()
p.add_argument("--out_file", default="viz.png")
args = p.parse_args()

log = json.load(sys.stdin)
train_log = log["train_log"]
params = log["params"]

fig, axs = plt.subplots(ncols=4, figsize=(24, 5), layout="constrained")

axs[0].plot([l["elbo"] for l in train_log])
axs[0].set_xlabel("step")
axs[0].set_title("elbo")

axs[1].plot(params["mu_gamma"])
axs[1].set_xlabel("pos")
axs[1].set_title("mu_gamma")

axs[2].plot(params["mu_d"])
axs[2].set_xlabel("pos")
axs[2].set_title("mu_d")

for t in ("t_max", "t_expect", "t_elbo", "t_step"):
    axs[3].plot([l[t] for l in train_log], label=t)
axs[3].set_xlabel("step")
axs[3].set_ylabel("ms")
axs[3].set_title("timing")
axs[3].legend(loc=0)

plt.savefig(args.out_file)
plt.close()

