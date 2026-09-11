# DFCP

Discrete Fragmentation Coagulation Processes (DFCP) fits a sequence of related
haplotype partitions across genomic loci. The active implementation is C++20;
the main executable trains a reference-panel model and can optionally impute
masked target alleles.

Build from the repository root:

```bash
./build.sh
```

Train a model on a prepared binary reference panel:

```bash
./build/impute REF.bin --init pbwt --pbwt_match_len 10
```

Add `--output_r_assign FILE` to write the final reference R assignments in the
`DFRA` format used by `eval_clusters`. To impute, provide the target sequence,
observed-loci file, and probability output as one positional group:

```bash
./build/impute REF.bin TARGET.bin observed_loci.txt probs.bin
```

The binary input format, command-line options, evaluation formats, and model
invariants are documented in [AGENTS.md](AGENTS.md). The chromosome 20 data
pipeline is documented in [scripts/1000g_phase3_v5b/README.md](scripts/1000g_phase3_v5b/README.md).

For the fastsimcoal cluster-evaluation fixture:

```bash
scripts/fsc_sim/run.sh
./build/impute data/fsc/prepared/ref.bin \
  --output_r_assign output/fsc_clusters/r_assign.bin
./build/eval_clusters \
  data/fsc/prepared/ref.bin \
  output/fsc_clusters/r_assign.bin \
  data/fsc/prepared/variant_pos.txt \
  data/fsc/ex_0_pop_1/ex_0_pop_1_1_true_trees.trees
```

The optional Beagle 4 comparison is described in `AGENTS.md`; it requires the
local source under `beagle/src` and a JDK. It is included in the FSC benchmark:

```bash
python3 scripts/fsc_sim/benchmark.py --javac javac --java java
```
