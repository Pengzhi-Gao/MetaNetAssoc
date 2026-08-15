# MetaNetAssoc scoring specification

This document fixes the behavior of version 0.1.0 so results can be audited and
reproduced.

## 1. Evidence preparation

For metabolite `m`, duplicate database evidence for the same MPI gene `i` is
combined with a bounded noisy-OR rule:

```text
w_MPI(m,i) = 1 - product_s(1 - w_s)
```

PPI edges are undirected. Reciprocal or repeated rows for the same gene pair are
collapsed, using the maximum confidence and retaining all source labels.

## 2. Supported paths

For each metabolite `m` and term `t`, the algorithm tests all pairs between MPI
genes `A_m` and term-member genes `B_t`. A path is retained when:

1. a direct PPI edge connects the two genes; or
2. `include_overlap = TRUE` and the MPI gene is itself a term member.

The second case receives PPI confidence 1 and path type `gene_overlap`.

## 3. Edge weights

```text
edge_weight(m,i,j,t) = w_MPI(m,i) * w_PPI(i,j) * hub_factor(j)
```

The default hub factor is one. `favor` maps the term-gene log degree into the
interval 0.5-1, so connected hubs receive relatively more weight without allowing
an edge weight above one. `penalize` uses `1 / sqrt(max(1, degree))` to suppress
promiscuous PPI hubs. Direct membership always has hub factor one.

## 4. Metabolite-term scores

```text
edge_count       = number of retained paths
weighted_sum     = sum(edge_weight)
possible_pairs   = |A_m| * |B_t|
normalized_score = weighted_sum / possible_pairs
```

Two coverage diagnostics are also returned:

```text
relation_coverage = connected MPI genes / |A_m|
term_coverage     = connected term genes / |B_t|
```

`normalized_score` is bounded between zero and one. It corrects the most direct
gene-set-size bias, but it is not a calibrated probability and depends on MPI/PPI
coverage. Raw counts, normalized values, coverage, and the complete path table
should therefore be reviewed together.

## 5. Interpretation boundary

The method measures network-supported functional coupling. It does not establish
direction of effect, metabolic flux, temporal order, or causal regulation. Those
claims require expression/abundance direction, experiments, or a separately
validated causal model.

