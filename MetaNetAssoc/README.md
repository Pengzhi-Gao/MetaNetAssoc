# MetaNetAssoc

## Version 0.2.0: DEG/HMDB workflow

For a step-by-step tutorial after installation, run
`vignette("deg-hmdb-workflow", package = "MetaNetAssoc")`. A fully offline
script is also provided at `system.file("examples", "deg_hmdb_workflow.R",
package = "MetaNetAssoc")`.

## Version 0.3.0: packaged MPI references

Human and mouse MPI references are now distributed with the package:

```r
human_mpi <- load_mpi_reference("human")
mouse_mpi <- load_mpi_reference("mouse")
```

They use `metabolite_kegg_id` (KEGG compound IDs). The human reference contains
the project-integrated KEGG, Reactome, Human-GEM, and BRENDA evidence. The mouse
reference is an Ensembl Compara orthology projection of that evidence. To use
these tables in the HMDB-first workflow, first add a documented HMDB--KEGG
cross-reference and expose an `hmdb_id` column.

The intended workflow is: (1) `build_ppi_network()` queries STRING from
differential genes (or accepts a downloaded STRING table); (2)
`build_mpi_network()` filters an evidence-traceable HMDB--gene table by the
differential genes and HMDB IDs; (3) `run_deg_enrichment()` calls
clusterProfiler and can write its raw result table with `export_file`; and (4)
`run_metanet()` scores the integrated networks. Human and mouse are supported.

```r
# deg: character vector or data.frame with a gene column; hmdb: HMDB IDs
ppi <- build_ppi_network(deg, species = "human")
mpi <- build_mpi_network(deg, hmdb, species = "human", mpi_database = my_mpi_db)
enr <- run_deg_enrichment(deg, species = "human", export_file = "go_results.csv")
fit <- run_metanet(hmdb, mpi, ppi, enr$enrichment, focus_genes = deg)
```

For a one-call variant, provide `deg` to `run_metanet()`; it will build missing
PPI/MPI/enrichment inputs. MPI evidence is never silently invented: provide an
HMDB--gene reference table for real analyses.

If no MPI rows are found, inspect `attr(mpi, "matching_diagnostics")`. It
reports how many DEG symbols and HMDB IDs occur in the reference table; a
non-zero result requires at least one *same-metabolite--same-DEG* edge.

`MetaNetAssoc` 是一个将候选代谢物连接到特征性生物学过程的 R 包。它使用
MPI（代谢物—基因）关系和 PPI（蛋白—蛋白）网络，构建可追溯的四层路径：

```text
Metabolite -> Relation_Gene -> Enrich_Gene -> Term
```

结果表示网络证据支持的潜在关联，而不是代谢物调控生物学过程的因果证明。

## 安装

```powershell
R CMD INSTALL D:/codex/MPI_codeBase/MetaNetAssoc
```

也可以在 R 会话中使用：

```r
# install.packages("devtools")
devtools::install("D:/codex/MPI_codeBase/MetaNetAssoc")
```

## 五分钟示例

```r
library(MetaNetAssoc)

example <- load_example_data()

fit <- run_metanet(
  metabolites = example$metabolites$metabolite_id,
  focus_genes = example$focus_genes,
  mpi = example$mpi,
  ppi = example$ppi,
  enrichment = example$enrichment
)

fit
head(result_scores(fit))
head(result_edges(fit))
result_qc(fit)

# 可选依赖 ggplot2
plot_term_ranking(fit, "GO:SIM0003")
plot_score_heatmap(fit)

# 可选依赖 igraph
plot_association_network(
  fit,
  metabolites = "Glutathione",
  terms = "Response to oxidative stress"
)
```

示例数据模拟生殖细胞阶段的代谢组—转录组联合分析。所有 `SIM...` 代谢物 ID、
`GO:SIM...` Term ID、MPI/PPI 关系和置信度均为教程用合成数据，不可作为生物学证据。

## 输入格式

建议先使用三个 `prepare_*()` 函数，把不同来源的列名转换为内部标准格式。

### MPI

| 标准列 | 含义 |
|---|---|
| `metabolite_id` | 推荐使用 HMDB ID；示例中为合成 ID |
| `metabolite_name` | 代谢物名称 |
| `relation_gene` | 代谢物关联基因符号 |
| `mpi_source` | 证据数据库 |
| `mpi_weight` | 0–1 关联置信度 |

```r
mpi <- prepare_mpi(
  my_mpi,
  metabolite_id = "HMDB_ID",
  metabolite_name = "Met_Name",
  gene = "Relation_Gene",
  source = "From_Database",
  weight = "Evidence_Score"
)
```

### PPI

标准列为 `gene_a`, `gene_b`, `ppi_source`, `ppi_score`。PPI 被视为无向网络，
互为反向的重复行会自动合并，避免重复计数。`build_ppi_network()` 会自动识别
STRING API 的 0–1 分值或下载文件常见的 0–1000 分值；其
`score_threshold` 始终按 0–1000 输入。

### 富集结果

内部格式为每行一个 `term_id`—`gene` 组合，并包含 `term_name`, `group`,
`p_value`, `p_adjust`。既支持长表，也支持 MetaScape/clusterProfiler 常见的
分隔基因列表：

```r
enrichment <- prepare_enrichment(
  metascape_table,
  term_id = "GO",
  term_name = "Description",
  gene = NULL,
  gene_list = "Genes",
  group = "Group",
  p_value = "PValue",
  p_adjust = "FDR"
)
```

也可以通过可选的 Bioconductor 依赖直接运行 GO 富集：

```r
enrichment <- run_enrichment(
  genes = my_genes,
  organism = "human",
  ontology = "BP",
  background = expressed_genes
)
```

## 分数解释

每个代谢物—Term 组合均返回：

- `edge_count`：原始 a2 路径数；
- `weighted_sum`：MPI 权重 × PPI 权重 × hub 系数的总和；
- `normalized_score`：`weighted_sum / (MPI基因数 × Term基因数)`；
- `relation_coverage`：有网络支持的 MPI 基因比例；
- `term_coverage`：被连接到的 Term 基因比例。

默认 `hub_mode = "none"`。`"favor"` 会提高高连接度 Term 基因的相对权重，
`"penalize"` 会抑制通用 hub；无论选择哪种模式，建议同时检查覆盖率和逐边明细。

## 推荐分析原则

1. 使用与测序实验一致的物种、基因 ID 和 PPI 版本。
2. 富集时提供“实验中可检测的基因”作为背景集。
3. 明确记录 MPI/PPI 数据库版本和筛选阈值。
4. 不要仅按原始连接数排序；Term 大小和通用 hub 都会造成偏差。
5. 将高分结果作为实验验证或文献核对的候选，而不是因果结论。
