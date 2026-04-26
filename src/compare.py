import argparse
import pandas as pd
import numpy as np
import os
import sys
import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from itertools import combinations
from upsetplot import UpSet, from_memberships


#visual metrics
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         11,
    "axes.titlesize":    13,
    "axes.titleweight":  "bold",
    "axes.labelsize":    11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    0.8,
    "xtick.labelsize":   10,
    "ytick.labelsize":   10,
    "legend.fontsize":   10,
    "legend.frameon":    False,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "savefig.facecolor": "white",
})
TOOL_PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7",
    "#56B4E9", "#D55E00", "#F0E442", "#000000",
]


#Custom functions

def jaccard(a, b):
    """Jaccard coefficient calculation"""
    if not a and not b:
        return np.nan
    u = len(a | b)
    return len(a & b) / u if u > 0 else 0.0


def get_tool_colors(tools):
    '''Function returns tool colours'''
    return {t: TOOL_PALETTE[i % len(TOOL_PALETTE)] for i, t in enumerate(sorted(tools))}


def wrap(text, width=30):
    return "\n".join(textwrap.wrap(str(text), width))


def gene_sets_for_sample(data, sample, tools):
    return {
        t: set(data[(data["sample"] == sample) & (data["tool"] == t)]["gene"].dropna())
        for t in tools
    }


def classify_confidence(n_tools, total):
    if n_tools == total:
        return "high_confidence"
    elif n_tools >= max(2, total // 2):
        return "moderate_confidence"
    else:
        return "low_confidence_single_tool"


def explain_discordance(gene, data, sample, tools):
    rows = data[(data["sample"] == sample) & (data["gene"] == gene)]
    detecting_tools = sorted(rows["tool"].unique())
    missing_tools = [t for t in tools if t not in detecting_tools]

    parts = []
    for _, row in rows.iterrows():
        ident = pd.to_numeric(row["identity"], errors="coerce")
        cov   = pd.to_numeric(row["coverage"], errors="coerce")
        ident_str = f"{ident:.1f}%" if pd.notna(ident) else "N/A"
        cov_str   = f"{cov:.1f}%"   if pd.notna(cov)   else "N/A"
        parts.append(
            f"{row['tool']} detected at {ident_str} identity, {cov_str} coverage"
        )

    explanation = "; ".join(parts)
    if missing_tools:
        explanation += f". Not found by: {', '.join(missing_tools)}."

    max_id = pd.to_numeric(rows["identity"], errors="coerce").max()
    if pd.notna(max_id) and max_id < 90:
        explanation += " [LOW IDENTITY — possible false positive or novel variant]"
    elif len(detecting_tools) == 1:
        explanation += " [SINGLE TOOL ONLY — treat with caution]"

    return explanation


def write_sample_tables(data, sample, tools, out_dir):
    sd = data[data["sample"] == sample].copy()
    tc = get_tool_colors(tools)

    sd.to_csv(f"{out_dir}/{sample}_all_hits.tsv", sep="\t", index=False)

    records = []
    for gene, gdf in sd.groupby("gene"):
        detecting = sorted(gdf["tool"].unique())
        n = len(detecting)
        missing = [t for t in tools if t not in detecting]

        # per-tool identity and coverage
        tool_identity = {}
        tool_coverage = {}
        for t in tools:
            tdf = gdf[gdf["tool"] == t]
            tool_identity[t] = round(pd.to_numeric(tdf["identity"], errors="coerce").mean(), 2) if not tdf.empty else None
            tool_coverage[t] = round(pd.to_numeric(tdf["coverage"], errors="coerce").mean(), 2) if not tdf.empty else None

        records.append({
            "gene":               gene,
            "drug_class":         gdf["drug_class"].dropna().iloc[0] if gdf["drug_class"].notna().any() else "unknown",
            "n_tools_detecting":  n,
            "tools_detecting":    "|".join(detecting),
            "tools_missing":      "|".join(missing) if missing else "none",
            "confidence":         classify_confidence(n, len(tools)),
            **{f"identity_{t}":  tool_identity[t] for t in tools},
            **{f"coverage_{t}":  tool_coverage[t] for t in tools},
            "note":               explain_discordance(gene, data, sample, tools),
        })

    summary = pd.DataFrame(records).sort_values(
        ["n_tools_detecting", "gene"], ascending=[False, True]
    )
    summary.to_csv(f"{out_dir}/{sample}_gene_detection_summary.tsv", sep="\t", index=False)

    # ── 3. High-confidence genes (all tools agree)
    high = summary[summary["confidence"] == "high_confidence"]
    high.to_csv(f"{out_dir}/{sample}_high_confidence_genes.tsv", sep="\t", index=False)

    # ── 4. Discordant genes — single tool only, with full explanation
    discordant = summary[summary["confidence"] == "low_confidence_single_tool"].copy()
    discordant.to_csv(f"{out_dir}/{sample}_discordant_single_tool.tsv", sep="\t", index=False)

    # ── 5. Pairwise overlap — lists actual gene names in each region
    gene_sets = gene_sets_for_sample(data, sample, tools)
    pair_rows = []
    for t1, t2 in combinations(tools, 2):
        s1, s2 = gene_sets.get(t1, set()), gene_sets.get(t2, set())
        shared   = sorted(s1 & s2)
        only_t1  = sorted(s1 - s2)
        only_t2  = sorted(s2 - s1)
        pair_rows.append({
            "tool_1":                t1,
            "tool_2":                t2,
            "jaccard":               round(jaccard(s1, s2), 4),
            "n_shared":              len(shared),
            "shared_genes":          "|".join(shared)  if shared   else "none",
            "n_only_tool_1":         len(only_t1),
            "genes_only_tool_1":     "|".join(only_t1) if only_t1  else "none",
            "n_only_tool_2":         len(only_t2),
            "genes_only_tool_2":     "|".join(only_t2) if only_t2  else "none",
        })
    pd.DataFrame(pair_rows).to_csv(
        f"{out_dir}/{sample}_pairwise_overlap.tsv", sep="\t", index=False
    )

    # ── 6. Human-readable discordance report
    write_discordance_report(summary, sample, tools, out_dir)

    return summary, gene_sets


def write_discordance_report(summary, sample, tools, out_dir):
    """Plain-text report listing exactly which genes differ and why."""
    lines = [
        f"AMR Tool Discordance Report — {sample}",
        "=" * 60,
        f"Total unique genes detected: {len(summary)}",
        f"Tools compared: {', '.join(tools)}",
        "",
    ]

    # high confidence
    high = summary[summary["confidence"] == "high_confidence"]
    lines += [
        f"HIGH CONFIDENCE — detected by all {len(tools)} tools ({len(high)} genes)",
        "-" * 40,
    ]
    for _, row in high.iterrows():
        lines.append(f"  {row['gene']}  [{row['drug_class']}]")
    lines.append("")

    # moderate
    mod = summary[summary["confidence"] == "moderate_confidence"]
    if not mod.empty:
        lines += [
            f"MODERATE CONFIDENCE — detected by majority of tools ({len(mod)} genes)",
            "-" * 40,
        ]
        for _, row in mod.iterrows():
            lines.append(
                f"  {row['gene']}  [{row['drug_class']}]"
                f"  detected by: {row['tools_detecting']}"
                f"  missing from: {row['tools_missing']}"
            )
        lines.append("")

    # discordant
    disc = summary[summary["confidence"] == "low_confidence_single_tool"]
    if not disc.empty:
        lines += [
            f"DISCORDANT — single tool only ({len(disc)} genes) — REQUIRES ATTENTION",
            "-" * 40,
        ]
        for _, row in disc.iterrows():
            lines.append(f"  {row['gene']}  [{row['drug_class']}]")
            lines.append(f"    {row['note']}")
        lines.append("")

    with open(f"{out_dir}/{sample}_discordance_report.txt", "w") as f:
        f.write("\n".join(lines))


def fig_upset_with_genes(gene_sets, summary, sample, tools, out_dir):
    """
    UpSet plot showing intersections, with gene names listed
    in a companion table beside each intersection region.
    """
    active = {t: gene_sets[t] for t in tools if gene_sets.get(t)}
    if len(active) < 2:
        return

    all_genes = sorted(set.union(*active.values()))
    memberships = [
        tuple(t for t in active if gene in active[t])
        for gene in all_genes
    ]

    try:
        upset_data = from_memberships(memberships, data=all_genes)

        fig = plt.figure(figsize=(12, 6))
        upset = UpSet(
            upset_data,
            subset_size="count",
            show_counts=True,
            sort_by="cardinality",
            totals_plot_elements=3,
        )
        upset.plot(fig)
        fig.suptitle(
            f"Gene Set Intersections Across Tools\n{sample}",
            y=1.02, fontsize=13, fontweight="bold"
        )
        fig.savefig(f"{out_dir}/{sample}_upset.png")
        plt.close(fig)
    except Exception as e:
        print(f"  WARNING: UpSet failed for {sample}: {e}", file=sys.stderr)

    # companion table: for each intersection combination, list gene names
    combo_records = []
    for combo in sorted(set(memberships), key=lambda x: (-len(x), x)):
        genes_in_combo = [g for g, m in zip(all_genes, memberships) if m == combo]
        combo_records.append({
            "intersection":  " & ".join(combo),
            "n_tools":       len(combo),
            "n_genes":       len(genes_in_combo),
            "genes":         "|".join(sorted(genes_in_combo)),
        })
    pd.DataFrame(combo_records).sort_values(
        ["n_tools", "n_genes"], ascending=[False, False]
    ).to_csv(f"{out_dir}/{sample}_intersection_gene_lists.tsv", sep="\t", index=False)


def fig_discordance_dotplot(summary, sample, tools, out_dir):
    """
    Dot plot: one row per gene, one column per tool.
    Dot present = detected. Color = identity. Size = coverage.
    Genes sorted by number of tools detecting them.
    Immediately shows which genes are missed by which tools.
    """
    tc = get_tool_colors(tools)
    genes = summary["gene"].tolist()
    if not genes:
        return

    # build presence matrix with identity values
    matrix = pd.DataFrame(index=genes, columns=tools, dtype=float)
    for gene in genes:
        for tool in tools:
            col = f"identity_{tool}"
            val = summary.loc[summary["gene"] == gene, col].values
            matrix.loc[gene, tool] = val[0] if len(val) > 0 and val[0] is not None else np.nan

    n_genes = len(genes)
    fig_h = max(6, n_genes * 0.35 + 2)
    fig, ax = plt.subplots(figsize=(max(6, len(tools) * 1.8), fig_h))

    for i, tool in enumerate(tools):
        for j, gene in enumerate(genes):
            val = matrix.loc[gene, tool]
            if pd.notna(val):
                color = plt.cm.RdYlGn((val - 60) / 40)  # 60-100% range
                ax.scatter(i, j, s=180, color=color,
                           edgecolors="#333333", linewidths=0.5, zorder=3)
            else:
                # not detected — show X
                ax.scatter(i, j, s=80, marker="x",
                           color="#cccccc", linewidths=1.2, zorder=2)

    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels(tools, rotation=30, ha="right")
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(genes, fontsize=8)
    ax.set_xlim(-0.5, len(tools) - 0.5)
    ax.set_ylim(-0.5, n_genes - 0.5)
    ax.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.5)

    # confidence band — shade high confidence rows
    high_genes = summary[summary["confidence"] == "high_confidence"]["gene"].tolist()
    for j, gene in enumerate(genes):
        if gene in high_genes:
            ax.axhspan(j - 0.4, j + 0.4, color="#e8f5e9", zorder=1)

    # colorbar for identity
    sm = plt.cm.ScalarMappable(
        cmap="RdYlGn",
        norm=plt.Normalize(vmin=60, vmax=100)
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.02)
    cbar.set_label("% Identity", fontsize=9)

    # legend
    detected_patch = mpatches.Patch(color="#2CA02C", label="Detected (dot = identity)")
    missed_patch   = mpatches.Patch(color="#cccccc", label="Not detected (×)")
    high_patch     = mpatches.Patch(color="#e8f5e9", label="High confidence (all tools)")
    ax.legend(
        handles=[detected_patch, missed_patch, high_patch],
        loc="lower right", fontsize=8
    )

    ax.set_title(
        f"Gene Detection Map — {sample}\n"
        f"Rows = genes | Columns = tools | Colour = % identity",
        fontsize=11
    )
    fig.savefig(f"{out_dir}/{sample}_gene_detection_map.png")
    plt.close(fig)


def fig_drug_class_discordance(data, summary, sample, tools, out_dir):
    """
    For each drug class: stacked bar showing how many genes are
    high/moderate/low confidence. Makes clear which drug classes
    have tool disagreement.
    """
    sd = data[data["sample"] == sample].copy()
    if sd.empty:
        return

    sd["drug_class_clean"] = (
        sd["drug_class"].astype(str)
        .str.split(" - ").str[0]
        .str.strip()
        .str.title()
    )

    merged = sd[["gene", "drug_class_clean"]].drop_duplicates().merge(
        summary[["gene", "confidence"]], on="gene", how="left"
    )

    conf_order = ["high_confidence", "moderate_confidence", "low_confidence_single_tool"]
    conf_labels = {
        "high_confidence":             "All tools agree",
        "moderate_confidence":         "Majority agree",
        "low_confidence_single_tool":  "Single tool only",
    }
    conf_colors = {
        "high_confidence":             "#2CA02C",
        "moderate_confidence":         "#FFBB00",
        "low_confidence_single_tool":  "#D62728",
    }

    grouped = (
        merged.groupby(["drug_class_clean", "confidence"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=conf_order, fill_value=0)
    )
    grouped = grouped.loc[grouped.sum(axis=1).nlargest(15).index]
    if grouped.empty:
        return

    fig, ax = plt.subplots(figsize=(9, max(4, len(grouped) * 0.5 + 1)))
    bottom = np.zeros(len(grouped))

    for conf in conf_order:
        if conf not in grouped.columns:
            continue
        vals = grouped[conf].values
        bars = ax.barh(
            grouped.index, vals, left=bottom,
            color=conf_colors[conf],
            label=conf_labels[conf],
            edgecolor="white", linewidth=0.4,
        )
        # annotate non-zero bars with count
        for bar, val in zip(bars, vals):
            if val > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + bar.get_height() / 2,
                    str(int(val)),
                    ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold"
                )
        bottom += vals

    ax.set_xlabel("Number of genes")
    ax.set_title(
        f"Drug Class Confidence — {sample}\n"
        f"Red = only one tool detected these genes"
    )
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.invert_yaxis()

    fig.savefig(f"{out_dir}/{sample}_drug_class_confidence.png")
    plt.close(fig)


def fig_identity_scatter(data, summary, sample, tools, out_dir):
    sd = data[data["sample"] == sample]
    multi = summary[summary["n_tools_detecting"] > 1]["gene"].tolist()
    if not multi:
        return

    tc = get_tool_colors(tools)
    fig, ax = plt.subplots(figsize=(max(7, len(multi) * 0.5 + 2), 5))

    for tool in tools:
        tool_data = sd[(sd["tool"] == tool) & (sd["gene"].isin(multi))].copy()
        if tool_data.empty:
            continue
        # align by gene position in multi list
        tool_data = tool_data.drop_duplicates(subset="gene")
        x = [multi.index(g) for g in tool_data["gene"] if g in multi]
        y = [
            tool_data.loc[tool_data["gene"] == multi[i], "identity"].values[0]
            for i in x
            if not tool_data.loc[tool_data["gene"] == multi[i], "identity"].empty
        ]
        # only plot where both x and y are valid
        xy = [(xi, yi) for xi, yi in zip(x, y) if pd.notna(yi)]
        if not xy:
            continue
        xs, ys = zip(*xy)
        ax.scatter(xs, ys, label=tool, color=tc[tool],
                   s=60, alpha=0.85, edgecolors="#333333",
                   linewidths=0.4, zorder=3)

    ax.set_xticks(range(len(multi)))
    ax.set_xticklabels(multi, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("% Identity to reference")
    ax.set_ylim(50, 105)
    ax.set_title(
        f"Identity Scores for Shared Genes — {sample}\n"
        f"Divergence between tools suggests database differences"
    )
    ax.legend(title="Tool")
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5)

    fig.savefig(f"{out_dir}/{sample}_identity_scatter_shared_genes.png")
    plt.close(fig)

def fig_global_jaccard(all_data, tools, out_dir):
    gene_sets = {
        t: set(all_data[all_data["tool"] == t]["gene"].dropna())
        for t in tools
    }
    matrix = pd.DataFrame(index=tools, columns=tools, dtype=float)
    for t1 in tools:
        for t2 in tools:
            matrix.loc[t1, t2] = jaccard(gene_sets[t1], gene_sets[t2])

    matrix.to_csv(f"{out_dir}/global_jaccard_matrix.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(6, 5))
    mask = np.eye(len(tools), dtype=bool)
    sns.heatmap(
        matrix.astype(float),
        annot=True, fmt=".2f",
        cmap="RdYlGn",
        vmin=0, vmax=1,
        mask=mask,
        linewidths=0.5,
        square=True,
        cbar_kws={"label": "Jaccard similarity"},
        ax=ax,
    )
    for i in range(len(tools)):
        ax.add_patch(plt.Rectangle((i, i), 1, 1, fill=True, color="#f5f5f5", lw=0))
        ax.text(i + 0.5, i + 0.5, "—", ha="center", va="center",
                fontsize=11, color="#aaaaaa")

    ax.set_title(
        "Global Jaccard Similarity Between Tools\n"
        "(all samples combined — higher = more agreement)"
    )
    ax.tick_params(axis="x", rotation=30)
    ax.tick_params(axis="y", rotation=0)
    fig.savefig(f"{out_dir}/global_jaccard_heatmap.png")
    plt.close(fig)


def fig_global_discordance(all_data, tools, samples, out_dir):
    rows = []
    for sample in samples:
        sd = all_data[all_data["sample"] == sample]
        n_per_gene = sd.groupby("gene")["tool"].nunique()
        rows.append({
            "sample":           sample,
            "All tools agree":  (n_per_gene == len(tools)).sum(),
            "Majority agree":   ((n_per_gene >= 2) & (n_per_gene < len(tools))).sum(),
            "Single tool only": (n_per_gene == 1).sum(),
        })

    df = pd.DataFrame(rows).set_index("sample")
    colors = ["#2CA02C", "#FFBB00", "#D62728"]

    fig, ax = plt.subplots(figsize=(max(8, len(samples)), 5))
    df.plot(kind="bar", stacked=True, color=colors,
            edgecolor="white", linewidth=0.4, ax=ax, width=0.65)

    ax.set_ylabel("Number of unique genes")
    ax.set_title(
        "Tool Consensus vs. Discordance per Sample\n"
        "Red = genes only one tool detected — potential false positives or novel variants"
    )
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(title="Confidence level", bbox_to_anchor=(1.01, 1), loc="upper left")

    fig.savefig(f"{out_dir}/global_consensus_discordance.png")
    plt.close(fig)
    df.reset_index().to_csv(f"{out_dir}/global_consensus_discordance.tsv", sep="\t", index=False)


def fig_gene_count_heatmap(all_data, tools, samples, out_dir):
    counts = (
        all_data.groupby(["sample", "tool"])
        .size()
        .reset_index(name="count")
    )
    pivot = (
        counts.pivot(index="sample", columns="tool", values="count")
        .fillna(0).astype(int)
        .reindex(columns=[t for t in tools if t in counts["tool"].unique()])
    )

    fig, ax = plt.subplots(figsize=(max(6, len(tools) * 1.8), max(5, len(samples) * 0.6 + 1)))
    sns.heatmap(
        pivot,
        annot=True, fmt="d",
        cmap="Blues",
        linewidths=0.4,
        linecolor="#dddddd",
        cbar_kws={"label": "Genes detected"},
        ax=ax,
    )
    ax.set_title("Total Genes Detected per Tool per Sample")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=30)
    ax.tick_params(axis="y", rotation=0)

    fig.savefig(f"{out_dir}/gene_count_heatmap.png")
    plt.close(fig)
    pivot.to_csv(f"{out_dir}/gene_count_per_tool_per_sample.tsv", sep="\t")


def write_global_tables(all_data, tools, samples, out_dir):
    # full data
    all_data.to_csv(f"{out_dir}/all_hits.tsv", sep="\t", index=False)

    # per-gene summary across all samples
    records = []
    for (sample, gene), gdf in all_data.groupby(["sample", "gene"]):
        detecting = sorted(gdf["tool"].unique())
        missing   = [t for t in tools if t not in detecting]
        records.append({
            "sample":            sample,
            "gene":              gene,
            "drug_class":        gdf["drug_class"].dropna().iloc[0] if gdf["drug_class"].notna().any() else "unknown",
            "n_tools_detecting": len(detecting),
            "tools_detecting":   "|".join(detecting),
            "tools_missing":     "|".join(missing) if missing else "none",
            "confidence":        classify_confidence(len(detecting), len(tools)),
            "mean_identity":     round(gdf["identity"].mean(), 2) if gdf["identity"].notna().any() else None,
            "mean_coverage":     round(gdf["coverage"].mean(), 2) if gdf["coverage"].notna().any() else None,
        })

    global_summary = pd.DataFrame(records).sort_values(
        ["sample", "n_tools_detecting"], ascending=[True, False]
    )
    global_summary.to_csv(f"{out_dir}/global_gene_summary.tsv", sep="\t", index=False)

    # genes that are discordant across ALL samples — systematic bias
    systematic = (
        global_summary[global_summary["confidence"] == "low_confidence_single_tool"]
        .groupby("gene")
        .agg(
            n_samples_discordant=("sample", "nunique"),
            always_detected_by=("tools_detecting", lambda x: "|".join(sorted(set("|".join(x).split("|"))))),
            drug_class=("drug_class", "first"),
        )
        .reset_index()
        .sort_values("n_samples_discordant", ascending=False)
    )
    systematic.to_csv(f"{out_dir}/systematically_discordant_genes.tsv", sep="\t", index=False)

    return global_summary



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",      nargs="+", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    frames = []
    for f in args.input:
        try:
            df = pd.read_csv(f, sep="\t")
            frames.append(df)
        except Exception as e:
            print(f"WARNING: could not read {f}: {e}", file=sys.stderr)

    if not frames:
        print("No data to compare.", file=sys.stderr)
        sys.exit(1)

    all_data = pd.concat(frames, ignore_index=True)
    for col in ["identity", "coverage", "start", "stop"]:
        all_data[col] = pd.to_numeric(all_data[col], errors="coerce")

    samples  = sorted(all_data["sample"].unique())
    tools    = sorted(all_data["tool"].unique())

    print(f"\nLoaded {len(all_data)} hits | {len(samples)} samples | {len(tools)} tools")
    print(f"Tools:   {', '.join(tools)}")
    print(f"Samples: {', '.join(samples)}\n")

    for sample in samples:
        sample_dir = os.path.join(args.output_dir, sample)
        fig_dir    = os.path.join(sample_dir, "figures")
        os.makedirs(fig_dir, exist_ok=True)
        print(f"  {sample}...")

        summary, gene_sets = write_sample_tables(
            all_data, sample, tools, sample_dir
        )

        fig_upset_with_genes(gene_sets, summary, sample, tools, fig_dir)
        fig_discordance_dotplot(summary, sample, tools, fig_dir)
        fig_drug_class_discordance(all_data, summary, sample, tools, fig_dir)
        fig_identity_scatter(all_data, summary, sample, tools, fig_dir)

    global_dir = os.path.join(args.output_dir, "global")
    global_fig = os.path.join(global_dir, "figures")
    os.makedirs(global_fig, exist_ok=True)
    print("\n  Global outputs...")

    write_global_tables(all_data, tools, samples, global_dir)
    fig_global_jaccard(all_data, tools, global_fig)
    fig_global_discordance(all_data, tools, samples, global_fig)
    fig_gene_count_heatmap(all_data, tools, samples, global_fig)

    print(f"""
Done.
  Per-sample folders : {args.output_dir}/{{sample}}/
  Global folder      : {global_dir}/
  Total hits         : {len(all_data)}
  Unique genes       : {all_data['gene'].nunique()}
  Samples            : {len(samples)}
  Tools              : {len(tools)}
""")


if __name__ == "__main__":
    main()