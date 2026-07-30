from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import chi2_contingency


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "output" / "eda_outputs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_CSV = ROOT / "output" / "Master_Phytochemical_Database_Clean.csv"

sns.set_theme(style="whitegrid")


def load_data() -> pd.DataFrame:
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Cleaned master dataset not found at {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)
    df = df.copy()
    df.columns = [c.strip().lower() for c in df.columns]
    for col in ["family", "plant_species", "plant_part", "phytochemical", "source_database"]:
        if col in df.columns:
            df[col] = df[col].fillna("Unknown").astype(str).str.strip()
    if "phytochemical" in df.columns:
        df["phytochemical"] = df["phytochemical"].str.lower()
    if "plant_species" in df.columns:
        df["plant_species"] = df["plant_species"].str.lower()
    return df


def split_list_field(series: pd.Series) -> list[str]:
    values: list[str] = []
    for value in series.fillna(""):
        if not isinstance(value, str):
            continue
        parts = re.split(r"[,;|/]+", value)
        for part in parts:
            cleaned = part.strip()
            if cleaned:
                values.append(cleaned)
    return values


def make_descriptive_stats(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    compounds_per_plant = (
        df.groupby("plant_species")["phytochemical"]
        .nunique()
        .reset_index(name="n_compounds")
        .sort_values("n_compounds", ascending=False)
    )

    if "phytochemical_class" in df.columns:
        class_counts = pd.Series(split_list_field(df["phytochemical_class"])).value_counts().head(20)
    else:
        class_counts = pd.Series([], dtype="int64")
    class_counts_df = class_counts.reset_index()
    class_counts_df.columns = ["therapeutic_class", "count"]

    plant_source_summary = (
        df.groupby("plant_species")["source_database"]
        .apply(lambda s: sorted(set(s)) if len(set(s)) > 0 else [])
        .reset_index(name="sources")
    )
    overlap_count = int((plant_source_summary["sources"].apply(len) > 1).sum())

    table1 = pd.DataFrame(
        [{
            "n_records": len(df),
            "n_plants": df["plant_species"].nunique(),
            "n_compounds": df["phytochemical"].nunique(),
            "n_families": df["family"].nunique(),
            "n_sources": df["source_database"].nunique(),
            "plants_with_multiple_sources": overlap_count,
        }]
    )
    return compounds_per_plant, class_counts_df, table1


def save_compound_frequency_plot(compounds_per_plant: pd.DataFrame) -> None:
    top = compounds_per_plant.head(20).copy()
    plt.figure(figsize=(12, 7), dpi=300)
    sns.barplot(data=top, x="n_compounds", y="plant_species", palette="viridis")
    plt.title("Compounds per Plant (Top 20)", fontsize=14)
    plt.xlabel("Number of Unique Compounds")
    plt.ylabel("Plant Species")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "figure_1_compounds_per_plant.png", dpi=300, bbox_inches="tight")
    plt.close()


def save_class_distribution_plot(class_counts_df: pd.DataFrame) -> None:
    plt.figure(figsize=(12, 7), dpi=300)
    sns.barplot(data=class_counts_df.head(20), x="count", y="therapeutic_class", palette="magma")
    plt.title("Distribution of Therapeutic Classes", fontsize=14)
    plt.xlabel("Count")
    plt.ylabel("Therapeutic Class")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "figure_2_therapeutic_class_distribution.png", dpi=300, bbox_inches="tight")
    plt.close()


def build_network(df: pd.DataFrame) -> nx.Graph:
    g = nx.Graph()
    plants = sorted(df["plant_species"].unique())
    compounds = sorted(df["phytochemical"].unique())
    g.add_nodes_from(plants, bipartite="plant")
    g.add_nodes_from(compounds, bipartite="compound")

    edge_df = df[["plant_species", "phytochemical"]].drop_duplicates()
    for _, row in edge_df.iterrows():
        g.add_edge(row["plant_species"], row["phytochemical"])
    return g


def save_network_plot(g: nx.Graph) -> None:
    top_plants = [p for p, deg in sorted(g.degree(), key=lambda x: x[1], reverse=True) if p in {n for n, d in g.nodes(data=True) if d.get("bipartite") == "plant"}][:20]
    if not top_plants:
        return
    relevant_nodes = set(top_plants)
    for plant in top_plants:
        relevant_nodes.update([nbr for nbr in g.neighbors(plant)])
    subg = g.subgraph(relevant_nodes).copy()
    pos = nx.bipartite_layout(subg, top_plants)

    plt.figure(figsize=(12, 10), dpi=300)
    nx.draw(subg, pos, with_labels=True, node_color=["#1f77b4" if subg.nodes[n].get("bipartite") == "plant" else "#ff7f0e" for n in subg.nodes], node_size=700, edge_color="gray", alpha=0.8)
    plt.title("Top 20 Plant–Phytochemical Network Clusters", fontsize=14)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "figure_3_network_top20.png", dpi=300, bbox_inches="tight")
    plt.close()

    nx.write_gexf(g, OUTPUT_DIR / "plant_phytochemical_network.gexf")


def compute_jaccard(df: pd.DataFrame) -> pd.DataFrame:
    kew_path_candidates = [ROOT / "kew_plant_list.csv", ROOT / "kew_plant_list.txt", OUTPUT_DIR / "kew_plant_list.csv"]
    kew_plant_list = None
    for path in kew_path_candidates:
        if path.exists():
            try:
                kew_plant_list = set(pd.read_csv(path, header=None).iloc[:, 0].astype(str).str.strip().str.lower())
            except Exception:
                kew_plant_list = set(path.read_text(encoding="utf-8").splitlines())
            break

    if kew_plant_list is None:
        # Fallback proxy: use the integrated plant list as a stand-in when no external Kew list is provided.
        kew_plant_list = set(df["plant_species"].unique())
        note = "No external Kew plant list was found; the integrated plant list was used as a proxy."
    else:
        note = "Kew plant list loaded from file."

    nmpdb_plants = set(df.loc[df["source_database"].str.lower() == "nmpdb", "plant_species"].unique())
    kew_set = set(kew_plant_list)
    intersection = kew_set & nmpdb_plants
    union = kew_set | nmpdb_plants
    jaccard = len(intersection) / len(union) if union else np.nan

    table2 = pd.DataFrame(
        [{
            "kew_plants": len(kew_set),
            "nmpdb_plants": len(nmpdb_plants),
            "intersection": len(intersection),
            "union": len(union),
            "jaccard_index": jaccard,
            "note": note,
        }]
    )
    return table2


def run_chi_square(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    df2 = df.copy()
    if "pharmacological_activities_compound" not in df2.columns:
        raise KeyError("pharmacological_activities_compound column is required for the chi-square analysis.")

    target_term = "Antiinflammatory"
    df2["target_therapeutic"] = df2["pharmacological_activities_compound"].fillna("").str.lower().str.contains(target_term.lower(), na=False)
    family_counts = df2["family"].value_counts().head(8)
    selected_families = family_counts.index.tolist()
    subset = df2[df2["family"].isin(selected_families)].copy()
    contingency = pd.crosstab(subset["family"], subset["target_therapeutic"])

    chi2, p_value, dof, expected = chi2_contingency(contingency)
    n = contingency.values.sum()
    cramers_v = np.sqrt(chi2 / (n * (min(contingency.shape) - 1))) if min(contingency.shape) > 1 else np.nan

    result = pd.DataFrame(
        [{
            "target_therapeutic_term": target_term,
            "chi_square_statistic": chi2,
            "degrees_of_freedom": dof,
            "p_value": p_value,
            "cramers_v": cramers_v,
            "n_rows": n,
        }]
    )
    return contingency, {"result": result, "contingency": contingency}


def save_chi_square_plot(contingency: pd.DataFrame) -> None:
    plt.figure(figsize=(10, 7), dpi=300)
    sns.heatmap(contingency, annot=True, fmt="d", cmap="coolwarm")
    plt.title("Contingency Table: Family vs Antiinflammatory Use")
    plt.ylabel("Family")
    plt.xlabel("Antiinflammatory Present")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "figure_4_chi_square_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()


def main() -> None:
    df = load_data()

    compounds_per_plant, class_counts_df, table1 = make_descriptive_stats(df)
    compounds_per_plant.to_csv(OUTPUT_DIR / "table_compounds_per_plant.csv", index=False)
    class_counts_df.to_csv(OUTPUT_DIR / "table_therapeutic_class_distribution.csv", index=False)
    table1.to_csv(OUTPUT_DIR / "table_1_overall_summary.csv", index=False)

    save_compound_frequency_plot(compounds_per_plant)
    save_class_distribution_plot(class_counts_df)

    g = build_network(df)
    save_network_plot(g)

    table2 = compute_jaccard(df)
    table2.to_csv(OUTPUT_DIR / "table_2_jaccard_similarity.csv", index=False)

    contingency, chi_result = run_chi_square(df)
    chi_result["result"].to_csv(OUTPUT_DIR / "table_chi_square_results.csv", index=False)
    contingency.to_csv(OUTPUT_DIR / "table_chi_square_contingency.csv")
    save_chi_square_plot(contingency)

    print("EDA completed successfully.")
    print(f"Figures and tables saved to: {OUTPUT_DIR}")
    print("\nTable 1 overall summary:")
    print(table1.to_string(index=False))
    print("\nTable 2 Jaccard summary:")
    print(table2.to_string(index=False))
    print("\nChi-square result:")
    print(chi_result["result"].to_string(index=False))


if __name__ == "__main__":
    main()
