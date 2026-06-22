"""
Perioperative A-V Proteomics Explorer
Streamlit version — reads pre-computed CSVs from outputs/tables/
Deploy: push repo to GitHub → Streamlit Cloud → set Main file path to streamlit_app/app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path
from scipy.stats import kruskal

st.set_page_config(
    page_title="A-V Proteomics Explorer",
    page_icon=None,
    layout="wide",
)

# ── Data paths ─────────────────────────────────────────────────────────────────
_here = Path(__file__).parent
# Works both locally (run from repo root or streamlit_app/) and on Streamlit Cloud
_candidates = [
    _here.parent / "outputs" / "tables",
    Path("outputs") / "tables",
    _here / ".." / "outputs" / "tables",
]
DATA_DIR = next((p.resolve() for p in _candidates if p.exists()), _candidates[0].resolve())


@st.cache_data
def load_data():
    def rd(name):
        p = DATA_DIR / name
        return pd.read_csv(p) if p.exists() else pd.DataFrame()

    av        = rd("03_av_results_full.csv")
    loadings  = rd("02_pca_loadings.csv")
    pca_var   = rd("02_pca_variance.csv")
    pca_sc    = rd("02_pca_scores.csv")
    gsea_h    = rd("04_gsea_hallmark.csv")
    gsea_k    = rd("04_gsea_kegg.csv")
    ecs       = rd("08_ecs_summary_hits.csv")
    eos8      = rd("08_eos_summary_hits.csv")
    eos9      = rd("09_eos_summary_hits.csv")
    cov       = rd("05_covariate_results.csv")
    skin_av   = rd("10_skin_av.csv")
    skin_kw   = rd("10_skin_surgery_kw.csv")
    return av, loadings, pca_var, pca_sc, gsea_h, gsea_k, ecs, eos8, eos9, cov, skin_av, skin_kw


(av, loadings, pca_var, pca_sc,
 gsea_h, gsea_k, ecs_hits, eos8_hits, eos9_hits,
 cov_res, skin_av, skin_kw) = load_data()

# ── Pre-process A-V table ──────────────────────────────────────────────────────
if not av.empty:
    av["neg_log10p"]    = -np.log10(av["P.Value"].clip(lower=1e-300))
    av["EntrezGeneID"]  = av["EntrezGeneID"].astype(str)
    av["sig_label"]     = np.where(av["adj.P.Val"] < 0.05, "FDR < 0.05",
                          np.where(av["P.Value"]   < 0.05, "p < 0.05", "NS"))

# ── Combine ECS / EOS hits ─────────────────────────────────────────────────────
def _norm_ecs(df):
    if df.empty:
        return pd.DataFrame()
    return df.assign(system="ECS", analysis="A-V").rename(
        columns={"ecs_role": "role"}
    )[["system", "analysis", "platform", "gene_symbol",
       "category", "role", "logFC", "p_value", "direction"]]

def _norm_eos8(df):
    if df.empty:
        return pd.DataFrame()
    return df.assign(system="EOS").rename(
        columns={"eos_category": "category", "eos_role": "role"}
    )[["system", "analysis", "platform", "gene_symbol",
       "category", "role", "logFC", "p_value", "direction"]]

def _norm_eos9(df):
    if df.empty:
        return pd.DataFrame()
    d = df.copy()
    d["system"]      = "EOS"
    d["platform"]    = "SomaScan"
    d["analysis"]    = d["analysis"].astype(str) + ": " + d.get("covariate", "").astype(str)
    d = d.rename(columns={"EntrezGeneSymbol": "gene_symbol",
                           "eos_category": "category", "eos_role": "role"})
    return d[["system", "analysis", "platform", "gene_symbol",
              "category", "role", "logFC", "p_value", "direction"]]

ecs_eos = pd.concat([_norm_ecs(ecs_hits), _norm_eos8(eos8_hits), _norm_eos9(eos9_hits)],
                    ignore_index=True).drop_duplicates(
                        subset=["system", "analysis", "platform", "gene_symbol"])

# ── GSEA combined ──────────────────────────────────────────────────────────────
gsea_h["collection"] = "Hallmark"
gsea_k["collection"] = "KEGG"
gsea_all = pd.concat([gsea_h, gsea_k], ignore_index=True).sort_values("padj")

def clean_pathway(s):
    return s.replace("HALLMARK_", "").replace("KEGG_", "").replace("_", " ").lower()

def parse_le(s):
    if pd.isna(s) or s == "":
        return set()
    return set(str(s).split(";"))

# ── Colours ────────────────────────────────────────────────────────────────────
COL_ART  = "#c0392b"
COL_VEN  = "#2980b9"
COL_NS   = "#cccccc"
COL_NOM  = "#1f78b4"
COL_SRC  = "#e31a1c"
COL_LE   = "#ff7f00"
COL_SRCH = "#9400d3"
SURG_COL = {"Head/Neck": "#66c2a5", "Laparoscopic": "#fc8d62", "Spine": "#8da0cb"}

# ── Helpers ────────────────────────────────────────────────────────────────────
def direction_color(logfc):
    return COL_VEN if logfc < 0 else COL_ART

def av_volcano_fig(df, highlight_genes=None, le_genes=None, title="A-V Volcano"):
    """Shared volcano builder used by multiple tabs."""
    df = df.copy()
    df["_hi"]  = df["EntrezGeneSymbol"].str.upper().isin(highlight_genes or set())
    df["_le"]  = df["EntrezGeneID"].isin(le_genes or set())
    df["_col"] = np.where(df["_hi"], "Searched",
                 np.where(df["_le"], "Leading edge",
                 np.where(df["sig_label"] == "FDR < 0.05", "FDR < 0.05",
                 np.where(df["sig_label"] == "p < 0.05",   "p < 0.05", "NS"))))
    df["_sz"]  = np.where(df["_hi"] | df["_le"], 8,
                 np.where(df["sig_label"] != "NS", 5, 3))
    df["_op"]  = np.where(df["_col"] == "NS", 0.22, 0.80)
    df["_tip"] = (
        "<b>" + df["Target"].fillna("") + "</b> (" + df["EntrezGeneSymbol"].fillna("") + ")<br>"
        + "logFC: " + df["logFC"].round(3).astype(str) + "<br>"
        + "p = " + df["P.Value"].map(lambda x: f"{x:.3g}") + "<br>"
        + df["direction"].fillna("")
    )
    col_map = {"Searched": COL_SRC, "Leading edge": COL_LE,
               "FDR < 0.05": "#ff7f00", "p < 0.05": COL_NOM, "NS": COL_NS}
    order   = ["NS", "p < 0.05", "FDR < 0.05", "Leading edge", "Searched"]
    fig = go.Figure()
    for grp in order:
        sub = df[df["_col"] == grp]
        if sub.empty:
            continue
        fig.add_trace(go.Scatter(
            x=sub["logFC"], y=sub["neg_log10p"],
            mode="markers",
            name=grp,
            marker=dict(color=col_map[grp], size=sub["_sz"],
                        opacity=sub["_op"].mean()),
            text=sub["_tip"], hoverinfo="text",
            customdata=sub["EntrezGeneSymbol"],
        ))
    # Reference lines
    max_y = df["neg_log10p"].max()
    fig.add_hline(y=-np.log10(0.05), line_dash="dot", line_color=COL_NOM, line_width=1)
    fig.add_vline(x=0, line_dash="dot", line_color="gray", line_width=1)
    # Annotate hits
    hits = df[df["_hi"] | (df["_le"] & (le_genes is not None))]
    for _, r in hits.iterrows():
        fig.add_annotation(x=r["logFC"], y=r["neg_log10p"],
                           text=r["EntrezGeneSymbol"], showarrow=True,
                           arrowhead=2, arrowsize=0.6,
                           font=dict(size=11, color=COL_SRC if r["_hi"] else COL_LE),
                           bgcolor="white",
                           bordercolor=COL_SRC if r["_hi"] else COL_LE)
    fig.update_layout(
        title=title,
        xaxis_title="log2FC (Arterial / Venous)",
        yaxis_title="−log10(p-value)",
        hovermode="closest",
        legend=dict(title="Significance"),
        height=500,
        margin=dict(t=50, b=40),
    )
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Layout
# ══════════════════════════════════════════════════════════════════════════════

st.title("Perioperative A-V Proteomics Explorer")
st.caption("7,481 proteins · 12 patients · simultaneous arterial + venous blood draws intraoperatively")

# Global gene search in sidebar
with st.sidebar:
    st.markdown("### Gene search")
    gene_raw   = st.text_input("Gene symbol(s)", placeholder="e.g. SPINK9  KLK7  IL6",
                               help="Space-, comma-, or semicolon-separated gene symbols")
    if st.button("Clear", use_container_width=True):
        gene_raw = ""
    searched = {g.strip().upper() for g in gene_raw.replace(",", " ").replace(";", " ").split()
                if g.strip()}

    if searched and not av.empty:
        st.markdown("---")
        st.markdown("**A-V results**")
        hits = av[av["EntrezGeneSymbol"].str.upper().isin(searched)].sort_values("P.Value")
        for _, r in hits.iterrows():
            fc_str = f"{r['logFC']:+.3f}"
            color  = COL_VEN if r["logFC"] < 0 else COL_ART
            st.markdown(
                f"**{r['EntrezGeneSymbol']}** "
                f"<span style='color:{color}'>{fc_str}</span> "
                f"p={r['P.Value']:.3g}",
                unsafe_allow_html=True,
            )

tab_vol, tab_pca, tab_gsea, tab_ecs, tab_skin = st.tabs([
    "A-V Volcano", "PCA", "GSEA", "ECS / EOS", "Skin Proteases"
])

# ══════════════════════════════════════════════════════════════════════════════
# Tab 1 — A-V Volcano
# ══════════════════════════════════════════════════════════════════════════════
with tab_vol:
    c1, c2, c3 = st.columns([2, 2, 1])
    sig_choice = c1.selectbox("Show proteins",
                              ["All (7,481)", "Nominal p < 0.05", "FDR < 0.05"],
                              key="v_sig")
    dir_choice = c2.selectbox("Direction",
                              ["Both", "Arterial > Venous", "Venous > Arterial"],
                              key="v_dir")

    if av.empty:
        st.warning("Run `R_scripts/03_av_limma.R` to generate `03_av_results_full.csv`.")
    else:
        df_vol = av.copy()
        if sig_choice == "Nominal p < 0.05":
            df_vol = df_vol[df_vol["sig_nom"]]
        elif sig_choice == "FDR < 0.05":
            df_vol = df_vol[df_vol["sig_fdr"]]
        if dir_choice != "Both":
            df_vol = df_vol[df_vol["direction"] == dir_choice]

        st.plotly_chart(
            av_volcano_fig(df_vol, highlight_genes=searched,
                           title=f"A-V Volcano — {len(df_vol):,} proteins"),
            use_container_width=True,
        )

        st.markdown("**Nominally significant proteins (p < 0.05)** — or gene search results")
        tbl = av[av["EntrezGeneSymbol"].str.upper().isin(searched)] if searched else \
              av[av["sig_nom"]].sort_values("P.Value").head(300)
        if not tbl.empty:
            st.dataframe(
                tbl[["EntrezGeneSymbol", "Target", "logFC", "P.Value", "adj.P.Val",
                      "direction", "sig_nom", "sig_fdr"]]
                .rename(columns={"EntrezGeneSymbol": "Gene", "Target": "Protein",
                                 "P.Value": "p", "adj.P.Val": "adj.p",
                                 "sig_nom": "Nominal", "sig_fdr": "FDR"})
                .assign(logFC=lambda d: d["logFC"].round(3),
                        p=lambda d: d["p"].map(lambda x: f"{x:.3g}"),
                        **{"adj.p": lambda d: d["adj.p"].map(lambda x: f"{x:.3g}")}),
                use_container_width=True, hide_index=True,
            )


# ══════════════════════════════════════════════════════════════════════════════
# Tab 2 — PCA
# ══════════════════════════════════════════════════════════════════════════════
with tab_pca:
    # ── Axis selectors ────────────────────────────────────────────────────────
    _pc_opts = [c for c in ["PC1", "PC2", "PC3", "PC4"] if c in loadings.columns] or ["PC1", "PC2", "PC3", "PC4"]
    _pct_dict = ({f"PC{i+1}": v for i, v in enumerate(pca_var["var_pct"].tolist())}
                 if not pca_var.empty else {})

    pa1, pa2, _ = st.columns([1, 1, 2])
    x_pc = pa1.selectbox("X axis", _pc_opts, index=0, key="pca_x")
    y_pc = pa2.selectbox("Y axis", _pc_opts, index=min(1, len(_pc_opts) - 1), key="pca_y")

    col_l, col_r = st.columns(2)

    with col_l:
        st.markdown(f"#### Protein loadings — {x_pc} vs {y_pc}")
        if loadings.empty:
            st.warning("Run `R_scripts/02_pca.R` to generate PCA outputs.")
        else:
            df_ld = loadings.copy()
            df_ld["_hit"] = df_ld["EntrezGeneSymbol"].str.upper().isin(searched)
            df_ld["_col"] = np.where(df_ld["_hit"], "Searched", "Other")
            df_ld["_sz"]  = np.where(df_ld["_hit"], 9, 4)
            df_ld["_op"]  = np.where(df_ld["_hit"], 1.0, 0.35)
            df_ld["_tip"] = ("<b>" + df_ld["Target"].fillna("") + "</b> ("
                             + df_ld["EntrezGeneSymbol"].fillna("") + ")<br>"
                             + f"{x_pc}: " + df_ld[x_pc].round(4).astype(str)
                             + f"  {y_pc}: " + df_ld[y_pc].round(4).astype(str))
            fig_ld = go.Figure()
            for grp, col in [("Other", COL_NS), ("Searched", COL_SRC)]:
                sub = df_ld[df_ld["_col"] == grp]
                if sub.empty:
                    continue
                fig_ld.add_trace(go.Scatter(
                    x=sub[x_pc], y=sub[y_pc], mode="markers", name=grp,
                    text=sub["_tip"], hoverinfo="text",
                    marker=dict(color=col, size=sub["_sz"], opacity=sub["_op"].mean())
                ))
            for _, r in df_ld[df_ld["_hit"]].iterrows():
                fig_ld.add_annotation(x=r[x_pc], y=r[y_pc], text=r["EntrezGeneSymbol"],
                                      showarrow=True, arrowhead=2,
                                      font=dict(size=10, color=COL_SRC),
                                      bgcolor="white", bordercolor=COL_SRC)
            fig_ld.update_layout(
                xaxis_title=f"{x_pc} ({_pct_dict.get(x_pc, 0):.1f}%)",
                yaxis_title=f"{y_pc} ({_pct_dict.get(y_pc, 0):.1f}%)",
                height=420, margin=dict(t=30))
            st.plotly_chart(fig_ld, use_container_width=True)

        st.markdown(f"#### Top {x_pc} contributors")
        if not loadings.empty:
            top_ld = (loadings.assign(**{f"_abs": loadings[x_pc].abs()})
                      .sort_values("_abs", ascending=False)
                      .head(100)
                      [["EntrezGeneSymbol", "Target", "PC1", "PC2", "PC3", "PC4"]]
                      .rename(columns={"EntrezGeneSymbol": "Gene", "Target": "Protein"}))
            st.dataframe(top_ld.style.format({"PC1": "{:.4f}", "PC2": "{:.4f}",
                                               "PC3": "{:.4f}", "PC4": "{:.4f}"}),
                         use_container_width=True, hide_index=True)

    with col_r:
        st.markdown(f"#### Sample scores — {x_pc} vs {y_pc}")
        if pca_sc.empty:
            st.info("Re-run `R_scripts/02_pca.R` to enable this plot.")
        else:
            art    = pca_sc[pca_sc["draw"] == "Arterial"].set_index("patient")
            ven    = pca_sc[pca_sc["draw"] == "Venous"].set_index("patient")
            paired = art.join(ven, lsuffix="_a", rsuffix="_v").dropna()
            fig_sc = go.Figure()
            for sg, col in SURG_COL.items():
                for draw, sym in [("Arterial", "circle"), ("Venous", "triangle-up")]:
                    sub = pca_sc[(pca_sc["surgery_group"] == sg) & (pca_sc["draw"] == draw)]
                    fig_sc.add_trace(go.Scatter(
                        x=sub[x_pc], y=sub[y_pc], mode="markers+text",
                        name=f"{sg} – {draw}",
                        text=sub["patient"], textposition="top center",
                        marker=dict(color=col, symbol=sym, size=11,
                                    line=dict(color="white", width=1)),
                    ))
            for pt in paired.index:
                r = paired.loc[pt]
                fig_sc.add_shape(type="line",
                    x0=r[f"{x_pc}_a"], y0=r[f"{y_pc}_a"],
                    x1=r[f"{x_pc}_v"], y1=r[f"{y_pc}_v"],
                    line=dict(color="gray", width=1, dash="dot"))
            fig_sc.update_layout(
                xaxis_title=f"{x_pc} ({_pct_dict.get(x_pc, 0):.1f}%)",
                yaxis_title=f"{y_pc} ({_pct_dict.get(y_pc, 0):.1f}%)",
                height=420, margin=dict(t=30),
                legend=dict(font=dict(size=10)))
            st.plotly_chart(fig_sc, use_container_width=True)
            st.caption("Lines connect arterial–venous pairs per patient.")

    # ── Covariate ANOVA section ───────────────────────────────────────────────
    if not pca_sc.empty:
        st.markdown("---")
        st.markdown("#### PC score differences by covariate (Kruskal-Wallis)")

        _cov_map = {
            "Draw (A vs V)":         "draw",
            "Surgery type":          "surgery_group",
            "Sex (male)":            "male",
            "Obese (BMI > 30)":      "obese",
            "Hypothermia (< 36°C)":  "hypothermia",
            "Hyperglycemia (> 150)": "hyperglycemia",
        }
        _avail_covs = {k: v for k, v in _cov_map.items() if v in pca_sc.columns}

        if len(_avail_covs) < len(_cov_map):
            st.caption("Re-run `R_scripts/02_pca.R` and commit `02_pca_scores.csv` "
                       "to unlock sex/obesity/hypothermia/hyperglycemia covariates.")

        test_label = st.selectbox("Test covariate", list(_avail_covs.keys()), key="pca_cov")
        test_col   = _avail_covs[test_label]

        df_test = pca_sc.copy()
        if df_test[test_col].dtype == object or df_test[test_col].dtype.name == "bool":
            df_test[test_col] = df_test[test_col].map(
                lambda v: ("Yes" if v else "No") if isinstance(v, bool) else v
            )

        grp_obj    = df_test.dropna(subset=[test_col]).groupby(test_col)
        grp_names  = list(grp_obj.groups.keys())
        _pal       = [COL_ART, COL_VEN, "#27ae60", "#f39c12", "#8e44ad"]

        ac1, ac2 = st.columns([1, 2])

        with ac1:
            st.markdown("**p-values across PCs**")
            kw_rows = []
            for pc in _pc_opts:
                vals = [g[pc].values for _, g in grp_obj]
                if len(vals) >= 2 and all(len(v) > 0 for v in vals):
                    stat, p = kruskal(*vals)
                    kw_rows.append({"PC": pc, "H": f"{stat:.2f}",
                                    "p": f"{p:.3g}", "*": "✓" if p < 0.05 else ""})
                else:
                    kw_rows.append({"PC": pc, "H": "—", "p": "—", "*": ""})
            st.dataframe(pd.DataFrame(kw_rows), use_container_width=True, hide_index=True)
            st.caption(f"Groups: {', '.join(str(g) for g in grp_names)}  |  "
                       f"N = {len(df_test.dropna(subset=[test_col]))} samples")

        with ac2:
            st.markdown(f"**{x_pc} scores by {test_label}**")
            rng     = np.random.default_rng(42)
            fig_str = go.Figure()
            for i, (grp_name, grp_df) in enumerate(grp_obj):
                col    = _pal[i % len(_pal)]
                jitter = rng.uniform(-0.12, 0.12, len(grp_df))
                fig_str.add_trace(go.Scatter(
                    x=np.full(len(grp_df), i) + jitter,
                    y=grp_df[x_pc].values,
                    mode="markers",
                    name=str(grp_name),
                    text=grp_df["patient"].astype(str) + " (" + grp_df["draw"] + ")",
                    hoverinfo="text+y",
                    marker=dict(color=col, size=9, opacity=0.8),
                ))
                mean_val = grp_df[x_pc].mean()
                fig_str.add_shape(type="line",
                    x0=i - 0.3, x1=i + 0.3, y0=mean_val, y1=mean_val,
                    line=dict(color=col, width=2.5))
            fig_str.update_layout(
                xaxis=dict(tickmode="array", tickvals=list(range(len(grp_names))),
                           ticktext=[str(g) for g in grp_names], title=test_label),
                yaxis_title=f"{x_pc} score",
                height=340, margin=dict(t=10, b=40),
                showlegend=False,
            )
            st.plotly_chart(fig_str, use_container_width=True)


# ══════════════════════════════════════════════════════════════════════════════
# Tab 3 — GSEA
# ══════════════════════════════════════════════════════════════════════════════
with tab_gsea:
    if gsea_all.empty:
        st.warning("Run `R_scripts/04_gsea.R` to generate GSEA outputs.")
    else:
        gc1, gc2, gc3 = st.columns([2, 2, 2])
        g_coll = gc1.selectbox("Collection", ["All", "Hallmark", "KEGG"])
        g_dir  = gc2.selectbox("Direction",
                               ["Both", "Venous enriched (NES < 0)", "Arterial enriched (NES > 0)"])
        g_padj = gc3.slider("Max adj.p", 0.001, 1.0, 0.25, 0.005)

        df_g = gsea_all.copy()
        if g_coll != "All":
            df_g = df_g[df_g["collection"] == g_coll]
        if "Venous" in g_dir:
            df_g = df_g[df_g["NES"] < 0]
        elif "Arterial" in g_dir:
            df_g = df_g[df_g["NES"] > 0]
        df_g = df_g[df_g["padj"] <= g_padj].sort_values("padj")
        df_g["pathway_clean"] = df_g["pathway"].map(clean_pathway)

        st.markdown(f"**{len(df_g)} pathways** match filters — click a row to inspect leading edge")
        sel = st.dataframe(
            df_g[["collection", "pathway_clean", "NES", "pval", "padj", "size"]]
            .rename(columns={"collection": "Collection", "pathway_clean": "Pathway",
                             "pval": "p-val", "padj": "adj.p", "size": "Size"})
            .assign(NES=lambda d: d["NES"].round(3),
                    **{"p-val": lambda d: d["p-val"].map(lambda x: f"{x:.3g}"),
                       "adj.p": lambda d: d["adj.p"].map(lambda x: f"{x:.3g}")}),
            use_container_width=True, hide_index=True,
            on_select="rerun", selection_mode="single-row",
        )

        le_genes = set()
        sel_title = "Select a pathway above"
        if sel and sel.selection.rows:
            row = df_g.iloc[sel.selection.rows[0]]
            sel_title  = clean_pathway(row["pathway"])
            le_ids     = parse_le(row.get("leadingEdge", ""))
            le_genes   = le_ids
            le_syms    = av[av["EntrezGeneID"].isin(le_ids)]["EntrezGeneSymbol"].tolist()
            st.info(f"**Leading edge ({len(le_syms)} genes):** {', '.join(le_syms[:30])}"
                    + (" …" if len(le_syms) > 30 else ""))

        st.plotly_chart(
            av_volcano_fig(av, highlight_genes=searched, le_genes=le_genes, title=sel_title),
            use_container_width=True,
        )


# ══════════════════════════════════════════════════════════════════════════════
# Tab 4 — ECS / EOS
# ══════════════════════════════════════════════════════════════════════════════
with tab_ecs:
    if ecs_eos.empty:
        st.warning("Run `R_scripts/08_endocannabinoid_analysis.R` and "
                   "`R_scripts/09_opioid_analysis.R` to generate ECS/EOS outputs.")
    else:
        ec1, ec2 = st.columns(2)
        e_sys  = ec1.selectbox("System", ["Both", "ECS", "EOS"])
        e_plat = ec2.selectbox("Platform", ["All", "Luminex", "SomaScan"])

        df_ecs = ecs_eos.copy()
        if e_sys != "Both":
            df_ecs = df_ecs[df_ecs["system"] == e_sys]
        if e_plat != "All":
            df_ecs = df_ecs[df_ecs["platform"] == e_plat]
        df_ecs = df_ecs.sort_values("p_value")

        st.markdown(f"**{len(df_ecs)} hits** — click a row to highlight on volcano")
        sel_e = st.dataframe(
            df_ecs[["system", "analysis", "platform", "gene_symbol",
                    "category", "logFC", "p_value", "direction"]]
            .rename(columns={"system": "System", "analysis": "Analysis",
                             "platform": "Platform", "gene_symbol": "Gene",
                             "category": "Category", "p_value": "p",
                             "direction": "Direction"})
            .assign(logFC=lambda d: d["logFC"].round(3),
                    p=lambda d: d["p"].map(lambda x: f"{x:.3g}")),
            use_container_width=True, hide_index=True,
            on_select="rerun", selection_mode="single-row",
        )

        sel_gene = set()
        if sel_e and sel_e.selection.rows:
            g = df_ecs.iloc[sel_e.selection.rows[0]]["gene_symbol"]
            sel_gene = {str(g).upper()}

        all_ecs_genes = set(ecs_eos["gene_symbol"].str.upper())
        df_ev = av.copy()
        df_ev["_grp"] = np.where(
            df_ev["EntrezGeneSymbol"].str.upper().isin(searched | sel_gene),
                "Searched/Selected",
            np.where(df_ev["EntrezGeneSymbol"].str.upper().isin(all_ecs_genes),
                "ECS/EOS", df_ev["sig_label"]))

        col_map_e = {"Searched/Selected": COL_SRC, "ECS/EOS": COL_LE,
                     "FDR < 0.05": "#ff7f00", "p < 0.05": COL_NOM, "NS": COL_NS}
        fig_ev = go.Figure()
        for grp in ["NS", "p < 0.05", "FDR < 0.05", "ECS/EOS", "Searched/Selected"]:
            sub = df_ev[df_ev["_grp"] == grp]
            if sub.empty:
                continue
            tip = ("<b>" + sub["Target"].fillna("") + "</b> ("
                   + sub["EntrezGeneSymbol"].fillna("") + ")<br>"
                   + "logFC: " + sub["logFC"].round(3).astype(str)
                   + "  p=" + sub["P.Value"].map(lambda x: f"{x:.3g}"))
            fig_ev.add_trace(go.Scatter(
                x=sub["logFC"], y=sub["neg_log10p"], mode="markers", name=grp,
                text=tip, hoverinfo="text",
                marker=dict(color=col_map_e[grp],
                            size=np.where(grp in ("ECS/EOS", "Searched/Selected"), 7, 3),
                            opacity=0.2 if grp == "NS" else 0.8)
            ))
        fig_ev.add_hline(y=-np.log10(0.05), line_dash="dot", line_color=COL_NOM, line_width=1)
        fig_ev.add_vline(x=0, line_dash="dot", line_color="gray", line_width=1)
        fig_ev.update_layout(
            title="ECS / EOS gene positions on A-V volcano",
            xaxis_title="log2FC (Arterial / Venous)",
            yaxis_title="−log10(p-value)",
            height=480, margin=dict(t=50),
        )
        st.plotly_chart(fig_ev, use_container_width=True)


# ══════════════════════════════════════════════════════════════════════════════
# Tab 5 — Skin Proteases
# ══════════════════════════════════════════════════════════════════════════════
with tab_skin:
    if skin_av.empty:
        st.warning("Run `R_scripts/10_skin_protease_analysis.R` to generate skin protease outputs.")
    else:
        sk1, sk2 = st.columns([2, 4])
        s_fam = sk1.selectbox("Gene family", ["All", "SPINK", "KLK", "SERPIN"])

        with sk2:
            st.info(
                "**SPINK9** (rank 1, p = 1.7×10⁻⁴) and **KLK7** (rank 3, p = 1.9×10⁻³) are the "
                "top two individual protein hits in the entire dataset. Both are skin-specific. "
                "All nominally significant family members show **venous > arterial**, consistent "
                "with surgical incision releasing epidermal proteins into the venous return."
            )

        df_sk = skin_av.copy()
        if s_fam != "All":
            df_sk = df_sk[df_sk["family"] == s_fam]
        df_sk = df_sk.sort_values("logFC")
        df_sk["label"] = (df_sk["EntrezGeneSymbol"] + " (" + df_sk["Target"].fillna("") + ")")
        df_sk["label"] = df_sk.apply(
            lambda r: r["EntrezGeneSymbol"]
            if (pd.isna(r["Target"]) or r["Target"] == "" or r["Target"] == r["EntrezGeneSymbol"])
            else f"{r['EntrezGeneSymbol']} ({r['Target']})", axis=1)
        df_sk["sig"]   = df_sk["P.Value"] < 0.05
        df_sk["color"] = df_sk["logFC"].map(lambda x: COL_VEN if x < 0 else COL_ART)

        SHAPE_MAP = {"SPINK": "circle", "KLK": "square", "SERPIN": "diamond"}

        fig_sk = go.Figure()
        for fam in ["SPINK", "KLK", "SERPIN"]:
            sub = df_sk[df_sk["family"] == fam]
            if sub.empty:
                continue
            tip = ("<b>" + sub["EntrezGeneSymbol"] + "</b> (" + sub["family"] + ")<br>"
                   + sub["Target"].fillna("") + "<br>"
                   + "logFC = " + sub["logFC"].round(3).astype(str) + "<br>"
                   + "p = " + sub["P.Value"].map(lambda x: f"{x:.3g}") + "<br>"
                   + sub["direction"])
            fig_sk.add_trace(go.Scatter(
                x=sub["logFC"], y=sub["label"],
                mode="markers", name=fam,
                text=tip, hoverinfo="text",
                marker=dict(
                    symbol=SHAPE_MAP.get(fam, "circle"),
                    color=sub["color"].tolist(),
                    size=sub["sig"].map(lambda s: 14 if s else 8).tolist(),
                    opacity=sub["sig"].map(lambda s: 0.9 if s else 0.45).tolist(),
                    line=dict(color="white", width=0.5),
                ),
            ))
        fig_sk.add_vline(x=0, line_dash="dot", line_color="gray", line_width=1)
        fig_sk.update_layout(
            title="Skin serine proteases & inhibitors: A-V gradient (larger dot = p < 0.05)",
            xaxis_title="log2FC (Arterial / Venous)   ← venous higher | arterial higher →",
            yaxis=dict(title="", tickfont=dict(size=10)),
            height=max(420, len(df_sk) * 22 + 80),
            margin=dict(l=200, t=50),
            legend=dict(title="Family (shape)"),
        )
        st.plotly_chart(fig_sk, use_container_width=True)

        col_t, col_kw = st.columns([3, 2])

        with col_t:
            st.markdown("**All detected family members**")
            disp_sk = (df_sk[["family", "EntrezGeneSymbol", "Target", "logFC",
                               "P.Value", "direction", "sig"]]
                       .rename(columns={"family": "Family", "EntrezGeneSymbol": "Gene",
                                        "Target": "Protein", "P.Value": "p",
                                        "direction": "Direction", "sig": "p<0.05"})
                       .sort_values("p")
                       .assign(logFC=lambda d: d["logFC"].round(3),
                               p=lambda d: d["p"].map(lambda x: f"{x:.3g}"),
                               **{"p<0.05": lambda d: d["p<0.05"].map(lambda v: "✓" if v else "")}))
            st.dataframe(disp_sk, use_container_width=True, hide_index=True)

        with col_kw:
            st.markdown("**Surgery-type gradient (exploratory, n = 3–5/group)**")
            if not skin_kw.empty:
                disp_kw = (skin_kw
                           .rename(columns={"gene_symbol": "Gene", "av_p": "A-V p",
                                            "kw_p": "KW p", "mean_HN": "Δ Head/Neck",
                                            "mean_Lap": "Δ Lap", "mean_Spine": "Δ Spine"})
                           .sort_values("A-V p")
                           .assign(**{"A-V p": lambda d: d["A-V p"].map(lambda x: f"{x:.3g}"),
                                      "KW p":  lambda d: d["KW p"].map(lambda x: f"{x:.3g}"),
                                      "Δ Head/Neck": lambda d: d["Δ Head/Neck"].round(3),
                                      "Δ Lap":       lambda d: d["Δ Lap"].round(3),
                                      "Δ Spine":     lambda d: d["Δ Spine"].round(3)}))
                st.dataframe(disp_kw, use_container_width=True, hide_index=True)
                st.caption("Δ = mean log2(arterial − venous). Negative = venous > arterial.")
            else:
                st.info("Re-run `10_skin_protease_analysis.R` to generate this table.")

        # Highlight skin genes on the main volcano
        st.markdown("---")
        st.markdown("**Skin protease family: positions on A-V volcano**")
        skin_genes = set(skin_av["EntrezGeneSymbol"].str.upper())
        df_sv = av.copy()
        df_sv["_grp"] = np.where(
            df_sv["EntrezGeneSymbol"].str.upper().isin(searched), "Searched",
            np.where(df_sv["EntrezGeneSymbol"].str.upper().isin(skin_genes), "Skin protease",
            df_sv["sig_label"]))
        skin_col_map = {"Searched": COL_SRC, "Skin protease": "#27ae60",
                        "FDR < 0.05": "#ff7f00", "p < 0.05": COL_NOM, "NS": COL_NS}
        fig_sv = go.Figure()
        for grp in ["NS", "p < 0.05", "FDR < 0.05", "Skin protease", "Searched"]:
            sub = df_sv[df_sv["_grp"] == grp]
            if sub.empty:
                continue
            tip = ("<b>" + sub["Target"].fillna("") + "</b> (" + sub["EntrezGeneSymbol"].fillna("")
                   + ")<br>logFC: " + sub["logFC"].round(3).astype(str)
                   + "  p=" + sub["P.Value"].map(lambda x: f"{x:.3g}"))
            fig_sv.add_trace(go.Scatter(
                x=sub["logFC"], y=sub["neg_log10p"], mode="markers", name=grp,
                text=tip, hoverinfo="text",
                marker=dict(color=skin_col_map[grp],
                            size=np.where(grp == "Skin protease", 7, 3),
                            opacity=0.2 if grp == "NS" else 0.85)
            ))
        sig_skin = df_sv[(df_sv["_grp"] == "Skin protease") & df_sv["sig_nom"]]
        for _, r in sig_skin.iterrows():
            fig_sv.add_annotation(x=r["logFC"], y=r["neg_log10p"],
                                  text=r["EntrezGeneSymbol"], showarrow=True, arrowhead=2,
                                  font=dict(size=10, color="#1a5e30"),
                                  bgcolor="white", bordercolor="#27ae60")
        fig_sv.add_hline(y=-np.log10(0.05), line_dash="dot", line_color=COL_NOM, line_width=1)
        fig_sv.add_vline(x=0, line_dash="dot", line_color="gray", line_width=1)
        fig_sv.update_layout(
            xaxis_title="log2FC (Arterial / Venous)",
            yaxis_title="−log10(p-value)",
            height=460, margin=dict(t=30),
        )
        st.plotly_chart(fig_sv, use_container_width=True)


# ── Footer ─────────────────────────────────────────────────────────────────────
st.markdown("---")
st.caption(
    "Prakash · Kurien · Chinn · UCSF Anesthesia & Perioperative Care · "
    "Data: SomaScan v4.1 (7,481 proteins) · N = 12 paired intraoperative draws"
)
