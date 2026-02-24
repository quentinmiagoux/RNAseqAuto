configfile: "config.yaml"

import csv
from pathlib import Path

PROJECT_DIR = Path(workflow.basedir).resolve()


def abs_path(p):
    p = Path(p)
    return str(p if p.is_absolute() else (PROJECT_DIR / p).resolve())


def r_escape(s):
    return (s or "").replace("\\", "\\\\").replace('"', '\\"')


def r_bool(value):
    return "TRUE" if value else "FALSE"


def normalize_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return [str(item).strip() for item in value if str(item).strip()]


def r_str_vector(values):
    if not values:
        return "character(0)"
    escaped = [f"\\\"{r_escape(item)}\\\"" for item in values]
    return f"c({', '.join(escaped)})"


def tag_value(value):
    text = str(value)
    return text.replace(" ", "").replace(".", "p").replace("-", "m")


def make_run_tag(cfg):
    return "_".join(
        [
            f"lfc{tag_value(cfg.get('lfc', 0.0))}",
            f"padj{tag_value(cfg.get('padj', 0.05))}",
            f"bm{tag_value(cfg.get('bm', 50))}",
            f"lfct{tag_value(cfg.get('lfc_threshold', 0.0))}",
            f"lfcse{tag_value(cfg.get('lfcse', 1))}",
        ]
    )


def disease_folder_name(name):
    cleaned = (name or "unknown").strip()
    return "_".join(cleaned.split())


COMPARISON_FIELDS = {"comparison", "group_a", "group_b", "rmd"}


def load_comparisons(path):
    comparisons = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = COMPARISON_FIELDS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"comparisons_file is missing required columns: {sorted(missing)}")
        for row in reader:
            comparisons.append({k: (v or "").strip() for k, v in row.items()})
    return comparisons


comparisons = load_comparisons(config["comparisons_file"])
comparisons_by_name = {row["comparison"]: row for row in comparisons}
disease_by_comparison = {
    row["comparison"]: disease_folder_name(row.get("group_a", "")) for row in comparisons
}
merge_replicates_enabled = bool(config.get("merge_replicates", False))
comparison_suffix = "_merged" if merge_replicates_enabled else "_unmerged"
active_comparisons = sorted(
    name for name in comparisons_by_name.keys() if name.endswith(comparison_suffix)
)
active_diseases = [disease_by_comparison[comparison] for comparison in active_comparisons]
samples_to_exclude = normalize_list(config.get("samples_to_exclude", []))
RUN_TAG = make_run_tag(config)
comparison_html_outputs = [
    f"results/{RUN_TAG}/{disease_by_comparison[comparison]}/{comparison}_{RUN_TAG}.html"
    for comparison in active_comparisons
]
comparison_deg_outputs = [
    f"results/{RUN_TAG}/{disease_by_comparison[comparison]}/DEGs_{comparison}_{RUN_TAG}.xlsx"
    for comparison in active_comparisons
]

dep_diseases = ["DMD", "DN", "GSDII", "LMNA", "DNM2"]
dep_xlsx_outputs = [
    f"results/{RUN_TAG}/{disease}/DEPs_{disease}_{RUN_TAG}.xlsx"
    for disease in dep_diseases
]

rna_prot_overlap_html = f"results/{RUN_TAG}/rnaseq_proteomics_overlap_{RUN_TAG}.html"
rna_prot_overlap_xlsx = f"results/{RUN_TAG}/rnaseq_proteomics_overlap_{RUN_TAG}.xlsx"




rule all:
    input:
        comparison_html_outputs,
        comparison_deg_outputs,
        dep_xlsx_outputs,
        "results/" + RUN_TAG + "/rnaseq_deg_comparison_" + RUN_TAG + ".html",
        "results/" + RUN_TAG + "/qc_all_individuals.html",
        "results/" + RUN_TAG + "/DEP.html",
        rna_prot_overlap_html,
        rna_prot_overlap_xlsx

rule render_rmd:
    input:
        rmd=lambda wc: abs_path(comparisons_by_name[wc.comparison]["rmd"]),
        expr_xlsx=lambda wc: abs_path(config["gene_expression_xlsx"]),
        coldata_xlsx=lambda wc: abs_path(config["coldata_xlsx"]),
    output:
        html=f"results/{RUN_TAG}" + "/{disease}/{comparison}_" + RUN_TAG + ".html",
        deg_xlsx=f"results/{RUN_TAG}" + "/{disease}/DEGs_{comparison}_" + RUN_TAG + ".xlsx",
    log:
        "logs/render_rmd/{disease}/{comparison}." + RUN_TAG + ".log"
    params:
        group_a=lambda wc: comparisons_by_name[wc.comparison]["group_a"],
        group_b=lambda wc: comparisons_by_name[wc.comparison]["group_b"],
        design_formula=lambda wc: r_escape(config.get("design_formula", "~ condition2")),
        merge_replicates=lambda wc: r_bool(config.get("merge_replicates", False)),
        lfc_threshold=lambda wc: config.get("lfc_threshold", 0.0),
        lfc=lambda wc: config.get("lfc", 0.0),
        padj=lambda wc: config.get("padj", 0.05),
        bm=lambda wc: config.get("bm", 50),
        lfcse=lambda wc: config.get("lfcse", 1),
        samples_to_exclude=lambda wc: r_str_vector(samples_to_exclude),
        output_root=abs_path(f"results/{RUN_TAG}"),
        disease_folder=lambda wc: wc.disease,
        run_tag=RUN_TAG,
    shell:
        r"""
        mkdir -p "results/{params.run_tag}/{params.disease_folder}" "results/{params.run_tag}/.knit/{wildcards.comparison}_{params.run_tag}" "logs/render_rmd/{wildcards.disease}"
        Rscript -e "rmarkdown::render(
          \"{input.rmd}\",
          output_file=\"{wildcards.comparison}_{params.run_tag}.html\",
          output_dir=\"results/{params.run_tag}/{params.disease_folder}\",
          intermediates_dir=\"results/{params.run_tag}/.knit/{wildcards.comparison}_{params.run_tag}\",
          params=list(
            group_a=\"{params.group_a}\",
            group_b=\"{params.group_b}\",
            comparison=\"{wildcards.comparison}\",
            run_tag=\"{params.run_tag}\",
            output_root=\"{params.output_root}\",
            disease_folder=\"{params.disease_folder}\",
            design_formula=\"{params.design_formula}\",
            merge_replicates={params.merge_replicates},
            lfc_threshold={params.lfc_threshold},
            lfc={params.lfc},
            padj={params.padj},
            bm={params.bm},
            lfcse={params.lfcse},
            samples_to_exclude={params.samples_to_exclude},
            expr_xlsx=\"{input.expr_xlsx}\",
            coldata_xlsx=\"{input.coldata_xlsx}\"
          )
        )" &> "{log}"
        """


rule render_qc_all_individuals:
    input:
        rmd=abs_path("rmd/RNAseq_QC_All_individuals_gold_standard.Rmd"),
        expr_xlsx=abs_path(config["gene_expression_xlsx"]),
        coldata_xlsx=abs_path(config["coldata_xlsx"]),
    output:
        "results/" + RUN_TAG + "/qc_all_individuals.html"
    log:
        "logs/render_rmd/qc_all_individuals.log"
    params:
        top_variable_genes=lambda wc: config.get("top_variable_genes", 500),
        samples_to_exclude=lambda wc: r_str_vector(samples_to_exclude),
    shell:
        r"""
        mkdir -p "results/{RUN_TAG}" "results/{RUN_TAG}/.knit/qc_all_individuals"
        Rscript -e "rmarkdown::render(
          \"{input.rmd}\",
          output_file=\"qc_all_individuals.html\",
          output_dir=\"results/{RUN_TAG}\",
          intermediates_dir=\"results/{RUN_TAG}/.knit/qc_all_individuals\",
          params=list(
            expr_xlsx=\"{input.expr_xlsx}\",
            coldata_xlsx=\"{input.coldata_xlsx}\",
            top_variable_genes={params.top_variable_genes},
            samples_to_exclude={params.samples_to_exclude}
          )
        )" &> "{log}"
        """


rule render_rnaseq_deg_comparison:
    input:
        rmd=abs_path("rmd/RNAseq.Rmd"),
        deg_files=comparison_deg_outputs,
    output:
        "results/" + RUN_TAG + "/rnaseq_deg_comparison_" + RUN_TAG + ".html"
    log:
        "logs/render_rmd/rnaseq_deg_comparison." + RUN_TAG + ".log"
    params:
        deg_root_dir=abs_path(f"results/{RUN_TAG}"),
        run_tag=RUN_TAG,
        lfc=lambda wc: config.get("lfc", 2),
        padj=lambda wc: config.get("padj", 0.05),
        lfc_threshold=lambda wc: config.get("lfc_threshold", 0.0),
        bm=lambda wc: config.get("bm", 50),
        lfcse=lambda wc: config.get("lfcse", 1),
    shell:
        r"""
        mkdir -p "results/{params.run_tag}" "results/{params.run_tag}/.knit/rnaseq_deg_comparison_{params.run_tag}"
        Rscript -e "rmarkdown::render(
          \"{input.rmd}\",
          output_file=\"rnaseq_deg_comparison_{params.run_tag}.html\",
          output_dir=\"results/{params.run_tag}\",
          intermediates_dir=\"results/{params.run_tag}/.knit/rnaseq_deg_comparison_{params.run_tag}\",
          params=list(
            deg_root_dir=\"{params.deg_root_dir}\",
            run_tag=\"{params.run_tag}\",
            lfc={params.lfc},
            padj={params.padj},
            lfc_threshold={params.lfc_threshold},
            bm={params.bm},
            lfcse={params.lfcse}
          )
        )" &> "{log}"
        """



rule render_dep:
    input:
        rmd=abs_path("rmd/DEP.Rmd"),
        msstats_quant=abs_path("data/msstats_quant_result.xls"),
    output:
        html="results/" + RUN_TAG + "/DEP.html",
        dep_xlsx=dep_xlsx_outputs,
    log:
        "logs/render_rmd/dep." + RUN_TAG + ".log"
    params:
        run_tag=RUN_TAG,
        output_root=abs_path(f"results/{RUN_TAG}"),
        lfc=lambda wc: config.get("lfc", 2),
        padj=lambda wc: config.get("padj", 0.05),
    shell:
        r"""
        mkdir -p "results/{params.run_tag}" "results/{params.run_tag}/.knit/dep_{params.run_tag}"
        Rscript -e "rmarkdown::render(
          \"{input.rmd}\",
          output_file=\"DEP.html\",
          output_dir=\"results/{params.run_tag}\",
          intermediates_dir=\"results/{params.run_tag}/.knit/dep_{params.run_tag}\",
          params=list(
            msstats_quant=\"{input.msstats_quant}\",
            output_root=\"{params.output_root}\",
            run_tag=\"{params.run_tag}\",
            lfc={params.lfc},
            padj={params.padj}
          )
        )" &> "{log}"
        """


rule render_rnaseq_proteomics_overlap:
    input:
        rmd=abs_path("rmd/RNAseq_Proteomics_Comparison.Rmd"),
        rna_deg_files=comparison_deg_outputs,
        dep_files=dep_xlsx_outputs,
    output:
        html=rna_prot_overlap_html,
        xlsx=rna_prot_overlap_xlsx,
    log:
        "logs/render_rmd/rnaseq_proteomics_overlap." + RUN_TAG + ".log"
    params:
        output_root=abs_path(f"results/{RUN_TAG}"),
        run_tag=RUN_TAG,
        min_count=3,
        lfc=lambda wc: config.get("lfc", 2),
        padj=lambda wc: config.get("padj", 0.05),
    shell:
        r"""
        mkdir -p "results/{params.run_tag}" "results/{params.run_tag}/.knit/rnaseq_proteomics_overlap_{params.run_tag}"
        Rscript -e "rmarkdown::render(
          \"{input.rmd}\",
          output_file=\"rnaseq_proteomics_overlap_{params.run_tag}.html\",
          output_dir=\"results/{params.run_tag}\",
          intermediates_dir=\"results/{params.run_tag}/.knit/rnaseq_proteomics_overlap_{params.run_tag}\",
          params=list(
            output_root=\"{params.output_root}\",
            run_tag=\"{params.run_tag}\",
            min_count={params.min_count},
            lfc={params.lfc},
            padj={params.padj}
          )
        )" &> "{log}"
        """
