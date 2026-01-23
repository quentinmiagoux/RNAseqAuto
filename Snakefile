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
merge_replicates_enabled = bool(config.get("merge_replicates", False))
comparison_suffix = "_merged" if merge_replicates_enabled else "_unmerged"
active_comparisons = sorted(
    name for name in comparisons_by_name.keys() if name.endswith(comparison_suffix)
)

rule all:
    input:
        expand("results/{comparison}.html", comparison=active_comparisons),
        "results/qc_all_individuals.html"

rule render_rmd:
    input:
        rmd=lambda wc: abs_path(comparisons_by_name[wc.comparison]["rmd"]),
        expr_xlsx=lambda wc: abs_path(config["gene_expression_xlsx"]),
        coldata_xlsx=lambda wc: abs_path(config["coldata_xlsx"]),
    output:
        "results/{comparison}.html"
    log:
        "logs/render_rmd/{comparison}.log"
    params:
        group_a=lambda wc: comparisons_by_name[wc.comparison]["group_a"],
        group_b=lambda wc: comparisons_by_name[wc.comparison]["group_b"],
        design_formula=lambda wc: r_escape(
            config.get("design_formula", "~ condition2")
        ),
        merge_replicates=lambda wc: r_bool(config.get("merge_replicates", False)),
        lfc_threshold=lambda wc: config.get("lfc_threshold", 0.0),
    shell:
        r"""
        mkdir -p results results/.knit/{wildcards.comparison}
        Rscript -e "rmarkdown::render(
          \"{input.rmd}\",
          output_file=\"{wildcards.comparison}.html\",
          output_dir=\"results\",
          intermediates_dir=\"results/.knit/{wildcards.comparison}\",
          params=list(
            group_a=\"{params.group_a}\",
            group_b=\"{params.group_b}\",
            comparison=\"{wildcards.comparison}\",
            design_formula=\"{params.design_formula}\",
            merge_replicates={params.merge_replicates},
            lfc_threshold={params.lfc_threshold},
            expr_xlsx=\"{input.expr_xlsx}\",
            coldata_xlsx=\"{input.coldata_xlsx}\"
          )
        )" &> "{log}"
        """

rule render_qc_all_individuals:
    input:
        rmd=abs_path("rmd/qc_all_individuals.Rmd"),
        expr_xlsx=abs_path(config["gene_expression_xlsx"]),
        coldata_xlsx=abs_path(config["coldata_xlsx"]),
    output:
        "results/qc_all_individuals.html"
    log:
        "logs/render_rmd/qc_all_individuals.log"
    params:
        top_variable_genes=lambda wc: config.get("top_variable_genes", 500),
    shell:
        r"""
        mkdir -p results results/.knit/qc_all_individuals
        Rscript -e "rmarkdown::render(
          \"{input.rmd}\",
          output_file=\"qc_all_individuals.html\",
          output_dir=\"results\",
          intermediates_dir=\"results/.knit/qc_all_individuals\",
          params=list(
            expr_xlsx=\"{input.expr_xlsx}\",
            coldata_xlsx=\"{input.coldata_xlsx}\",
            top_variable_genes={params.top_variable_genes}
          )
        )" &> "{log}"
        """
