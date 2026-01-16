configfile: "config.yaml"

import csv
from pathlib import Path

PROJECT_DIR = Path(workflow.basedir).resolve()

def abs_path(p):
    p = Path(p)
    return str(p if p.is_absolute() else (PROJECT_DIR / p).resolve())

def r_escape(s):
    return (s or "").replace("\\", "\\\\").replace('"', '\\"')

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

rule all:
    input:
        expand("results/{comparison}.html", comparison=sorted(comparisons_by_name.keys()))

rule render_rmd:
    input:
        rmd=lambda wc: abs_path(comparisons_by_name[wc.comparison]["rmd"]),
        expr_xlsx=lambda wc: abs_path(config["gene_expression_xlsx"]),
        coldata_xlsx=lambda wc: abs_path(config["coldata_xlsx"]),
    output:
        "results/{comparison}.html"
    params:
        group_a=lambda wc: comparisons_by_name[wc.comparison]["group_a"],
        group_b=lambda wc: comparisons_by_name[wc.comparison]["group_b"],
        design_formula=lambda wc: r_escape(
            config.get("design_formula", "~ Reprogrammation + Age + sex + condition2")
        ),
    shell:
        r"""
        mkdir -p results
        Rscript -e "rmarkdown::render(
          \"{input.rmd}\",
          output_file=\"{wildcards.comparison}.html\",
          output_dir=\"results\",
          params=list(
            group_a=\"{params.group_a}\",
            group_b=\"{params.group_b}\",
            comparison=\"{wildcards.comparison}\",
            design_formula=\"{params.design_formula}\",
            expr_xlsx=\"{input.expr_xlsx}\",
            coldata_xlsx=\"{input.coldata_xlsx}\"
          )
        )"
        """
