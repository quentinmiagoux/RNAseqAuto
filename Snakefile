configfile: "config.yaml"

import csv
from pathlib import Path

COMPARISON_FIELDS = {"comparison", "group_a", "group_b", "rmd"}


def load_comparisons(path):
    comparisons = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = COMPARISON_FIELDS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"comparisons_file is missing required columns: {sorted(missing)}"
            )
        for row in reader:
            comparisons.append({key: (value or "").strip() for key, value in row.items()})
    return comparisons


comparisons = load_comparisons(config["comparisons_file"])
comparisons_by_name = {row["comparison"]: row for row in comparisons}

rule all:
    input:
        expand("results/{comparison}.html", comparison=comparisons_by_name.keys())


rule init_results:
    output:
        touch("results/.initialized")
    shell:
        """
        mkdir -p results
        """


rule render_rmd:
    input:
        rmd=lambda wc: comparisons_by_name[wc.comparison]["rmd"],
        init="results/.initialized",
    output:
        "results/{comparison}.html"
    params:
        group_a=lambda wc: comparisons_by_name[wc.comparison]["group_a"],
        group_b=lambda wc: comparisons_by_name[wc.comparison]["group_b"],
        design_formula=lambda wc: config.get(
            "design_formula", "~ Reprogrammation + Age + sex + condition2"
        ),
    shell:
        r"""
        Rscript -e "rmarkdown::render('{input.rmd}', output_file='{output}', params=list(group_a='{params.group_a}', group_b='{params.group_b}', comparison='{wildcards.comparison}', design_formula='{params.design_formula}'))"
        """
