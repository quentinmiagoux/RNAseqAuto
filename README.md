# RNAseqAuto

Ce dépôt contient un exemple de pipeline **Snakemake** pour lancer automatiquement des compilations **R Markdown** (Rmd) à partir d'une liste de comparaisons.

## Structure

- `Snakefile` : définition du workflow Snakemake.
- `config.yaml` : fichier de configuration pointant vers la table des comparaisons.
- `comparisons.tsv` : table TSV listant les comparaisons et le Rmd à compiler.
- `rmd/comparison_report.Rmd` : exemple de template R Markdown.
- `results/` : sorties HTML générées.

## Pré-requis

- Snakemake (ex: `pip install snakemake`)
- R + le package `rmarkdown`

## Utilisation

1. Modifiez `comparisons.tsv` pour ajouter vos comparaisons :

```tsv
comparison	group_a	group_b	rmd
control_vs_treatment	control	treatment	rmd/comparison_report.Rmd
```

2. Lancez Snakemake :

```bash
snakemake --cores 1
```

Les rapports HTML seront générés dans `results/`.

## Personnalisation

- Modifiez `rmd/comparison_report.Rmd` pour intégrer vos analyses.
- Ajoutez d'autres colonnes dans `comparisons.tsv` et adaptez le `Snakefile` si besoin.
