# Pipeline de RNA-Seq Automatizado com Snakemake (GSE161218)

Este repositório contém o pipeline completo e automatizado para a análise de expressão gênica diferencial (RNA-Seq) utilizando dados públicos do dataset **GSE161218**. O projeto foi desenvolvido como requisito para a disciplina de **RIB0305:Laboratório de Bioinformática**.

## 🧬 Contexto Biológico do Dataset
O estudo analisa o perfil transcricional de linhagens celulares de Leucemia Mieloide Aguda (LMA) submetidas a tratamentos farmacológicos com dois inibidores específicos:
* MEKi: Inibidor da via MEK/MAPK.
* FLT3i: Inibidor da tirosina quinase FLT3.
* ctrl: Amostras de controle (não tratadas).

O objetivo do pipeline é identificar genes diferencialmente expressos (DEGs) e as vias biológicas alteradas por esses tratamentos.


## 🛠️ Ferramentas Utilizadas e Pipeline
O fluxo de trabalho foi automatizado utilizando o gerenciador de fluxo **Snakemake**, integrando as seguintes ferramentas:

1. Download de Dados: sra-tools (prefetch + fasterq-dump)
2. Controle de Qualidade: FastQC + fastp (para remoção de adaptadores e bases de baixa qualidade)
3. Relatório de Qualidade: MultiQC (agrupamento dos relatórios em um único HTML)
4. Indexação e Quantificação: Salmon (quasi-alignment contra o transcriptoma de referência do GENCODE)
5. Análise Estatística: R + tximport + DESeq2 (normalização e expressão diferencial)
6.Enriquecimento Funcional: R + clusterProfiler (Ontologia Gênica - GO e KEGG)


## 📂 Estrutura do Repositório

```text
P8-GSE161218-STALINE/
├── .gitignore               # Arquivos pesados ignorados (FASTQ, referências)
├── README.md                # Manual de execução (Este arquivo)
├── RELATORIO.md             # Relatório científico final com os resultados
├── Snakefile                # Fluxo de trabalho automatizado do Snakemake
├── CONFIG/
│   ├── config.yaml          # Metadados e caminhos do projeto
│   └── samples.tsv          # Desenho experimental (amostras vs condições)
├── envs/
│   └── environment.yaml     # Ambiente Conda com todas as dependências
├── scripts/
│   ├── deseq2_analysis.R    # Script R para expressão diferencial
│   └── enrichment_analysis.R# Script R para GO/KEGG com clusterProfiler
└── results/                 # Diretório de saídas leve (tabelas e gráficos)
    ├── multiqc_report.html
    └── deseq2/
        ├── pca_plot.png
        ├── volcano_plot.png
        ├── heatmap_plot.png
        ├── ma_plot.png
        └── functional_enrichment.png