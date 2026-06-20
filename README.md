# Pipeline de RNA-Seq Automatizado com Snakemake (GSE161218)

Este repositório contém o pipeline completo e automatizado para a análise de expressão gênica diferencial (RNA-Seq) utilizando dados públicos do dataset GSE161218. O projeto foi desenvolvido como requisito prático para a disciplina RIB0305: Laboratório de Bioinformática.

## 🧬 Contexto Biológico do Dataset
O estudo analisa o perfil transcriptômico da linhagem celular H1792 de Câncer de Pulmão de Células Não Pequenas (CPCNP), portadora de mutação ativadora no oncogene KRAS. A investigação avalia a resposta celular frente a tratamentos farmacológicos baseados em monoterapia e estratégias de reposicionamento de fármacos:

* MEKi: Inibidor da via MEK/MAPK (trametinib), alvo central downstream do KRAS.
* FLT3i: Inibidor da tirosina quinase FLT3 (quizartinib), avaliado como bloqueador de rotas acessórias de escape tumoral.
* ctrl: Amostras de controle basal (tratadas apenas com o veículo DMSO).

O objetivo do pipeline é mapear os rearranjos transcricionais globais, identificar Genes Diferencialmente Expressos (DEGs) e desvendar as vias biológicas (KEGG) envolvidas no sinergismo dessas drogas para superação da resistência terapêutica.

## 🛠️ Ferramentas Utilizadas e Pipeline
O fluxo de trabalho foi totalmente automatizado utilizando o gerenciador de fluxo Snakemake, integrando as seguintes ferramentas bioinformáticas:

1. Download de Dados: `sra-tools` (`prefetch` + `fasterq-dump`) para obtenção dos arquivos FastQ públicos.
2. Controle de Qualidade: `FastQC` para análise diagnóstica inicial das leituras brutas.
3. Trimming e Pré-processamento: `fastp` para remoção automatizada de adaptadores e filtragem de bases de baixa qualidade (Q < 20).
4. Relatório Consolidado: `MultiQC` para agregação dos logs de qualidade em um único relatório HTML interativo.
5. Indexação e Quantificação: `Salmon` para o quase-alinhamento de alta velocidade contra o transcriptoma de referência humano (`GENCODE` Release v50).
6. Expressão Gênica Diferencial: Ambiente `R` com os pacotes `tximport` (mapeamento transcrito-gene) e `DESeq2` (normalização por distribuição Binomial Negativa e testes de Wald).
7. Análise de Enriquecimento Funcional: Ambiente `R` utilizando o pacote `WebGestaltR` para a execução do algoritmo GSEA (Gene Set Enrichment Analysis) contra o banco de dados de vias do `KEGG`.

## 📂 Estrutura do Repositório

P8-GSE161218-STALINE/
├── .gitignore                # Arquivos pesados ignorados (FASTQ, índices, referências)
├── README.md                 # Manual de execução e introdução (Este arquivo)
├── RELATORIO.md              # Relatório científico final formatado com a discussão biológica
├── Snakefile                 # Script de automação do fluxo de trabalho do Snakemake
├── CONFIG/
│   ├── config.yaml          # Parâmetros, metadados e caminhos globais do projeto
│   └── samples.tsv          # Desenho experimental e metadados (amostras vs condições)
├── envs/
│   └── environment.yaml     # Ambiente Conda com todas as dependências e versões fixadas
├── scripts/
│   ├── deseq2_analysis.R    # Script R para controle de qualidade, PCA e expressão diferencial
│   └── enrichment_analysis.R# Script R para análise de enriquecimento GSEA com WebGestaltR
└── results/                  # Diretório de saídas e produtos gráficos da análise
    ├── multiqc_report.html   # Relatório consolidado de qualidade do sequenciamento
    └── deseq2/
        ├── pca_plot.png      # Gráfico de Análise de Componentes Principais
        ├── volcano_plot.png  # Gráfico Volcano destacando a significância dos DEGs
        ├── heatmap_plot.png  # Clusterização hierárquica dos top genes de maior variância
        ├── ma_plot.png       # Gráfico de dispersão M-A (Log-fold change vs abundância média)
        └── functional_enrichment.png  # Gráfico do GSEA com as principais vias KEGG alteradas
