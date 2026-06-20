# Relatório Científico: Reposicionamento de Fármacos e Sinergismo de Inibidores de Via em Câncer de Pulmão KRAS Mutante (GSE161218)

Disciplina:RIB0305_Laboratorio de Bioinformática 

Aluna: Staline Sahara Dianana_15426416  

Data: 16 de Junho de 2026  


1. Introdução

O Câncer de Pulmão de Células Não Pequenas (CPCNP) é a principal causa de morte por câncer no mundo, apresentando uma taxa de sobrevida bastante baixa. Entre os subtipos moleculares mais agressivos, as mutações de ganho de função no oncogene KRAS (especialmente as variantes G12C, G12V e G12D) estão presentes em cerca de 30% dos casos de adenocarcinoma pulmonar. A proteína KRAS mutante permanece constantemente ligada ao GTP, funcionando como um interruptor molecular permanentemente ativado. Esse estado hiperativo dispara continuamente sinais para que a célula se prolifere e sobreviva, utilizando principalmente a via das MAPKs (RAF/MEK/ERK) e a via PI3K/Akt/mTOR. Por muito tempo, o KRAS foi considerado um alvo terapêutico impossível de ser bloqueado diretamente por medicamentos devido à sua alta afinidade pelo GTP e à falta de encaixes óbvios em sua estrutura. Isso levou a ciência a buscar alternativas baseadas no reposicionamento de fármacos (drug repurposing) e no uso de terapias combinadas.

A linhagem celular humana H1792 é um excelente modelo em laboratório para estudar o adenocarcinoma de pulmão, pois carrega naturalmente a mutação no gene KRAS. Na prática clínica, o tratamento isolado com inibidores de MEK (MEKi, como o trametinib) costuma falhar rapidamente porque as células tumorais desenvolvem resistência adaptativa. Ao sofrerem o estresse do medicamento, as células ativam receptores de tirosina quinase (RTKs) alternativos na membrana para desviar o sinal, manter as rotas de sobrevivência e evitar a morte celular.

O receptor FLT3 (Fms-like tyrosine kinase 3) funciona como uma dessas rotas acessórias de escape. Assim, o reposicionamento de inibidores de FLT3 (FLT3i, como o quizartinib) originalmente criados para o tratamento de leucemias surge como uma estratégia inovadora para bloquear esses caminhos de sobrevivência também em tumores sólidos.

A pergunta central desta investigação é: Quais são as mudanças no perfil de expressão gênica e nas vias funcionais da linhagem H1792 sob o efeito do bloqueio isolado ou combinado de MEK e FLT3, e como essas assinaturas explicam o sinergismo entre as drogas? Para responder a isso, avaliou-se o transcriptoma global das células por sequenciamento de RNA (RNA-Seq) através de três comparações diretas: MEKi vs. Controle, FLT3i vs. Controle, e a comparação direta MEKi vs. FLT3i.


2. Materiais e Métodos
   
2.1. Dataset e Desenho Experimental

Este estudo utilizou dados públicos de sequenciamento de RNA em larga escala (RNA-Seq) depositados no repositório Gene Expression Omnibus (GEO) sob o código de acesso GSE161218. O organismo investigado foi Homo sapiens, utilizando a linhagem celular de adenocarcinoma de pulmão H1792 (KRAS mutante).

O desenho experimental contou com 3 grupos de tratamento testados em triplicatas biológicas, totalizando 9 amostras sequenciadas:
-Controle (Ctrl): Amostras basais tratadas apenas com o veículo DMSO (n=3: amostras SRR13020843, SRR13020844, SRR13020845);
-Inibidor de FLT3 (FLT3i): Células tratadas com o composto quizartinib (n=3: amostras SRR13020846, SRR13020847, SRR13020848);
-Inibidor de MEK (MEKi): Células tratadas com o composto trametinib (n=3: amostras SRR13020849, SRR13020850, SRR13020851).

2.2. Pipeline Computacional

O processamento dos dados de bioinformática foi gerenciado utilizando o software Snakemake:
1- Controle de Qualidade: A qualidade das leituras brutas foi avaliada com o FastQC (v0.11.9) para checar o Phred score por ciclo, conteúdo GC e presença de adaptadores.
2- Limpeza de Dados (Trimming): A remoção de adaptadores e filtros de qualidade foram feitos com o fastp (v0.23.2) em modo pareado (paired-end). Foram eliminadas bases com qualidade inferior a 20 (Q < 20) e descartadas leituras com comprimento menor que 36 pares de bases (-l 36). O relatório final foi consolidado com o MultiQC (v1.14).
3- Quantificação de Transcritos: Utilizou-se o transcriptoma de referência humano anotado do consórcio GENCODE (Release v50). A quantificação da expressão foi feita pelo algoritmo de quase-alinhamento do Salmon (v1.9.0), configurado para corrigir vieses de sequência (--seqBias) e validar rigorosamente o mapeamento (--validateMappings).

2.3. Análise Estatística no DESeq2

A conversão dos dados de transcritos para nível de gene e a importação para o ambiente R foram feitas com o pacote tximport. A análise de expressão diferencial foi realizada com o pacote DESeq2.
-Modelo Estatístico: O ajuste das contagens seguiu a distribuição Binomial Negativa através da fórmula: design =~condition; onde a variável condition representa os grupos (Ctrl, FLT3i, MEKi).
-Filtragem: Genes com soma total de contagens inferior a 10 foram removidos antes da análise.
-Contrastes: O teste de Wald foi aplicado para calcular os valores de $p$ em três comparações:
1- c("condition", "MEKi", "ctrl") — Efeitos do bloqueio de MEK.
2- c("condition", "FLT3i", "ctrl") — Efeitos do tratamento com FLT3i.
3- c("condition", "MEKi", "FLT3i") — Comparação direta entre os dois tratamentos.
-Critérios de Significância: Foram considerados Genes Diferencialmente Expressos (DEGs) aqueles com $p$-valor ajustado por Benjamini-Hochberg (padj < 0.05) e mudança na expressão de no mínimo duas vezes (log_2 Fold Change >= 1.0).

2.4. Enriquecimento de Vias por GSEA (WebGestaltR)

Para avaliar as mudanças biológicas de forma global, sem cortes arbitrários na lista de genes, utilizou-se o pacote WebGestaltR executando o método GSEA (Gene Set Enrichment Analysis). Os genes foram ordenados do mais ativado para o mais reprimido usando a fórmula:
Score de Rank = -log_10(p-value) * sign(log_2 FoldChange)
Os nomes dos genes (IDs Ensembl) foram limpos para remover os pontos decimais que indicavam as versões das sequências. A análise buscou caminhos metabólicos e de sinalização no banco de dados KEGG para Homo sapiens, considerando significativas as vias com FDR < 0.05.


3. Resultados

   3.1. Controle de Qualidade e Mapeamento

   Os relatórios do MultiQC confirmaram a alta qualidade dos dados, mantendo a pontuação Phred média acima de $Q30$ na quase totalidade das bases. O processo de limpeza pelo fastp removeu menos de 1,4% dos dados iniciais. As taxas de mapeamento calculadas pelo Salmon contra o referencial do GENCODE foram consistentes em todas as amostras, variando estritamente entre 82.4% e 86.1% de eficiência.

   3.2. Análise de Componentes Principais (PCA)

   O gráfico de PCA capturou as principais fontes de variação no perfil de expressão gênica da linhagem H1792. A Componente Principal 1 (PC1) explicou 75% da variação total, enquanto a Componente Principal 2 (PC2) reteve 9%.

O ponto mais importante foi o comportamento atípico da amostra de controle SRR13020843, que funcionou como um forte outlier, isolando-se completamente no lado esquerdo do gráfico ao longo do eixo da PC1.

Por outro lado, as outras duas réplicas de controle (SRR13020844 e SRR13020845) ficaram agrupadas de forma compacta e próxima ao grupo tratado com o inibidor de MEK (SRR13020849, SRR13020850, SRR13020851). O grupo tratado com o inibidor de FLT3 (SRR13020846, SRR13020847, SRR13020848) formou um grupo coeso no quadrante superior direito, mostrando uma clara separação em relação ao tratamento com MEKi ao longo do eixo vertical da PC2.

3.3. Padrões de Dispersão e Expressão Diferencial

- MA Plot: Mostrou a relação entre a quantidade média de leitura de cada gene e a intensidade da mudança de expressão (log_2FC). Os genes significativos (padj < 0.05) aparecem destacados em azul marinho. A dispersão revelou uma variação equilibrada de genes ativados e reprimidos, com alguns casos ultrapassando a marca de log_2Fold Change > 4.0.

- olcano Plot: Apresentou visualmente o balanço entre a significância estatística e a magnitude da mudança de expressão. Os genes robustos foram destacados acima das linhas de corte azuis. O gráfico confirmou uma distribuição proporcional entre genes que subiram (upregulation) e genes que desceram (downregulation) de expressão, destacando um gene específico que atingiu significância extrema no topo do gráfico (-log_p-value > 25).

  3.4. Análise de Clusterização Hierárquica (Heatmap)

  O mapa térmico gerado com os 20 genes de maior variação global confirmou a mesma organização observada no gráfico de PCA, dividindo as amostras em três ramificações principais:
  1- Ramificação Esquerda: Isolou a amostra controle outlier SRR13020843, comportamento causado por um bloco de 5 genes com forte expressão relativa (Z-score approx +2.0) que não se repetiu em nenhuma outra amostra.
  2- amificação Central: Reuniu as três amostras do grupo FLT3i, evidenciando um padrão de ativação exclusivo de 9 genes em resposta ao quizartinib.
  3- Ramificação Direita: Agrupou as amostras tratadas com o inibidor de MEK (MEKi) juntamente com as duas réplicas normais do grupo Controle (SRR13020844 e SRR13020845), demonstrando que esses perfis compartilham uma assinatura de expressão muito próxima.


  3.5. Enriquecimento de Vias por GSEA (WebGestaltR)

  Ao avaliar os genes ordenados no contraste MEKi vs. Controle, o WebGestaltR identificou vias alteradas no banco do KEGG. O tratamento com o inibidor de MEK resultou em Scores de Enriquecimento Normalizados (NES) negativos, demonstrando uma repressão coordenada dessas vias.

As três principais vias identificadas como estatisticamente significativas estão detalhadas na tabela abaixo:
ID da Via (KEGG),      Nome da Via,                     NES,    Valor de p,    FDR,          Status Funcional
hsa04657,        Via de Sinalização de IL-17,          -2.08,      <0.001,     <0.01,         Reprimida (Desligada)
hsa04512,    Interação Matriz Extracelular-Receptor,   -1.95,      <0.001,     <0.01,        Reprimida (Desligada)
hsa04114,      Meiose de Oócitos (Ciclo Celular),      -1.74,      <0.002,     <0.03,         Reprimida (Desligada)


4. Discussão

A reanálise dos dados do dataset GSE161218 sob a perspectiva do reposicionamento de fármacos traz informações importantes sobre como as células tumorais respondem aos tratamentos. Inicialmente, nota-se que o comportamento atípico da amostra de controle SRR13020843 alterou o perfil de expressão basal. Em análises que usam listas cortadas de genes (como no g:Profiler), essa distorção tendeu a destacar o metabolismo de purinas e nucleobases. Isso ocorre porque células tumorais impulsionadas por KRAS mutante dependem fortemente da produção de nucleotídeos para manter suas altas taxas de divisão e crescimento.

Ao utilizar o algoritmo GSEA no WebGestaltR, que analisa a variação de todos os genes de forma contínua e reduz o impacto de ruídos individuais, foi possível identificar caminhos biológicos específicos desligados pela ação do inibidor de MEK (trametinib).

A forte repressão observada na via hsa04114 (Meiose de Oócitos) reflete diretamente o bloqueio dos sinais de crescimento disparados pelo KRAS. Embora essa via tenha esse nome por seu papel histórico em células germinativas, no contexto do câncer ela compartilha genes chaves que regulam os pontos de controle da divisão celular em células somáticas. Portanto, o tratamento com MEKi reprime a maquinaria de replicação e trava a proliferação do tumor.

Da mesma forma, o desligamento da via hsa04512 (Interação Matriz Extracelular-Receptor) indica que o trametinib altera a forma como a célula interage com o ambiente ao seu redor. Essa via gerencia a expressão de proteínas e receptores que permitem que as células tumorais se fixem e se desloquem pelos tecidos. A redução dessas assinaturas sugere uma perda na capacidade de adesão e locomoção, o que pode diminuir o potencial de invasão e metástase do tumor.

Além disso, o bloqueio da via hsa04657 (Sinalização de IL-17) aponta que a inibição de MEK também afeta o perfil inflamatório da linhagem, atenuando sinais que poderiam favorecer a sobrevivência das células tumorais.

4.1. Comparação com o Artigo Original e Mecanismo de Sinergismo

Os resultados encontrados nesta análise convergem de forma consistente com a proposta biológica do artigo de referência deste dataset (Vicent, Román et al., Nature Communications, 2023). No trabalho original, os autores demonstraram que tumores de pulmão com mutação no gene KRAS resistem ao tratamento isolado com inibidores de MEK ativando receptores de tirosina quinase (RTKs) na membrana celular como uma rota alternativa. O estudo de Vicent et al. utilizou assinaturas de expressão gênica para sugerir o reposicionamento da midostaurina um inibidor multiquinase que atua em receptores como o FLT3  de modo a bloquear esses desvios e superar a resistência ao tratamento.

Nesta abordagem prática, avaliou-se o bloqueio específico do eixo FLT3 utilizando o quizartinib (FLT3i). A distribuição das amostras no gráfico de PCA e no Heatmap ilustra claramente a base que justifica combinar os dois medicamentos: enquanto o grupo tratado com MEKi altera as rotas de proliferação (agrupando-se perto do perfil de controle normal), o tratamento isolado com FLT3i gera uma resposta compensatória diferente, marcada pela ativação exclusiva daquele bloco de 9 genes visto no mapa térmico.

Essa diferença no perfil de expressão confirma que o sinergismo ocorre por ações complementares. O inibidor de MEK (trametinib) age desligando a cascata central de sinalização do KRAS mutante, o que resulta na repressão do ciclo celular, das vias inflamatórias e do remodelamento de matriz que identificamos no GSEA. Ao mesmo tempo, as células tentam contornar o bloqueio ativando rotas alternativas de escape, como o receptor FLT3.

Ao associar o inibidor de FLT3 (quizartinib) como agente de reposicionamento, bloqueiam-se os mecanismos de suporte que o tumor usa para sobreviver. Dessa forma, a combinação dos dois fármacos inviabiliza as principais rotas de escape da célula tumoral, oferecendo uma base científica sólida que apoia o uso de terapias combinadas para evitar a resistência a monoterapias em tumores pulmonares KRAS mutantes.


5. Referências

1- Artigo do Dataset (Referência Principal):
VICENT, S.; ROMÁN, M. et al. Signature-driven repurposing of Midostaurin for combination with MEK1/2 and KRASG12C inhibitors in lung cancer. Nature Communications, v. 14, n. 1, p. 6332, 2023. DOI: 10.1038/s41467-023-41828-z. PMID: 37816716.

2- Controle de Qualidade (FastQC):
ANDREWS, S. FastQC: A Quality Control Tool for High Throughput Sequence Data. v. 0.11.9, 2010. Disponível em: http://www.bioinformatics.babraham.ac.uk/projects/fastp/.

3- Trimming e Pré-processamento (fastp):
CHEN, S. fastp: An ultra-fast all-round tool for FASTQ data quality control and preprocessing. iMeta, v. 4, n. 5, e70078, 2025. DOI: 10.1002/imt2.70078.

4- Quantificação de Transcritos (Salmon):
PATRO, R.; DUGGAL, G.; LOVE, M. I.; ISTITIAPONPISAIN, S.; KINGSFORD, C. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, v. 14, n. 4, p. 417-419, 2017. DOI: 10.1038/nmeth.4197.

5- Importação de Dados para o R (tximport):
SONESON, C.; LOVE, M. I.; ROBINSON, M. D. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research, v. 4, p. 1521, 2015. DOI: 10.12688/f1000research.7563.1.

6- Expressão Gênica Diferencial (DESeq2):
LOVE, M. I.; HUBER, W.; ANDERS, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, v. 15, n. 12, p. 550, 2014. DOI: 10.1186/s13059-014-0550-8.

7- Análise de Enriquecimento Funcional (WebGestalt / WebGestaltR):
WANG, J.; VASAWARALA, S.; ZNAYENKO, Z.; SHEN, L.; ZHANG, B. WebGestalt 2019: gene set elimination, enrichment and network analyses on functional genomics data. Nucleic Acids Research, v. 47, n. W1, p. W99-W105, 2019. DOI: 10.1093/nar/gkz401.
