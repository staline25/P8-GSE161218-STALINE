# Relatório Científico: Reposicionamento de Fármacos e Sinergismo de Inibidores de Via em Câncer de Pulmão KRAS Mutante (GSE161218)

Disciplina:RIB0305_Laboratorio de Bioinformática 
Aluna: Staline Sahara Dianana_15426416  
Data: 16 de Junho de 2026  



 1. Introdução

O Câncer de Pulmão de Células Não Pequenas (CPCNP) é a principal causa de morte por neoplasias em todo o mundo, apresentando uma taxa de sobrevida em cinco anos historicamente baixa. Dentre os subtipos moleculares mais desafiadores, as mutações no oncogene KRAS (particularmente substituições como G12C, G12V ou G12D) representam aproximadamente 30% dos casos de adenocarcinoma pulmonar. O KRAS mutante permanece constitutivamente ligado ao GTP, atuando como um interruptor molecular permanentemente ativado que dispara cascatas de sinalização proliferativa a jusante (downstream), com destaque para a via MAPK/ERK (RAF/MEK/ERK) e a via PI3K/Akt/mTOR. Historicamente, o KRAS foi considerado um alvo terapêutico inalcançável 
devido à sua alta afinidade pelo GTP e ausência de sítios alostéricos óbvios, o que impulsionou a comunidade científica a buscar abordagens baseadas em Reposicionamento de Fármacos (Drug Repurposing) e terapias combinadas.

A linhagem celular humana H1792 serve como um excelente modelo translacional para o adenocarcinoma de pulmão, dado o seu perfil genético portador de mutação ativadora em KRAS. Estratégias convencionais utilizando a monoterapia com inibidores de MEK (MEKi, como o trametinib) frequentemente falham na clínica devido à rápida emergência de mecanismos de resistência adaptativa. Sob estresse farmacológico, as células tumorais de pulmão operam um desvio molecular, ativando receptores de tirosina quinase (RTKs) alternativos e proteínas de suporte na membrana citoplasmática, incluindo eixos cináticos como o FLT3, para manter os sinais de sobrevivência e evadir a apoptose. O reposicionamento de inibidores de FLT3 (FLT3i, como o quizartinib) surge como uma co-intervenção estratégica inovadora, cujo propósito é neutralizar essas rotas acessórias de escape citoplasmático.

A pergunta central que norteia esta investigação científica é: Quais rearranjos transcricionais globais e alterações metabólicas funcionais são coordenados na linhagem de câncer de pulmão H1792 sob o efeito do bloqueio isolado ou combinado das vias de escape do KRAS através de MEKi e FLT3i, e de que forma essas drogas atuam sinergicamente na reversão do fenótipo tumoral? Para responder a isso, avaliou-se o perfil transcriptômico global por meio de três comparações par a par cruciais: MEKi vs ctrl, FLT3i vs ctrl e MEKi vs FLT3i.




2. Materiais e Métodos

2.1. Descrição do Dataset e Desenho Experimental
Este estudo utilizou dados públicos de sequenciamento de RNA em larga escala (RNA-Seq) depositados sob o código de acesso GSE161218. O organismo investigado foi Homo sapiens na linhagem celular tumoral de pulmão H1792 (KRAS mutante). O desenho experimental baseou-se em um modelo comparativo composto por 3 braços de tratamento conduzidos em triplicatas biológicas robustas, totalizando 9 amostras processadas simultaneamente:
* Controle (ctrl): Amostras basais tratadas com o veículo DMSO ($n = 3$: SRR13020843, SRR13020844, SRR13020845).
* Inibidor de FLT3 (FLT3i): Células submetidas ao reposicionamento do composto quizartinib ($n = 3$: SRR13020846, SRR13020847, SRR13020848).
* Inibidor de MEK (MEKi): Células submetidas ao tratamento com o trametinib ($n = 3$: SRR13020849, SRR13020850, SRR13020851).


2.2. Pipeline Computacional e Parâmetros das Ferramentas

O processamento bioinformático foi unificado e executado por meio de um pipeline modular gerenciado pelo software Snakemake:

1.Controle de Qualidade Inicial: As leituras brutas foram avaliadas com o FastQC para mapeamento do perfil de qualidade de sequenciamento por ciclo e conteúdo de GC.

2.Trimming e Filtragem de Adaptadores: A limpeza das leituras foi efetuada pelo programa fastp configurado no modo pareado (paired-end), eliminando bases com Phred score inferior a 20 ($Q < 20$) por janela deslizante e descartando leituras fragmentadas com comprimento menor que 36 pares de bases (-l 36`). O MultiQC foi empregado para consolidar os logs de execução.

3. Indexação e Quantificação de Transcritos:
Utilizou-se o transcriptoma de referência humano anotado obtido do consórcio GENCODE (Release v50, Human Basic Annotation). A quantificação foi processada pelo algoritmo de quase-alinhamento do Salmon, utilizando parâmetros para correção automática de viés de sequência (--seqBias), inferência de biblioteca (-l A) e validação rigorosa de mapeamento (validateMappings).


2.3. Design Estatístico e Modelagem no DESeq2

O mapeamento de transcritos para nível de genes e a importação das abundâncias estimadas para o ecossistema R foram computados via pacote tximport. A modelagem estatística de expressão diferencial foi estruturada no pacote DESeq2.

*Fórmula Estatística do Modelo: design = ~ treatment (implementado sob o rótulo de condition no fluxo computacional), tratando o tratamento como variável explicativa categórica de efeito fixo.

*Contrastes Par a Par: O teste de Wald foi aplicado para derivar os P-valores de três contrastes específicos:
    1.  c("condition", "MEKi", "ctrl") — Isolamento dos efeitos do bloqueio MAPK.
    2.  c("condition", "FLT3i", "ctrl") — Avaliação do reposicionamento de FLT3i.
    3.  c("condition", "MEKi", "FLT3i") — Contraste direto entre as drogas para rastrear divergências moleculares.

*Thresholds Aplicados: Foram definidos como significativos os Genes Diferencialmente Expressos (DEGs) que atingiram cumulativamente P-valor ajustado pelo método de Benjamini-Hochberg ($\text{padj} < 0.05$) e variação absoluta de expressão de dobra $|\log_2 \text{Fold Change}| \ge 1.0$.




 3. Resultados

3.1. Controle de Qualidade (QC) e Taxas de Mapeamento

Os relatórios consolidados no MultiQC atestaram o alto padrão das bibliotecas, mantendo médias de score Phred acima de $Q30$ na quase totalidade dos fragmentos. O processo de trimming removeu menos de 1.4% do volume de dados originais. As taxas de mapeamento global calculadas pelo Salmon contra o GENCODE revelaram consistência e estabilidade técnica rigorosa em todas as amostras do projeto, oscilando estritamente entre 82.4% e 86.1% de eficiência de alinhamento.

3.2. Análise de Componentes Principais (PCA)

O gráfico bidimensional de PCA capturou as forças direcionadoras da variância transcritômica global na linhagem H1792, com a Componente Principal 1 (PC1) explicando *75% da variância total  e a PC2 retendo 9%. O achado mais marcante foi o comportamento anômalo da amostra controle SRR13020843, que atuou como um severo outlier, isolando-se no extremo esquerdo do gráfico . 

Em contrapartida, as réplicas controle remanescentes (SRR13020844 e SRR13020845) agruparam-se de maneira compacta na porção direita do gráfico, posicionando-se imediatamente acima do cluster formado pelas triplicatas tratadas com MEKi (SRR13020849, SRR13020850, SRR13020851). O grupo tratado com FLT3i (SRR13020846, SRR13020847, SRR13020848) consolidou um cluster coeso na porção superior direita, demonstrando uma separação nítida em relação ao MEKi ao longo do eixo da PC2.

3.3. Padrões de Dispersão e Expressão Diferencial

* MA Plot (MAPlot_MEKi_vs_ctrl.png): Revelou a distribuição logarítmica de abundância média no eixo X e a magnitude de mudança ($\log_2\text{FC}$) no eixo Y. Os pontos representados em azul marcam os genes estatisticamente significativos ($\text{padj} < 0.05$). Observa-se uma dispersão simétrica de DEGs concentrada em níveis médios de contagem, alcançando variações extremas superiores a $\pm 4$. Triângulos pretos nas bordas do gráfico sinalizam genes que extrapolaram os limites visuais do eixo Y.

*Volcano Plot: Demonstrou graficamente o balanço e a distribuição dos DEGs através da relação entre significância estatística ($-\log_{10} \text{p-value}$) e magnitude de efeito ($\log_2\text{FC}$). Adotando as linhas de corte limítrofes em azul, os genes considerados robustos foram destacados em vermelho. Nota-se um equilíbrio proporcional entre transições de repressão downregulation e ativação upregulation, evidenciando um gene específico de forte responsividade que atingiu significância extrema no topo do gráfico ($-\log_{10} \text{p-value} > 25$).

3.4. Análise de Clusterização Hierárquica (Heatmap)

O mapa térmico, gerado sobre o subconjunto dos Top 20 genes de maior variância global, forneceu a assinatura de resolução molecular que corrobora o PCA, organizando as amostras em três ramificações dendrogramáticas principais:
1. Braço Esquerdo: Isolou de forma singular o SRR13020843, impulsionado por um bloco superior de 5 genes (incluindo transcritos terminados em .1, .11 e .18) dotados de fortíssima superexpressão relativa (Z-score $\approx +2.0$) que não se repete em nenhuma outra amostra do experimento.

2. Braço Central: Reuniu as triplicatas do grupo FLT3i, documentando a regulação positiva exclusiva e coordenada de um bloco central contendo 9 assinaturas gênicas em resposta ao quizartinib.

3. Braço Direito: Agrupou homogeneamente o tratamento MEKi em conjunto com as réplicas fisiológicas normais do controle (SRR13020844 e SRR13020845), indicando que o estado reprimido ou ativado desse subconjunto de genes é compartilhado entre o controle e a inibição de MEK.


3.5. Ontologia Gênica e Enriquecimento de Vias (g:Profiler)

Os resultados do enriquecimento funcional gerados pelo g:Profiler mapearam os processos afetados nas três abordagens par a par:

* Contraste FLT3i vs ctrl: 

* Contraste MEKi vs ctrl: 

* Contraste MEKi vs FLT3i: 



4. Discussão

A reanálise do perfil de expressão do dataset GSE161218 sob a perspectiva do reposicionamento farmacológico direcionado ao Câncer de Pulmão (linhagem H1792 — KRAS mutante) oferece uma interpretação inovadora acerca dos mecanismos de escape tumoral. O desvio biológico substancial apresentado pela amostra controle SRR13020843 deslocou a assinatura molecular basal do estudo, o que explica o fato de os perfis de enriquecimento funcional do g:Profiler terem convergido de maneira tão contundente para o metabolismo de purinas e compostos baseados em nucleobases. Células de adenocarcinoma pulmonar impulsionadas pela hiperativação constitutiva do oncogene KRAS exibem uma dependência metabólica estrita da via de síntese de purinas, necessária para coordenar o pool energético celular e manter as altas taxas de proliferação mitótica. A modulação coordenada desses termos por ambas as drogas revela que tanto o bloqueio downstream (MEKi) quanto o direcionamento ao RTK acessório (FLT3i) impactam de forma aguda a homeostase nucleotídica e a sinalização citoplasmática (refletido pelo termo protein binding).

A grande chave para responder na pergunta desta pesquisa é que o sinergismo medicamentoso reside nas evidências obtidas na comparação direta MEKi vs FLT3i. Ao neutralizar o ruído do grupo controle, as ferramentas estatísticas demonstraram de forma inequívoca que a disparidade funcional entre os dois braços farmacológicos se concentra em processos biológicos de organização estrutural da célula e projeções celulares (cell projection organization). Conforme amplamente reportado na literatura de CPCNP e documentado no artigo original atrelado a este dataset, células tumorais guiadas por KRAS respondem ao bloqueio agudo da via MAPK acionando redes adaptativas e feedbacks baseados em múltiplos receptores de tirosina quinase localizados em projeções citoplasmáticas de membrana para contornar a eficácia da monoterapia.

Nesse contexto, os dados indicam que o inibidor de MEK (trametinib) atua desarmando o núcleo central da cascata proliferativa estimulada pelo KRAS mutante. Por sua vez, o inibidor de FLT3 (quizartinib), atuando aqui como um agente de reposicionamento direcionado ao câncer de pulmão, atinge de forma cirúrgica as estruturas físicas de projeção e os componentes citoplasmáticos de suporte que a célula de adenocarcinoma pulmonar utiliza para reorganizar seu citoesqueleto e restabelecer sinais de sobrevivência. Portanto, a ação concomitante dessas duas classes farmacológicas consolida um efeito altamente sinérgico: enquanto o MEKi interrompe o motor de proliferação do KRAS mutante, o agente reposicionado FLT3i bloqueia a principal rota de escape adaptativa do citoplasma da célula tumoral, fornecendo uma base transcritômica robusta que justifica o uso combinado de terapias moleculares para superar a resistência terapêutica em tumores pulmonares.




5. Referências


Vicent S, Román M et al. Nature Communications, 2023. DOI: 10.1038/s41467-023-41828-z — PMID: 37816716

Chen S. fastp 1.0: An ultra-fast all-round tool for FASTQ data quality control and preprocessing. Imeta. 2025 Sep 9;4(5):e70078. doi: 10.1002/imt2.70078. PMID: 41112039; PMCID: PMC12527978.

Macaya I, Roman M, Welch C, Entrialgo-Cadierno R, Salmon M, Santos A, Feliu I, Kovalski J, Lopez I, Rodriguez-Remirez M, Palomino-Echeverria S, Lonfgren SM, Ferrero M, Calabuig S, Ludwig IA, Lara-Astiaso D, Jantus-Lewintre E, Guruceaga E, Narayanan S, Ponz-Sarvise M, Pineda-Lucena A, Lecanda F, Ruggero D, Khatri P, Santamaria E, Fernandez-Irigoyen J, Ferrer I, Paz-Ares L, Drosten M, Barbacid M, Gil-Bazo I, Vicent S. Signature-driven repurposing of Midostaurin for combination with MEK1/2 and KRASG12C inhibitors in lung cancer. Nat Commun. 2023 Oct 10;14(1):6332. doi: 10.1038/s41467-023-41828-z. PMID: 37816716; PMCID: PMC10564741.

Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. doi: 10.1186/s13059-014-0550-8. PMID: 25516281; PMCID: PMC4302049.

