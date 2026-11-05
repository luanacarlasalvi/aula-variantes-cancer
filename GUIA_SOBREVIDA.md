# Guia da aula — Análise de sobrevida em câncer colorretal

**Bioinformática aplicada ao estudo do câncer · IEC – PUC Minas · Prof. Dr. Fábio Ribeiro Queiroz**

Nesta parte da aula, vamos responder a uma pergunta clássica da oncologia: **quais características do paciente e do tumor estão associadas a viver mais ou menos tempo?** Usaremos uma coorte **fictícia** de 500 pacientes com câncer colorretal, com dados clínicos e mutações nos mesmos 13 genes do painel da aula anterior.

> **Sobre os dados.** Os pacientes são simulados, mas as mutações são reais (coordenadas GRCh38, anotadas no Ensembl VEP) e as frequências e os efeitos foram calibrados pela literatura. Os arquivos estão em `sobrevida/dados/`:
> - `coorte_colorretal_clinico.tsv` — uma linha por paciente: idade, sexo, estágio (I a IV), lado do tumor, status de instabilidade de microssatélites (MSI-H ou MSS), **tempo de seguimento em meses** e **óbito** (1 = faleceu; 0 = vivo no último contato, ou seja, **censurado**).
> - `coorte_colorretal.maf` — as mutações de cada tumor, no mesmo formato MAF que vocês geraram na aula anterior.

**Tempo previsto: 1 hora.** Cada bloco abaixo leva de 3 a 8 minutos.

---

## Antes de começar

**1. Atualize o seu codespace** (se ele é da aula anterior, os arquivos de sobrevida ainda não estão nele). No terminal:

```bash
git pull
ls sobrevida/dados
```

**2. Abra o R** digitando `R` no terminal e teclando Enter. O prompt muda para `>`.

```bash
R
```

**3. Cole um bloco por vez**, na ordem. Os gráficos são salvos em `resultados/sobrevida/` — abra-os clicando no arquivo `.png` no Explorador (à esquerda). Para sair do R ao final: `q()` e responda `n`.

> Copie só o conteúdo dos quadros cinza, nunca as linhas com ```` ``` ````. Se você se perder, saia do R e rode a análise inteira de uma vez: `Rscript sobrevida/analise_sobrevida.R`.

---

## Bloco 1 — Carregar e conhecer a coorte

Carregamos as funções de gráfico da aula e a tabela clínica. Transformamos estágio e MSI em **fatores** (categorias com uma ordem e uma categoria de referência: estágio I e MSS) e criamos `idade10`, a idade dividida por 10, para que o hazard ratio da idade seja lido "a cada 10 anos".

```r
source("sobrevida/funcoes_sobrevida.R")
coorte <- read.delim("sobrevida/dados/coorte_colorretal_clinico.tsv", stringsAsFactors = FALSE)
coorte$estadio    <- factor(coorte$estadio, levels = c("I", "II", "III", "IV"))
coorte$status_msi <- factor(coorte$status_msi, levels = c("MSS", "MSI-H"))
coorte$idade10    <- coorte$idade / 10          # HR por 10 anos de idade
head(coorte)
table(coorte$estadio)
table(coorte$status_msi)
table(coorte$obito)                             # 1 = óbito, 0 = vivo (censurado)
```

**Observe:** quantos pacientes há em cada estágio? Quantos morreram e quantos estão censurados?

---

## Bloco 2 — Quanto tempo acompanhamos esses pacientes?

Antes de olhar a sobrevida, precisamos saber se o seguimento é longo o bastante. A **mediana de seguimento** é calculada com o "Kaplan-Meier invertido": o evento passa a ser a censura.

```r
# Mediana de seguimento pelo método de Kaplan-Meier "invertido":
# o evento passa a ser a censura (quem saiu do estudo vivo).
seguimento <- survfit(Surv(tempo_meses, 1 - obito) ~ 1, data = coorte)
seguimento
```

**Observe:** a mediana (coluna `median`) está em meses. É suficiente para falar de sobrevida em 5 anos?

---

## Bloco 3 — Juntar as mutações do MAF à tabela clínica

Para cada gene, criamos uma coluna com 1 (tumor com mutação não sinônima) ou 0 (selvagem). Depois olhamos como as mutações se distribuem entre os grupos clínicos.

```r
genes <- c("APC", "TP53", "KRAS", "BRAF", "PIK3CA")
mutacoes <- genes_mutados("sobrevida/dados/coorte_colorretal.maf", coorte$Tumor_Sample_Barcode, genes)
coorte <- merge(coorte, mutacoes, by = "Tumor_Sample_Barcode")
colMeans(coorte[, genes])                       # frequência de cada gene na coorte
round(prop.table(table(coorte$estadio, coorte$TP53), margin = 1), 2)   # TP53 por estágio
round(prop.table(table(coorte$status_msi, coorte$BRAF), margin = 1), 2) # BRAF por MSI
```

**Observe:** a frequência de TP53 é igual nos quatro estágios? E o BRAF V600E é igual em MSS e MSI-H? Guarde essas duas tabelas: elas explicam os blocos 6 a 8.

---

## Bloco 4 — Curvas de Kaplan-Meier: a coorte inteira e por estágio

`Surv(tempo_meses, obito)` cria o objeto de sobrevida: o tempo de cada paciente e se ele terminou em óbito (1) ou censura (0). `survfit()` estima a curva de **Kaplan-Meier**; cada degrau é um óbito e cada tracinho (+) é um paciente censurado. O `survdiff()` faz o **teste de log-rank**, que compara as curvas.

```r
km_todos <- grafico_km(Surv(tempo_meses, obito) ~ 1, coorte, "01_km_coorte.png",
                       "Sobrevida global da coorte")
summary(km_todos, times = c(12, 36, 60))        # sobrevida em 1, 3 e 5 anos

km_estadio <- grafico_km(Surv(tempo_meses, obito) ~ estadio, coorte, "02_km_estadio.png",
                         "Sobrevida global por estágio", rotulos = paste("Estágio", levels(coorte$estadio)))
summary(km_estadio, times = 60)                 # sobrevida em 5 anos por estágio
survdiff(Surv(tempo_meses, obito) ~ estadio, data = coorte)   # teste de log-rank
```

**Abra** `01_km_coorte.png` e `02_km_estadio.png`. **Observe:** qual a sobrevida em 5 anos (60 meses) em cada estágio? Compare com os dados do SEER: 91% (doença localizada), 74% (regional) e 13% (à distância). Repare na tabela "Pacientes em risco" embaixo do gráfico: por que ela diminui com o tempo?

---

## Bloco 5 — Kaplan-Meier por característica molecular

Agora comparamos pacientes com e sem cada alteração molecular, **um fator de cada vez**.

```r
grafico_km(Surv(tempo_meses, obito) ~ status_msi, coorte, "03_km_msi.png", "MSI-H x MSS")
grafico_km(Surv(tempo_meses, obito) ~ KRAS, coorte, "04_km_kras.png", "KRAS mutado x selvagem",
           rotulos = c("KRAS selvagem", "KRAS mutado"))
grafico_km(Surv(tempo_meses, obito) ~ TP53, coorte, "05_km_tp53.png", "TP53 mutado x selvagem",
           rotulos = c("TP53 selvagem", "TP53 mutado"))
grafico_km(Surv(tempo_meses, obito) ~ BRAF, coorte, "06_km_braf.png", "BRAF V600E x selvagem",
           rotulos = c("BRAF selvagem", "BRAF V600E"))
```

**Abra** as figuras `03` a `06`. **Observe:** quais diferenças são significativas pelo log-rank? Anote a sua conclusão para MSI-H, KRAS, TP53 e BRAF — vamos revisá-la no bloco 7.

---

## Bloco 6 — Modelo de Cox univariado

O **modelo de riscos proporcionais de Cox** estima o **hazard ratio (HR)**: quantas vezes o risco instantâneo de óbito é maior (HR > 1) ou menor (HR < 1) em um grupo em relação ao grupo de referência. Aqui ajustamos um modelo separado para cada variável.

```r
univariado <- cox_univariado(c("estadio", "idade10", "sexo", "lado", "status_msi",
                               "APC", "TP53", "KRAS", "BRAF", "PIK3CA"), coorte)
univariado
```

**Observe:** as colunas `IC95_inf` e `IC95_sup` são o intervalo de confiança de 95%. Quando ele cruza o 1, o efeito não é significativo. Quais variáveis parecem associadas à sobrevida?

---

## Bloco 7 — Modelo de Cox multivariado

Agora colocamos as variáveis **juntas** no mesmo modelo. Cada HR passa a ser o efeito daquela variável **mantendo as outras constantes** (ajustado). O forest plot compara o HR de cada fator antes (círculo vazio) e depois (círculo cheio) do ajuste.

```r
modelo <- coxph(Surv(tempo_meses, obito) ~ estadio + idade10 + status_msi + BRAF + KRAS + TP53,
                data = coorte)
multivariado <- tabela_hr(modelo)
multivariado
forest_plot(univariado, multivariado,
            termos  = c("TP53", "BRAF", "status_msiMSI-H", "KRAS"),
            rotulos = c("TP53 mutado", "BRAF V600E", "MSI-H", "KRAS mutado"),
            arquivo = "07_forest_uni_x_multi.png", titulo = "Antes e depois do ajuste")
```

**Abra** `07_forest_uni_x_multi.png`. **Observe:** o que aconteceu com TP53? E com BRAF? Por quê? (Dica: volte às duas tabelas do bloco 3.)

---

## Bloco 8 — O que estava escondendo o BRAF?

O BRAF V600E é muito mais comum nos tumores MSI-H, que têm prognóstico melhor. Vamos olhar só os tumores MSS.

```r
table(BRAF = coorte$BRAF, MSI = coorte$status_msi)
mss <- subset(coorte, status_msi == "MSS")
grafico_km(Surv(tempo_meses, obito) ~ BRAF, mss, "08_km_braf_so_mss.png",
           "BRAF V600E entre os tumores MSS", rotulos = c("BRAF selvagem", "BRAF V600E"))
```

**Abra** `08_km_braf_so_mss.png` e compare com `06_km_braf.png`.

---

## Bloco 9 — A mesma pergunta pelo maftools

O maftools junta o MAF e os dados clínicos e faz a análise de sobrevida por gene diretamente. O oncoplot mostra a coorte inteira com o estágio e o MSI anotados embaixo.

```r
suppressPackageStartupMessages(library(maftools))
maf <- read.maf("sobrevida/dados/coorte_colorretal.maf",
                clinicalData = "sobrevida/dados/coorte_colorretal_clinico.tsv", verbose = FALSE)
png(file.path(PASTA_FIG, "09_oncoplot_coorte.png"), width = 1900, height = 1150, res = 150)
oncoplot(maf, top = 10, clinicalFeatures = c("status_msi", "estadio"), sortByAnnotation = TRUE,
         annotationOrder = c("MSS", "MSI-H"))
dev.off()
png(file.path(PASTA_FIG, "10_mafsurvival_braf.png"), width = 1500, height = 1100, res = 150)
mafSurvival(maf, genes = "BRAF", time = "tempo_meses", Status = "obito", isTCGA = FALSE)
dev.off()
cat("Figuras salvas em", PASTA_FIG, "\n")
```

**Abra** `09_oncoplot_coorte.png` e `10_mafsurvival_braf.png`. **Observe:** o barplot de carga mutacional (TMB) no topo — o que distingue os tumores MSI-H? O resultado do `mafSurvival` para BRAF concorda com o Kaplan-Meier univariado ou com o Cox multivariado?

---

## Para fechar a discussão

1. Por que um gene pode parecer prognóstico quando analisado sozinho e deixar de ser depois do ajuste?
2. Por que um gene pode parecer **sem** efeito sozinho e se revelar importante depois do ajuste?
3. O MSI-H continua significativo no modelo multivariado? Isso quer dizer que ele não tem efeito?
4. Quais variáveis você incluiria em um modelo para orientar decisões clínicas? E o que faltou nesta coorte (tratamento, por exemplo)?

## Problemas comuns

| Sintoma | Solução |
|---|---|
| `cannot open file 'sobrevida/...'` | O R não está na pasta do projeto. Saia (`q()`, `n`), rode `cd /workspaces/aula-variantes-cancer` e abra o R de novo |
| `ls: cannot access 'sobrevida/dados'` | Faltou o `git pull` |
| `object 'coorte' not found` (ou `mutacoes`, `modelo`...) | Um bloco anterior não foi rodado. Volte e rode os blocos em ordem |
| `could not find function "grafico_km"` | Faltou a primeira linha do Bloco 1 (`source(...)`) |
| O R mostra `+` e não responde | Ficou um parêntese ou aspa sem fechar. Tecle `Esc` (ou `Ctrl` + `C`) e cole o bloco de novo |
| Quero recomeçar do zero | Saia do R e rode `Rscript sobrevida/analise_sobrevida.R` |
