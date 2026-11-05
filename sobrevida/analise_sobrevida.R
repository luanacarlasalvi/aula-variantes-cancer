# =============================================================================
#  Análise de sobrevida — coorte fictícia de 500 pacientes com câncer colorretal
#  Bioinformática aplicada ao estudo do câncer · IEC – PUC Minas
#
#  Rodar tudo de uma vez (resgate):   Rscript sobrevida/analise_sobrevida.R
#  Ou abrir o R (digite R no terminal) e colar um bloco de cada vez,
#  seguindo o GUIA_SOBREVIDA.md. As figuras vão para resultados/sobrevida/.
# =============================================================================

# ---- Bloco 1: carregar a coorte ----------------------------------------------
source("sobrevida/funcoes_sobrevida.R")
coorte <- read.delim("sobrevida/dados/coorte_colorretal_clinico.tsv", stringsAsFactors = FALSE)
coorte$estadio    <- factor(coorte$estadio, levels = c("I", "II", "III", "IV"))
coorte$status_msi <- factor(coorte$status_msi, levels = c("MSS", "MSI-H"))
coorte$idade10    <- coorte$idade / 10          # HR por 10 anos de idade
head(coorte)
table(coorte$estadio)
table(coorte$status_msi)
table(coorte$obito)                             # 1 = óbito, 0 = vivo (censurado)

# ---- Bloco 2: tempo de seguimento --------------------------------------------
# Mediana de seguimento pelo método de Kaplan-Meier "invertido":
# o evento passa a ser a censura (quem saiu do estudo vivo).
seguimento <- survfit(Surv(tempo_meses, 1 - obito) ~ 1, data = coorte)
seguimento

# ---- Bloco 3: juntar as mutações do MAF ---------------------------------------
genes <- c("APC", "TP53", "KRAS", "BRAF", "PIK3CA")
mutacoes <- genes_mutados("sobrevida/dados/coorte_colorretal.maf", coorte$Tumor_Sample_Barcode, genes)
coorte <- merge(coorte, mutacoes, by = "Tumor_Sample_Barcode")
colMeans(coorte[, genes])                       # frequência de cada gene na coorte
round(prop.table(table(coorte$estadio, coorte$TP53), margin = 1), 2)   # TP53 por estágio
round(prop.table(table(coorte$status_msi, coorte$BRAF), margin = 1), 2) # BRAF por MSI

# ---- Bloco 4: Kaplan-Meier da coorte inteira e por estágio --------------------
km_todos <- grafico_km(Surv(tempo_meses, obito) ~ 1, coorte, "01_km_coorte.png",
                       "Sobrevida global da coorte")
summary(km_todos, times = c(12, 36, 60))        # sobrevida em 1, 3 e 5 anos

km_estadio <- grafico_km(Surv(tempo_meses, obito) ~ estadio, coorte, "02_km_estadio.png",
                         "Sobrevida global por estágio", rotulos = paste("Estágio", levels(coorte$estadio)))
summary(km_estadio, times = 60)                 # sobrevida em 5 anos por estágio
survdiff(Surv(tempo_meses, obito) ~ estadio, data = coorte)   # teste de log-rank

# ---- Bloco 5: Kaplan-Meier por variável molecular ------------------------------
grafico_km(Surv(tempo_meses, obito) ~ status_msi, coorte, "03_km_msi.png", "MSI-H x MSS")
grafico_km(Surv(tempo_meses, obito) ~ KRAS, coorte, "04_km_kras.png", "KRAS mutado x selvagem",
           rotulos = c("KRAS selvagem", "KRAS mutado"))
grafico_km(Surv(tempo_meses, obito) ~ TP53, coorte, "05_km_tp53.png", "TP53 mutado x selvagem",
           rotulos = c("TP53 selvagem", "TP53 mutado"))
grafico_km(Surv(tempo_meses, obito) ~ BRAF, coorte, "06_km_braf.png", "BRAF V600E x selvagem",
           rotulos = c("BRAF selvagem", "BRAF V600E"))

# ---- Bloco 6: Cox univariado (uma variável por vez) ----------------------------
univariado <- cox_univariado(c("estadio", "idade10", "sexo", "lado", "status_msi",
                               "APC", "TP53", "KRAS", "BRAF", "PIK3CA"), coorte)
univariado

# ---- Bloco 7: Cox multivariado (todas juntas) ----------------------------------
modelo <- coxph(Surv(tempo_meses, obito) ~ estadio + idade10 + status_msi + BRAF + KRAS + TP53,
                data = coorte)
multivariado <- tabela_hr(modelo)
multivariado
forest_plot(univariado, multivariado,
            termos  = c("TP53", "BRAF", "status_msiMSI-H", "KRAS"),
            rotulos = c("TP53 mutado", "BRAF V600E", "MSI-H", "KRAS mutado"),
            arquivo = "07_forest_uni_x_multi.png", titulo = "Antes e depois do ajuste")

# ---- Bloco 8: o que estava escondendo o BRAF? -----------------------------------
table(BRAF = coorte$BRAF, MSI = coorte$status_msi)
mss <- subset(coorte, status_msi == "MSS")
grafico_km(Surv(tempo_meses, obito) ~ BRAF, mss, "08_km_braf_so_mss.png",
           "BRAF V600E entre os tumores MSS", rotulos = c("BRAF selvagem", "BRAF V600E"))

# ---- Bloco 9: a mesma pergunta pelo maftools ------------------------------------
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
