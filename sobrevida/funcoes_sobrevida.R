# =============================================================================
# Funções auxiliares da aula de análise de sobrevida
#
# Servem apenas para DESENHAR os gráficos (curva de Kaplan-Meier com tabela de
# pacientes em risco e forest plot) e para mostrar as tabelas de forma legível.
# Os cálculos estatísticos ficam no guia: Surv(), survfit(), survdiff(), coxph().
# Usam só o R básico + o pacote survival (que já vem com o maftools).
# =============================================================================
suppressPackageStartupMessages(library(survival))
# Garante acentos corretos nos gráficos, mesmo que o terminal não esteja em UTF-8
invisible(suppressWarnings(
  if (!l10n_info()$`UTF-8`) for (loc in c("C.UTF-8", "en_US.UTF-8", "pt_BR.UTF-8"))
    if (Sys.setlocale("LC_CTYPE", loc) != "") break
))

PASTA_FIG <- "resultados/sobrevida"
dir.create(PASTA_FIG, recursive = TRUE, showWarnings = FALSE)

CORES <- c("#1b6f8a", "#d1495b", "#edae49", "#5b8e3e", "#7b5ea7", "#555555")

# Lê o MAF e cria uma coluna 0/1 por gene (1 = paciente com mutação não sinônima)
genes_mutados <- function(arquivo_maf, ids, genes) {
  maf <- read.delim(arquivo_maf, stringsAsFactors = FALSE)
  nao_sinonimas <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                     "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site")
  maf <- maf[maf$Variant_Classification %in% nao_sinonimas, ]
  out <- data.frame(Tumor_Sample_Barcode = ids)
  for (g in genes) {
    portadores <- unique(maf$Tumor_Sample_Barcode[maf$Hugo_Symbol == g])
    out[[g]] <- as.integer(ids %in% portadores)
  }
  out
}

# Curva de Kaplan-Meier com tabela de pacientes em risco e p do log-rank
grafico_km <- function(formula, dados, arquivo, titulo, rotulos = NULL,
                       tempos = seq(0, 120, by = 24), xlab = "Meses desde o diagnóstico") {
  fit <- survfit(formula, data = dados)
  n_grupos <- if (is.null(fit$strata)) 1 else length(fit$strata)
  if (is.null(rotulos)) rotulos <- if (n_grupos == 1) "Todos" else sub(".*=", "", names(fit$strata))
  cores <- CORES[seq_len(n_grupos)]
  p_txt <- ""
  if (n_grupos > 1) {
    lr <- survdiff(formula, data = dados)
    p <- pchisq(lr$chisq, df = n_grupos - 1, lower.tail = FALSE)
    p_txt <- paste0("log-rank p ", ifelse(p < 0.001, "< 0,001", paste("=", format(round(p, 3), decimal.mark = ","))))
  }
  png(file.path(PASTA_FIG, arquivo), width = 1700, height = 1250, res = 170)
  layout(matrix(1:2, ncol = 1), heights = c(3.2, 1.05 + 0.18 * n_grupos))
  par(mar = c(4.2, 8.5, 3, 1.5), las = 1)
  plot(fit, col = cores, lwd = 2.4, mark.time = TRUE, cex = 0.6, xlim = range(tempos), xaxt = "n",
       xlab = xlab, ylab = "Probabilidade de sobrevida", main = titulo, yaxt = "n", font.main = 1)
  axis(1, at = tempos); axis(2, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), "%"))
  abline(h = seq(0, 1, 0.2), col = "grey92", lty = 1)
  lines(fit, col = cores, lwd = 2.4, mark.time = TRUE, cex = 0.6)
  legend("bottomleft", legend = rotulos, col = cores, lwd = 3, bty = "n", inset = 0.02)
  if (p_txt != "") mtext(p_txt, side = 3, adj = 1, line = 0.2, cex = 0.95)
  # tabela de pacientes em risco
  s <- summary(fit, times = tempos, extend = TRUE)
  par(mar = c(1, 8.5, 1.8, 1.5))
  plot(NA, xlim = range(tempos), ylim = c(0.5, n_grupos + 0.5), axes = FALSE, xlab = "", ylab = "")
  title("Pacientes em risco", adj = 0, font.main = 1, cex.main = 0.95, line = 0.6)
  estratos <- if (n_grupos == 1) rep(1, length(s$time)) else as.integer(s$strata)
  for (g in seq_len(n_grupos)) {
    y <- n_grupos - g + 1
    text(s$time[estratos == g], y, s$n.risk[estratos == g], col = cores[g], cex = 0.95)
    mtext(rotulos[g], side = 2, at = y, las = 1, line = 0.3, col = cores[g], cex = 0.75)
  }
  dev.off()
  cat("Figura salva:", file.path(PASTA_FIG, arquivo), "\n")
  if (p_txt != "") cat(p_txt, "\n")
  invisible(fit)
}

# Tabela de hazard ratios de um modelo de Cox
tabela_hr <- function(modelo) {
  s <- summary(modelo)
  ci <- s$conf.int; co <- s$coefficients
  t <- data.frame(HR = ci[, "exp(coef)"], IC95_inf = ci[, "lower .95"],
                  IC95_sup = ci[, "upper .95"], p = co[, "Pr(>|z|)"],
                  row.names = rownames(co))
  t <- round(t, 3)
  t$p <- ifelse(t$p < 0.001, "<0,001", format(t$p, decimal.mark = ","))
  t
}

# Cox univariado para várias variáveis, uma de cada vez
cox_univariado <- function(variaveis, dados) {
  do.call(rbind, lapply(variaveis, function(v) {
    tabela_hr(coxph(as.formula(paste("Surv(tempo_meses, obito) ~", v)), data = dados))
  }))
}

# Forest plot comparando dois conjuntos de HR (ex.: univariado x multivariado)
forest_plot <- function(uni, multi, termos, rotulos, arquivo, titulo) {
  png(file.path(PASTA_FIG, arquivo), width = 1700, height = 250 + 120 * length(termos), res = 170)
  par(mar = c(4.5, 13, 3, 1.5), las = 1)
  y <- rev(seq_along(termos))
  lim <- range(c(uni[termos, 2:3], multi[termos, 2:3]), 0.2, 5)
  plot(NA, xlim = lim, ylim = c(0.4, length(termos) + 0.6), log = "x", yaxt = "n",
       xlab = "Hazard ratio (escala log)", ylab = "", main = titulo, font.main = 1)
  abline(v = 1, lty = 2, col = "grey40")
  axis(2, at = y, labels = rotulos)
  for (i in seq_along(termos)) {
    tt <- termos[i]
    segments(uni[tt, 2], y[i] + 0.15, uni[tt, 3], y[i] + 0.15, col = CORES[6], lwd = 2)
    points(uni[tt, 1], y[i] + 0.15, pch = 21, bg = "white", col = CORES[6], cex = 1.4)
    segments(multi[tt, 2], y[i] - 0.15, multi[tt, 3], y[i] - 0.15, col = CORES[1], lwd = 2)
    points(multi[tt, 1], y[i] - 0.15, pch = 19, col = CORES[1], cex = 1.4)
  }
  legend("topright", c("Univariado", "Multivariado (ajustado)"), pch = c(21, 19), col = CORES[c(6, 1)],
         pt.bg = "white", bty = "n", cex = 0.9)
  dev.off()
  cat("Figura salva:", file.path(PASTA_FIG, arquivo), "\n")
}

cat("Funções da aula de sobrevida carregadas.\n")
