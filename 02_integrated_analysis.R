# 広島大学 駐車料金政策SP調査：統合分析コード ----
# 架空データ専用。実測結果ではない。
# RStudioでこのファイルと analysis_long.csv、design.csv のあるフォルダを
# 作業ディレクトリに設定してから Source を押す。
# 初回のみ install.packages(c("survival", "ordinal", "nnet"))
# 元の02〜12の順序・モデル式を保持。Pythonによる生成は別工程。
# results_R 内の既存モデルを一部再利用する元コードの仕様があるため、
# 新規実行は results_R がない専用フォルダで行うこと。
# CLMMの収束診断は未解決の場合がある。推論は暫定扱い。

if (!all(file.exists(c("analysis_long.csv", "design.csv")))) {
  stop("analysis_long.csv と design.csv のある作業フォルダを設定してください。")
}
for (pkg in c("survival", "ordinal", "nnet")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package: ", pkg)
}
library(survival) # strata() と cluster() を最初から使用可能にする
if (dir.exists("results_R") && length(list.files("results_R")) > 0L) {
  stop("既存結果との混在を避けるため、results_R のない専用フォルダで実行してください。")
}


# 02 基礎集計・政策選好MNL・政策受容性CLMMの推定 ----
# 元ファイル: 02_analyze_R.R
# Hiroshima University SP dummy-data demonstration.
# Run in RStudio: setwd("<folder containing this script and CSVs>"); source("02_analyze_R.R")
# These are fictional data; any p-value is NOT empirical evidence.
# Install first (only once): install.packages(c("survival", "ordinal"))

needed <- c("survival", "ordinal")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)>0) {
  stop(paste("必要なRパッケージをインストールしてください: install.packages(c(",
             paste(sprintf('"%s"', missing), collapse=", "), "))"))
}
if (!file.exists("analysis_long.csv")) stop("分析フォルダを作業ディレクトリに設定してください。analysis_long.csv が見つかりません。")

x <- read.csv("analysis_long.csv", fileEncoding="UTF-8-BOM", stringsAsFactors=FALSE)
x$choice_set_id <- factor(x$choice_set_id)
x$respondent_id <- factor(x$respondent_id)
x$info_group <- factor(x$info_group, levels = paste0("I",1:5))
x$alternative <- factor(x$alternative, levels=c("A","B","SQ"))

# Checks before estimation.
stopifnot(nrow(x)==18000L,
          length(unique(x$respondent_id))==1000L,
          length(unique(x$choice_set_id))==6000L,
          all(table(x$choice_set_id)==3L),
          all(tapply(x$chosen, x$choice_set_id, sum)==1L))
if (!dir.exists("results_R")) dir.create("results_R")

cat("====================\n")
cat("架空データの基本集計\n")
cat("====================\n")
cat("政策選択（A, B, SQ）:\n"); print(table(x$alternative[x$chosen==1]))
cat("情報群×回答者の人数:\n"); print(table(x$info_group[!duplicated(x$respondent_id)]))
cat("賛否（1=反対, 2=中立, 3=賛成）:\n"); print(table(x$support, useNA="ifany"))

# MNL expressed as a conditional logit. Choice attributes are alternative specific.
# Info group/main respondent variables are constant within a choice set, so a main
# effect alone CANNOT be identified. We use reform interactions and user interactions.
# Each individual contributes six tasks -> clustered standard errors by respondent.
formula_choice <- chosen ~ price_1000 + gate_ic + gate_lpr +
  safety_basic + safety_plus + parking_live + parking_reserve +
  bus_disc50 + bus_disc150 + bus_plus1 + bus_plus2 +
  is_reform + reform_info2 + reform_info3 + reform_info4 + reform_info5 +
  price_parking_user + parking_live_user + parking_reserve_user +
  bus_disc50_nonuser + bus_disc150_nonuser +
  bus_plus1_nonuser + bus_plus2_nonuser +
  strata(choice_set_id) + cluster(respondent_id)

mnl <- survival::clogit(formula_choice, data=x, method="approximate")
cat("\n====================\n")
cat("条件付きロジット (MNL): クラスタロバストSE\n")
cat("====================\n")
print(summary(mnl))
writeLines(capture.output(summary(mnl)), "results_R/01_mnl_summary.txt", useBytes=TRUE)

# Ordered logit for reform alternatives only, as no SQ support was asked.
sup <- x[x$is_reform==1 & !is.na(x$support), , drop=FALSE]
sup$support_ordered <- ordered(sup$support, levels=c(1,2,3), labels=c("反対","中立","賛成"))
stopifnot(nrow(sup)==12000L)
formula_ordered <- support_ordered ~ price_1000 + gate_ic + gate_lpr +
  safety_basic + safety_plus + parking_live + parking_reserve +
  bus_disc50 + bus_disc150 + bus_plus1 + bus_plus2 +
  info_group + parking_user + price_parking_user +
  parking_live_user + parking_reserve_user +
  bus_disc50_nonuser + bus_disc150_nonuser +
  bus_plus1_nonuser + bus_plus2_nonuser

# Simple model first to verify the ordering and signs; standard errors here ignore
# within-person dependence and should not be the final inferential result.
ordered_simple <- ordinal::clm(formula_ordered, data=sup, link="logit")
cat("\n====================\n")
cat("順序ロジット (CLM): 動作確認用\n")
cat("====================\n")
print(summary(ordered_simple))
writeLines(capture.output(summary(ordered_simple)), "results_R/02_ordered_clm_check.txt", useBytes=TRUE)

# Main ordered model: person-specific random intercept for repeated judgments.
# May take time for 1000 individuals (nAGQ=1 is Laplace approximation).
formula_clmm <- update(formula_ordered, . ~ . + (1 | respondent_id))
cat("\n混合効果順序ロジット (CLMM) の推定を開始します。\n")
ordered_mixed <- tryCatch(
  ordinal::clmm(formula_clmm, data=sup, link="logit", nAGQ=1, Hess=TRUE),
  error=function(e) {message("CLMM推定に失敗: ", conditionMessage(e)); NULL}
)
if (!is.null(ordered_mixed)) {
  print(summary(ordered_mixed))
  writeLines(capture.output(summary(ordered_mixed)), "results_R/03_ordered_clmm_summary.txt", useBytes=TRUE)
}

# Intention categories are not ordinal; only basic cross-tab here.
cat("\n====================\n")
cat("利用意向（1=利用,2=非利用,3=分からない）×現利用者\n")
cat("====================\n")
print(with(sup, table(parking_user, parking_intention)))

write.csv(data.frame(alternative=names(table(x$alternative[x$chosen==1])),
                     count=as.integer(table(x$alternative[x$chosen==1]))),
          "results_R/04_choice_counts.csv", row.names=FALSE, fileEncoding="UTF-8")
write.csv(data.frame(parking_user=rep(c(0,1),each=3),
                     intention=rep(1:3,2),
                     count=as.vector(t(table(sup$parking_user,sup$parking_intention)))),
          "results_R/05_intention_counts.csv", row.names=FALSE, fileEncoding="UTF-8")
cat("\n完了。結果は results_R フォルダに保存しました。\n")


# 統合時の追加：後続処理が必要とする推定済みモデルを明示的に保存。
if (is.null(ordered_mixed)) stop("CLMMの推定に失敗したため、後続処理を停止します。")
saveRDS(mnl, "results_R/mnl_choice.rds")
saveRDS(ordered_mixed, "results_R/ordered_mixed.rds")


# 03 料金別の賛成確率 ----
# 元ファイル: 03_predict_support.R
# ==========================================
# 駐車料金と政策受容性の分析
# CLMMによる予測確率の計算
# ==========================================

# 1. 推定済みのCLMMから係数を取得
b <- ordered_mixed$beta
a <- ordered_mixed$alpha

# 2. 比較する料金と利用者区分を設定
scenario <- expand.grid(
  price_yen = c(700, 1500, 3500),
  parking_user = c(0, 1)
)

# 3. 料金を1000円単位に変換
scenario$price_1000 <- scenario$price_yen / 1000

# 4. 線形予測値を計算
eta <- 
  b["price_1000"] * scenario$price_1000 +
  b["parking_user"] * scenario$parking_user +
  b["price_parking_user"] *
  scenario$price_1000 * scenario$parking_user

# 5. 賛成・中立・反対の確率を計算
scenario$oppose <- plogis(a[1] - eta)

scenario$support <- 1 - plogis(a[2] - eta)

scenario$neutral <- 1 -
  scenario$oppose - scenario$support

# 6. パーセント表示に変換
scenario$oppose <- scenario$oppose * 100
scenario$neutral <- scenario$neutral * 100
scenario$support <- scenario$support * 100

# 7. 利用者区分を表示
scenario$user_type <- ifelse(
  scenario$parking_user == 1,
  "利用者",
  "非利用者"
)

# 8. 計算結果を表示
print(
  scenario[, c(
    "user_type",
    "price_yen",
    "oppose",
    "neutral",
    "support"
  )],
  row.names = FALSE
)

# 9. CSVとして保存
write.csv(
  scenario,
  "results_R/06_support_probabilities.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)# ==========================================
# 駐車料金と賛成確率のグラフ
# ==========================================

# 利用者・非利用者のデータを分ける
nonuser <- scenario[scenario$parking_user == 0, ]
user <- scenario[scenario$parking_user == 1, ]

# グラフをPNGとして保存
png(
  "results_R/07_support_by_price.png",
  width = 1000,
  height = 700
)

# 非利用者
plot(
  nonuser$price_yen,
  nonuser$support,
  type = "b",
  pch = 16,
  lwd = 2,
  ylim = c(0, 40),
  xlab = "Monthly parking fee (JPY)",
  ylab = "Support probability (%)",
  main = "Parking fee and policy support"
)

# 利用者
lines(
  user$price_yen,
  user$support,
  type = "b",
  pch = 17,
  lwd = 2,
  lty = 2
)

# 凡例
legend(
  "topright",
  legend = c("Non-user", "Parking user"),
  pch = c(16, 17),
  lty = c(1, 2),
  bty = "n"
)

# グラフを保存
dev.off()


# 04 改善施策別の賛成確率 ----
# 元ファイル: 04_compare_measures(1).R
# ============================================================
# Parking SP simulation (1,000 respondents)
# Compare improvements at 1,500 JPY/month using fitted CLMM
# Scenario prediction: information group I1 and random effect = 0
# NOTE: This is for simulation practice, not real survey evidence.
# ============================================================

# Run this script from within the parking_sp_simulation_1000 folder.
model_path <- file.path("results_R", "ordered_mixed.rds")
if (!file.exists(model_path)) {
  stop("Model not found: results_R/ordered_mixed.rds. Check getwd() and your saved model.")
}

ordered_mixed <- readRDS(model_path)
b <- ordered_mixed$beta
a <- ordered_mixed$alpha

policies <- data.frame(
  policy = c("No improvement", "Number plate gate", "Basic facility repair",
             "Enhanced safety", "Bus discount 150 JPY", "Bus +2/hour",
             "Live parking information", "Parking reservation"),
  gate_lpr       = c(0, 1, 0, 0, 0, 0, 0, 0),
  safety_basic   = c(0, 0, 1, 0, 0, 0, 0, 0),
  safety_plus    = c(0, 0, 0, 1, 0, 0, 0, 0),
  bus_disc150    = c(0, 0, 0, 0, 1, 0, 0, 0),
  bus_plus2      = c(0, 0, 0, 0, 0, 1, 0, 0),
  parking_live   = c(0, 0, 0, 0, 0, 0, 1, 0),
  parking_reserve= c(0, 0, 0, 0, 0, 0, 0, 1)
)

scenario <- expand.grid(policy_id = seq_len(nrow(policies)),
                        parking_user = c(0L, 1L))
scenario$policy <- policies$policy[scenario$policy_id]
scenario$price_yen <- 1500
scenario$price_1000 <- scenario$price_yen / 1000
scenario$user_type <- ifelse(scenario$parking_user == 1L,
                             "Parking user", "Non-user")

# Build a design matrix matching the fitted model's coefficients.
X <- matrix(0, nrow = nrow(scenario), ncol = length(b),
            dimnames = list(NULL, names(b)))
required <- c("price_1000", "parking_user", "price_parking_user",
              names(policies)[-1], "parking_live_user", "parking_reserve_user",
              "bus_disc150_nonuser", "bus_plus2_nonuser")
missing_vars <- setdiff(required, colnames(X))
if (length(missing_vars) > 0) {
  stop("Coefficient(s) not found in fitted model: ",
       paste(missing_vars, collapse = ", "))
}

X[, "price_1000"] <- scenario$price_1000
X[, "parking_user"] <- scenario$parking_user
X[, "price_parking_user"] <- scenario$price_1000 * scenario$parking_user
for (v in names(policies)[-1]) {
  X[, v] <- policies[[v]][scenario$policy_id]
}
X[, "parking_live_user"] <- X[, "parking_live"] * scenario$parking_user
X[, "parking_reserve_user"] <- X[, "parking_reserve"] * scenario$parking_user
X[, "bus_disc150_nonuser"] <- X[, "bus_disc150"] * (1 - scenario$parking_user)
X[, "bus_plus2_nonuser"] <- X[, "bus_plus2"] * (1 - scenario$parking_user)

# Cumulative logit model probabilities, conditional on random effect = 0.
eta <- as.vector(X %*% b)
scenario$oppose <- 100 * plogis(a[1] - eta)
scenario$support <- 100 * (1 - plogis(a[2] - eta))
scenario$neutral <- 100 - scenario$oppose - scenario$support
baseline <- scenario$support[scenario$policy_id == 1]
scenario$gain_vs_no_improvement_pp <- scenario$support -
  baseline[match(scenario$parking_user, c(0L, 1L))]

result <- scenario[, c("policy", "user_type", "price_yen", "oppose",
                       "neutral", "support", "gain_vs_no_improvement_pp")]
cat("\nPolicy support predictions (CLMM; simulated data):\n")
print(result, row.names = FALSE, digits = 4)

dir.create("results_R", showWarnings = FALSE)
write.csv(result, file = "results_R/08_policy_support.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

# Two bars per policy; English labels avoid font issues on Windows.
nonuser <- scenario$support[scenario$parking_user == 0]
user <- scenario$support[scenario$parking_user == 1]
bar_data <- rbind("Non-user" = nonuser, "Parking user" = user)
png("results_R/09_policy_comparison.png", width = 1400, height = 900)
par(mar = c(5, 16, 4, 2))
barplot(bar_data, beside = TRUE, horiz = TRUE,
        names.arg = policies$policy, las = 1,
        col = c("gray35", "gray75"),
        xlim = c(0, max(bar_data) * 1.2),
        xlab = "Predicted support probability (%)",
        main = "Improvements at 1,500 JPY/month (simulated CLMM)")
legend("topright", legend = rownames(bar_data),
       fill = c("gray35", "gray75"), bty = "n")
dev.off()

cat("\nSaved:\n",
    "results_R/08_policy_support.csv\n",
    "results_R/09_policy_comparison.png\n", sep = "")


# 05 18課題の政策受容性 ----
# 元ファイル: 05_compare_task18(1).R
# ===============================================================
# SP study: actual 18-task A/B packages (simulated data, CLMM)
# Predictions are conditional on random intercept = 0 and I1 info.
# IMPORTANT: CLMM support model was fitted only for A/B; SQ is
# excluded because the questionnaire did not ask SQ support.
# These are support probabilities, NOT MNL choice probabilities.
# ===============================================================

# 0. Check files. Run from the UTPL Lab parking_sp_simulation_1000 folder.
stopifnot(file.exists("design.csv"))
stopifnot(file.exists("results_R/ordered_mixed.rds"))

# 1. Load actual experiment design and fitted mixed ordinal model.
d <- read.csv("design.csv", stringsAsFactors = FALSE)
model <- readRDS("results_R/ordered_mixed.rds")
b <- model$beta
a <- model$alpha

# 2. Restrict to A/B: 18 tasks x 2 alternatives = 36 packages.
d <- d[d$alternative %in% c("A", "B"), ]
d <- d[order(d$task_id, match(d$alternative, c("A", "B"))), ]
stopifnot(nrow(d) == 36L, length(unique(d$task_id)) == 18L)

# 3. Evaluate each A/B package for current nonusers and users.
sc <- rbind(transform(d, parking_user = 0L),
            transform(d, parking_user = 1L))
sc <- sc[order(sc$parking_user, sc$task_id,
               match(sc$alternative, c("A", "B"))), ]
row.names(sc) <- NULL
sc$price_1000 <- sc$price_yen / 1000
sc$user_type <- ifelse(sc$parking_user == 1L, "Parking user", "Non-user")

# 4. Create indicator variables following the original design codes.
sc$gate_ic         <- as.integer(sc$gate_code == 1L)
sc$gate_lpr        <- as.integer(sc$gate_code == 2L)
sc$safety_basic    <- as.integer(sc$safety_code == 1L)
sc$safety_plus     <- as.integer(sc$safety_code == 2L)
sc$parking_live    <- as.integer(sc$parking_service_code == 1L)
sc$parking_reserve <- as.integer(sc$parking_service_code == 2L)
sc$bus_disc50      <- as.integer(sc$bus_discount_code == 1L)
sc$bus_disc150     <- as.integer(sc$bus_discount_code == 2L)
sc$bus_plus1       <- as.integer(sc$bus_frequency_code == 1L)
sc$bus_plus2       <- as.integer(sc$bus_frequency_code == 2L)

sc$price_parking_user   <- sc$price_1000 * sc$parking_user
sc$parking_live_user    <- sc$parking_live * sc$parking_user
sc$parking_reserve_user <- sc$parking_reserve * sc$parking_user
sc$bus_disc50_nonuser   <- sc$bus_disc50 * (1L - sc$parking_user)
sc$bus_disc150_nonuser  <- sc$bus_disc150 * (1L - sc$parking_user)
sc$bus_plus1_nonuser    <- sc$bus_plus1 * (1L - sc$parking_user)
sc$bus_plus2_nonuser    <- sc$bus_plus2 * (1L - sc$parking_user)

# 5. Construct the design matrix in exactly the fitted order.
X <- matrix(0, nrow(sc), length(b), dimnames = list(NULL, names(b)))
vars <- c("price_1000", "gate_ic", "gate_lpr", "safety_basic",
          "safety_plus", "parking_live", "parking_reserve",
          "bus_disc50", "bus_disc150", "bus_plus1", "bus_plus2",
          "parking_user", "price_parking_user", "parking_live_user",
          "parking_reserve_user", "bus_disc50_nonuser",
          "bus_disc150_nonuser", "bus_plus1_nonuser", "bus_plus2_nonuser")
missing_coef <- setdiff(vars, colnames(X))
if (length(missing_coef)) stop("Missing model coefficients: ",
                               paste(missing_coef, collapse = ", "))
X[, vars] <- as.matrix(sc[, vars])
# The info group is I1 (no additional information): info_groupI2..I5 = 0.

# 6. Predicted CLMM probabilities, random intercept fixed at zero.
eta <- drop(X %*% b)
sc$oppose  <- 100 * plogis(a[1] - eta)
sc$support <- 100 * (1 - plogis(a[2] - eta))
sc$neutral <- 100 - sc$oppose - sc$support
stopifnot(all(abs(sc$oppose + sc$neutral + sc$support - 100) < 1e-8))

# 7. Save the 72 package/group predictions.
dir.create("results_R", showWarnings = FALSE)
output_columns <- c("task_id", "block", "alternative", "user_type",
                    "price_yen", "gate_code", "safety_code",
                    "parking_service_code", "bus_discount_code",
                    "bus_frequency_code", "oppose", "neutral", "support")
pred <- sc[, output_columns]
write.csv(pred, "results_R/10_task18_predictions.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# 8. A vs B differences within each task, separately by user group.
AA <- sc[sc$alternative == "A", c("task_id", "parking_user", "support")]
BB <- sc[sc$alternative == "B", c("task_id", "parking_user", "support")]
comparison <- merge(AA, BB, by = c("task_id", "parking_user"),
                    suffixes = c("_A", "_B"))
comparison <- comparison[order(comparison$parking_user, comparison$task_id), ]
comparison$user_type <- ifelse(comparison$parking_user == 1L,
                               "Parking user", "Non-user")
comparison$B_minus_A_pp <- comparison$support_B - comparison$support_A
comparison$higher_support <- ifelse(comparison$B_minus_A_pp > 0, "B",
                                    ifelse(comparison$B_minus_A_pp < 0, "A", "Equal"))
comparison <- comparison[, c("task_id", "user_type", "support_A",
                             "support_B", "B_minus_A_pp", "higher_support")]
write.csv(comparison, "results_R/11_task18_A_vs_B.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# 9. Generate easy-to-read, separate figures for the two groups.
plot_comparison <- function(group_name, filename) {
  z <- comparison[comparison$user_type == group_name, ]
  z <- z[order(z$task_id), ]
  v <- rbind(A = z$support_A, B = z$support_B)
  png(filename, width = 1700, height = 900, res = 150)
  op <- par(mar = c(5.2, 5, 4, 1))
  barplot(v, beside = TRUE, col = c("gray30", "gray75"),
          names.arg = sprintf("%02d", z$task_id),
          ylim = c(0, max(35, max(v) * 1.15)),
          xlab = "Task ID (01-18)", ylab = "Predicted support probability (%)",
          main = paste("18 SP tasks: A vs B -", group_name))
  legend("topright", legend = c("Policy A", "Policy B"),
         fill = c("gray30", "gray75"), bty = "n")
  par(op)
  dev.off()
}
plot_comparison("Non-user", "results_R/12_task18_support_nonusers.png")
plot_comparison("Parking user", "results_R/13_task18_support_users.png")

cat("\nDone: 18 tasks x A/B x 2 user groups =", nrow(pred), "predictions.\n")
cat("Files saved under results_R:\n")
cat(" - 10_task18_predictions.csv\n")
cat(" - 11_task18_A_vs_B.csv\n")
cat(" - 12_task18_support_nonusers.png\n")
cat(" - 13_task18_support_users.png\n\n")
cat("First six A/B comparisons:\n")
print(head(comparison, 6), row.names = FALSE)


# 06 18課題の政策選択確率 ----
# 元ファイル: 06_compare_MNL_task18(1).R
# ================================================================
# Parking SP simulation, 1,000 respondents
# MNL / conditional logit: P(choose A, B, SQ) for all 18 SP tasks.
# Evaluated separately for current parking users / non-users.
# Information scenario: I1 (no additional information).
# Uses the ORIGINAL 18-task experiment design (design.csv).
# Simulated-data exercise: NOT observed behavior or real policy effects.
# ================================================================

# Run from: UTPL Lab / parking_sp_simulation_1000
if (!file.exists("design.csv") || !file.exists("analysis_long.csv")) {
  stop("design.csv or analysis_long.csv not found. Check getwd().")
}
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Please install the survival package: install.packages('survival')")
}
library(survival)   # Makes strata() and cluster() available in the formula.
dir.create("results_R", showWarnings = FALSE)

# 1. Load previously fitted MNL, or fit MNL ONLY (not slow CLMM).
model_path <- "results_R/mnl_choice.rds"
if (file.exists(model_path)) {
  cat("Loading saved MNL: ", model_path, "\n", sep = "")
  mnl <- readRDS(model_path)
} else {
  cat("Saved MNL not found. Fitting MNL ONLY (CLMM is NOT re-estimated).\n")
  x <- read.csv("analysis_long.csv", fileEncoding = "UTF-8-BOM",
                stringsAsFactors = FALSE)
  x$choice_set_id <- factor(x$choice_set_id)
  x$respondent_id <- factor(x$respondent_id)
  stopifnot(nrow(x) == 18000L,
            length(unique(x$choice_set_id)) == 6000L,
            all(table(x$choice_set_id) == 3L),
            all(tapply(x$chosen, x$choice_set_id, sum) == 1L))

  # Same formula used in 02_analyze_R.R.
  formula_choice <- chosen ~ price_1000 + gate_ic + gate_lpr +
    safety_basic + safety_plus + parking_live + parking_reserve +
    bus_disc50 + bus_disc150 + bus_plus1 + bus_plus2 +
    is_reform + reform_info2 + reform_info3 + reform_info4 + reform_info5 +
    price_parking_user + parking_live_user + parking_reserve_user +
    bus_disc50_nonuser + bus_disc150_nonuser +
    bus_plus1_nonuser + bus_plus2_nonuser +
    strata(choice_set_id) + cluster(respondent_id)

  mnl <- survival::clogit(formula_choice, data = x, method = "approximate")
  saveRDS(mnl, model_path)
  cat("MNL saved: ", model_path, "\n", sep = "")
}

# 2. 18 tasks x 3 alternatives (A, B, SQ) x 2 user types = 108 rows.
d <- read.csv("design.csv", fileEncoding = "UTF-8-BOM",
              stringsAsFactors = FALSE)
d$task_id <- as.integer(d$task_id)
d$alternative <- as.character(d$alternative)
stopifnot(nrow(d) == 54L, length(unique(d$task_id)) == 18L,
          all(table(d$task_id) == 3L),
          all(vapply(split(d$alternative, d$task_id),
                     function(z) setequal(z, c("A", "B", "SQ")), logical(1))))
sc <- rbind(transform(d, parking_user = 0L),
            transform(d, parking_user = 1L))
sc <- sc[order(sc$parking_user, sc$task_id,
               match(sc$alternative, c("A", "B", "SQ"))), ]
row.names(sc) <- NULL
sc$user_type <- ifelse(sc$parking_user == 1L, "Parking user", "Non-user")
sc$price_1000 <- sc$price_yen / 1000

# 3. Same indicator coding as original analysis_long.csv.
sc$gate_ic         <- as.integer(sc$gate_code == 1L)
sc$gate_lpr        <- as.integer(sc$gate_code == 2L)
sc$safety_basic    <- as.integer(sc$safety_code == 1L)
sc$safety_plus     <- as.integer(sc$safety_code == 2L)
sc$parking_live    <- as.integer(sc$parking_service_code == 1L)
sc$parking_reserve <- as.integer(sc$parking_service_code == 2L)
sc$bus_disc50      <- as.integer(sc$bus_discount_code == 1L)
sc$bus_disc150     <- as.integer(sc$bus_discount_code == 2L)
sc$bus_plus1       <- as.integer(sc$bus_frequency_code == 1L)
sc$bus_plus2       <- as.integer(sc$bus_frequency_code == 2L)
sc$is_reform       <- as.integer(sc$alternative != "SQ")
sc$price_parking_user   <- sc$price_1000 * sc$parking_user
sc$parking_live_user    <- sc$parking_live * sc$parking_user
sc$parking_reserve_user <- sc$parking_reserve * sc$parking_user
sc$bus_disc50_nonuser   <- sc$bus_disc50 * (1L - sc$parking_user)
sc$bus_disc150_nonuser  <- sc$bus_disc150 * (1L - sc$parking_user)
sc$bus_plus1_nonuser    <- sc$bus_plus1 * (1L - sc$parking_user)
sc$bus_plus2_nonuser    <- sc$bus_plus2 * (1L - sc$parking_user)
# I1 baseline information: reform_info2/3/4/5 are ALL zero.
for (v in paste0("reform_info", 2:5)) sc[[v]] <- 0L

# 4. Reconstruct MNL utility and apply softmax WITHIN each A/B/SQ set.
b <- stats::coef(mnl)
if (anyNA(b)) stop("MNL has missing coefficient(s). Check summary(mnl).")
missing_vars <- setdiff(names(b), names(sc))
if (length(missing_vars)) {
  stop("Predictor(s) missing from scenario: ", paste(missing_vars, collapse = ", "))
}
X <- as.matrix(sc[, names(b), drop = FALSE])
storage.mode(X) <- "double"
sc$utility <- drop(X %*% b)
sc$choice_probability <- NA_real_
for (idx in split(seq_len(nrow(sc)),
                  interaction(sc$parking_user, sc$task_id, drop = TRUE))) {
  v <- sc$utility[idx]
  e <- exp(v - max(v))  # stable softmax
  sc$choice_probability[idx] <- 100 * e / sum(e)
}
stopifnot(all(is.finite(sc$choice_probability)),
          all(abs(vapply(split(sc$choice_probability,
                               interaction(sc$parking_user, sc$task_id,
                                           drop = TRUE)), sum, numeric(1)) - 100) < 1e-8))

# 5. Save 108 predicted choice probabilities, including the SQ alternative.
output_columns <- c("task_id", "block", "alternative", "user_type",
                    "price_yen", "gate_code", "safety_code",
                    "parking_service_code", "bus_discount_code",
                    "bus_frequency_code", "choice_probability")
pred <- sc[, output_columns]
write.csv(pred, "results_R/14_mnl_task18_probabilities.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# 6. One row per task and user type: probabilities add up to 100%.
make_comparison <- function(group_name) {
  z <- sc[sc$user_type == group_name, ]
  z <- z[order(z$task_id, match(z$alternative, c("A", "B", "SQ"))), ]
  P <- matrix(z$choice_probability, ncol = 3, byrow = TRUE,
              dimnames = list(NULL, c("p_A", "p_B", "p_SQ")))
  ans <- data.frame(task_id = sort(unique(z$task_id)),
                    user_type = group_name, P, check.names = FALSE)
  ans$B_minus_A_pp <- ans$p_B - ans$p_A
  ans$most_likely_choice <- c("A", "B", "SQ")[max.col(P, ties.method = "first")]
  ans
}
comparison <- rbind(make_comparison("Non-user"),
                    make_comparison("Parking user"))
write.csv(comparison, "results_R/15_mnl_task18_comparison.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# 7. 18 task-wise grouped bars (A/B/SQ), separately for each group.
make_plot <- function(group_name, path) {
  z <- comparison[comparison$user_type == group_name, ]
  z <- z[order(z$task_id), ]
  P <- rbind(A = z$p_A, B = z$p_B, SQ = z$p_SQ)
  png(path, width = 1800, height = 900, res = 150)
  op <- par(mar = c(5.1, 5, 4, 1))
  barplot(P, beside = TRUE, col = c("gray25", "gray55", "gray85"),
          names.arg = sprintf("%02d", z$task_id), ylim = c(0, 100),
          xlab = "Task ID (01-18)", ylab = "Predicted choice probability (%)",
          main = paste("MNL: Choice A / B / Status quo -", group_name))
  legend("topright", legend = rownames(P),
         fill = c("gray25", "gray55", "gray85"), bty = "n")
  par(op)
  dev.off()
}
make_plot("Non-user", "results_R/16_mnl_choice_nonusers.png")
make_plot("Parking user", "results_R/17_mnl_choice_users.png")

cat("\nMNL 18-task choice probabilities completed (simulated data).\n")
cat("Information group = I1; model-based probabilities (NOT actual observed shares).\n")
cat("Rows:", nrow(pred), "= 18 tasks x 3 alternatives x 2 user types.\n")
cat("Saved files:\n")
cat(" results_R/14_mnl_task18_probabilities.csv\n")
cat(" results_R/15_mnl_task18_comparison.csv\n")
cat(" results_R/16_mnl_choice_nonusers.png\n")
cat(" results_R/17_mnl_choice_users.png\n\n")
cat("First six tasks (Non-user):\n")
print(head(comparison[comparison$user_type == "Non-user", ], 6),
      row.names = FALSE, digits = 4)


# 07 情報提供と政策選好MNL ----
# 元ファイル: 07_information_effects(1).R
# ============================================================================
# 07_information_effects.R
# Hiroshima University parking SP / simulated 1,000-person data
# Information groups I1-I5: descriptive statistics and MNL interactions
# Run from UTPL Lab/parking_sp_simulation_1000 in RStudio:
#   source("07_information_effects.R", encoding = "UTF-8")
# This script DOES NOT re-estimate the slow CLMM.
# All numbers are from hypothetical data; NOT real-world policy effects.
# ============================================================================

# 0. Check environment. Use the folder that contains analysis_long.csv.
if (!file.exists("analysis_long.csv") || !file.exists("design.csv")) {
  stop("analysis_long.csv / design.csv が見つかりません。getwd()を確認してください。")
}
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("survival が必要です: install.packages('survival')")
}
library(survival)  # strata() and cluster() in model formulas
dir.create("results_R", showWarnings = FALSE, recursive = TRUE)

# 1. Load and validate repeated SP responses.
x <- read.csv("analysis_long.csv", fileEncoding = "UTF-8-BOM",
              stringsAsFactors = FALSE)
x$info_group <- factor(x$info_group, levels = paste0("I", 1:5))
x$alternative <- factor(x$alternative, levels = c("A", "B", "SQ"))
x$choice_set_id <- factor(x$choice_set_id)
x$respondent_id <- factor(x$respondent_id)
stopifnot(nrow(x) == 18000L,
          length(unique(x$respondent_id)) == 1000L,
          length(unique(x$choice_set_id)) == 6000L,
          all(table(x$choice_set_id) == 3L),
          all(tapply(x$chosen, x$choice_set_id, sum) == 1L),
          !anyNA(x$info_group))

# 2. Descriptive summaries, using proper denominators.
# Respondent counts: one row per PERSON (not 18 alternative rows).
person <- x[!duplicated(x$respondent_id),
            c("respondent_id", "info_group", "parking_user")]
n_per_group <- as.data.frame(with(person, table(info_group, parking_user)))
names(n_per_group) <- c("info_group", "parking_user", "respondents")
n_per_group$user_type <- ifelse(n_per_group$parking_user == 1L,
                                 "Parking user", "Non-user")
n_per_group <- n_per_group[, c("info_group", "user_type", "respondents")]
write.csv(n_per_group, "results_R/18_info_group_respondents.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# Policy choice: 6 choice sets per person, 1 selected alternative per set.
obs_choice <- with(x[x$chosen == 1L, ], table(info_group, alternative))
choice_counts <- as.data.frame(obs_choice, stringsAsFactors = FALSE)
names(choice_counts) <- c("info_group", "alternative", "count")
choice_totals <- tapply(choice_counts$count, choice_counts$info_group, sum)
choice_counts$share_pct <- 100 * choice_counts$count /
  as.numeric(choice_totals[as.character(choice_counts$info_group)])
write.csv(choice_counts, "results_R/19_info_group_choices.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# Support: A/B only; 12 evaluations per respondent.
s <- x[x$alternative != "SQ" & !is.na(x$support), , drop = FALSE]
s$support_label <- factor(s$support, levels = 1:3,
                          labels = c("Oppose", "Neutral", "Support"))
obs_support <- with(s, table(info_group, support_label))
support_counts <- as.data.frame(obs_support, stringsAsFactors = FALSE)
names(support_counts) <- c("info_group", "support_label", "count")
support_totals <- tapply(support_counts$count, support_counts$info_group, sum)
support_counts$share_pct <- 100 * support_counts$count /
  as.numeric(support_totals[as.character(support_counts$info_group)])
write.csv(support_counts, "results_R/20_info_group_support.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n== I1-I5 descriptive statistics (SIMULATED DATA) ==\n")
cat("Respondents by information group and user status:\n")
print(n_per_group, row.names = FALSE)
cat("\nObserved choice share (%):\n")
print(choice_counts[, c("info_group", "alternative", "share_pct")],
      row.names = FALSE, digits = 4)
cat("\nObserved support share (%):\n")
print(support_counts[, c("info_group", "support_label", "share_pct")],
      row.names = FALSE, digits = 4)

# 3. MNL interaction model.
# Within a 3-option choice set, information-group dummies are constant;
# they cannot be identified as standalone main effects in conditional logit.
# Baseline model already uses reform_info2-5 = I-group x is_reform.
# Add information x price for I2-I5 to test DIFFERENCES IN FEE SENSITIVITY.
# Add I2 x both updated-gate levels and I5 x both safety levels as targeted
# information x policy-attribute interaction demonstrations.
for (g in 2:5) {
  indicator <- as.integer(x$info_group == paste0("I", g))
  x[[paste0("price_info", g)]] <- x$price_1000 * indicator
}
x$gate_ic_info2 <- x$gate_ic * as.integer(x$info_group == "I2")
x$gate_lpr_info2 <- x$gate_lpr * as.integer(x$info_group == "I2")
x$safety_basic_info5 <- x$safety_basic * as.integer(x$info_group == "I5")
x$safety_plus_info5 <- x$safety_plus * as.integer(x$info_group == "I5")

baseline_terms <- c(
  "price_1000", "gate_ic", "gate_lpr", "safety_basic", "safety_plus",
  "parking_live", "parking_reserve", "bus_disc50", "bus_disc150",
  "bus_plus1", "bus_plus2", "is_reform",
  paste0("reform_info", 2:5),
  "price_parking_user", "parking_live_user", "parking_reserve_user",
  "bus_disc50_nonuser", "bus_disc150_nonuser",
  "bus_plus1_nonuser", "bus_plus2_nonuser"
)
extra_terms <- c(paste0("price_info", 2:5),
                 "gate_ic_info2", "gate_lpr_info2",
                 "safety_basic_info5", "safety_plus_info5")
all_terms <- c(baseline_terms, extra_terms)
stopifnot(all(all_terms %in% names(x)))
model_formula <- as.formula(paste(
  "chosen ~", paste(all_terms, collapse = " + "),
  "+ strata(choice_set_id) + cluster(respondent_id)"
))

model_path <- "results_R/mnl_info_interactions_v1.rds"
if (file.exists(model_path)) {
  cat("\nLoading saved info-interaction MNL (no re-estimation).\n")
  model_info <- readRDS(model_path)
} else {
  cat("\nFitting MNL with information-group interactions...\n")
  model_info <- survival::clogit(model_formula, data = x,
                                method = "approximate")
  saveRDS(model_info, model_path)
}
ss <- summary(model_info)
writeLines(capture.output(ss), "results_R/21_mnl_info_interactions.txt",
           useBytes = TRUE)
cf <- ss$coefficients
se_column <- if ("robust se" %in% colnames(cf)) "robust se" else "se(coef)"
p_column <- grep("^Pr\\(", colnames(cf), value = TRUE)[1]
coef_table <- data.frame(
  term = rownames(cf), coefficient = as.numeric(cf[, "coef"]),
  robust_se = as.numeric(cf[, se_column]),
  z = as.numeric(cf[, "z"]), p_value = as.numeric(cf[, p_column]),
  row.names = NULL
)
write.csv(coef_table, "results_R/22_mnl_info_coefficients.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nMNL: information-related coefficients (cluster-robust SE):\n")
info_coefficients <- coef_table[coef_table$term %in%
  c(paste0("reform_info", 2:5), extra_terms), ]
print(info_coefficients, row.names = FALSE, digits = 4)

# 4. Predicted choice probabilities under SAME 18 tasks, separately
# for each information group and parking-user status.
# We change the I-group as a MODEL SCENARIO, not an estimated causal effect.
d <- read.csv("design.csv", fileEncoding = "UTF-8-BOM",
              stringsAsFactors = FALSE)
stopifnot(nrow(d) == 54L, all(table(d$task_id) == 3L))
sc <- do.call(rbind, lapply(paste0("I", 1:5), function(g) {
  do.call(rbind, lapply(0:1, function(u) {
    tmp <- d
    tmp$info_group <- g
    tmp$parking_user <- u
    tmp
  }))
}))
row.names(sc) <- NULL
sc$info_group <- factor(sc$info_group, levels = paste0("I", 1:5))
sc$user_type <- ifelse(sc$parking_user == 1L, "Parking user", "Non-user")
sc$price_1000 <- sc$price_yen / 1000
sc$gate_ic <- as.integer(sc$gate_code == 1L)
sc$gate_lpr <- as.integer(sc$gate_code == 2L)
sc$safety_basic <- as.integer(sc$safety_code == 1L)
sc$safety_plus <- as.integer(sc$safety_code == 2L)
sc$parking_live <- as.integer(sc$parking_service_code == 1L)
sc$parking_reserve <- as.integer(sc$parking_service_code == 2L)
sc$bus_disc50 <- as.integer(sc$bus_discount_code == 1L)
sc$bus_disc150 <- as.integer(sc$bus_discount_code == 2L)
sc$bus_plus1 <- as.integer(sc$bus_frequency_code == 1L)
sc$bus_plus2 <- as.integer(sc$bus_frequency_code == 2L)
sc$is_reform <- as.integer(sc$alternative != "SQ")
sc$price_parking_user <- sc$price_1000 * sc$parking_user
sc$parking_live_user <- sc$parking_live * sc$parking_user
sc$parking_reserve_user <- sc$parking_reserve * sc$parking_user
sc$bus_disc50_nonuser <- sc$bus_disc50 * (1L - sc$parking_user)
sc$bus_disc150_nonuser <- sc$bus_disc150 * (1L - sc$parking_user)
sc$bus_plus1_nonuser <- sc$bus_plus1 * (1L - sc$parking_user)
sc$bus_plus2_nonuser <- sc$bus_plus2 * (1L - sc$parking_user)
for (g in 2:5) {
  indicator <- as.integer(sc$info_group == paste0("I", g))
  sc[[paste0("reform_info", g)]] <- sc$is_reform * indicator
  sc[[paste0("price_info", g)]] <- sc$price_1000 * indicator
}
sc$gate_ic_info2 <- sc$gate_ic * as.integer(sc$info_group == "I2")
sc$gate_lpr_info2 <- sc$gate_lpr * as.integer(sc$info_group == "I2")
sc$safety_basic_info5 <- sc$safety_basic * as.integer(sc$info_group == "I5")
sc$safety_plus_info5 <- sc$safety_plus * as.integer(sc$info_group == "I5")

b <- coef(model_info)
if (anyNA(b)) stop("MNL has NA coefficients. Check 21_mnl_info_interactions.txt")
if (length(setdiff(names(b), names(sc))) > 0L) {
  stop("Missing scenario variables: ",
       paste(setdiff(names(b), names(sc)), collapse = ", "))
}
X <- as.matrix(sc[, names(b), drop = FALSE])
storage.mode(X) <- "double"
sc$utility <- drop(X %*% b)
sc$choice_probability <- NA_real_
for (idx in split(seq_len(nrow(sc)),
                  interaction(sc$info_group, sc$parking_user, sc$task_id,
                              drop = TRUE))) {
  utilities <- sc$utility[idx]
  weights <- exp(utilities - max(utilities))
  sc$choice_probability[idx] <- 100 * weights / sum(weights)
}
stopifnot(nrow(sc) == 540L, !anyNA(sc$choice_probability),
          all(abs(vapply(split(sc$choice_probability,
                                interaction(sc$info_group, sc$parking_user,
                                            sc$task_id, drop = TRUE)),
                         sum, numeric(1)) - 100) < 1e-8))

pred <- sc[, c("info_group", "user_type", "task_id", "alternative",
               "price_yen", "choice_probability")]
write.csv(pred, "results_R/23_mnl_info_18tasks_predictions.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# Mean predicted SQ choice probability over 18 tasks, not observed SQ rate.
sq <- sc[sc$alternative == "SQ", ]
sq_summary <- aggregate(choice_probability ~ info_group + user_type,
                        data = sq, FUN = mean)
names(sq_summary)[3] <- "mean_predicted_SQ_pct"
sq_summary <- sq_summary[order(sq_summary$user_type,
                               sq_summary$info_group), ]
row.names(sq_summary) <- NULL
write.csv(sq_summary, "results_R/24_mnl_info_sq_summary.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nMean predicted SQ choice probability over 18 tasks (%):\n")
print(sq_summary, row.names = FALSE, digits = 4)

# 5. Graph: same 18 tasks, differing I1-I5 information conditions.
# Construct a row-aligned matrix, without relying on aggregate row order.
get_sq <- function(user_label) {
  z <- sq_summary[sq_summary$user_type == user_label, ]
  z$mean_predicted_SQ_pct[match(paste0("I", 1:5), as.character(z$info_group))]
}
mat <- rbind("Non-user" = get_sq("Non-user"),
             "Parking user" = get_sq("Parking user"))
colnames(mat) <- paste0("I", 1:5)
png("results_R/25_mnl_info_sq_comparison.png", width = 1300,
    height = 820, res = 140)
barplot(mat, beside = TRUE, col = c("gray30", "gray75"),
        ylim = c(0, max(mat) * 1.3),
        main = "Predicted status-quo choice by information group (simulated)",
        xlab = "Information group", ylab = "Predicted SQ choice (%)")
legend("topright", legend = rownames(mat),
       fill = c("gray30", "gray75"), bty = "n")
dev.off()

cat("\nDONE: Information-group descriptive stats + MNL interactions.\n")
cat("Outputs under results_R: 18_info_group_respondents.csv through ")
cat("25_mnl_info_sq_comparison.png\n")
cat("Saved model: results_R/mnl_info_interactions_v1.rds\n")
cat("NOTE: This is a simulation. Differences are not demonstrated causal effects.\n")
cat("NOTE: CLMM support-by-information interactions are a SEPARATE, later step.\n")


# 08 情報提供と政策受容性CLMM ----
# 元ファイル: 08_clmm_information_effects(1).R
# ============================================================================
# 08_clmm_information_effects.R
# Hiroshima University parking SP / 1,000-person SIMULATED sample
# Information I1-I5 and policy acceptance (A and B ratings)
# Run from folder with analysis_long.csv and results_R/ordered_mixed.rds:
#   source("08_clmm_information_effects.R", encoding="UTF-8")
# This script fits one NEW interaction CLMM (or loads its saved RDS).
# It does not re-estimate the earlier baseline model.
# All significance results are from GENERATED data, NOT evidence about people.
# ============================================================================

if (!requireNamespace("ordinal", quietly = TRUE)) {
  stop("R package ordinal is required: install.packages('ordinal')")
}
if (!file.exists("analysis_long.csv")) {
  stop("analysis_long.csv not found. Use getwd() to check the working folder.")
}
if (!file.exists("results_R/ordered_mixed.rds")) {
  stop("Baseline CLMM missing: results_R/ordered_mixed.rds")
}
dir.create("results_R", recursive = TRUE, showWarnings = FALSE)

x <- read.csv("analysis_long.csv", fileEncoding = "UTF-8-BOM",
              stringsAsFactors = FALSE)
sup <- x[x$is_reform == 1L & !is.na(x$support), , drop = FALSE]
sup$respondent_id <- factor(sup$respondent_id)
sup$info_group <- factor(sup$info_group, levels = paste0("I", 1:5))
sup$support_ordered <- ordered(sup$support, levels = 1:3,
                               labels = c("Oppose", "Neutral", "Support"))
stopifnot(nrow(sup) == 12000L, length(unique(sup$respondent_id)) == 1000L,
          all(table(sup$respondent_id) == 12L), !anyNA(sup$info_group))

# Interaction questions were prespecified in the preceding MNL exercise:
#   all I2-I5 x parking fee (whether fee sensitivity differs by information),
#   I2 x updated-gate features, I5 x safety improvements.
for (g in 2:5) {
  sup[[paste0("price_info", g)]] <-
    sup$price_1000 * as.integer(sup$info_group == paste0("I", g))
}
sup$gate_ic_info2 <- sup$gate_ic * as.integer(sup$info_group == "I2")
sup$gate_lpr_info2 <- sup$gate_lpr * as.integer(sup$info_group == "I2")
sup$safety_basic_info5 <- sup$safety_basic * as.integer(sup$info_group == "I5")
sup$safety_plus_info5 <- sup$safety_plus * as.integer(sup$info_group == "I5")

# Match 02_analyze_R.R's original CLMM fixed-effect specification.
fixed_terms <- c(
  "price_1000", "gate_ic", "gate_lpr", "safety_basic", "safety_plus",
  "parking_live", "parking_reserve", "bus_disc50", "bus_disc150",
  "bus_plus1", "bus_plus2", "info_group", "parking_user",
  "price_parking_user", "parking_live_user", "parking_reserve_user",
  "bus_disc50_nonuser", "bus_disc150_nonuser",
  "bus_plus1_nonuser", "bus_plus2_nonuser"
)
interactions <- c(paste0("price_info", 2:5), "gate_ic_info2",
                  "gate_lpr_info2", "safety_basic_info5", "safety_plus_info5")
stopifnot(all(c(fixed_terms, interactions) %in% names(sup)))
formula_full <- as.formula(paste(
  "support_ordered ~", paste(c(fixed_terms, interactions), collapse = " + "),
  "+ (1 | respondent_id)"
))

model_file <- "results_R/clmm_info_interactions_v1.rds"
if (file.exists(model_file)) {
  cat("Loading previously fitted information-interaction CLMM...\n")
  fit <- readRDS(model_file)
} else {
  cat("Fitting information-interaction CLMM (not the original CLMM)...\n")
  fit <- ordinal::clmm(formula_full, data = sup,
                       link = "logit", nAGQ = 1, Hess = TRUE)
  saveRDS(fit, model_file)
  cat("Saved new model: ", model_file, "\n", sep = "")
}
writeLines(capture.output(summary(fit)),
           "results_R/26_clmm_info_model_summary.txt", useBytes = TRUE)

# Coefficients, Wald SE and false-discovery-rate adjustment for the 8
# interaction contrasts. Adjustment accounts for the exploratory 8 tests.
beta <- fit$beta
V <- vcov(fit)
missing_in_vcov <- setdiff(names(beta), rownames(V))
if (length(missing_in_vcov)) {
  stop("Could not locate these fitted coefficients in vcov(): ",
       paste(missing_in_vcov, collapse = ", "))
}
se <- sqrt(diag(V)[names(beta)])
if (anyNA(se) || any(!is.finite(se))) {
  stop("Non-finite SE detected; inspect results_R/26_clmm_info_model_summary.txt")
}
coef_table <- data.frame(term = names(beta), estimate = unname(beta),
                         std_error = unname(se),
                         z = unname(beta / se),
                         p_value = unname(2 * pnorm(-abs(beta / se))),
                         stringsAsFactors = FALSE)
coef_table$p_BH_8_interactions <- NA_real_
idx <- which(coef_table$term %in% interactions)
stopifnot(length(idx) == 8L)
coef_table$p_BH_8_interactions[idx] <- p.adjust(coef_table$p_value[idx],
                                                method = "BH")
write.csv(coef_table, "results_R/27_clmm_info_coefficients.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n== CLMM interactions (SIMULATED DATA) ==\n")
print(coef_table[idx, ], row.names = FALSE, digits = 4)

# Likelihood ratio test of all 8 interaction coefficients as a pre-defined block.
# If optimizer gradients are not sufficiently small, treat p-values cautiously.
base <- readRDS("results_R/ordered_mixed.rds")
ll_old <- as.numeric(logLik(base))
ll_new <- as.numeric(logLik(fit))
LR <- max(0, 2 * (ll_new - ll_old))
df_LR <- length(interactions)
p_LR <- pchisq(LR, df = df_LR, lower.tail = FALSE)
max_gradient <- if (!is.null(fit$gradient)) max(abs(fit$gradient)) else NA_real_
comparison <- data.frame(baseline_logLik = ll_old,
                         interaction_logLik = ll_new,
                         LR_statistic = LR, df = df_LR, p_value = p_LR,
                         max_abs_gradient = max_gradient)
write.csv(comparison, "results_R/28_clmm_info_omnibus_test.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nJoint test of 8 added terms: LR=", round(LR, 3),
    ", df=", df_LR, ", p=", signif(p_LR, 4), "\n", sep = "")
cat("Max absolute optimizer gradient: ", max_gradient, "\n", sep = "")
if (is.finite(max_gradient) && max_gradient > 0.001) {
  warning("CLMM optimizer gradient exceeds 0.001. Treat inference cautiously.")
}

# Predicted support for no-other-improvements policies at 700/1500/3500 JPY.
# Information groups and current use differ. Person random intercept set to 0;
# these are conditional predictions, NOT population-averaged probabilities.
sce <- expand.grid(info_group = paste0("I", 1:5),
                   parking_user = c(0L, 1L),
                   price_yen = c(700, 1500, 3500),
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
sce$user_type <- ifelse(sce$parking_user == 1L,
                        "Parking user", "Non-user")
sce$price_1000 <- sce$price_yen / 1000
X <- matrix(0, nrow = nrow(sce), ncol = length(beta),
            dimnames = list(NULL, names(beta)))
required <- c("price_1000", "parking_user", "price_parking_user",
              paste0("info_groupI", 2:5), interactions)
if (!all(required %in% colnames(X))) {
  stop("Unexpected CLMM coefficient names: ",
       paste(setdiff(required, colnames(X)), collapse = ", "))
}
X[, "price_1000"] <- sce$price_1000
X[, "parking_user"] <- sce$parking_user
X[, "price_parking_user"] <- sce$price_1000 * sce$parking_user
for (g in 2:5) {
  flag <- as.integer(sce$info_group == paste0("I", g))
  X[, paste0("info_groupI", g)] <- flag
  X[, paste0("price_info", g)] <- sce$price_1000 * flag
}
# Gate and safety improvement x information interactions = 0 in these scenarios.
eta <- drop(X %*% beta)
alpha <- fit$alpha
stopifnot(length(alpha) == 2L)
sce$oppose_pct <- 100 * plogis(alpha[1] - eta)
sce$support_pct <- 100 * (1 - plogis(alpha[2] - eta))
sce$neutral_pct <- 100 - sce$oppose_pct - sce$support_pct
stopifnot(all(abs(sce$oppose_pct + sce$neutral_pct +
                  sce$support_pct - 100) < 1e-9))
comparison_ids <- paste(sce$parking_user, sce$price_yen)
ref <- sce[sce$info_group == "I1", c("parking_user", "price_yen", "support_pct")]
sce$vs_I1_pp <- sce$support_pct - ref$support_pct[
  match(comparison_ids, paste(ref$parking_user, ref$price_yen))
]
write.csv(sce[, c("info_group", "user_type", "price_yen", "oppose_pct",
                  "neutral_pct", "support_pct", "vs_I1_pp")],
          "results_R/29_clmm_info_support_by_price.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# Two unambiguous plots: same policy except information condition and fee.
for (u in c(0L, 1L)) {
  dat <- sce[sce$parking_user == u, ]
  mat <- sapply(paste0("I", 1:5), function(g) {
    tmp <- dat[dat$info_group == g, ]
    tmp$support_pct[match(c(700, 1500, 3500), tmp$price_yen)]
  })
  png_file <- if (u == 0L)
    "results_R/30_clmm_info_nonusers.png" else
    "results_R/31_clmm_info_users.png"
  png(png_file, width = 1350, height = 850, res = 135)
  matplot(c(700, 1500, 3500), mat, type = "b", lty = 1:5, pch = 1:5,
          col = 1:5, lwd = 2, xlab = "Monthly parking fee (JPY)",
          ylab = "Predicted support probability (%)",
          ylim = c(0, max(mat) * 1.18), xaxt = "n",
          main = paste0("CLMM: Fee x information (",
                        if (u == 0L) "Non-user" else "Parking user",
                        "; simulated data)"))
  axis(1, at = c(700, 1500, 3500))
  legend("topright", legend = paste0("I", 1:5), col = 1:5,
         lty = 1:5, pch = 1:5, bty = "n")
  dev.off()
}

cat("\nCLMM information-group support predictions (at 1500 JPY):\n")
print(sce[sce$price_yen == 1500,
          c("info_group", "user_type", "support_pct", "vs_I1_pp")],
      row.names = FALSE, digits = 4)
cat("\nDONE. Files in results_R:\n")
cat(" 26_clmm_info_model_summary.txt\n",
    " 27_clmm_info_coefficients.csv\n",
    " 28_clmm_info_omnibus_test.csv\n",
    " 29_clmm_info_support_by_price.csv\n",
    " 30_clmm_info_nonusers.png\n",
    " 31_clmm_info_users.png\n",
    " clmm_info_interactions_v1.rds\n", sep = "")
cat("NOTE: All results are model outputs from SIMULATED respondents.\n")


# 09 駐車場利用意向MNL ----
# 元ファイル: 09_parking_intention_mnl(1).R
# ============================================================================
# 09_parking_intention_mnl.R
# Hiroshima University parking SP - SIMULATED 1,000-respondent dataset
# Parking-use intention analysis: Use / Do not use / Unsure (MNL)
# Run in RStudio from parking_sp_simulation_1000:
#   source("09_parking_intention_mnl.R", encoding = "UTF-8")
#
# IMPORTANT:
#   - Intention data exist only for A/B (not status quo, SQ).
#   - Each respondent supplies 12 intention responses (6 tasks x A/B).
#   - nnet::multinom gives a pooled exploratory MNL fit. The repeated answers
#     are not independent; naive standard errors must NOT be used as formal
#     population inference. Cluster-robust/participant-level analysis is a
#     separate validation step.
#   - ALL DATA ARE SYNTHETIC. Predictions are not actual survey evidence.
# ============================================================================

if (!requireNamespace("nnet", quietly = TRUE)) {
  stop("R package 'nnet' is needed. Run: install.packages('nnet')")
}
if (!file.exists("analysis_long.csv") || !file.exists("design.csv")) {
  stop("analysis_long.csv or design.csv missing. Check getwd() in RStudio.")
}
dir.create("results_R", showWarnings = FALSE, recursive = TRUE)

# 1. Input checks and intention category definitions.
x <- read.csv("analysis_long.csv", fileEncoding = "UTF-8-BOM",
              stringsAsFactors = FALSE)
x <- x[x$alternative %in% c("A", "B"), , drop = FALSE]
needed <- c("respondent_id", "alternative", "task_id", "parking_intention",
            "parking_user", "info_group", "price_yen", "price_1000",
            "price_parking_user", "gate_ic", "gate_lpr", "safety_basic",
            "safety_plus", "parking_live", "parking_reserve", "bus_disc50",
            "bus_disc150", "bus_plus1", "bus_plus2", "parking_live_user",
            "parking_reserve_user", "bus_disc50_nonuser", "bus_disc150_nonuser",
            "bus_plus1_nonuser", "bus_plus2_nonuser")
if (!all(needed %in% names(x))) {
  stop("Missing columns: ", paste(setdiff(needed, names(x)), collapse = ", "))
}
stopifnot(nrow(x) == 12000L, length(unique(x$respondent_id)) == 1000L,
          !anyNA(x$parking_intention), all(x$parking_intention %in% 1:3),
          all(table(x$respondent_id) == 12L),
          all(x$parking_user %in% 0:1))

# Codes in original simulated SP files: 1=Use, 2=Do not use, 3=Unsure.
# Choose "Do not use" as the reference outcome of multinomial logit.
x$intention_category <- factor(x$parking_intention,
  levels = c(2, 1, 3), labels = c("Do not use", "Use", "Unsure"))
x$info_group <- factor(x$info_group, levels = paste0("I", 1:5))
x$respondent_id <- factor(x$respondent_id)
x$user_type <- ifelse(x$parking_user == 1L, "Parking user", "Non-user")

# 2. Observed DESCRIPTIVE rates, including proper within-group denominator.
obs <- as.data.frame(with(x, table(user_type, intention_category)),
                     stringsAsFactors = FALSE)
names(obs) <- c("user_type", "intention_category", "n_responses")
obs$share_pct <- ave(obs$n_responses, obs$user_type,
                     FUN = function(z) 100 * z / sum(z))
write.csv(obs, "results_R/32_intention_observed.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n== Parking intention observed shares (SIMULATED DATA) ==\n")
cat("A/B ONLY: 12,000 intention answers from 1,000 respondents.\n")
print(obs, row.names = FALSE, digits = 5)

# 3. Pooled multinomial-logit model; reuse fitted model if it already exists.
# Price x current-use and selected parking/bus interactions allow heterogeneous
# associations by current parking use. Information group is included as a
# control, not interpreted as a demonstrated causal treatment effect.
terms <- c(
  "price_1000", "parking_user", "price_parking_user",
  "gate_ic", "gate_lpr", "safety_basic", "safety_plus",
  "parking_live", "parking_reserve", "bus_disc50", "bus_disc150",
  "bus_plus1", "bus_plus2", "parking_live_user", "parking_reserve_user",
  "bus_disc50_nonuser", "bus_disc150_nonuser",
  "bus_plus1_nonuser", "bus_plus2_nonuser", "info_group"
)
fit_formula <- as.formula(paste("intention_category ~",
                                paste(terms, collapse = " + ")))
model_path <- "results_R/parking_intention_mnl_v1.rds"
if (file.exists(model_path)) {
  cat("\nLoading saved parking-intention MNL (no re-estimation)...\n")
  fit <- readRDS(model_path)
} else {
  cat("\nFitting parking-intention MNL (nnet::multinom)...\n")
  fit <- nnet::multinom(fit_formula, data = x,
                       trace = FALSE, maxit = 500, Hess = TRUE)
  saveRDS(fit, model_path)
  cat("Saved model: ", model_path, "\n", sep = "")
}
if (!is.null(fit$convergence) && fit$convergence != 0) {
  warning("MNL optimizer convergence != 0. Coefficients may be unreliable.")
}
writeLines(capture.output(summary(fit)),
           "results_R/33_intention_model_summary.txt", useBytes = TRUE)

# Estimate and ordinary SE for exploration only. NOT cluster-robust.
co <- coef(fit)
if (is.null(dim(co))) co <- matrix(co, nrow = 1L)
s <- summary(fit)$standard.errors
if (is.null(dim(s))) s <- matrix(s, nrow = 1L)
coef_table <- data.frame(
  outcome_vs_reference = rep(rownames(co), each = ncol(co)),
  term = rep(colnames(co), times = nrow(co)),
  coefficient = as.vector(t(co)),
  naive_se_not_cluster_robust = as.vector(t(s)),
  stringsAsFactors = FALSE
)
write.csv(coef_table, "results_R/34_intention_coefficients.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# 4. Construct data frames compatible with the original attributes.
build_attributes <- function(z) {
  z$info_group <- factor(z$info_group, levels = paste0("I", 1:5))
  z$price_1000 <- z$price_yen / 1000
  z$gate_ic <- as.integer(z$gate_code == 1L)
  z$gate_lpr <- as.integer(z$gate_code == 2L)
  z$safety_basic <- as.integer(z$safety_code == 1L)
  z$safety_plus <- as.integer(z$safety_code == 2L)
  z$parking_live <- as.integer(z$parking_service_code == 1L)
  z$parking_reserve <- as.integer(z$parking_service_code == 2L)
  z$bus_disc50 <- as.integer(z$bus_discount_code == 1L)
  z$bus_disc150 <- as.integer(z$bus_discount_code == 2L)
  z$bus_plus1 <- as.integer(z$bus_frequency_code == 1L)
  z$bus_plus2 <- as.integer(z$bus_frequency_code == 2L)
  z$price_parking_user <- z$price_1000 * z$parking_user
  z$parking_live_user <- z$parking_live * z$parking_user
  z$parking_reserve_user <- z$parking_reserve * z$parking_user
  z$bus_disc50_nonuser <- z$bus_disc50 * (1L - z$parking_user)
  z$bus_disc150_nonuser <- z$bus_disc150 * (1L - z$parking_user)
  z$bus_plus1_nonuser <- z$bus_plus1 * (1L - z$parking_user)
  z$bus_plus2_nonuser <- z$bus_plus2 * (1L - z$parking_user)
  z$user_type <- ifelse(z$parking_user == 1L, "Parking user", "Non-user")
  z
}

predict_3 <- function(z) {
  p <- as.matrix(predict(fit, newdata = z, type = "probs"))
  expected <- c("Do not use", "Use", "Unsure")
  if (!all(expected %in% colnames(p))) {
    stop("Unexpected multinom outcome names: ",
         paste(colnames(p), collapse = ", "))
  }
  stopifnot(all(is.finite(p)), all(abs(rowSums(p) - 1) < 1e-6))
  z$p_do_not_use_pct <- 100 * p[, "Do not use"]
  z$p_use_pct <- 100 * p[, "Use"]
  z$p_unsure_pct <- 100 * p[, "Unsure"]
  z
}

# 5. Predict use intentions at 700 / 1500 / 3500 JPY with no improvements,
#    holding all other policy attributes at the status quo, info group I1.
price_only <- expand.grid(price_yen = c(700L, 1500L, 3500L),
                          parking_user = c(0L, 1L),
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
price_only$info_group <- "I1"
for (a in c("gate_code", "safety_code", "parking_service_code",
            "bus_discount_code", "bus_frequency_code")) price_only[[a]] <- 0L
price_pred <- predict_3(build_attributes(price_only))
cols_price <- c("user_type", "parking_user", "info_group", "price_yen",
                "p_use_pct", "p_do_not_use_pct", "p_unsure_pct")
write.csv(price_pred[, cols_price],
          "results_R/35_intention_by_price.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
cat("\n== Predicted parking-use intention by monthly fee ==\n")
cat("No improvements; I1; model predictions from simulated responses.\n")
print(price_pred[, c("user_type", "price_yen", "p_use_pct",
                      "p_do_not_use_pct", "p_unsure_pct")],
      row.names = FALSE, digits = 4)

# 6. Predict intention for the ORIGINAL 18 survey tasks, policy A and B.
d <- read.csv("design.csv", fileEncoding = "UTF-8-BOM",
              stringsAsFactors = FALSE)
d <- d[d$alternative %in% c("A", "B"), , drop = FALSE]
stopifnot(nrow(d) == 36L, all(table(d$task_id) == 2L))
sc <- do.call(rbind, lapply(0:1, function(u) {
  q <- d
  q$parking_user <- u
  q$info_group <- "I1"
  q
}))
row.names(sc) <- NULL
task_pred <- predict_3(build_attributes(sc))
stopifnot(nrow(task_pred) == 72L)
cols_tasks <- c("task_id", "alternative", "user_type", "parking_user",
                "info_group", "price_yen", "gate_code", "safety_code",
                "parking_service_code", "bus_discount_code",
                "bus_frequency_code", "p_use_pct", "p_do_not_use_pct",
                "p_unsure_pct")
write.csv(task_pred[, cols_tasks],
          "results_R/36_intention_task18_predictions.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# 7. Two graphs comparing predicted use intention by monthly fee and user type.
for (u in c(0L, 1L)) {
  pp <- price_pred[price_pred$parking_user == u, ]
  pp <- pp[order(pp$price_yen), ]
  png_path <- if (u == 0L)
    "results_R/37_intention_price_nonusers.png" else
    "results_R/38_intention_price_users.png"
  png(png_path, width = 1250, height = 800, res = 130)
  y <- rbind(pp$p_use_pct, pp$p_do_not_use_pct, pp$p_unsure_pct)
  rownames(y) <- c("Use", "Do not use", "Unsure")
  barplot(y, beside = TRUE,
          col = c("#3175A3", "#DB9E4B", "#8A95A0"),
          names.arg = paste0(pp$price_yen, " JPY"),
          ylim = c(0, 100),
          ylab = "Predicted probability (%)",
          main = paste0("Parking intention by fee - ",
                        if (u == 0L) "Non-user" else "Parking user",
                        " (simulated)"),
          legend.text = rownames(y),
          args.legend = list(x = "topright", bty = "n"))
  dev.off()
}

cat("\nDONE: Parking-use intention MNL (SIMULATED DATA).\n")
cat("Outcomes: 1=Use; 2=Do not use; 3=Unsure. A/B only (no SQ).\n")
cat("Rows predicted for the 18 tasks: ", nrow(task_pred),
    " = 18 tasks x A/B x 2 user groups.\n", sep = "")
cat("Files in results_R:\n")
cat(" 32_intention_observed.csv\n",
    " 33_intention_model_summary.txt\n",
    " 34_intention_coefficients.csv\n",
    " 35_intention_by_price.csv\n",
    " 36_intention_task18_predictions.csv\n",
    " 37_intention_price_nonusers.png\n",
    " 38_intention_price_users.png\n",
    " parking_intention_mnl_v1.rds\n", sep = "")
cat("CAUTION: Standard errors are NOT respondent-cluster-adjusted.\n")
cat("CAUTION: All results come from SYNTHETIC data, not real participants.\n")


# 10 データ・収束診断とクラスタ標準誤差 ----
# 元ファイル: 10_validate_models(1).R
# ============================================================================
# 10_validate_models.R
# Hiroshima University parking SP: checks for SYNTHETIC-data analysis.
# Loads saved model objects, no model is re-estimated.
# Validates source structure, documents CLMM optimizer gradients, computes
# respondent-clustered sandwich SE for nnet::multinom parking intention MNL.
# Run: source("10_validate_models.R", encoding = "UTF-8")
# ============================================================================

if (!file.exists("analysis_long.csv")) {
  stop("analysis_long.csv not found. Verify getwd() is the UTPL Lab analysis directory.")
}
dir.create("results_R", showWarnings = FALSE, recursive = TRUE)
model_path <- "results_R/parking_intention_mnl_v1.rds"
if (!file.exists(model_path)) {
  stop("Missing saved intention model: ", model_path,
       ". First run 09_parking_intention_mnl.R.")
}

x <- read.csv("analysis_long.csv", fileEncoding = "UTF-8-BOM",
              stringsAsFactors = FALSE)
sp <- x[x$alternative %in% c("A", "B"), , drop = FALSE]
sp$intention_category <- factor(sp$parking_intention,
  levels = c(2L, 1L, 3L), labels = c("Do not use", "Use", "Unsure"))
sp$info_group <- factor(sp$info_group, levels = paste0("I", 1:5))
sp$respondent_id <- factor(sp$respondent_id)

checks <- data.frame(
  item = c("all_alternative_rows", "choice_sets", "respondents",
           "reform_attitude_intention_rows", "intention_responses_per_person",
           "three_alternatives_per_choice_set", "one_chosen_per_choice_set",
           "no_missing_intention", "intention_codes_1_2_3"),
  value = c(nrow(x), length(unique(x$choice_set_id)),
            length(unique(x$respondent_id)), nrow(sp),
            if (all(table(sp$respondent_id) == 12L)) "YES" else "NO",
            if (all(table(x$choice_set_id) == 3L)) "YES" else "NO",
            if (all(tapply(x$chosen, x$choice_set_id, sum) == 1L)) "YES" else "NO",
            if (!anyNA(sp$parking_intention)) "YES" else "NO",
            if (all(sp$parking_intention %in% 1:3)) "YES" else "NO"),
  stringsAsFactors = FALSE)
write.csv(checks, "results_R/39_data_validation.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
cat("\n== Data structure checks ==\n")
print(checks, row.names = FALSE)
if (nrow(x) != 18000L || nrow(sp) != 12000L ||
    length(unique(x$respondent_id)) != 1000L ||
    !all(table(sp$respondent_id) == 12L) ||
    !all(table(x$choice_set_id) == 3L) ||
    !all(tapply(x$chosen, x$choice_set_id, sum) == 1L) ||
    anyNA(sp$parking_intention) || !all(sp$parking_intention %in% 1:3)) {
  stop("Data structure inconsistency. Check results_R/39_data_validation.csv.")
}

# CLMM diagnostics from SAVED fits, no extra optimization here.
paths <- c("results_R/ordered_mixed.rds",
           "results_R/clmm_info_interactions_v1.rds")
model_names <- c("Baseline CLMM", "Information-interaction CLMM")
diag_clmm <- do.call(rbind, lapply(seq_along(paths), function(i) {
  if (!file.exists(paths[i])) {
    return(data.frame(model = model_names[i], saved = FALSE,
                      optimizer_code = NA_character_,
                      max_abs_gradient = NA_real_,
                      gradient_over_0_001 = NA,
                      stringsAsFactors = FALSE))
  }
  m <- readRDS(paths[i])
  opt_code <- tryCatch(m$optRes$convergence, error = function(e) NULL)
  if (is.null(opt_code)) opt_code <- tryCatch(m$optRes$code,
                                              error = function(e) NULL)
  if (is.null(opt_code) || !length(opt_code)) opt_code <- NA_character_
  gr <- tryCatch(m$gradient, error = function(e) NULL)
  g <- if (length(gr) && all(is.finite(gr))) max(abs(gr)) else NA_real_
  data.frame(model = model_names[i], saved = TRUE,
             optimizer_code = as.character(opt_code[1]),
             max_abs_gradient = g,
             gradient_over_0_001 = if (is.finite(g)) g > 0.001 else NA,
             stringsAsFactors = FALSE)
}))
write.csv(diag_clmm, "results_R/40_model_convergence.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n== Saved CLMM optimizer checks ==\n")
print(diag_clmm, row.names = FALSE, digits = 6)

# Compute clustered SE for the saved nnet multinomial logistic fit, without
# re-estimating. Score per respondent is sum of 12 row-level score vectors.
# Sandwich = H^{-1} (sum_g score_g score_g') H^{-1}, with finite-cluster factor.
fit <- readRDS(model_path)
beta <- as.matrix(stats::coef(fit))
if (is.null(rownames(beta)) || nrow(beta) != 2L) {
  stop("Expected 2 nonreference outcome rows in multinomial coefficients.")
}
if (!all(rownames(beta) %in% levels(sp$intention_category))) {
  stop("Outcome labels in fitted model do not match the source data.")
}
X <- stats::model.matrix(stats::delete.response(stats::terms(fit)), data = sp)
if (!all(colnames(beta) %in% colnames(X))) {
  stop("Model-matrix columns do not match saved model coefficients: ",
       paste(setdiff(colnames(beta), colnames(X)), collapse = ", "))
}
X <- X[, colnames(beta), drop = FALSE]
p <- as.matrix(stats::predict(fit, newdata = sp, type = "probs"))
if (!all(rownames(beta) %in% colnames(p)) ||
    nrow(p) != nrow(sp) || any(!is.finite(p))) {
  stop("Predicted multinomial probabilities are malformed.")
}
if (any(abs(rowSums(p) - 1) > 1e-7)) {
  stop("Multinomial probabilities do not add to one.")
}

V <- stats::vcov(fit)
N <- nrow(sp)
P <- length(as.vector(t(beta)))
if (!all(dim(V) == c(P, P))) stop("Saved MNL vcov dimension mismatch.")
naive_se <- as.vector(t(summary(fit)$standard.errors))
# This guards against silently using the wrong coefficient ordering.
if (length(naive_se) != P ||
    max(abs(sqrt(diag(V)) - naive_se)) > 1e-5) {
  stop("Mismatch in coefficient order between vcov() and coef().")
}

score_blocks <- lapply(rownames(beta), function(k) {
  resid_k <- as.numeric(sp$intention_category == k) - p[, k]
  sweep(X, 1L, resid_k, `*`)
})
S <- do.call(cbind, score_blocks)
G <- nlevels(sp$respondent_id)
S_person <- rowsum(S, group = as.integer(sp$respondent_id), reorder = TRUE)
if (nrow(S_person) != G) stop("Unexpected number of respondent clusters.")
meat <- crossprod(S_person)
V_cluster <- V %*% meat %*% V
# Small-sample scaling, conventional CR1 adjustment.
V_cluster <- V_cluster * (G / (G - 1)) * ((N - 1) / (N - P))
se_cluster <- sqrt(pmax(0, diag(V_cluster)))
est <- as.vector(t(beta))
term <- rep(colnames(beta), times = nrow(beta))
outcome <- rep(rownames(beta), each = ncol(beta))
valid <- is.finite(se_cluster) & se_cluster > 0
if (!all(valid)) stop("Nonpositive/nonfinite respondent-clustered standard error.")
z_cluster <- est / se_cluster
out <- data.frame(outcome_vs_do_not_use = outcome, term = term,
                  estimate = est, naive_se = naive_se,
                  respondent_cluster_se = se_cluster,
                  z_cluster = z_cluster,
                  p_cluster = 2 * stats::pnorm(-abs(z_cluster)),
                  stringsAsFactors = FALSE)
write.csv(out, "results_R/41_intention_cluster_robust.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n== Parking-intention MNL: person-clustered standard errors ==\n")
cat("Clusters:", G, "; responses:", N,
    "; estimated coefficients:", P, "\n")
print(out[out$term %in% c("price_1000", "parking_user", "price_parking_user"), ],
      row.names = FALSE, digits = 4)

notes <- c(
  "MODEL VALIDATION (SIMULATED DATA ONLY)",
  "The dataset comprises 1,000 persons, 6 choice tasks per person,",
  "and 12 parking-intention responses per person (A/B in each task).",
  "Intention MNL: respondent-clustered sandwich SE computed by score aggregation.",
  "Clustered SE handle within-person correlation for inference; they do NOT",
  "model unobserved individual heterogeneity or validate a causal claim.",
  "CLMM optimization gradients are diagnostic only. Values > 0.001 flag",
  "the need for optimization-sensitivity / refitting checks before inference.",
  "Model fit, identifiability, calibration and real-data analysis remain pending.",
  "Synthetic-data significance is NOT evidence of real-world effects."
)
writeLines(notes, "results_R/42_validation_notes.txt", useBytes = TRUE)
cat("\nDONE: model validation (no model was re-estimated).\n")
cat("Created under results_R:\n")
cat(" 39_data_validation.csv\n 40_model_convergence.csv\n")
cat(" 41_intention_cluster_robust.csv\n 42_validation_notes.txt\n")
if (any(diag_clmm$gradient_over_0_001, na.rm = TRUE)) {
  cat("WARNING: At least one CLMM gradient is above 0.001.\n")
}


# 11 モデル適合度 ----
# 元ファイル: 11_fit_assessment(1).R
# ============================================================================
# 11_fit_assessment.R
# Hiroshima University parking SP -- SYNTHETIC 1,000-person exercise
# Model fit assessment using ALREADY SAVED models. NO re-estimation.
# Run from the directory containing analysis_long.csv and results_R/:
#   source("11_fit_assessment.R", encoding = "UTF-8")
# ============================================================================

required_files <- c(
  "analysis_long.csv",
  "results_R/mnl_choice.rds",
  "results_R/mnl_info_interactions_v1.rds",
  "results_R/ordered_mixed.rds",
  "results_R/clmm_info_interactions_v1.rds",
  "results_R/parking_intention_mnl_v1.rds"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required file(s) not found: ", paste(missing_files, collapse = ", "),
       ". Check getwd() and the existing R models; this script does NOT refit.")
}
for (p in c("survival", "ordinal", "nnet")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Package not installed: ", p)
  }
}
# Loading namespaces makes the saved models' print/logLik/predict methods available.
dir.create("results_R", showWarnings = FALSE, recursive = TRUE)

models <- list(
  "Choice MNL (baseline)" = readRDS("results_R/mnl_choice.rds"),
  "Choice MNL (info interactions)" = readRDS("results_R/mnl_info_interactions_v1.rds"),
  "Attitude CLMM (baseline)" = readRDS("results_R/ordered_mixed.rds"),
  "Attitude CLMM (info interactions)" = readRDS("results_R/clmm_info_interactions_v1.rds"),
  "Intention multinom" = readRDS("results_R/parking_intention_mnl_v1.rds")
)

safe_metric <- function(expr) {
  tryCatch(as.numeric(expr)[1L], error = function(e) NA_real_)
}
get_optimizer_code <- function(m) {
  if (inherits(m, "clmm")) {
    code <- m$optRes$convergence
    if (is.null(code)) code <- m$optRes$code
  } else if (inherits(m, "multinom")) {
    code <- m$convergence
  } else {
    code <- NA_integer_
  }
  if (is.null(code) || length(code) == 0L) return(NA_character_)
  as.character(code[1L])
}

# 1. Summary of fit statistics, each based on the model's own outcome.
outcome_group <- c("A/B/SQ choice", "A/B/SQ choice",
                   "A/B policy attitude", "A/B policy attitude",
                   "A/B parking intention")
fit_summary <- do.call(rbind, lapply(seq_along(models), function(i) {
  m <- models[[i]]
  ll <- stats::logLik(m)
  gr <- if (inherits(m, "clmm") && !is.null(m$gradient))
    max(abs(m$gradient)) else NA_real_
  data.frame(
    model = names(models)[i], outcome = outcome_group[i],
    observations = safe_metric(stats::nobs(m)),
    parameter_count = as.integer(attr(ll, "df")),
    logLik = as.numeric(ll), AIC = safe_metric(stats::AIC(m)),
    BIC = safe_metric(stats::BIC(m)),
    optimizer_code = get_optimizer_code(m),
    max_abs_gradient = gr, stringsAsFactors = FALSE
  )
}))
write.csv(fit_summary, "results_R/43_fit_summary.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n== Existing model fit summary (SIMULATED DATA) ==\n")
print(fit_summary, row.names = FALSE, digits = 5)

# 2. Nested model comparisons; ONLY compare within the SAME outcome family.
# An LR test assumes both fits used the same records, likelihood definition,
# and were optimized sufficiently well. CLMM gradients are separately flagged.
compare_pair <- function(base_name, full_name, outcome) {
  a <- models[[base_name]]
  b <- models[[full_name]]
  ll_a <- stats::logLik(a)
  ll_b <- stats::logLik(b)
  df_added <- as.integer(attr(ll_b, "df") - attr(ll_a, "df"))
  lr <- 2 * (as.numeric(ll_b) - as.numeric(ll_a))
  if (!is.finite(lr) || df_added < 1) {
    p <- NA_real_
  } else {
    p <- stats::pchisq(max(0, lr), df = df_added, lower.tail = FALSE)
  }
  data.frame(
    outcome = outcome, baseline_model = base_name, extended_model = full_name,
    added_parameters = df_added,
    LR_statistic = lr, LR_p_value = p,
    delta_AIC_extended_minus_baseline = stats::AIC(b) - stats::AIC(a),
    delta_BIC_extended_minus_baseline = safe_metric(stats::BIC(b)) -
      safe_metric(stats::BIC(a)),
    stringsAsFactors = FALSE
  )
}
model_comparison <- rbind(
  compare_pair("Choice MNL (baseline)", "Choice MNL (info interactions)",
               "A/B/SQ choice"),
  compare_pair("Attitude CLMM (baseline)", "Attitude CLMM (info interactions)",
               "A/B policy attitude")
)
write.csv(model_comparison, "results_R/44_fit_model_comparison.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n== Baseline vs interaction fit (same outcome only) ==\n")
print(model_comparison, row.names = FALSE, digits = 5)

# 3. In-sample predictive diagnostics for choice and intention.
# These are NOT out-of-sample tests and cannot establish generalization.
x <- read.csv("analysis_long.csv", fileEncoding = "UTF-8-BOM",
              stringsAsFactors = FALSE)
stopifnot(nrow(x) == 18000L, sum(x$chosen) == 6000L,
          length(unique(x$choice_set_id)) == 6000L)
x$info_group <- factor(x$info_group, levels = paste0("I", 1:5))
for (g in 2:5) {
  flag <- as.integer(x$info_group == paste0("I", g))
  x[[paste0("price_info", g)]] <- x$price_1000 * flag
}
x$gate_ic_info2 <- x$gate_ic * as.integer(x$info_group == "I2")
x$gate_lpr_info2 <- x$gate_lpr * as.integer(x$info_group == "I2")
x$safety_basic_info5 <- x$safety_basic * as.integer(x$info_group == "I5")
x$safety_plus_info5 <- x$safety_plus * as.integer(x$info_group == "I5")

predict_choice_fit <- function(m, label) {
  b <- stats::coef(m)
  if (anyNA(b) || length(setdiff(names(b), names(x)))) {
    stop("Missing/NA choice-model coefficient predictors: ", label,
         "; missing=", paste(setdiff(names(b), names(x)), collapse = ", "))
  }
  X <- as.matrix(x[, names(b), drop = FALSE])
  storage.mode(X) <- "double"
  eta <- drop(X %*% b)
  group <- x$choice_set_id
  # Normalize utilities within each of 6,000 choice sets.
  centered <- eta - ave(eta, group, FUN = max)
  e <- exp(centered)
  p <- e / ave(e, group, FUN = sum)
  if (any(!is.finite(p)) || any(p <= 0)) stop("Invalid MNL predictions")
  split_idx <- split(seq_along(p), group)
  acc <- mean(vapply(split_idx,
    function(idx) x$chosen[idx[which.max(p[idx])]] == 1L, logical(1)))
  logloss <- -mean(log(p[x$chosen == 1L]))
  data.frame(model = label, outcome = "A/B/SQ choice",
             n_answer_units = length(split_idx),
             in_sample_log_loss = logloss,
             in_sample_top_choice_accuracy = acc,
             reference_log_loss = log(3),
             reference_accuracy = 1/3,
             stringsAsFactors = FALSE)
}

choice_diagnostics <- rbind(
  predict_choice_fit(models[["Choice MNL (baseline)"]],
                     "Choice MNL (baseline)"),
  predict_choice_fit(models[["Choice MNL (info interactions)"]],
                     "Choice MNL (info interactions)")
)

sp <- x[x$alternative %in% c("A", "B"), , drop = FALSE]
sp$intention_category <- factor(sp$parking_intention,
  levels = c(2L, 1L, 3L), labels = c("Do not use", "Use", "Unsure"))
intent_model <- models[["Intention multinom"]]
p_int <- as.matrix(stats::predict(intent_model, newdata = sp, type = "probs"))
truth <- as.character(sp$intention_category)
idx_truth <- match(truth, colnames(p_int))
if (nrow(p_int) != nrow(sp) || anyNA(idx_truth) ||
    any(!is.finite(p_int)) || any(p_int <= 0) ||
    any(abs(rowSums(p_int) - 1) > 1e-6)) {
  stop("Intention predictions and classes are inconsistent.")
}
classes <- colnames(p_int)[max.col(p_int, ties.method = "first")]
shares <- prop.table(table(sp$intention_category))
intention_diagnostics <- data.frame(
  model = "Intention multinom", outcome = "A/B parking intention",
  n_answer_units = nrow(sp),
  in_sample_log_loss = -mean(log(p_int[cbind(seq_len(nrow(sp)), idx_truth)])),
  in_sample_top_choice_accuracy = mean(classes == truth),
  reference_log_loss = -sum(as.numeric(shares) * log(as.numeric(shares))),
  reference_accuracy = max(as.numeric(shares)),
  stringsAsFactors = FALSE
)

predictive_summary <- rbind(choice_diagnostics, intention_diagnostics)
write.csv(predictive_summary, "results_R/45_in_sample_prediction.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n== IN-SAMPLE prediction (NOT a validation/test-set score) ==\n")
print(predictive_summary, row.names = FALSE, digits = 5)

notes <- c(
  "MODEL FIT ASSESSMENT: SYNTHETIC RESPONDENTS ONLY",
  "All 5 saved models were READ from disk. No model was re-estimated.",
  "Do NOT compare AIC/BIC directly across DIFFERENT response variables:",
  "policy choice, policy attitude, and parking intention are distinct outcomes.",
  "44_fit_model_comparison.csv compares nested models within the same outcome.",
  "AIC lower is better within the comparable model pair, all else equal.",
  "45_in_sample_prediction.csv uses the same data as model estimation.",
  "It checks calculation and in-sample fit, NOT out-of-sample forecasting ability.",
  "The equal-probability reference for A/B/SQ choice is 1/3 per alternative.",
  "The intention reference uses overall empirical simulated class frequencies.",
  "CLMM optimizer gradients in 43_fit_summary.csv should be reviewed:",
  "earlier runs showed 0.00943 (baseline) and 0.02006 (interaction).",
  "Therefore CLMM LR p-values and coefficient p-values are PROVISIONAL.",
  "Respondent-level repeat answers and simulated nature limit interpretation.",
  "Real survey effects, causal effects, and population prediction are NOT established."
)
writeLines(notes, "results_R/46_fit_notes.txt", useBytes = TRUE)

cat("\nDONE: model-fit diagnostics; NO refitting.\n")
cat("Created under results_R:\n")
cat(" 43_fit_summary.csv\n")
cat(" 44_fit_model_comparison.csv\n")
cat(" 45_in_sample_prediction.csv\n")
cat(" 46_fit_notes.txt\n")
cat("NOTE: CLMM inference remains provisional due to optimizer gradients.\n")


# 12 卒論用推定結果表 ----
# 元ファイル: 12_export_thesis_tables(1).R
# ==========================================================================
# 12_export_thesis_tables.R
# Export all estimated coefficients as thesis-ready CSV tables.
# Does NOT re-estimate or alter any model.
# Requires previously saved results_R/*.rds files in current working folder.
# For simulated data only; replace with real estimates after actual survey.
# ==========================================================================
files <- c(
  choice = "results_R/mnl_choice.rds",
  attitude = "results_R/ordered_mixed.rds",
  intention = "results_R/parking_intention_mnl_v1.rds"
)
missing <- files[!file.exists(files)]
if (length(missing)) stop("Missing saved models: ", paste(missing, collapse=", "))
if (!requireNamespace("survival", quietly=TRUE) ||
    !requireNamespace("ordinal", quietly=TRUE) ||
    !requireNamespace("nnet", quietly=TRUE)) {
  stop("Install required packages: survival, ordinal, nnet")
}
dir.create("results_R", showWarnings=FALSE)

# Policy choice MNL (survival::clogit); variance incorporates
# respondent clustering if the original model specified cluster(id).
choice <- readRDS(files[["choice"]])
b <- stats::coef(choice)
V <- stats::vcov(choice)
s <- sqrt(diag(V))[names(b)]
if (anyNA(s)) stop("Choice model SE mismatch: inspect vcov(choice)")
z <- as.numeric(b/s)
choice_table <- data.frame(
  variable=names(b), estimate=as.numeric(b), clustered_se=as.numeric(s),
  z_value=z, p_value=2*stats::pnorm(-abs(z)),
  odds_ratio=exp(as.numeric(b)), check.names=FALSE)
write.csv(choice_table, "results_R/47_thesis_choice_MNL_coefficients.csv",
          row.names=FALSE, fileEncoding="UTF-8")

# Policy attitude CLMM (ordinal::clmm). Wald statistics are provisional
# until optimizer gradient / numerical convergence is resolved.
att <- readRDS(files[["attitude"]])
b <- att$beta
V <- stats::vcov(att)
s <- sqrt(diag(V))[names(b)]
if (anyNA(s)) stop("CLMM SE mismatch: inspect vcov(att)")
z <- as.numeric(b/s)
att_table <- data.frame(
  variable=names(b), estimate=as.numeric(b), model_based_se=as.numeric(s),
  z_value=z, p_value=2*stats::pnorm(-abs(z)), check.names=FALSE)
write.csv(att_table, "results_R/48_thesis_attitude_CLMM_coefficients.csv",
          row.names=FALSE, fileEncoding="UTF-8")

# Parking intention MNL: use existing RESPONDENT-CLUSTERED SE, not
# naive SE. 10_validate_models.R must have been run.
robust_path <- "results_R/41_intention_cluster_robust.csv"
if (file.exists(robust_path)) {
  intention <- read.csv(robust_path, fileEncoding="UTF-8-BOM",
                        check.names=FALSE)
  write.csv(intention, "results_R/49_thesis_intention_MNL_coefficients.csv",
            row.names=FALSE, fileEncoding="UTF-8")
} else {
  warning("Intention cluster-SE table not found: run 10_validate_models.R",
          " to create results_R/41_intention_cluster_robust.csv")
}

# Fit index table that belongs in the same model family comparisons
fit_path <- "results_R/43_fit_summary.csv"
if (file.exists(fit_path)) {
  fit <- read.csv(fit_path, fileEncoding="UTF-8-BOM", check.names=FALSE)
  write.csv(fit, "results_R/50_thesis_model_fit.csv",
            row.names=FALSE, fileEncoding="UTF-8")
}
cat("DONE: Thesis tables exported (no model fitting).\n")
cat(" 47_thesis_choice_MNL_coefficients.csv\n")
cat(" 48_thesis_attitude_CLMM_coefficients.csv\n")
cat(" 49_thesis_intention_MNL_coefficients.csv (if cluster table exists)\n")
cat(" 50_thesis_model_fit.csv (if summary exists)\n")
cat("NOTICE: SIMULATED data only; CLMM inference is provisional.\n")
