# 선택한 모형을 Python 코드로. statsmodels에 같은 모형이 없으면 rpy2로 R 모형을 그대로 부른다.
py_for <- function(i, file, sep = ",", enc = "UTF-8") {
  k <- S$L[[i]]; n <- k$name; m <- S$fits[[as.character(i)]]
  Q <- function(x) as.character(toJSON(x, auto_unbox = TRUE))
  col <- function(v) sprintf("d[%s]", Q(lab(v)))
  bt <- function(s) { for (j in seq_along(S$labels)) s <- gsub(sprintf("vv%d__", j), sprintf("Q(%s)", Q(S$labels[j])), s, fixed = TRUE); gsub(")^2", ")**2", s, fixed = TRUE) }
  has_re <- length(S$re) > 0; off <- if (!is.null(S$off)) sprintf("np.log(%s)", col(S$off))
  offa <- if (is.null(off)) "" else paste0(", offset=", off)
  yv <- col(S$yv); body <- NULL
  zero <- c("# {{0이냐 아니냐 부분 (상수만)||zero vs non-zero part (constant only)}}", "z = (d0[ycol] == 0).astype(int)", "m_zero = sm.Logit(z, np.ones(len(z))).fit(disp=0)", "print(m_zero.summary())")
  pos <- c("d0, d = d, d[d[ycol] > 0].copy()  # {{양수 부분만||positive part only}}")
  X <- sprintf("X = dmatrix(\"\"\"%s\"\"\", d, return_type=\"dataframe\")", bt(S$rhs))
  glmc <- function(fam) c(X, sprintf("m = sm.GLM(%s, X, family=sm.families.%s%s).fit()", yv, fam, offa))
  dsc <- function(cls, extra = "") c(X, sprintf("m = %s(%s, X%s%s).fit(method=\"bfgs\", maxiter=2000, disp=0)", cls, yv, offa, extra))
  gl <- "# {{참고: GLM의 산포 추정 방식이 달라 표준오차와 AIC가 앱과 조금 다를 수 있음 (계수는 같음)||Note: dispersion is estimated differently, so standard errors and AIC can differ from the app (coefficients match)}}"
  lnm <- sprintf("m = sm.OLS(np.log(%s)%s, X).fit()", yv, if (is.null(off)) "" else paste(" -", off))
  if (!has_re) body <- switch(n,
    "Gaussian" = c(X, sprintf("m = sm.OLS(%s, X).fit()", yv)),
    "Gaussian (log)" = c(gl, glmc("Gaussian(sm.families.links.Log())")),
    "Gamma (log)" = c(gl, glmc("Gamma(sm.families.links.Log())")),
    "Gamma (identity)" = c(gl, glmc("Gamma(sm.families.links.Identity())")),
    "Gamma (inverse)" = c(gl, glmc("Gamma(sm.families.links.InversePower())")),
    "Inverse Gaussian (log)" = c(gl, glmc("InverseGaussian(sm.families.links.Log())")),
    "Lognormal" = c(X, lnm, sprintf("print(\"{{원척도 AIC:||Raw-scale AIC:}}\", m.aic + 2 + 2 * np.log(%s).sum())", yv)),
    "Hurdle Lognormal" = c(pos, X, lnm, zero),
    "Hurdle Gamma" = c(pos, gl, glmc("Gamma(sm.families.links.Log())"), zero),
    "Tweedie (log)" = c("# {{var_power는 앱이 추정한 값||var_power is the value estimated by the app}}", gl, glmc(sprintf("Tweedie(link=sm.families.links.Log(), var_power=%.6f)", family_params(m)[1]))),
    "Poisson" = glmc("Poisson()"),
    "NB2" = dsc("sm.NegativeBinomial", ", loglike_method=\"nb2\""),
    "NB1" = dsc("sm.NegativeBinomial", ", loglike_method=\"nb1\""),
    "Generalized Poisson" = dsc("sm.GeneralizedPoisson", ", p=1"),
    "ZI Poisson" = dsc("sm.ZeroInflatedPoisson"),
    "Hurdle Poisson" = c(pos, dsc("TruncatedLFPoisson"), zero),
    "Zero-truncated Poisson" = dsc("TruncatedLFPoisson"),
    "Beta" = c(X, sprintf("m = BetaModel(%s, X).fit(disp=0)", yv)),
    "Zero-inflated Beta" = c(pos, X, sprintf("m = BetaModel(%s, X).fit(disp=0)", yv), zero),
    "Binomial (logit)" = glmc("Binomial()"),
    "Binomial (probit)" = glmc("Binomial(sm.families.links.Probit())"),
    "Binomial (cloglog)" = glmc("Binomial(sm.families.links.CLogLog())"),
    "Binomial" = c(X, "# {{[성공, 실패]||[successes, failures]}}", if (!is.null(S$fl)) sprintf("Y = np.column_stack([%s, %s])", yv, col(S$fl)) else sprintf("Y = np.column_stack([%1$s, %2$s - %1$s])", yv, if (!is.null(S$tr)) col(S$tr) else ""), "m = sm.GLM(Y, X, family=sm.families.Binomial()).fit()"),
    "Cumulative (logit)" = , "Cumulative (probit)" = if (S$rhs != "1") c(X, "X = X.drop(columns=\"Intercept\")",
      sprintf("y = pd.Series(pd.Categorical(%s, ordered=True), index=d.index)", yv),
      sprintf("m = OrderedModel(y, X, distr=\"%s\").fit(method=\"bfgs\", maxiter=2000, disp=0)", k$fam)))
  else if (length(S$re) == 1 && n %in% c("Gaussian", "Lognormal") && is.null(off)) {
    lhs <- if (n == "Gaussian") sprintf("Q(%s)", Q(lab(S$yv))) else sprintf("np.log(Q(%s))", Q(lab(S$yv)))
    body <- c(sprintf("m = smf.mixedlm(\"\"\"%s ~ %s\"\"\", d, groups=%s%s).fit(reml=False, method=[\"powell\", \"lbfgs\"])", lhs, bt(S$rhs), col(S$re),
      if (is.null(S$slope)) "" else sprintf(", re_formula=\"\"\"~%s\"\"\"", bt(S$slope))),
      if (n == "Lognormal") sprintf("print(\"{{원척도 AIC:||Raw-scale AIC:}}\", -2 * m.llf + 2 * (len(m.params) + 1) + 2 * np.log(%s).sum())", yv))
  }
  # 따옴표/역슬래시가 있는 열 이름은 수식 문자열에서 깨질 수 있어 rpy2로
  if (any(grepl("[\"'\\\\]", lab(names(S$d))))) body <- NULL
  if (is.null(body)) {
    r <- strsplit(code_for(i, file, sep, enc), "\n")[[1]]
    r <- r[!grepl("^plot\\(|^testDispersion|^summary\\(m\\)$", r)]; r <- sub("; library(DHARMa)", "", r, fixed = TRUE)
    if (any(grepl("'''", r, fixed = TRUE))) return(js(list(code = "# {{열 이름에 따옴표가 있어 Python 코드를 만들기 어렵습니다. 열 이름을 바꿔 주십시오.||Cannot build Python code because a column name contains quotes. Please rename it.}}", native = FALSE)))
    out <- c(sprintf("# %s", n), "# {{Python(statsmodels)에는 같은 모형이 없어서 rpy2로 R의 모형을 그대로 부릅니다.||statsmodels has no reliable equivalent, so this calls the same R model through rpy2.}}",
      "# {{준비: R 설치, R에서 install.packages(c(\"glmmTMB\", \"ordinal\")), 터미널에서 pip install rpy2||Setup: install R, run install.packages(c(\"glmmTMB\", \"ordinal\")) in R, then pip install rpy2 in a terminal}}",
      "import rpy2.robjects as ro", "", "ro.r(r'''", r, "''')", "print(ro.r(\"summary(m)\"))")
    return(js(list(code = paste(out, collapse = "\n"), native = FALSE)))
  }
  used <- names(S$d); fac <- used[sapply(S$d, is.factor)]
  imp <- c("import numpy as np", "import pandas as pd", "import statsmodels.api as sm",
    if (has_re) "import statsmodels.formula.api as smf" else "from patsy import dmatrix",
    if (any(grepl("TruncatedLF", body))) "from statsmodels.discrete.truncated_model import TruncatedLFPoisson",
    if (any(grepl("BetaModel", body))) "from statsmodels.othermod.betareg import BetaModel",
    if (any(grepl("OrderedModel", body))) "from statsmodels.miscmodels.ordinal_model import OrderedModel")
  o <- c(sprintf("# %s", n), "# {{준비: pip install numpy pandas statsmodels||Setup: pip install numpy pandas statsmodels}}", imp, "",
    sprintf("d = pd.read_csv(%s, sep=%s, encoding=\"%s\", skipinitialspace=True,", Q(file), Q(sep), if (enc == "UTF-8") "utf-8-sig" else "cp949"),
    "                keep_default_na=False, na_values=[\"\", \"NA\", \"NaN\", \"N/A\", \".\"])",
    sprintf("used = %s", Q(I(lab(used)))), "d = d.dropna(subset=used).copy()", sprintf("ycol = %s", Q(lab(S$yv))))
  for (f in setdiff(fac, S$re)) { lv <- levels(S$d[[f]])
    o <- c(o, sprintf("%1$s = pd.Categorical(%1$s, categories=%2$s)", col(f), if (is.numeric(S$raw[[f]])) paste0("[", paste(lv, collapse = ", "), "]") else Q(I(lv)))) }
  for (z in S$lx) o <- c(o, sprintf("%1$s = np.log(%1$s)", col(z)))
  for (z in c(S$zx, if (S$zy) S$yv)) o <- c(o, sprintf("%1$s = (%1$s - %1$s.mean()) / %1$s.std(){{  # 표준화||  # standardize}}", col(z)))
  if (S$kind == "binary") o <- c(o, sprintf("%1$s = (%1$s == %2$s).astype(int)", yv, Q(S$succ)))
  o <- c(o, "", body, "print(m.summary())")
  js(list(code = paste(o, collapse = "\n"), native = TRUE))
}
