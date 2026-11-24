# GLM 분포 비교 엔진. webR(브라우저)와 일반 R에서 모두 실행 가능.
suppressPackageStartupMessages({library(glmmTMB); library(jsonlite)})
S <- new.env()
js <- function(x) as.character(toJSON(x, auto_unbox = TRUE, digits = 8, na = "null", null = "null"))
NSIM <- 250

load_csv <- function(path, sep = ",") {
  d <- read.csv(path, sep = sep, check.names = FALSE, encoding = "UTF-8",
                na.strings = c("", "NA", "NaN", "N/A", "."), strip.white = TRUE)
  S$labels <- names(d)
  names(d) <- sprintf("vv%d__", seq_along(d))
  for (i in seq_along(d)) if (is.logical(d[[i]])) d[[i]] <- as.integer(d[[i]])
  S$raw <- d
  js(list(n = nrow(d), cols = I(lapply(seq_along(d), function(i) {
    x <- d[[i]]; v <- x[!is.na(x)]
    list(label = S$labels[i], numeric = is.numeric(x), uniq = length(unique(v)), int = is.numeric(x) && all(abs(v - round(v)) < 1e-9),
         missing = sum(is.na(x)), min = if (is.numeric(x) && length(v)) min(v) else NA, max = if (is.numeric(x) && length(v)) max(v) else NA, kind = detect(x),
         range = if (is.numeric(x) && length(v)) paste(signif(range(v), 4), collapse = " ~ ")
                 else paste(head(unique(as.character(v)), 4), collapse = ", "))
  }))))
}

detect <- function(x) {
  x <- x[!is.na(x)]; u <- unique(x)
  if (!length(u)) return("none")
  if (!is.numeric(x)) return(if (length(u) == 2) "binary" else "none")
  int <- all(abs(x - round(x)) < 1e-9)
  if (length(u) == 2 && all(u %in% 0:1)) return("binary")
  if (int && min(x) >= 1 && max(x) <= 11 && length(u) >= 3) return("ordinal")
  if (int && min(x) > 0 && max(x) >= 100 && length(u) > 20) return("pos")
  if (int && min(x) >= 0) return("count")
  if (min(x) >= 0 && max(x) <= 1) return("prop")
  if (min(x) > 0) return("pos")
  if (min(x) >= 0) return("zeropos")
  "real"
}

# wide -> long. cols: 문항 열 번호, id: 행마다 다른 ID 열 번호(0이면 새로 만듦)
to_long <- function(cols, id = 0) {
  d <- S$raw; n <- nrow(d); keep <- setdiff(seq_along(d), c(cols, id)); other <- d[rep(seq_len(n), length(cols)), keep, drop = FALSE]
  nm <- S$labels[keep]; nw <- c("participant", "item", "rating"); while (any(nw %in% nm)) nw <- paste0(nw, "_")
  pid <- if (id > 0) as.character(d[[id]]) else sprintf("P%03d", seq_len(n))
  o <- cbind(data.frame(rep(pid, length(cols)), rep(S$labels[cols], each = n), unlist(d[cols], use.names = FALSE)), other)
  names(o) <- c(nw, nm); o <- o[order(rep(seq_len(n), length(cols))), ]
  f <- tempfile(); write.csv(o, f, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  paste(readLines(f, encoding = "UTF-8"), collapse = "
")
}

# {{한국어||English}} 표기는 화면에서 언어에 맞게 풀린다. bf는 양쪽에 같은 값을 넣는 sprintf
bf <- function(f, ...) { p <- strsplit(substr(f, 3, nchar(f) - 2), "||", fixed = TRUE)[[1]]
  paste0("{{", sprintf(p[1], ...), "||", sprintf(p[2], ...), "}}") }

# cmp: AIC를 서로 비교할 수 있는 묶음인지 (FALSE면 참고용)
cand <- function(name, fam, zi = "~0", engine = "tmb", cmp = TRUE, link = "log", note = "")
  list(name = name, fam = fam, zi = zi, engine = engine, cmp = cmp, link = link, note = note)

candidates <- function(kind, y, has_re, links, offset) {
  z <- any(y == 0)
  L <- switch(kind,
    real = list(cand("Gaussian", "gaussian()", link = "identity"),
      cand("Student t", "t_family()", link = "identity", note = "{{꼬리가 두꺼울 때||for heavy tails}}"),
      cand("Skew normal", "skewnormal()", link = "identity", note = "{{비대칭일 때||for asymmetry}}")),
    pos = c(list(cand("Gaussian", "gaussian()", link = "identity"),
      cand("Gamma (log)", "Gamma(link='log')"),
      cand("Lognormal", "lognormal", engine = "ln", note = "{{log(Y)에 정규분포, AIC는 원척도로 보정함||normal on log(Y); AIC corrected to the raw scale}}")),
      if (!has_re) list(cand("Inverse Gaussian (log)", "inverse.gaussian(link='log')", engine = "glm")),
      if (links) list(cand("Gaussian (log)", "gaussian(link='log')"),
        cand("Gamma (identity)", "Gamma(link='identity')", link = "identity"),
        cand("Gamma (inverse)", "Gamma(link='inverse')", link = "inverse"))),
    zeropos = list(cand("Tweedie (log)", "tweedie(link='log')"),
      cand("Hurdle Gamma", "ziGamma(link='log')", zi = "~1", note = "{{0 확률은 상수로 가정||zero probability assumed constant}}"),
      cand("Hurdle Lognormal", "lognormal", engine = "ln", zi = "~1", note = "{{0 확률은 상수로 가정||zero probability assumed constant}}"),
      cand("Gaussian", "gaussian()", link = "identity", cmp = FALSE, note = "{{0의 점확률을 다루지 않아 AIC 비교가 어려움||ignores the point mass at 0, so AIC is not comparable}}")),
    count = c(list(cand("Poisson", "poisson()"), cand("NB2", "nbinom2()", note = "{{음이항, 분산 = μ + μ²/θ||negative binomial, variance = μ + μ²/θ}}"),
      cand("NB1", "nbinom1()", note = "{{음이항, 분산이 평균에 비례||negative binomial, variance proportional to mean}}"), cand("Generalized Poisson", "genpois()")),
      if (z) list(cand("ZI Poisson", "poisson()", zi = "~1"), cand("ZI NB2", "nbinom2()", zi = "~1"),
        cand("Hurdle Poisson", "truncated_poisson()", zi = "~1"),
        cand("Hurdle NB2", "truncated_nbinom2()", zi = "~1"))
      else list(cand("Zero-truncated Poisson", "truncated_poisson()", note = "{{0이 구조적으로 불가능할 때||when zeros are structurally impossible}}"),
        cand("Zero-truncated NB2", "truncated_nbinom2()", note = "{{0이 구조적으로 불가능할 때||when zeros are structurally impossible}}"))),
    prop = { edge <- any(y == 0 | y == 1)
      c(if (!edge) list(cand("Beta", "beta_family()", link = "logit")),
        if (edge) list(cand("Ordered Beta", "ordbeta()", link = "logit", note = "{{0과 1을 허용||allows 0 and 1}}")),
        if (z && !any(y == 1)) list(cand("Zero-inflated Beta", "beta_family()", zi = "~1", link = "logit")),
        list(cand("Gaussian", "gaussian()", link = "identity", cmp = !edge,
          note = if (edge) "{{경계값의 점확률 때문에 AIC 비교가 어려움||AIC not comparable because of point mass at the boundary}}" else "{{범위를 벗어난 예측 가능||can predict outside 0-1}}"))) },
    binary = list(cand("Binomial (logit)", "binomial('logit')", link = "logit"),
      cand("Binomial (probit)", "binomial('probit')", link = "probit"),
      cand("Binomial (cloglog)", "binomial('cloglog')", link = "cloglog")),
    trials = list(cand("Binomial", "binomial()", link = "logit"),
      cand("Beta-binomial", "betabinomial()", link = "logit", note = "{{과산포 허용||allows overdispersion}}")),
    ordinal = lapply(c("logit", "probit"), function(l)
      cand(paste0("Cumulative (", l, ")"), l, engine = "ord", link = l)))
  if (offset) L <- Filter(function(k) k$link == "log", L)
  L
}

lab <- function(k) S$labels[as.integer(gsub("\\D", "", k))]

setup <- function(cfgjson) {
  cfg <- fromJSON(cfgjson)
  v <- function(i) sprintf("vv%d__", as.integer(i))
  yv <- v(cfg$y); num <- v(cfg$num); fac <- v(cfg$fac); re <- v(cfg$re)
  off <- if (length(cfg$off)) v(cfg$off); tr <- if (length(cfg$tri)) v(cfg$tri)
  d <- S$raw[c(yv, num, fac, re, off, tr)]; n0 <- nrow(d); d <- d[complete.cases(d), , drop = FALSE]
  warn <- character(); kind <- cfg$kind
  if (nrow(d) < 10) stop("{{결측을 뺀 행이 10개 미만이라 계산이 어렵습니다.||Fewer than 10 rows remain after removing missing values.}}")
  for (k in c(num, off, tr)) if (!is.numeric(d[[k]])) stop(bf("{{'%s' 열이 숫자 형식인지 확인해 주십시오.||Column '%s' is not numeric.}}", lab(k)))
  if (kind == "ordinal") {
    if (!is.numeric(d[[yv]])) stop("{{순서형 Y는 1, 2, 3처럼 숫자로 코딩해 주십시오.||Please code ordinal Y as numbers such as 1, 2, 3.}}")
    if (length(unique(d[[yv]])) < 3) stop("{{순서형 Y는 단계가 3개 이상 필요합니다.||Ordinal Y needs at least 3 levels.}}")
  } else if (kind == "binary") {
    u <- sort(unique(d[[yv]])); if (length(u) != 2) stop("{{이분형 Y의 값이 두 종류인지 확인해 주십시오.||Binary Y must have exactly two values.}}")
    S$succ <- u[2]; d[[yv]] <- as.integer(d[[yv]] == u[2])
    if (!is.numeric(u) || any(u != c(0, 1))) warn <- c(warn, bf("{{Y의 '%s'를 1(성공)로 두었습니다.||'%s' in Y is treated as 1 (success).}}", u[2]))
  } else if (!is.numeric(d[[yv]])) stop("{{Y가 숫자 형식인지 확인해 주십시오.||Y is not numeric.}}")
  y <- d[[yv]]; int <- all(abs(y - round(y)) < 1e-9)
  if (any(!is.finite(y))) stop("{{Y에 무한대 값이 있는지 확인해 주십시오.||Y contains infinite values.}}")
  bad <- switch(kind, pos = any(y <= 0), zeropos = any(y < 0), count = any(y < 0) || !int,
    prop = any(y < 0 | y > 1), FALSE)
  if (bad) stop("{{Y 값이 고른 유형과 맞지 않을 수 있습니다 (음수, 0, 소수 여부를 확인해 주십시오).||Y does not match the selected type (check for negatives, zeros or decimals).}}")
  if (kind == "trials") {
    if (is.null(tr)) stop("{{시행 수 열을 지정해 주십시오.||Please assign a trials column.}}")
    if (!int || any(y < 0 | y > d[[tr]]) || any(d[[tr]] < 1 | d[[tr]] %% 1 != 0)) stop("{{성공 수와 시행 수 값을 확인해 주십시오.||Invalid successes/trials values.}}")
  }
  if (length(unique(y)) < 2) stop("{{Y 값이 모두 같습니다.||All Y values are identical.}}")
  if (!is.null(off)) {
    if (any(d[[off]] <= 0)) stop("{{노출량에 0 이하인 값이 있는지 확인해 주십시오.||Exposure must be greater than 0.}}")
    if (!kind %in% c("pos", "zeropos", "count")) stop("{{노출량은 로그 링크를 쓰는 Y 유형에서만 쓸 수 있습니다.||Exposure is only available for Y types that use a log link.}}")
  }
  for (k in c(fac, re)) d[[k]] <- droplevels(factor(d[[k]]))
  for (k in c(fac, re)) if (nlevels(d[[k]]) < 2) stop(bf("{{'%s'의 수준이 하나뿐입니다.||'%s' has only one level.}}", lab(k)))
  for (k in num) if (length(unique(d[[k]])) < 2) stop(bf("{{'%s' 값이 모두 같습니다.||All values of '%s' are identical.}}", lab(k)))
  S$zx <- if (isTRUE(cfg$zx)) num else character(); S$zy <- isTRUE(cfg$zy)
  for (k in S$zx) d[[k]] <- as.numeric(scale(d[[k]]))
  if (S$zy) {
    if (kind %in% c("binary", "trials", "ordinal")) stop("{{이 Y 유형은 표준화를 지원하지 않습니다.||This Y type cannot be standardized.}}")
    if (!is.null(off)) stop("{{Y를 표준화하면 노출량(offset)은 지원하지 않습니다.||Exposure (offset) cannot be used with standardized Y.}}")
    d[[yv]] <- y <- as.numeric(scale(y))
    warn <- c(warn, "{{Y를 표준화해서 Gaussian만 돌렸습니다. 표준화해도 분포 모양은 그대로라 치우친 Y는 여전히 치우쳐 있을 수 있습니다.||Y is standardized, so only Gaussian was run. Standardizing keeps the shape: a skewed Y does not become normal.}}")
    if (kind %in% c("count", "prop", "zeropos")) warn <- c(warn, "{{원래 Y가 횟수, 비율, 0이 섞인 값이라 Gaussian이 부적절할 수 있습니다. raw Y로도 비교해 보십시오.||The original Y is a count, proportion or zero-mixed value, so Gaussian may be inappropriate. Compare with raw Y too.}}")
  }
  if (length(S$zx)) warn <- c(warn, "{{수치 X를 분석에 쓴 행의 평균과 표준편차로 표준화했습니다. 계수는 X가 1 SD 늘 때의 변화입니다.||Numeric X were standardized with the mean and SD of the analyzed rows. Coefficients are changes per 1 SD of X.}}")
  rhsx <- c(num, fac)
  rhs <- if (!length(rhsx)) "1" else if (isTRUE(cfg$inter) && length(rhsx) > 1)
    sprintf("(%s)^2", paste(rhsx, collapse = " + ")) else paste(rhsx, collapse = " + ")
  slope <- if (length(cfg$slope)) v(cfg$slope)
  if (kind == "ordinal" && !is.null(slope)) { slope <- NULL; warn <- c(warn, "{{순서형에서는 랜덤 기울기를 지원하지 않아 랜덤 절편만 썼습니다.||Random slopes are not supported for ordinal Y; only random intercepts were used.}}") }
  if (!is.null(slope) && !(slope %in% rhsx)) stop("{{랜덤 기울기 변수를 X에도 넣어 주십시오.||The random slope variable must also be an X.}}")
  ret <- if (length(re)) paste(sprintf("(1%s | %s)", if (is.null(slope)) "" else paste(" +", slope), re), collapse = " + ")
  S$tr <- tr; S$slope <- slope; S$off <- off; S$rhs <- rhs; S$re <- re; S$yv <- yv; S$kind <- kind; S$d <- d; S$y <- y
  lhs <- if (kind == "trials") sprintf("cbind(%s, %s - %s)", yv, tr, yv) else yv
  S$f <- paste(lhs, "~", paste(c(rhs, ret, if (!is.null(off)) sprintf("offset(log(%s))", off)), collapse = " + "))
  S$L <- candidates(kind, y, length(re) > 0, isTRUE(cfg$links), !is.null(off)); S$fits <- list(); dis <- list()
  if (S$zy) { dis <- lapply(S$L, `[[`, "name"); dis <- dis[dis != "Gaussian"]
    S$L <- list(cand("Gaussian", "gaussian()", link = "identity")) }
  X <- model.matrix(as.formula(paste("~", rhs)), d); p <- ncol(X); n <- nrow(d)
  if (qr(X)$rank < p) stop("{{X끼리 완전히 겹치는 조합이 있어 계수를 추정하기 어렵습니다. 변수나 교호작용을 줄여 주십시오.||Some X are perfectly collinear, so coefficients cannot be estimated. Remove variables or interactions.}}")
  if (n0 > n) warn <- c(warn, bf("{{결측 때문에 %d행을 뺐습니다. 모든 분포를 같은 %d행으로 비교합니다.||%d rows were dropped for missing values. All distributions are compared on the same %d rows.}}", n0 - n, n))
  eff <- if (kind == "binary") min(table(y)) else n
  if (eff / p < 10) warn <- c(warn, bf("{{데이터가 적습니다 (유효 관측 %d개, 고정효과 계수 %d개). 분포 순위가 불안정할 수 있습니다.||Small data (%d effective observations, %d fixed coefficients). The ranking may be unstable.}}", as.integer(eff), p))
  if (n > 5000) warn <- c(warn, "{{표본이 커서 작은 어긋남도 유의하게 나옵니다. p값보다 그래프와 어긋난 크기를 보십시오.||With a large sample even tiny misfit is significant. Look at the plots and the size of the misfit rather than p-values.}}")
  for (k in re) {
    if (nlevels(d[[k]]) < 5) warn <- c(warn, bf("{{랜덤효과 '%s'의 집단이 %d개뿐이라 분산 추정이 불안정합니다.||Random effect '%s' has only %d groups, so its variance estimate is unstable.}}", lab(k), nlevels(d[[k]])))
    if (nlevels(d[[k]]) == n) warn <- c(warn, bf("{{'%s'는 행마다 값이 다릅니다. 관측 단위 랜덤효과라 잔차 분산과 구분되지 않을 수 있습니다.||'%s' differs on every row. As an observation-level random effect it may not be separable from residual variance.}}", lab(k)))
  }
  if (length(re) && kind == "pos" && !S$zy) warn <- c(warn, "{{Inverse Gaussian은 랜덤효과가 있으면 같은 기준의 우도를 얻기 어려워 뺐습니다.||Inverse Gaussian was left out: with random effects a likelihood on the same footing is hard to obtain.}}")
  if (!length(re)) warn <- c(warn, "{{랜덤효과가 없습니다. 한 사람/집단에서 여러 번 측정했다면 ID를 랜덤효과로 넣는 것을 권합니다.||No random effect. If the same person or group was measured repeatedly, add its ID as a random effect.}}")
  if (kind == "binary") warn <- c(warn, "{{0/1 자료의 분포는 이항분포로 정해져 있습니다. 여기서 비교하는 것은 링크 함수이고, 분포 진단으로는 문제가 잘 드러나지 않을 수 있습니다.||0/1 data are always binomial. What is compared here is the link function, and distribution diagnostics detect little.}}")
  if (kind == "ordinal") {
    gap <- setdiff(seq(min(y), max(y)), y)
    if (int && length(gap)) warn <- c(warn, bf("{{중간 단계 %s은(는) 응답이 없어서, 나온 단계끼리 이어서 계산합니다.||Level(s) %s were never chosen; the observed levels are treated as adjacent.}}", paste(gap, collapse = ", ")))
    for (k in fac) { one <- names(which(tapply(y, d[[k]], function(v) all(v == min(y)) || all(v == max(y)))))
      if (length(one)) warn <- c(warn, bf("{{'%s'의 %s은(는) 응답이 모두 척도 끝값이라 계수와 표준오차가 매우 커질 수 있습니다.||In '%s', %s has only end-of-scale responses, so its coefficient and standard error may blow up.}}", lab(k), paste(one, collapse = ", "))) }
    if (int && length(fac)) { part <- unique(unlist(lapply(fac, function(k) names(which(tapply(y, d[[k]], function(v) length(unique(v))) < length(unique(y)))))))
      if (length(part)) warn <- c(warn, "{{일부 문항이나 집단에서는 안 쓰인 단계가 있습니다. 정상이고, 전체 데이터의 단계를 함께 써서 계산합니다.||Some items or groups do not use every level. That is fine: the levels of the whole data are shared.}}") }
  }
  if (kind == "ordinal") warn <- c(warn, "{{누적 모형의 비례 가정은 이 앱이 검사하지 않습니다. ordinal::nominal_test 등으로 따로 확인하십시오.||This app does not test the proportionality assumption of cumulative models. Check it separately, e.g. with ordinal::nominal_test.}}")
  if (kind == "pos" && int) warn <- c(warn, "{{Y가 정수입니다. 횟수라면 유형을 카운트로 바꾸고, 밀리초처럼 반올림된 측정값이면 그대로 두십시오.||Y is integer-valued. If it is a count, switch the type to Count; if it is a rounded measure such as milliseconds, keep it.}}")
  js(list(n = n, warn = I(warn), formula = S$f, disabled = I(dis), cands = I(lapply(S$L, `[`, c("name", "note", "cmp")))))
}

catchw <- function(expr) {
  w <- character()
  v <- tryCatch(withCallingHandlers(expr, warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") }),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "fiterr"))
  list(v = v, w = unique(w))
}

fit_model <- function(k) {
  f <- as.formula(S$f); d <- S$d
  if (k$engine == "tmb")
  {
    a <- list(formula = f, data = d, family = eval(parse(text = k$fam)), ziformula = as.formula(k$zi))
    if (k$link %in% c("identity", "inverse") && grepl("Gamma", k$fam)) {
      p <- ncol(model.matrix(as.formula(paste("~", S$rhs)), d))
      a$start <- list(beta = c(if (k$link == "inverse") 1 / mean(S$y) else mean(S$y), rep(0, p - 1)))
    }
    do.call(glmmTMB, a)
  }
  else if (k$engine == "ln") {
    # 표준 로그정규. glmmTMB의 lognormal()은 원척도 SD가 상수라 다른 모형이라서 쓰지 않음
    d <- d[S$y > 0, , drop = FALSE]
    f2 <- as.formula(sub("^[^~]*~", sprintf("log(%s) ~", S$yv), S$f))
    do.call(glmmTMB, list(formula = f2, data = d, family = gaussian()))
  } else if (k$engine == "glm")
    glm(f, data = d, family = eval(parse(text = k$fam)), control = glm.control(maxit = 100))
  else {
    d[[S$yv]] <- ordered(d[[S$yv]])
    if (length(S$re)) do.call(ordinal::clmm, list(formula = f, data = d, link = k$fam))
    else do.call(ordinal::clm, list(formula = f, data = d, link = k$fam))
  }
}

sim_model <- function(m, k) {
  n <- nrow(S$d)
  if (k$engine == "tmb") {
    s <- simulate(m, nsim = NSIM, seed = 1)
    sim <- sapply(s, function(x) if (is.matrix(x)) x[, 1] else as.numeric(x))
    list(sim = sim, mu = as.numeric(predict(m, type = "response", re.form = NA)))
  } else if (k$engine == "ln") {
    set.seed(1); d <- S$d; X <- model.matrix(as.formula(paste("~", S$rhs)), d)
    e0 <- as.numeric(X %*% fixef(m)$cond); if (!is.null(S$off)) e0 <- e0 + log(d[[S$off]])
    Z <- model.matrix(as.formula(paste("~ 1", if (!is.null(S$slope)) paste("+", S$slope))), d)
    sg <- sigma(m); pz <- mean(S$y == 0); vc <- VarCorr(m)$cond
    sim <- sapply(seq_len(NSIM), function(i) {
      e <- e0
      for (g in S$re) { V <- matrix(as.numeric(vc[[g]]), ncol(Z)); R <- chol(V + diag(1e-10, ncol(Z))); G <- d[[g]]
        b <- matrix(rnorm(nlevels(G) * ncol(Z)), nlevels(G)) %*% R; e <- e + rowSums(Z * b[as.integer(G), , drop = FALSE]) }
      exp(rnorm(n, e, sg)) * (runif(n) >= pz)
    })
    list(sim = sim, mu = exp(e0))
  } else if (k$engine == "glm") {
    set.seed(1); mu <- as.numeric(fitted(m)); phi <- deviance(m) / n
    list(sim = matrix(statmod::rinvgauss(n * NSIM, mean = rep(mu, NSIM), dispersion = phi), n), mu = mu)
  } else {
    set.seed(1); X <- model.matrix(as.formula(paste("~", S$rhs)), S$d)[, -1, drop = FALSE]
    b <- m$beta; eta0 <- if (length(b)) as.numeric(X[, names(b), drop = FALSE] %*% b) else numeric(n)
    th <- m$alpha; lv <- sort(unique(S$y))
    P <- switch(k$fam, logit = plogis, probit = pnorm, cloglog = function(x) 1 - exp(-exp(x)))
    sim <- sapply(seq_len(NSIM), function(i) {
      e <- eta0
      for (g in S$re) { G <- S$d[[g]]; sdv <- sqrt(ordinal::VarCorr(m)[[g]][1]); e <- e + rnorm(nlevels(G), 0, sdv)[G] }
      u <- runif(n); lv[1 + rowSums(u > P(outer(-e, th, "+")))]
    })
    list(sim = sim, mu = eta0)
  }
}

diagnose <- function(y, sim, mu) {
  n <- length(y); set.seed(2)
  if (any(!is.finite(sim))) stop("{{시뮬레이션 값이 유한하지 않음||simulated values are not finite}}")
  lo <- rowSums(sim < y); eq <- rowSums(sim == y)
  u <- (lo + runif(n) * (eq + 1)) / (NSIM + 1)
  ks <- suppressWarnings(ks.test(u, "punif"))
  two <- function(s, o) min(1, 2 * min(mean(c(s, o) >= o), mean(c(s, o) <= o)))
  c0 <- rowMeans(sim); v0 <- sum((y - c0)^2); vs <- colSums((sim - c0)^2)
  r <- v0 / mean(vs); pd <- two(vs, v0)
  zero <- if (any(y == 0)) { z0 <- sum(y == 0); zs <- colSums(sim == 0)
    list(obs = z0, exp = mean(zs), p = two(zs, z0)) }
  vary <- length(unique(round(mu, 10))) > 1
  hs <- if (vary) suppressWarnings(cor.test(abs(u - .5), mu, method = "spearman", exact = FALSE))
  tr <- if (vary) suppressWarnings(cor.test(u, mu, method = "spearman", exact = FALSE))
  D <- unname(ks$statistic); iss <- character(); lev <- 0
  flag <- function(l, t) { lev <<- max(lev, l); iss <<- c(iss, t) }
  if (ks$p.value < .05 && D > .04) flag(if (D > .1) 2 else 1, bf("{{분포 모양 어긋남 (D=%.3f)||shape misfit (D=%.3f)}}", D))
  if (pd < .05 && (r > 1.2 || r < .8)) flag(if (r > 1.5 || r < .67) 2 else 1, bf(if (r > 1) "{{과산포 (관측/기대 분산 %.2f)||overdispersion (observed/expected variance %.2f)}}" else "{{과소산포 (관측/기대 분산 %.2f)||underdispersion (observed/expected variance %.2f)}}", r))
  if (!is.null(zero) && zero$p < .05 && abs(zero$obs - zero$exp) > .03 * n)
    flag(if (abs(zero$obs - zero$exp) > .1 * n) 2 else 1, bf("{{0의 개수 %d, 모델 기대 %.0f||%d zeros, model expects %.0f}}", as.integer(zero$obs), zero$exp))
  if (!is.null(hs) && hs$p.value < .01 && abs(hs$estimate) > .1) flag(1, "{{예측값에 따라 퍼짐이 달라짐||spread changes with predicted value}}")
  if (!is.null(tr) && tr$p.value < .01 && abs(tr$estimate) > .15) flag(1, "{{잔차에 추세가 남음||trend left in residuals}}")
  # 정규 척도 Q-Q: 무작위화 분위수 잔차. 양 끝은 남기고 가운데만 솎음
  qi <- if (n > 200) sort(unique(c(1:25, round(seq(26, n - 25, length.out = 150)), (n - 24):n))) else seq_len(n)
  id <- if (n > 400) sort(sample(n, 400)) else seq_len(n)
  disc <- all(abs(y - round(y)) < 1e-9) && diff(range(y)) <= 25
  int <- all(abs(y - round(y)) < 1e-9)
  br <- if (disc) seq(min(y) - .5, max(y) + .5, 1) else if (int) {
    q99 <- unname(quantile(y, .99)); bw <- max(1, ceiling((q99 - min(y) + 1) / 24))
    seq(min(y) - .5, by = bw, length.out = max(2, ceiling((q99 - min(y) + 1) / bw)) + 1) } else {
    a <- unname(quantile(y, c(0, .99))); if (a[1] == a[2]) a <- range(y); seq(a[1], a[2], length.out = 25) }
  nb <- length(br) - 1
  h <- function(x) tabulate(pmin(nb, pmax(1, findInterval(x, br, rightmost.closed = TRUE))), nb) / n
  H <- apply(sim, 2, h)
  list(level = lev, issues = I(iss), D = D, ksp = ks$p.value, ratio = r, pd = pd, zero = zero,
    pp = list(x = qnorm((qi - .5) / n), y = qnorm(sort(u)[qi]), lo = qnorm(qbeta(.025, qi, n + 1 - qi)), hi = qnorm(qbeta(.975, qi, n + 1 - qi))),
    rv = list(x = rank(mu)[id] / n, y = u[id]),
    hist = list(cut = max(y) > max(br), br = br, obs = I(h(y)), mid = I(apply(H, 1, median)), lo = I(apply(H, 1, quantile, .05)), hi = I(apply(H, 1, quantile, .95))))
}

fit_one <- function(i) {
  k <- S$L[[i]]; a <- catchw(fit_model(k)); m <- a$v
  base <- list(i = i, name = k$name, note = k$note, cmp = k$cmp, link = k$link)
  if (inherits(m, "fiterr")) return(js(c(base, list(ok = FALSE, msg = m$msg))))
  ok <- TRUE; msg <- character()
  if (k$engine %in% c("tmb", "ln")) {
    if (m$fit$convergence != 0) { ok <- FALSE; msg <- c(msg, "{{수렴하지 않음||did not converge}}") }
    if (!isTRUE(m$sdr$pdHess)) { ok <- FALSE; msg <- c(msg, "{{헤시안이 양의 정부호가 아님 (추정 불안정)||Hessian not positive definite (unstable estimates)}}") }
    if (length(S$re)) {
      vc <- VarCorr(m)$cond
      sds <- unlist(lapply(vc, function(x) attr(x, "stddev")))
      cr <- unlist(lapply(vc, function(x) { q <- attr(x, "correlation"); q[upper.tri(q)] }))
      sc <- if (k$link == "identity" && k$engine == "tmb") sd(S$y) else 1
      if (any(sds < 1e-3 * sc) || any(abs(cr) > .99)) msg <- c(msg, "{{singular 가능성: 랜덤효과 분산이 0에 가깝거나 상관이 ±1||possibly singular: random-effect variance near 0 or correlation near ±1}}")
    }
  } else if (k$engine == "glm") { if (!m$converged) { ok <- FALSE; msg <- c(msg, "{{수렴하지 않음||did not converge}}") } }
  else {
    if (inherits(m, "clmm")) { if (any(!is.finite(suppressWarnings(sqrt(diag(vcov(m))))))) { ok <- FALSE; msg <- c(msg, "{{표준오차 계산 실패 (추정 불안정)||standard errors could not be computed (unstable)}}") } }
    else if (m$convergence$code > 0) { ok <- FALSE; msg <- c(msg, "{{수렴하지 않음||did not converge}}") }
  }
  l0 <- logLik(m); df <- attr(l0, "df"); n <- nrow(S$d); ll <- as.numeric(l0)
  if (k$engine == "ln") {
    ll <- ll - sum(log(S$y[S$y > 0]))
    if (k$zi == "~1") { nz <- sum(S$y == 0); ll <- ll + nz * log(nz / n) + (n - nz) * log(1 - nz / n); df <- df + 1 }
  }
  if (!is.finite(ll)) { ok <- FALSE; msg <- c(msg, "{{우도를 계산하지 못함||likelihood could not be computed}}") }
  S$fits[[as.character(i)]] <- m
  out <- c(base, list(ok = ok, msg = paste(msg, collapse = " / "), rwarn = paste(a$w, collapse = " / "),
    ll = ll, df = df, aic = -2 * ll + 2 * df,
    aicc = if (n > df + 1) -2 * ll + 2 * df + 2 * df * (df + 1) / (n - df - 1) else NA,
    bic = -2 * ll + log(n) * df))
  if (ok) {
    g <- catchw({ s <- sim_model(m, k); diagnose(S$y, s$sim, s$mu) })
    if (inherits(g$v, "fiterr")) out$diagerr <- g$v$msg else out$diag <- g$v
  }
  js(out)
}

coefs <- function(i) {
  m <- S$fits[[as.character(i)]]; k <- S$L[[i]]
  co <- if (k$engine %in% c("tmb", "ln")) summary(m)$coefficients else list(cond = summary(m)$coefficients)
  rows <- do.call(c, lapply(c("cond", "zi"), function(w) {
    A <- co[[w]]; if (is.null(A)) return(NULL)
    cz <- if (k$engine == "glm") qt(.975, m$df.residual) else 1.96
    lapply(seq_len(nrow(A)), function(j) list(part = w, term = rownames(A)[j], est = A[j, 1], se = A[j, 2],
      lo = A[j, 1] - cz * A[j, 2], hi = A[j, 1] + cz * A[j, 2], p = A[j, 4]))
  }))
  js(list(rows = I(rows), link = k$link, name = k$name))
}

code_for <- function(i, file, sep = ",", enc = "UTF-8") {
  k <- S$L[[i]]
  q <- function(v) paste0("`", gsub("`", "\\\\`", lab(v)), "`")
  bt <- function(s) { for (j in seq_along(S$labels)) s <- gsub(sprintf("vv%d__", j), q(sprintf("vv%d__", j)), s, fixed = TRUE); s }
  used <- names(S$d); fac <- used[sapply(S$d, is.factor)]
  o <- c(sprintf("# %s: {{분포 비교 앱에서 선택한 모델||model chosen in GLM Distribution Explorer}}", k$name), "# {{준비(처음 한 번)||setup (first time only)}}: install.packages(c(\"glmmTMB\", \"DHARMa\", \"ordinal\"))", "library(glmmTMB); library(DHARMa)",
    sprintf("d <- read.csv(%s, sep = %s, check.names = FALSE, fileEncoding = \"%s\",", encodeString(file, quote = "\""), encodeString(sep, quote = "\""), if (enc == "UTF-8") "UTF-8-BOM" else "CP949"),
    "              na.strings = c(\"\", \"NA\", \"NaN\", \"N/A\", \".\"), strip.white = TRUE)",
    sprintf("used <- c(%s)", paste(encodeString(lab(used), quote = "\""), collapse = ", ")),
    "d <- d[complete.cases(d[used]), ]")
  if (length(fac)) o <- c(o, sprintf("d$%1$s <- factor(d$%1$s)", sapply(fac, q)))
  for (z in c(S$zx, if (S$zy) S$yv)) o <- c(o, sprintf("d$%1$s <- as.numeric(scale(d$%1$s))", q(z)))
  if (S$kind == "binary") o <- c(o, sprintf("d$%1$s <- as.integer(d$%1$s == sort(unique(d$%1$s))[2])", q(S$yv)))
  if (k$engine == "tmb") o <- c(o, sprintf("m <- glmmTMB(%s, data = d,\n             family = %s, ziformula = %s)", bt(S$f), k$fam, paste0(k$zi, if (k$link %in% c("identity", "inverse") && grepl("Gamma", k$fam))
    sprintf(",\n             start = list(beta = c(%s, rep(0, %d)))", if (k$link == "inverse") sprintf("1 / mean(d$%s)", q(S$yv)) else sprintf("mean(d$%s)", q(S$yv)),
      ncol(model.matrix(as.formula(paste("~", S$rhs)), S$d)) - 1L))),
    "summary(m)", "plot(simulateResiduals(m, n = 1000))", "testDispersion(m)")
  else if (k$engine == "ln") o <- c(o, if (k$zi == "~1") c("d0 <- d", sprintf("d <- d[d$%s > 0, ]  # {{양수 부분만||positive part only}}", q(S$yv)), sprintf("m_zero <- glm((%s == 0) ~ 1, data = d0, family = binomial())  # {{0이냐 아니냐 부분||zero vs non-zero part}}", q(S$yv))),
    sprintf("m <- glmmTMB(%s, data = d)", bt(sub("^[^~]*~", sprintf("log(%s) ~", S$yv), S$f))), "summary(m)", "plot(simulateResiduals(m, n = 1000))",
    sprintf("AIC(m) + 2 * sum(log(d$%s))%s  # {{원척도로 보정한 AIC||AIC corrected to the raw scale}}", q(S$yv), if (k$zi == "~1") " + AIC(m_zero)" else ""))
  else if (k$engine == "glm") o <- c(o, sprintf("m <- glm(%s, data = d, family = %s)", bt(S$f), k$fam), "summary(m)")
  else o <- c(o, "library(ordinal)", sprintf("d$%1$s <- ordered(d$%1$s)", q(S$yv)),
    sprintf("m <- %s(%s, data = d, link = \"%s\")", if (length(S$re)) "clmm" else "clm", bt(S$f), k$fam), "summary(m)")
  paste(o, collapse = "\n")
}

demo_csv <- function(w) {
  set.seed(42)
  if (w == "rt") {
    d <- expand.grid(trial = 1:12, condition = c("A", "B", "C"), participant = sprintf("P%02d", 1:24))
    b <- rnorm(24, 0, .2); mu <- exp(6.3 + c(0, .12, .25)[d$condition] + b[d$participant])
    d$rt <- round(rgamma(nrow(d), 6, 6 / mu), 2)
  }
  if (w == "count") {
    d <- data.frame(group = rep(c("control", "treat"), 150), age = round(runif(300, 18, 65)))
    d$errors <- rnbinom(300, size = 1.2, mu = exp(1 + .5 * (d$group == "treat") + .01 * (d$age - 40)))
  }
  if (w == "likert") {
    d <- expand.grid(ui = c("old", "new"), item = 1:3, participant = sprintf("P%02d", 1:40))
    e <- .9 * (d$ui == "new") + rnorm(40, 0, 1)[d$participant] + rlogis(nrow(d))
    d$score <- findInterval(e, c(-2, -.7, .6, 2)) + 1
  }
  paste(capture.output(write.csv(d, row.names = FALSE, quote = FALSE)), collapse = "\n")
}
