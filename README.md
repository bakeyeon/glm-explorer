# GLM Distribution Explorer

**Open the app: https://bakeyeon.github.io/glm-explorer/**

Fit every plausible GLM/GLMM distribution to your data at once, compare them visually, and copy ready-to-run R and Python code. English and Korean.

> **This is a quick first look, not a final analysis.** It was made to run a rough pass over your data early on. Use it to see which distributions look plausible, then write analysis code that fits your own data and design. The exported code is a starting point, not something to paste into a paper as is.

Everything runs in your browser with [webR](https://docs.r-wasm.org/webr/latest/). Your data file is never uploaded. The first load takes about a minute because R and its packages are downloaded from the webR CDN.

## What it does

- Y types: real, positive, non-negative with zeros, count, proportion, binary, successes/trials, ordinal
- Suggests variable roles, and converts wide survey data to long format
- Random intercepts and slopes, interactions, exposure offsets, optional standardization
- ΔAIC/AICc/BIC ranking, Q-Q plots side by side, residual plots, simulated vs observed Y
- Coefficients and p-values appear only after you choose a distribution
- Exports R (`glmmTMB`, `ordinal`) and Python (`statsmodels`, or `rpy2` when there is no reliable equivalent)

## Limits

- This is for exploration. Decide from your design and theory first, and recheck the final model in R before publishing.
- The diagnostic grades are rules of thumb chosen for this app.
- AIC ranks candidates relative to each other; the best one can still fit badly.
- `Lognormal` means a normal model on log(Y) with AIC corrected to the raw scale. It is not `glmmTMB::lognormal()`, which has a constant SD on the raw scale.
- Fits were checked against simulated data and base R/lme4/MASS. The generated `rpy2` wrappers and DHARMa lines have not been run.

## Development

Edit `engine.R` (fitting), `py.R` (Python export) or `app.src.html` (page), then run `python build.py` to rebuild `index.html`. Preview with `python -m http.server`.

## License

Code in this repository is MIT. R and R packages are not bundled; they are downloaded at run time and keep their own licenses (including GPL and AGPL).

---

## 한국어

여러 GLM/GLMM 분포를 한 번에 돌려 그래프로 비교하고, R과 Python 코드를 바로 복사할 수 있는 웹앱입니다. 화면 오른쪽 위에서 한국어로 바꿀 수 있습니다. 데이터 파일은 브라우저 안에서만 계산되고 업로드되지 않습니다.

> **분석 초기에 대충 한번 돌려보는 용도로 만들었습니다.** 어떤 분포가 그럴듯한지 대강 보시고, 실제 분석은 데이터와 연구 설계에 맞게 코드를 다시 짜세요. 내보낸 코드는 출발점일 뿐, 그대로 논문에 쓰라고 만든 것이 아닙니다.
