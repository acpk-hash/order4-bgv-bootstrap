# 全 baseline 参数实验汇总

本报告由本地真实日志生成。`total/extract/linear1/linear2` 来自密文自举日志的 `time for ...` 行；`pass=True` 表示日志包含 `### bts finished, everything ok ###`。

## 1. Ma-large-p 同参数实验

这部分是与我们的实现最严格可比的结果：同一个 HElib 大 p 自举路径、同一个 `fatboot` 参数、同一个实际日志解析得到的溢出界 `B`。

| 参数 | p | h | B | type-B 可行 | Ma type-A total | Ma type-B total | ours sparse total | ours order4 total | 主要结论 |
|---|---:|---:|---:|---|---:|---:|---:|---:|---|
| A | 17 | 14 | 18 | False | 43.298 | -- | -- | -- | type-B support 不可注入，不能做同支持 cleaner 对比 |
| B | 127 | 12 | 17 | False | 89.127 | -- | -- | -- | type-B support 不可注入，不能做同支持 cleaner 对比 |
| C | 257 | 12 | 17 | False | 87.408 | -- | -- | -- | type-B support 不可注入，不能做同支持 cleaner 对比 |
| D | 8191 | 12 | 17 | True | 53.348 | 170.470 | 167.317 | -- | sparse 相对 type-B total 提升 1.85% |
| E | 65537 | 12 | 17 | True | 62.706 | 207.442 | 194.446 | 132.236 | order4 相对 type-B total 提升 36.25% |

| 参数 | Ma type-B extract | ours sparse extract | ours order4 extract | sparse extract speedup | order4 extract speedup |
|---|---:|---:|---:|---:|---:|
| A | -- | -- | -- | -- | -- |
| B | -- | -- | -- | -- | -- |
| C | -- | -- | -- | -- | -- |
| D | 129.200 | 126.193 | -- | 2.33 | -- |
| E | 156.520 | 143.937 | 81.391 | 8.04 | 48.00 |

| 参数 | generic A / deg / terms | our chosen A / deg / terms | decomposition |
|---|---|---|---|
| A | A=37, deg=None, terms=None | A=None, deg=None, terms=None | d=None, est_mults=None |
| B | A=35, deg=None, terms=None | A=None, deg=None, terms=None | d=None, est_mults=None |
| C | A=35, deg=None, terms=None | A=None, deg=None, terms=None | d=None, est_mults=None |
| D | A=35, deg=1223, terms=612 | A=45, deg=1223, terms=611 | d=2, est_mults=55 |
| E | A=35, deg=1223, terms=612 | A=256, deg=1223, terms=307 | d=4, est_mults=41 |

实际多项式系数和可读 TeX 形式已导出到 `docs/cleaner_polynomial_forms_p65537_B17.json` 和 `docs/cleaner_polynomial_forms_p65537_B17.tex`。

## 2. Homomorphic-NTT 文章参数

这部分跑的是 NTT 文章自己的 optimized/baseline 对照。它使用另一个 patched HElib executor，因此目前不能直接声称我们的 order-four cleaner 已经接入这条密文路径。

| 参数 | p | h | B | variant | status | key-independent done(s) | extract | total | pass |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| I | 65537 | 26 | 13 | ntt_article_baseline | timeout | -- | -- | -- | False |
| I | 65537 | 26 | 13 | ntt_article_optimized | ok | 712.098 | 5.472 | 20.564 | True |
| II | 8191 | 24 | 12 | ntt_article_baseline | ok | 371.710 | 5.670 | 60.543 | True |
| II | 8191 | 24 | 12 | ntt_article_optimized | ok | 79.552 | 5.872 | 21.567 | True |
| III | 131071 | 26 | 13 | ntt_article_baseline | ok | 1611.070 | 5.924 | 178.141 | True |
| III | 131071 | 26 | 13 | ntt_article_optimized | ok | 1770.870 | 5.220 | 22.584 | True |

| 参数 | NTT optimized total | NTT baseline total | total speedup | 备注 |
|---|---:|---:|---:|---|
| I | 20.564 | -- | -- | baseline timeout |
| II | 21.567 | 60.543 | 64.38 | same article comparison |
| III | 22.584 | 178.141 | 87.32 | same article comparison |

## 3. Bootstrapping_Polyfunctions 本地可运行项

已执行 `polyfunctions_toy_repeat1.log`。该程序是 toy 参数，不是 Ma-large-p 的同参数集合；日志中有 3 个 successful marker。生成新 digit polynomials 仍需要 Magma。

## 4. 可以写进论文的边界

- 我们的主 claim 应放在 Ma-large-p type-B exact cleaner 路径下：D/E 是同参数、同实现路径、同密文自举可比。
- A/B/C 的真实 `B` 使 `(2B+1)^2 < p` 不成立，不能把 toy 支持下搜到的 cleaner 拿来声称同参数可用。
- E 参数 `p=65537` 有四阶结构，`A=256` 使项数从 612 降到 307，并触发 order-four evaluator，这是当前最强实验证据。
- NTT 文章的 optimized 路线验证了同态线性变换优化的重要性，但它是另一个 executor；在未移植前只能作为外部 benchmark。
- Polyfunctions 库本地只能跑 toy；新参数多项式生成依赖 Magma，因此不能在本机完成同参数大 p 再生成。

