#!/usr/bin/env python3
"""Export concrete bounded-support cleaner polynomial forms.

The generated JSON keeps the complete coefficient list.  The generated TeX
keeps the readable part used in the short theory note.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent))

from sat_cleaning_circuit_search import (  # noqa: E402
    centered_mod,
    degree,
    newton_interpolation_coeffs,
    newton_to_monomial,
    polynomial_eval,
    ps_cost_for_degree,
    support_table,
)


def latex_signed_coeff(c: int) -> str:
    if c == 1:
        return "+"
    if c == -1:
        return "-"
    if c >= 0:
        return f"+{c}"
    return str(c)


def latex_monomial(term: dict, variable: str = "X") -> str:
    exp = term["exponent"]
    coeff = term["coeff_centered"]
    sign = latex_signed_coeff(coeff)
    if exp == 0:
        body = ""
    elif exp == 1:
        body = variable
    else:
        body = f"{variable}^{{{exp}}}"
    if abs(coeff) == 1 and body:
        return f"{sign}{body}"
    return f"{sign}{body}" if body else sign


def join_terms_for_latex(terms: list[dict], variable: str = "X") -> str:
    if not terms:
        return "0"
    text = " ".join(latex_monomial(t, variable=variable) for t in terms)
    return text[1:] if text.startswith("+") else text


def polynomial_for_aux(p: int, B: int, aux: int) -> dict:
    points, target, _support, collisions = support_table(p, B, aux)
    if collisions:
        return {
            "aux": aux,
            "valid": False,
            "reason": "support map is not injective",
            "collisions": collisions,
        }

    poly = newton_to_monomial(points, newton_interpolation_coeffs(points, target, p), p)
    terms = [
        {
            "exponent": i,
            "coeff_mod_p": c % p,
            "coeff_centered": centered_mod(c, p),
        }
        for i, c in enumerate(poly)
        if c % p
    ]
    verified = all(polynomial_eval(poly, x, p) == y % p for x, y in zip(points, target))
    residues_mod4 = {
        str(r): sum(1 for t in terms if t["exponent"] % 4 == r) for r in range(4)
    }

    result = {
        "aux": aux,
        "valid": True,
        "support_size": len(points),
        "degree": degree(poly),
        "nonzero_terms": len(terms),
        "leading_coeff_mod_p": terms[-1]["coeff_mod_p"],
        "leading_coeff_centered": terms[-1]["coeff_centered"],
        "terms": terms,
        "first_terms": terms[:12],
        "last_terms": terms[-12:],
        "residue_counts_mod4": residues_mod4,
        "verified_on_support": verified,
        "paterson_stockmeyer_cost_for_full_degree": ps_cost_for_degree(degree(poly)),
    }

    if (aux * aux + 1) % p == 0:
        order4_ok = True
        q_terms = []
        linear = 0
        for term in terms:
            exp = term["exponent"]
            if exp == 1:
                linear = term["coeff_mod_p"]
            elif exp >= 3 and exp % 4 == 3:
                q_terms.append(
                    {
                        "exponent": (exp - 3) // 4,
                        "coeff_mod_p": term["coeff_mod_p"],
                        "coeff_centered": term["coeff_centered"],
                    }
                )
            else:
                order4_ok = False
        identity_ok = all(
            (
                polynomial_eval(poly, (aux * x) % p, p)
                + aux * polynomial_eval(poly, x, p)
                - aux * x
            )
            % p
            == 0
            for x in points
        )
        result["order_four"] = {
            "aux_squared_is_minus_one": True,
            "linear_coeff_mod_p": linear,
            "linear_coeff_centered": centered_mod(linear, p),
            "Q_degree": max((t["exponent"] for t in q_terms), default=-1),
            "Q_nonzero_terms": len(q_terms),
            "Q_first_terms": q_terms[:12],
            "Q_last_terms": q_terms[-12:],
            "Q_paterson_stockmeyer_cost": ps_cost_for_degree(
                max((t["exponent"] for t in q_terms), default=-1)
            ),
            "shape_verified": order4_ok,
            "identity_P_AX_plus_A_P_X_equals_AX_verified": identity_ok,
            "estimated_ciphertext_multiplications": {
                "generic_PS_for_P_degree": ps_cost_for_degree(degree(poly))["mults"],
                "order_four_path": ps_cost_for_degree(
                    max((t["exponent"] for t in q_terms), default=-1)
                )["mults"]
                + 4,
                "notes": "x^2, x^3, x^4, Q(x^4), then multiply by x^3; linear term uses no ciphertext-ciphertext multiplication.",
            },
        }
    return result


def write_tex(result: dict, path: Path) -> None:
    aux35 = next(row for row in result["polynomials"] if row["aux"] == 35)
    aux256 = next(row for row in result["polynomials"] if row["aux"] == 256)
    order4 = aux256["order_four"]
    replacements = {
        "@@AUX35_DEGREE@@": str(aux35["degree"]),
        "@@AUX35_TERMS@@": str(aux35["nonzero_terms"]),
        "@@AUX256_DEGREE@@": str(aux256["degree"]),
        "@@AUX256_TERMS@@": str(aux256["nonzero_terms"]),
        "@@P35_FIRST@@": join_terms_for_latex(aux35["terms"][:4]),
        "@@P35_LAST@@": join_terms_for_latex(aux35["terms"][-4:]),
        "@@P256_FIRST@@": join_terms_for_latex(aux256["terms"][:4]),
        "@@P256_LAST@@": join_terms_for_latex(aux256["terms"][-4:]),
        "@@Q_FIRST@@": join_terms_for_latex(order4["Q_first_terms"][:4], variable="Y"),
        "@@Q_LAST@@": join_terms_for_latex(order4["Q_last_terms"][-4:], variable="Y"),
        "@@GENERIC_MULTS@@": str(
            order4["estimated_ciphertext_multiplications"]["generic_PS_for_P_degree"]
        ),
        "@@ORDER4_MULTS@@": str(
            order4["estimated_ciphertext_multiplications"]["order_four_path"]
        ),
    }
    tex = r"""\section{同一参数下的实际多项式对比}

本节固定实验参数
\[
  p=65537,\quad B=17,\quad
  S_A=\{hA+\ell:\ |h|,|\ell|\le 17\},
  \quad |S_A|=1225.
\]
Ma 等人的 large-\(p\) 实验库在 \(t=-1\) 的 non-power-of-\(p\) 路径中，
默认选择
\[
  A=2\lceil B_{\rm rec}\rceil+1=35,
\]
并在 \(\mathbb F_p\) 上插值得到满足
\[
  P_{35}(35h+\ell)=\ell,\qquad |h|,|\ell|\le 17
\]
的 exact cleaner。我们的 SAT/结构搜索保留同一个 bounded-support
cleaner 问题，但把 auxiliary radix 改成
\[
  A=256,\qquad 256^2\equiv -1\pmod {65537}.
\]
由此得到另一个 exact cleaner \(P_{256}\)，满足
\[
  P_{256}(256h+\ell)=\ell,\qquad |h|,|\ell|\le 17.
\]

\begin{table}[h]
\centering
\resizebox{\linewidth}{!}{%
\begin{tabular}{lrrrll}
\toprule
来源/实现 & aux & degree & nonzero & 结构 & 正确性 \\
\midrule
HElib 原生 digit polynomial & -- & \(p\) 或 \((e-1)(p-1)+1\) & 大 & 全局 \(p\)-进制抽取 & HElib 插值/CH 公式 \\
Ma large-\(p\) 默认库 & 35 & @@AUX35_DEGREE@@ & @@AUX35_TERMS@@ & odd powers & support 全验证 \\
本文 SAT/结构搜索 & 256 & @@AUX256_DEGREE@@ & @@AUX256_TERMS@@ & \(X,X^{4j+3}\) & support 全验证 \\
本文 order-four evaluator & 256 & @@AUX256_DEGREE@@ & @@AUX256_TERMS@@ & \(\frac12X+X^3Q(X^4)\) & 恒等式验证 \\
\bottomrule
\end{tabular}
}
\caption{同一 \(p=65537,B=17\) bounded-support 问题下的多项式形式。}
\end{table}

这里的 ``HElib 原生'' 指 \texttt{buildDigitPolynomial} 或 Chen--Han
\texttt{compute\_magic\_poly} 路径：它解决的是全局 \(p\)-进制 digit
extraction，而不是 bounded support 上的 two-coordinate cleaner。
因此它没有 \(A=35\) 或 \(A=256\) 的 support 几何结构；在大 \(p=65537\)
下直接使用会得到 degree 至少为 \(p\) 量级的多项式。我们的 baseline
实验实际比较的是 Ma large-\(p\) 库中的 bounded-support cleaner。
对应代码路径为
\begin{center}
\scriptsize
\begin{tabular}{ll}
HElib 原生: & \texttt{baselines/HElib/src/extractDigits.cpp}\\
Ma large-\(p\) 默认库: &
\texttt{baselines/BGV-Boot-for-Large-p}\\
Ma homomorphic NTT 库: &
\texttt{baselines/bgv-bootstrapping-with-homomorphic-NTT}\\
本文实现: &
\texttt{src/HElib\_auxradix\_opt/src/extractDigits.cpp}
\end{tabular}
\end{center}

另外，Polyfunctions baseline 对应 GIKV/polyfunction 路线，目录里已有的
\texttt{poly*.txt} 主要是小 \(p\) 和小 \(e\) 的 digit extraction/retain
多项式样例，不是 \(p=65537,B=17\) 的 bounded-support cleaner。Ma 的
homomorphic NTT baseline 中缓存的 \texttt{12289\_3\_ZZX.txt} 也对应另一组参数。
因而在本实验同一参数集下，可直接比较的 Ma bounded-support cleaner 是
默认 \(A=35\) 的 \(P_{35}\)；我们的 SAT/结构搜索把同一问题的 auxiliary
radix 改为 \(A=256\)。

\paragraph{Ma 默认 cleaner 的实际形式。}
记 \(P_{35}\) 为 Ma large-\(p\) 默认路径得到的多项式。它的非零指数为
所有奇数 \(1,3,\ldots,1223\)，共 \(612\) 项。前若干项和最高若干项为
\[
\begin{aligned}
P_{35}(X)
  &= @@P35_FIRST@@+\cdots\\
  &\quad @@P35_LAST@@ \pmod{65537}.
\end{aligned}
\]
完整 \(612\) 项系数见
\texttt{docs/cleaner\_polynomial\_forms\_p65537\_B17.json}。

\paragraph{SAT/结构搜索得到的 cleaner。}
本文搜索得到的 \(P_{256}\) 的实际形式是
\[
  P_{256}(X)=32769X+X^3Q(X^4),
\]
其中 \(32769=2^{-1}\bmod 65537\)，且 \(Q\) 有
degree \(305\)、\(306\) 个非零项。展开到 \(X\) 的多项式时，
\[
\begin{aligned}
P_{256}(X)
  &= @@P256_FIRST@@+\cdots\\
  &\quad @@P256_LAST@@ \pmod{65537}.
\end{aligned}
\]
对应的 \(Q\) 的前后若干项为
\[
\begin{aligned}
Q(Y)
  &= @@Q_FIRST@@+\cdots\\
  &\quad @@Q_LAST@@ \pmod{65537}.
\end{aligned}
\]

\paragraph{正确性。}
两个 cleaner 的正确性都不是抽样测试，而是在全部 \(1225\) 个 support
点上验证：
\[
  P_A(hA+\ell)=\ell,\qquad |h|,|\ell|\le 17.
\]
对于 \(A=256\)，还验证了 order-four 恒等式
\[
  P_{256}(256X)+256P_{256}(X)=256X
\]
在全部 support 上成立；因为两边次数都小于 \(|S_A|=1225\)，这个 support
恒等式也等价于低次数代表元中的多项式恒等式。

\paragraph{节省的 cost。}
从多项式本身看，\(P_{35}\) 到 \(P_{256}\) 把非零单项式数从 \(612\)
降到 \(307\)，减少
\[
  1-\frac{307}{612}=49.84\%.
\]
order-four evaluator 不再直接对 degree \(1223\) 的 \(P\) 做通用求值，
而是先算 \(X^2,X^3,X^4\)，再对 degree \(305\) 的 \(Q(X^4)\) 求值。
按同一 Paterson--Stockmeyer 估算，generic degree-\(1223\) 路径约
@@GENERIC_MULTS@@ 次 ciphertext multiplication，order-four 路径约
@@ORDER4_MULTS@@ 次，理论乘法数减少约
\[
  1-\frac{@@ORDER4_MULTS@@}{@@GENERIC_MULTS@@}
  =43.84\%.
\]
以 \texttt{/home/user/experiments/baselines/BGV-Boot-for-Large-p}
中的 Ma large-\(p\) baseline 为正式对比对象，repeat-3 实测中，
默认 \(A=35\) 的 extract 均值为 \(160.78\)s、total 均值为
\(213.10\)s。切到 \(A=256\) sparse cleaner 后，extract 为
\(147.91\)s、total 为 \(199.40\)s；再启用 order-four evaluator 后，
extract 为 \(80.77\)s、total 为 \(131.18\)s。相对该 baseline，
sparse cleaner 的 extract/total 加速分别为 \(8.01\%\) 和
\(6.43\%\)，order-four 完整路径的 extract/total 加速分别为
\(49.76\%\) 和 \(38.44\%\)。
"""
    for old, new in replacements.items():
        tex = tex.replace(old, new)
    path.write_text(tex + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=65537)
    parser.add_argument("--B", type=int, default=17)
    parser.add_argument("--aux", type=int, nargs="+", default=[35, 256])
    parser.add_argument("--json-output", type=Path, required=True)
    parser.add_argument("--tex-output", type=Path, required=True)
    args = parser.parse_args()

    result = {
        "p": args.p,
        "B": args.B,
        "support_formula": "S_A={hA+lo : |h|,|lo| <= B}",
        "polynomials": [polynomial_for_aux(args.p, args.B, aux) for aux in args.aux],
    }
    args.json_output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_tex(result, args.tex_output)


if __name__ == "__main__":
    main()
