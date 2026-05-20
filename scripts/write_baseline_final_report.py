#!/usr/bin/env python3
"""Write a concise Chinese report for the baseline matrix experiments."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any


ROOT = Path(__file__).resolve().parent.parent

DONE_RE = re.compile(r" done in ([0-9.]+) seconds")


def metric(variant: dict[str, Any], key: str) -> float | None:
    try:
        return float(variant["timing"]["stats"][key]["mean"])
    except Exception:
        return None


def fmt(x: Any, digits: int = 3) -> str:
    if x is None:
        return "--"
    if isinstance(x, float):
        return f"{x:.{digits}f}"
    return str(x)


def pct(base: float | None, opt: float | None) -> float | None:
    if base is None or opt is None or base == 0:
        return None
    return (base - opt) / base * 100.0


def done_times(log_path: str | None) -> list[float]:
    if not log_path:
        return []
    path = Path(log_path)
    if not path.exists():
        return []
    return [float(x) for x in DONE_RE.findall(path.read_text(errors="replace"))]


def variant_brief(v: dict[str, Any]) -> str:
    if v.get("status") in {"skipped", "not_integrated", "not_run"}:
        return v.get("reason", v.get("status", "--"))
    if v.get("status") == "timeout" or v.get("timeout_marker"):
        return "timeout"
    if v.get("status") == "failed":
        return "failed"
    if v.get("passed_marker") is False:
        return "no correctness marker"
    return "ok" if v.get("passed_marker") else str(v.get("status", "--"))


def add_ma_section(lines: list[str], rows: list[dict[str, Any]]) -> None:
    lines.append("## 1. Ma-large-p 同参数实验")
    lines.append("")
    lines.append(
        "这部分是与我们的实现最严格可比的结果：同一个 HElib 大 p 自举路径、同一个 "
        "`fatboot` 参数、同一个实际日志解析得到的溢出界 `B`。"
    )
    lines.append("")
    lines.append("| 参数 | p | h | B | type-B 可行 | Ma type-A total | Ma type-B total | ours sparse total | ours order4 total | 主要结论 |")
    lines.append("|---|---:|---:|---:|---|---:|---:|---:|---:|---|")
    for row in rows:
        variants = row.get("variants", {})
        type_b = variants.get("ma_typeB", {})
        sparse = variants.get("ours_sparse", {})
        order4 = variants.get("ours_order4", {})
        base_total = metric(type_b, "total")
        order_total = metric(order4, "total")
        sparse_total = metric(sparse, "total")
        conclusion = "type-B support 不可注入，不能做同支持 cleaner 对比"
        if base_total is not None:
            if order_total is not None:
                conclusion = f"order4 相对 type-B total 提升 {fmt(pct(base_total, order_total), 2)}%"
            elif sparse_total is not None:
                conclusion = f"sparse 相对 type-B total 提升 {fmt(pct(base_total, sparse_total), 2)}%"
        lines.append(
            f"| {row['label']} | {row['p']} | {row['h']} | {fmt(row.get('B'), 0)} | "
            f"{(row.get('cleaner_search') or {}).get('type_b_feasible')} | "
            f"{fmt(metric(variants.get('ma_typeA', {}), 'total'))} | "
            f"{fmt(base_total)} | {fmt(sparse_total)} | {fmt(order_total)} | {conclusion} |"
        )
    lines.append("")
    lines.append("| 参数 | Ma type-B extract | ours sparse extract | ours order4 extract | sparse extract speedup | order4 extract speedup |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for row in rows:
        variants = row.get("variants", {})
        base_extract = metric(variants.get("ma_typeB", {}), "extract")
        sparse_extract = metric(variants.get("ours_sparse", {}), "extract")
        order4_extract = metric(variants.get("ours_order4", {}), "extract")
        lines.append(
            f"| {row['label']} | {fmt(base_extract)} | {fmt(sparse_extract)} | "
            f"{fmt(order4_extract)} | {fmt(pct(base_extract, sparse_extract), 2)} | "
            f"{fmt(pct(base_extract, order4_extract), 2)} |"
        )
    lines.append("")
    lines.append("| 参数 | generic A / deg / terms | our chosen A / deg / terms | decomposition |")
    lines.append("|---|---|---|---|")
    for row in rows:
        search = row.get("cleaner_search") or {}
        generic = search.get("generic_aux_result") or {}
        chosen = search.get("chosen_practical") or {}
        decomp = chosen.get("best_decomposition") or {}
        lines.append(
            f"| {row['label']} | A={search.get('generic_aux')}, deg={generic.get('degree')}, "
            f"terms={generic.get('nonzero_terms')} | A={chosen.get('aux')}, "
            f"deg={chosen.get('degree')}, terms={chosen.get('nonzero_terms')} | "
            f"d={decomp.get('modulus')}, est_mults={decomp.get('estimated_mults')} |"
        )
    lines.append("")
    lines.append(
        "实际多项式系数和可读 TeX 形式已导出到 "
        "`docs/cleaner_polynomial_forms_p65537_B17.json` 和 "
        "`docs/cleaner_polynomial_forms_p65537_B17.tex`。"
    )
    lines.append("")


def add_ntt_section(lines: list[str], rows: list[dict[str, Any]]) -> None:
    lines.append("## 2. Homomorphic-NTT 文章参数")
    lines.append("")
    lines.append(
        "这部分跑的是 NTT 文章自己的 optimized/baseline 对照。它使用另一个 patched HElib "
        "executor，因此目前不能直接声称我们的 order-four cleaner 已经接入这条密文路径。"
    )
    lines.append("")
    lines.append("| 参数 | p | h | B | variant | status | key-independent done(s) | extract | total | pass |")
    lines.append("|---|---:|---:|---:|---|---|---:|---:|---:|---|")
    for row in rows:
        for name, var in row.get("variants", {}).items():
            if name == "ours_ciphertext":
                continue
            dones = done_times(var.get("log"))
            lines.append(
                f"| {row['label']} | {row['p']} | {row['h']} | {fmt(row.get('B'), 0)} | "
                f"{name} | {variant_brief(var)} | {fmt(dones[0] if dones else None)} | "
                f"{fmt(metric(var, 'extract'))} | {fmt(metric(var, 'total'))} | "
                f"{var.get('passed_marker')} |"
            )
    lines.append("")
    lines.append("| 参数 | NTT optimized total | NTT baseline total | total speedup | 备注 |")
    lines.append("|---|---:|---:|---:|---|")
    for row in rows:
        variants = row.get("variants", {})
        opt = variants.get("ntt_article_optimized", {})
        base = variants.get("ntt_article_baseline", {})
        opt_total = metric(opt, "total")
        base_total = metric(base, "total")
        note = "baseline timeout" if base.get("status") == "timeout" or base.get("timeout_marker") else "same article comparison"
        lines.append(
            f"| {row['label']} | {fmt(opt_total)} | {fmt(base_total)} | "
            f"{fmt(pct(base_total, opt_total), 2)} | {note} |"
        )
    lines.append("")


def add_polyfunctions_section(lines: list[str], log_path: Path) -> None:
    lines.append("## 3. Bootstrapping_Polyfunctions 本地可运行项")
    lines.append("")
    if not log_path.exists():
        lines.append("未找到 Polyfunctions toy 日志。该库生成新多项式需要 Magma，本环境没有 `magma`。")
        lines.append("")
        return
    text = log_path.read_text(errors="replace")
    ok_count = text.count("successful!")
    lines.append(
        f"已执行 `{log_path.name}`。该程序是 toy 参数，不是 Ma-large-p 的同参数集合；"
        f"日志中有 {ok_count} 个 successful marker。生成新 digit polynomials 仍需要 Magma。"
    )
    lines.append("")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=ROOT / "results/baseline_matrix/summary_repeat1.json")
    parser.add_argument("--poly-log", type=Path, default=ROOT / "results/baseline_matrix/polyfunctions_toy_repeat1.log")
    parser.add_argument("--output", type=Path, default=ROOT / "results/baseline_matrix/final_experiment_report.md")
    args = parser.parse_args()

    summary = json.loads(args.summary.read_text())
    ma_rows = [r for r in summary["parameters"] if r.get("article") == "Ma-large-p"]
    ntt_rows = [r for r in summary["parameters"] if r.get("article") == "Ma-NTT"]

    lines: list[str] = []
    lines.append("# 全 baseline 参数实验汇总")
    lines.append("")
    lines.append(
        "本报告由本地真实日志生成。`total/extract/linear1/linear2` 来自密文自举日志的 "
        "`time for ...` 行；`pass=True` 表示日志包含 `### bts finished, everything ok ###`。"
    )
    lines.append("")
    add_ma_section(lines, ma_rows)
    add_ntt_section(lines, ntt_rows)
    add_polyfunctions_section(lines, args.poly_log)
    lines.append("## 4. 可以写进论文的边界")
    lines.append("")
    lines.append("- 我们的主 claim 应放在 Ma-large-p type-B exact cleaner 路径下：D/E 是同参数、同实现路径、同密文自举可比。")
    lines.append("- A/B/C 的真实 `B` 使 `(2B+1)^2 < p` 不成立，不能把 toy 支持下搜到的 cleaner 拿来声称同参数可用。")
    lines.append("- E 参数 `p=65537` 有四阶结构，`A=256` 使项数从 612 降到 307，并触发 order-four evaluator，这是当前最强实验证据。")
    lines.append("- NTT 文章的 optimized 路线验证了同态线性变换优化的重要性，但它是另一个 executor；在未移植前只能作为外部 benchmark。")
    lines.append("- Polyfunctions 库本地只能跑 toy；新参数多项式生成依赖 Magma，因此不能在本机完成同参数大 p 再生成。")
    lines.append("")

    args.output.write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
