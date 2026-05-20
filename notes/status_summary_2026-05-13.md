# 二代 BGV 自举优化：当前状态与方案总结

Date: 2026-05-13

## 一、已完成的工作

### 1. Order-4 Character-Filtered Cleaner（已验证，核心贡献）

**数学原理：** 选择 A=256 ∈ F_{65537}^× 满足 A²≡-1 mod p，利用 Galois 协变性
P_A(AX) + A·P_A(X) = AX 强制 cleaner 的 monomial support 只含 k≡3 mod 4 的项。

**结构分解：** P_A(X) = (1/2)X + X³·Q(X⁴)，其中 deg(Q) = 305。

**实验结果：**
- extract: 160.8s → 80.8s（-49.8%）
- total: 213.1s → 131.2s（-38.4%）

### 2. Parallel CoeffToSlot（已验证，工程贡献）

**原理：** HElib newBTS pipeline 中两个独立的 coeffToSlot 变换用 std::async 并行执行。

**实验结果：**
- linear2: 45.3s → 23.9s（-47.2%）
- total: 213.1s → 182.4s（-14.4%）

### 3. 组合加速（已验证）

Order-4 + Parallel CoeffToSlot 组合：
- linear1=4.6s, linear2=24.7s, extract=86.4s, total=118.5s
- **vs baseline 213.1s: 1.80x 总加速**

## 二、理论贡献（已完成定理形式化）

### 统一的 Galois-Covariant Digit Extraction 定理

**定理：** 对任意素数 p，r|p-1，A∈F_p^× 且 ord(A)=r：
1. Cleaner 满足协变性 P_A(AX) + A·P_A(X) = AX
2. 精确项数 |T| = (d-1)(|S_A|-1)/r + 1（Schur 引理保证精确）
3. 结构分解 P_A(X) = c₁X + X^{r-1}·Q(X^r)
4. 求值加速 √r 倍

**普适性推论：** 对所有 p > (2B+1)² 且 p≡1 mod 4 的参数（覆盖所有 p≥65537 的实际场景），
order-4 cleaner 自动适用，无需参数调整。

### 定理文件位置
`/home/luck/xzy/0424project/github_order4_cleaner/paper/theorem_galois_covariant.tex`

## 三、线性变换优化探索（负面结果，但有理论价值）

### 尝试过的方案及结论

| 方案 | 数学正确性 | 密文验证 | Wall-clock 加速 | 失败原因 |
|------|-----------|---------|----------------|---------|
| Rader97 (13-stage) | ✓ | ✗ (噪声爆炸) | N/A | 13次key-switch累积噪声超出budget |
| Good-Thomas (2-stage) | ✓ (F_97上) | ✗ (多项式环不成立) | N/A | 标量分解不能直接用于多项式矩阵 |
| Frobenius-Orbit BSGS | ✓ (关系验证) | ✓ | ✗ (变慢) | 不减少非零对角线数，只影响常数 |
| Asymmetric BSGS | ✓ (论文方法) | ✓ | ✗ (变慢) | degree-305太小，优势需degree>32768 |
| BSGS参数调整 (g=16) | ✓ | ✓ | ✗ (变慢) | 增加baby-step数反而增加key-switch |

### 根本性结论

对于 D=96 的 dense Vandermonde 矩阵（非 power-of-two 分圆环），在 HElib 的
key-switching + hoisting 框架下，**BSGS 已是近最优算法**。原因：
- 矩阵是 dense 的（96条非零对角线），无法通过稀疏性减少计算
- 多 stage 方案的 hoisting precomputation overhead（~0.6s/stage）抵消 automorphism 减少的收益
- Frobenius 轨道关系只影响 plaintext 常数（便宜操作），不影响 key-switch 数量（贵操作）

## 四、论文定位建议

### 标题建议
"Galois-Covariant Digit Extraction for Large-Prime BGV Bootstrapping"

### 贡献列表
1. **统一理论框架**：Galois-covariant cleaner 定理（适用于任意 p, r|p-1）
2. **精确项数公式**：|T| = (d-1)(|S_A|-1)/r + 1（Schur 引理证明是精确值）
3. **结构化求值算法**：P(X) = c₁X + X^{r-1}·Q(X^r)，加速因子 √r
4. **实验验证**：p=65537 上 extract -49.8%，total -38.4%（单独），1.80x（组合）
5. **普适性证明**：对所有 p≥65537 且 p≡1 mod 4 自动适用

### 与 baseline 的区分
- Chen-Han: 降低 depth
- GIKV: 降低 degree (via null lattice CVP)
- Ma et al.: 降低 degree (via bounded support)
- **本文: 在固定 degree 下，降低 evaluation cost (via Galois structure)**

### 线性变换侧
- Parallel CoeffToSlot 作为工程优化（wall-clock 加速 47%）
- Good-Thomas/Frobenius 理论分析作为 negative evidence + future work
- 说明非 power-of-two 分圆环上的线性变换优化是一个 open problem

## 五、待完成事项

1. [ ] 在其他 large-p 参数上验证（如 p=40961, p=786433）
2. [ ] 完善论文 LaTeX 写作
3. [ ] 整理实验数据为论文表格格式
4. [ ] 补充 order-6 的理论分析（对 p≡3 mod 4 的情况）
5. [ ] 与 Ma et al. baseline 的详细对比表
