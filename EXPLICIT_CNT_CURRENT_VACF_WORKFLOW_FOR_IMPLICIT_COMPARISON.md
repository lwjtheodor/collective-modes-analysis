# 显式 CNT：电流模、VACF 与跨盒长外推工作流

**用途。** 本文把目前 `(8,8)` 显式原子 CNT 中水的纵向/横向电流模分析、轴向 VACF 的构造、tagged--collective 投影以及 10L → 20L 外推整理为一个可迁移的技术基线。它刻意把“已经由资产验证的步骤”与“供 implicit-CNT 比较时仍待检验的假设”分开；不能把本文中的有效阻尼或跨盒长规律直接当作热力学极限的流体动力学定律。

**工作目录。** `H:/gcmc_explore/translational_anomaly/02_isf_collective_modes/`

**符号。** (L_z) 为实际盒长；(R_{\rm eff}) 为水氧的有效圆柱半径；(N_O) 为氧数；(k_n=2\pi n/L_z)；(q_\theta=m/R_{\rm eff})；(q=\sqrt{k_n^2+q_\theta^2})。本文所有低频水动力学讨论均应以**物理波数** (k_n) 或 (q) 配对，而不可仅按离散整数 (n) 配对。

---

## 1. 可引用的资产、输入与脚本入口

### 1.1 轴向 `m=0` 的 10L--20L 主档案

| 目的 | 权威资产或脚本 | 输入/适用范围 | 备注 |
|---|---|---|---|
| 10L 的 tagged--collective 顶角、(F_s)、
  Φ_J 与轴向重构 | `results/collective_mode_response/tagged_collective_vertex_m0_88_L10_100fs_1ns_3rep/2026-08-27/` | 10L、100 fs、共同首 1 ns、3 条**同一构型的 velocity seeds** | `m0_tagged_collective_vertex_arrays.npz`、`m0_tagged_collective_vertex.csv`、`orthogonal_projection_summary.csv` 为首选读入物。 |
| 20L 的对应实测量及 matched-​(k) 对比 | `results/collective_mode_response/L20_axial_m0_current_vertex_100fs_10ns_4rep/2026-08-28/` | 20L、100 fs、固定首 10 ns、4 条同一构型的 velocity seeds | 这是当前 10L → 20L 外推的验证目标；先读其 `README.md`、QA 和紧凑 CSV/NPZ。 |
| 20L 远端原始分析根 | `/lustre/home/users/ewu/vb_gcmc/MD/stage_L20_axial_m0_mode_vertex_100fs_10ns_4rep_20260828/` | 每条轨迹为 `id z vz` 的氧位点 dump | 远端 Python 分析必须在 `HB_analysis` 环境中运行。 |
| 20L 原始 VACF 远端来源 | `/lustre/home/users/ewu/vb_gcmc/MD/stage_vacf_tail_8_8_L20_4rep_veryweaknh_nomom_20ns_100fs_20260826/` | `rep{1..4}/VACF88_L20_rep*_oxygen_id_z_vz_100fs_20ns.dump` | 实际电流/VACF比较固定首 10 ns；不要把更长轨迹的后段混入。 |
| 20L fixed-10-ns VACF 主档案 | `results/collective_mode_response/vacf_msd_alpha_88_L20_veryweakNH_nomom_100fs_fixed10ns/2026-08-28/` | (N_O=2660, L_z=2016.79996\ \AA) | `vacf_msd_alpha_10ns_mean_sem.csv`、`alpha_minima_5to200ps_10ns.csv` 是所用 ensemble 输入。 |
| 通用轴向 VACF 脚本 | `scripts/analyze_vacf_tail.py` | 需要 `vz`；full-water 需 `type`，氧类型默认 3；氧-only dump 无 `type` 亦可 | 同时输出实验室系与逐帧水氧 COM 去除后的 peculiar VACF，带 cadence 和 block 审计。 |
| 更严格的 VACF--MSD--α 每 replica 处理 | `stage_vacf_tail_8_8_L2L5_4rep_weaknh_nomom_10ns_100fs_20260821/analyze_rep_vacf_msd.py`；20L 归并脚本 `stage_L20_fixed10ns_vacf_msd_alpha_20260828/merge_fixed10ns_ensemble.py` | 需要可按 oxygen ID 对齐的 (z,v_z) 序列 | 跨盒长图必须携带各 case 的 thermostat、采样窗口和是否 no-momentum-removal。 |

### 1.2 圆柱/螺旋 ​((k_z,m)) 主档案

| 目的 | 权威资产或脚本 | 最低输入字段 | 结论边界 |
|---|---|---|---|
| (J_z,J_\theta,J_L,J_T) 的一般 ​((n,m)) 筛查 | `results/collective_mode_response/88_10L_helical_cylindrical_current_nm_screen/2026-08-27/` | `id mol type x y z ix iy iz vx vy vz` | 仅 10L、100 fs、共同首 1 ns、3 velocity seeds，(n=0..2,m=0..3)。对应生成脚本：`scripts/analyze_88_10L_helical_cylindrical_currents.py`。 |
| (m=0)、相对 (k_z) 的周向横电流 (J_T=J_\theta) | `results/collective_mode_response/CT_theta_m0_kn_88_L10_100fs_1ns_3rep/2026-08-27/` | full-water 字段 | (n=0..16)。零中心连续谱/1/e 时间是筛查，不等同可精确外推的横声阻尼。 |
| 2L/5L/10L 的 matched-​(k_z) (m=0) 周向横电流比较 | `results/collective_mode_response/CT_theta_m0_crosslength_matched_k_88_L5L10_available_assets/2026-08-28/` | 10 fs、相同物理 (k_z) | 2L 与 10L 的同协议比较支持相同物理 (k_z) 下的可比表示；5L 对比有协议差，不可仅解释为尺寸效应。 |
| (C_{LT},C_{TL}) 的完整性检验 | `results/collective_mode_response/88_10L_helical_CLT_cross_term_100fs_1ns_3seed/2026-08-28/` | full-water 字段 | (n=0..10,m=1..4)，验证有限时延的交叉项写法。 |
| (k_z\sim m/R) 粘弹性带的短时交叉项 | `results/collective_mode_response/88_10L_helical_CLT_mixed_viscoelastic_band_100fs_1ns_3seed/2026-08-28/` | full-water 字段 | 对轴向/周向 **0--2 ps** 残差可显著；尚未和该模的 (F_s)+vertex 一起闭合到 tagged VACF。 |
| (k_z=0,m=1,2,3) 的横向有限-​(q) 对照 | `results/collective_mode_response/88_10L_kz0_m_transverse_control/2026-08-28/` | 来自上述 helical 原数据 | 脚本：`scripts/analyze_kz0_m_transverse_control.py`。该组 (q=m/R_{\rm eff}\ne0)，不能测量 (q\to0) 截距。 |

### 1.3 需要明确排除的旧构造

不要用下列目录的曲线做定量外推：

- `results/collective_mode_response/empirical_10L_modes_kinterpolated_to_20L/2026-08-28/`：曾对已经投影后的 (C_n(t)) 插值，丢失了原始 (F_s)、Φ 和顶角的各自 (k) 依赖，时间尺度与相位均会失真。
- `results/collective_mode_response/sound_mode_10L_to_20L_constantW_k3half_prediction/2026-08-28/`：把 10L 的**离散逐模**权重原样放到 20L，在同一物理 cutoff 下漏掉 Δ(k\propto1/L_z) 的有限盒测度修正，强度过大。
- 仅以 (n) 对齐不同 (L_z) 的图：它们不是同一物理波数。必须用 `10L n ↔ 20L 2n` 或更一般的 matched-​(k) 对应关系。

---

## 2. 前处理合同：坐标、速度、采样与不确定度

### 2.1 最低原始输入

1. **仅轴向 (m=0) LA/vertex/VACF：** 每帧稳定氧 ID、`z`、`vz`，另有 (L_z,N_O)、integration timestep 与 dump stride。`z` 必须能以 PBC 最小映像连续展开；单纯 wrapped `z` 不足以安全算相位随时间的变化。
2. **周向、径向或 (m\ne0) 模：** 每帧 `id mol type x y z ix iy iz vx vy vz`，并记录 CNT 轴线/中心定义和氧 type。当前显式 CNT 数据采用 oxygen type 3。
3. **每个 case 的 protocol manifest：** CNT 是否显式且是否可动、water/CNT 力场、温度与热浴类型/阻尼、是否移除总动量、盒长/半径/水数、dump cadence、可用时窗、配置独立性和 velocity seed 关系。

### 2.2 必须统一的去漂移约定

对每一帧氧水速度先在**笛卡尔坐标**去掉瞬时水氧 COM 速度：

\[
\mathbf v'_i(t)=\mathbf v_i(t)-N_O^{-1}\sum_j\mathbf v_j(t).
\]

然后才投影到 (z,r,\theta)。这避免无 `fix momentum` 的残余整体漂移形成伪长尾。比较时保留原始 lab-frame 曲线作为诊断，但构造模式、peculiar VACF 与 tagged 投影时使用 (mathbf v')。

`m=0,k_z=0` 在上述操作后不是总平移守恒量：例如 (\sum_i v'_{z,i}=0) 恒为零；周向 (m=0,k_z=0) 则是连接的整体扭转/环流自由度，是否衰减由壁耦合、热浴和角动量交换决定，不能用它定义连续介质的 γ₀。

### 2.3 采样、窗口与 replica

- 高频形貌、短时 (C_{LT})：10 fs 更合适；100 fs 不能补回 <0.1 ps 的信息。
- 中等时间 VACF、低 (k) 声模：100 fs 可用，但必须固定每个 replica 的相同物理窗口（当前 20L 为首 10 ns）。
- velocity seeds 估计动力学散布，却不等同独立的初始构型。报告中应写 “seed-conditional SEM”，不能把它标成构型独立的总体置信区间。
- 每个结果包应有 `metadata.json`、per-replica 数据/日志、ensemble CSV/NPZ、`README.md`、QA 文档和 `FINISHED.txt`；提交成功或队列状态不构成完成证据。

---

## 3. 电流模的定义：以传播波矢决定 L/T，而非以速度分量命名

对水氧定义复相电流（规范化采用 (N_O^{-1/2})）：

\[
J_a(k_n,m;t)=N_O^{-1/2}\sum_{j=1}^{N_O}v'_{a,j}(t)
\exp\{i[k_nz_j(t)+m\theta_j(t)]\},\qquad a\in\{z,\theta\}.
\]

这里 `longitudinal/transverse` 是相对于二维切向波矢
(\mathbf q=(k_n,q_\theta)) 的名称。令 (a=k_n/q, b=q_\theta/q)，则

\[
J_L=aJ_z+bJ_\theta,\qquad
J_T=-bJ_z+aJ_\theta.
\]

因此四种常被混淆的极限为：

| 波矢 | 纵支 (J_L) | 横支 (J_T) | 速度物理含义 |
|---|---|---|---|
| (m=0,k_z\ne0) | (J_z) | (J_\theta) | 相对轴向传播；周向速度确实是横支。 |
| (m\ne0,k_z=0) | (J_\theta) | (-J_z) | 相对周向传播；轴向速度此时是横支。 |
| (m\ne0,k_z\ne0) | (aJ_z+bJ_\theta) | (-bJ_z+aJ_\theta) | 两个笛卡尔/圆柱分量均混合。 |
| (m=k_z=0) | 无定义 | 无定义 | 没有传播方向；另行定义总平移/扭转零模。 |

这也回答了“周向速度是否天然为纵向”这一点：**不是**。例如 (m=0,k_z\ne0) 时 (v_\theta) 是横向；而 (m\ne0,k_z=0) 时它才是纵向。

### 3.1 相关函数矩阵与交叉项

按同一 ((k_n,m)) 定义

\[
C_{AB}(t)=\langle J_A(t_0+t)J_B^*(t_0)\rangle_{t_0},\qquad A,B\in\{L,T\}.
\]

完整的周向电流核不是只取对角项：

\[
C_{\theta\theta}(t)=b^2C_{LL}(t)+a^2C_{TT}(t)+ab\,[C_{LT}(t)+C_{TL}(t)].
\]

对轴向核相应为

\[
C_{zz}(t)=a^2C_{LL}(t)+b^2C_{TT}(t)-ab\,[C_{LT}(t)+C_{TL}(t)].
\]

只有 (t=0) 且用平稳互易对称性时，才可将交叉部分写成 (2ab\operatorname{Re}C_{LT}(0))。此前的 10L audit 给出低 (k_z,m=1) 周向 0--2 ps 交叉修正 <0.14%，但在 (k_z\sim m/R) 粘弹性带可达 10--16% 的短时量级。因此：低-​(k_z) 周向中期尾可以近似对角；想修复轴向短时残差时不能省略 (C_{LT},C_{TL})。

---

## 4. 从轨迹到 (C_{JJ})、谱和阻尼的实际流程

1. **字段与帧审计。** 核对每帧原子数、氧 ID 恒定性、box、stride、units、总时长；从连续坐标得到 (z_i(t))，以 (x,y) 和 CNT 轴中心计算 (θ_i(t))。
2. **去水 COM 后构造 (J_z,J_\theta)。** 正/负 (k) 已由实场表示处理；不额外乘一个 “±k 简并因子”。
3. **all-origin ACF。** 逐 seed 计算 (C_{JJ}(t))，随后再 ensemble 平均；给出未归一化 (C_{JJ}(0)) 和 Φ：
   \[
   \Phi_J(k,m,t)=C_{JJ}(k,m,t)/C_{JJ}(k,m,0).
   \]
4. **频域。** 在明确的 window、detrend 与归一化下 FFT / Welch；(S_{JJ}) 或 current spectrum 的峰位置、半宽均应回到时域 ACF 核查。把 correlation spectrum 与绝对 PSD 区分记录。
5. **低 (k) 的 DHO 描述（仅作有效参数化）。** 常用拟合核为
   \[
   \Phi_J(t)=e^{-\Gamma t}\,[A\cos(\omega_dt)+B\sin(\omega_dt)].
   \]
   (B/A) 是有限窗口/投影下的有效相位表征，不能把它直接解释为“严格自相关在 (t=0) 的导数不为零”。对严格实、平稳自相关，零时斜率应为零；若复模、实部投影、数据离散化或模型窗口使该约束不显式，应以对称延拓/拟合残差另作检验。
6. **跨长度只匹配物理 (k)。** 10L `n` 对 20L `2n`。因 (J\propto N_O^{-1/2})，所存 (C_{JJ}(0)) 是 intensive；若改用未归一化的 Σ(v e^{ikz})，强度会平凡地随 (N_O) 成比例，不能被称为集体增强。

当前 10L--20L intensity audit 的归一化检查使用 (√{N_O}c_n)、√{N_O}a_n 与 (N_OW_n)，而不是直接比 (c_n,a_n,W_n)。这是 finite-box 变化下唯一可解释的比较方式。

---

## 5. VACF、MSD 与 α 的构造

### 5.1 直接轴向 peculiar VACF

\[
C_{vv,z}(t)=\left\langle N_O^{-1}\sum_i v'_{z,i}(t_0+t)v'_{z,i}(t_0)\right\rangle_{t_0}.
\]

应同时留存有量纲的 (C_{vv,z})（\(\AA^2\,{\rm ps}^{-2}\)）和若需要展示的 (C(t)/C(0))。跨长度做面积、负瓣深度或外推对比时绝不可暗中按零时方差再归一化。

对每条轨迹 all-origin 计算后再平均；block 或 seed 级曲线用来检查窗口稳健性。完整负瓣的面积、宽度、深度必须用相邻零点界定，不能将 0--100 ps 的截断段称为完整 lobe。

### 5.2 由 VACF 生成 MSD 与瞬时指数

\[
I(t)=\int_0^tC_{vv,z}(s)\,ds,\qquad
J(t)=\int_0^tI(u)\,du,\qquad
\mathrm{MSD}_z(t)=2J(t),\qquad
\alpha_z(t)=\frac{tI(t)}{J(t)}.
\]

数值实现使用一致的积分规则，并把 ODE 形式的 α 与直接 MSD 的 decade estimator 交叉核验。`\alpha_{\min}`、第一回正峰位置/高度、负瓣退出时刻等都应从同一条 ensemble 曲线和同一采样网格提取；per-seed 值只报告为 conditional spread。

### 5.3 已有尺度现象的正确表述

2L--10L 的 raw 100 fs 轴向 peculiar VACF 显示第一负瓣、α 最小值与第一回正特征会随 (L_z) 推迟。它是“中等时间、mode-resolved memory”的实测有限尺寸现象，尚不是已证明的渐近线性定律。20L 与 10L 在 thermostat、窗口与 seed 结构上不同，因此只能用 protocol-matched 控制或显式模型项剥离后再拟合长度律。

---

## 6. tagged--collective 顶角与无自由幅度 VACF 重构

这是目前最重要、也最适合搬到 implicit CNT 的闭合骨架。

### 6.1 每个 tagged 氧的实投影场

对 (m=0) 轴向正 (k_n) 实场，定义

\[
X_n^{(i)}(t)=\sqrt2\,\operatorname{Re}\left[e^{-ik_nz_i(t)}J_z(k_n,t)\right].
\]

在零时刻构造中心化的静态协方差和 tagged--velocity 顶角：

\[
K_{nm}=\langle\delta X_n^{(i)}\delta X_m^{(i)}\rangle,
\qquad c_n=\langle\delta X_n^{(i)}\delta v'_{z,i}\rangle.
\]

用伪逆（并报告条件数、截断规则）求

\[
\mathbf a=K^+\mathbf c,\qquad
W_n=\frac{c_na_n}{\langle(\delta v'_z)^2\rangle}.
\]

这是统计投影得到的权重，**不是**任意拟合的振幅。相邻 (k) 模若强相关，简单逐模 (c_n^2/K_{nn}) 会 double count，因此应以全矩阵 (K^+c) 为主，并将逐模量留作诊断。

### 6.2 重构式与必要的 (F_s)

\[
C^{\rm proj}_{vv,z}(t)=\sum_{n=1}^{n_{\max}}
W_n\,F_s(k_n,t)\,\Phi_{J,z}(k_n,t).
\]

其中

\[
F_s(k_n,t)=\left\langle e^{ik_n[z_i(t_0+t)-z_i(t_0)]}\right\rangle_{i,t_0}
\]

必须从同一轨迹、同一 (k)、同一去漂移/位置定义测量。将 (F_s\approx1)、用 total ISF 替代 (F_s)，或只插值已经相乘后的 (C_n=W_nF_s\Phi)，都会损坏中等时间的幅度或相位；这正是早期外推失败的主要方法学原因之一。

### 6.3 长度与粒子数标度

由于 (J^{(N)}=N_O^{-1/2}\sum_i\cdots)，若相同物理 (k) 的连续低-​(k) 权重谱在有限范围内近似常数，则离散模式权重应满足

\[
W_n\propto \Delta k=\frac{2\pi}{L_z}.
\]

故在同一物理 cutoff 下，20L 相对 10L 的**逐模**权重应近似减半，而模式数加倍；实正场已含 ±(k) 对，不能再加第二个两倍因子。更稳健的直接量是 (N_OW_n) 的 matched-​(k) 比较。

### 6.4 当前 20L 的自检结果

在 20L fixed-first-10-ns archive 中，`n=1..20` 的**本体系实测**重构在 10--200 ps 达到 (R^2=0.9795)、RMSE (3.288\times10^{-5}\ \AA^2\,\mathrm{ps}^{-2})。移除 n=1 后，形貌相关降至约 0.53、RMSE 约 (2.31\times10^{-4}\)，且 n=1 占第一负瓣有符号面积约 48%。这证明当前投影在该协议下能重构中期轴向 VACF；不证明只由声模已完全闭合短时和长时记忆。

---

## 7. 10L → 20L 的外推：可转移量、不可转移量与推荐顺序

### 7.1 必须逐项外推，不能插值成品曲线

对 10L 的 matched-​(k) 低波数测量，分别拟合/插值：

1. (N_OW(k)) 或连续谱密度；
2. (F_s(k,t)) 的 (k) 依赖；
3. 声载波 (omega(k))（低 (k) 区可由 (c_s\simeq16.27\ \AA\,\mathrm{ps}^{-1}) 作先验/交叉核验）；
4. 阻尼 γ((k)) 与相位/非 DHO 残差；
5. 再在目标 20L 的离散 (k_n) 网格上重建 (W_nF_s\Phi_n) 并求和。

推荐脚本与其用途：

- `results/collective_mode_response/L20_axial_m0_current_vertex_100fs_10ns_4rep/2026-08-28/scripts/compare_10L_20L_current_vertex.py`：matched-​(k) 的 (N_OW,F_s,\Phi_J) 与 no-free-amplitude 总重构对比。它会更新该 archive 的说明/QA，复制到新实验目录再运行。
- `.../scripts/fit_matched_phi_dho.py`：matched-​(k) DHO 参数化和相位诊断。
- `.../scripts/extrapolate_10L_to_20L_n1.py`：严格的 withheld 20L n=1 低波数 transfer test。
- `.../scripts/audit_l20_n1_damping_identifiability.py`：改变拟合起点/窗口与 seed 的 n=1 阻尼稳健性审计。
- `.../scripts/fit_l20_lowk_gamma_offset_models.py`：比较 (Ak^{3/2}\)、γ₀+​(Ak^{3/2}\)、自由幂律与带 offset 的自由幂律。
- `.../scripts/matched_k_current_intensity_audit.py`：检查 (N_O) 归一化及 (C_{JJ}(0)) 的 matched-​(k) 强度。

### 7.2 已验证到什么程度

- 10L `n=1..8` 对 20L **withheld n=1**：载波频率预测 (0.05016\) vs 实测 (0.04887\ \mathrm{ps}^{-1})；(N_OW_1) 预测 (2.034\) vs 实测 (1.967)；(F_s) 的 0--100 ps RMSE 为 (2.49\times10^{-5})。这些部分可转移。
- 同一盲测中 γ 预测 (0.00341\)--(0.00396\ \mathrm{ps}^{-1})，实测有效值约 (0.00882\ \mathrm{ps}^{-1})。故 10L 不能在没有额外输入时外推 20L 最低模耗散。
- 20L n=1 的 effective DHO γ 在 50--250 ps 变窗口时仍为 (0.00874\)--(0.00901\ \mathrm{ps}^{-1})，100 ps 处 seed-conditional SEM (0.00038\ \mathrm{ps}^{-1})。这不支持“只是拟合窗太短”，但仍不能将它宣布为真正的 (q=0) 水动力学阻尼。
- 20L 的低-​(k) 有效 γ 更偏好带 offset 的形式：n≤4 时 γ₀+​(Ak^{3/2}) 给 γ₀=0.00560 ps⁻¹，并以 ΔAICc=174.4 胜过无截距式；n≤8 时最佳形式为 γ₀+​(Ak^{1.793})，γ₀=0.00688 ps⁻¹。它们只是当前壁耦合、热浴与记忆核折叠后的 effective crossover 描述。

### 7.3 外推的最小成功判据

1. 用 10L 模型预测后，不调自由总振幅；
2. 同时报告 10--200 ps shape correlation、RMSE、第一负瓣面积/退出时间、第一回正峰位置/高度；
3. 拆开 (W)、(F_s)、Φ 的误差；
4. 若 γ 为唯一失配项，只能说 “需要目标体系的耗散/记忆核输入”，而不是归咎于权重或 ±(k)；
5. 禁止以拟合完的 target VACF 反调某个全局幅度后再称为“预测”。

---

## 8. 用于 implicit CNT 的比较合同

### 8.1 可以直接沿用的核心

若 implicit CNT 给出水氧的位置和速度，以下定义可完全照搬：

- peculiar (C_{vv,z})、MSD、α；
- (J_z(k_n,m=0))、(F_s)、Φ_J；
- tagged--collective (K,c,W) 与无自由幅度重构；
- matched physical (k)、(N_O^{-1/2}) 电流归一化、(N_OW) 强度审计；
- 完整零点界定的 VACF lobe 指标与 seed/configuration 分层不确定度。

### 8.2 implicit CNT 必须额外固定的协议

为让 explicit vs implicit 的差异可解释，至少逐项匹配或记录：

| 类别 | 必须一致/给出的量 | 不一致时的处理 |
|---|---|---|
| 几何 | (L_z)、可及半径或等效氧密度剖面、轴定义 | 用 matched (k) 与 matched accessible cross-section；不要仅引用名义 CNT 直径。 |
| 水的状态 | (N_O) 或线密度、温度、water potential、压力/化学势条件 | 分开建 protocol family；不合并曲线。 |
| 墙 | 显式 CNT 的刚柔性/声子/摩擦与 implicit wall potential、时间依赖 | 这是最可能改变 γ、短时 (C_{LT}) 和横向耗散的物理差别，应作为主比较对象。 |
| 动力学 | thermostat 类型、阻尼时间、是否移除总动量、采样 cadence、起始窗口 | 这些直接改变低频核；至少做一个同 protocol 的显式/implicit pair。 |
| 统计 | independent configurations 数、每个配置的 velocity seeds、统一可用时窗 | 不可把同一母构型的多 seed 当作配置 replica。 |

implicit wall 若没有原子 CNT，就不能直接分解 CNT 声子或 CNT--water 力的相关通道；但不妨碍水电流和 tagged 投影的比较。若希望解释差异，需要额外保存 implicit wall force/torque 与每个水分子的 wall force，或至少保存壁势的明确表达与固定轴线。对于 (m\ne0) 分析，仍必须有 (x,y) 和可重复定义的 θ。

### 8.3 建议的最小 pairwise 实验

1. 先做同 (L_z,N_O,R_{\rm eff},T)、同 100 fs cadence、同 10 ns 窗口的显式/implicit `m=0` 轴向对比。
2. 在 (k_1,k_2,k_3,k_4) 分别比较 (C_{JJ}(0))、Φ、(F_s)、(N_OW) 和 direct/projection (C_{vv,z})。这能判别差异首先来自静态顶角、自扩散去相干、传播相位还是耗散。
3. 再做 (m=1)、matched (q_\theta=1/R_{\rm eff}) 的 (k_z=0) 有限-​(q) 横向控制，以及少数 (k_z\sim m/R) 点的 (C_{LT},C_{TL})。
4. 最后才拟合 γ((k))。要判断是否真有零波数阻尼，需一族不同 (L_z) 的 (m=0,k_{\min}) 或一族不同 (R) 的 (m=1,k_z=0)；单个 ​(m=1) 有限-​(q) 点不能回答该问题。

---

## 9. 目前方法学尚未闭合的点

1. 轴向 (m=0) 的低-​(k) tagged--collective 重构在中等时间已很强，但它不是包含全部模式、全部 memory kernel 的严格 Mori--Zwanzig 闭合。
2. (m\ne0) 的 (C_{LT},C_{TL}) 已证实会修正粘弹性区间短时核；同模式的 (F_s(k,m,t)) 与静态 vertex 尚未完整纳入，因而尚不能定量归因到 tagged (C_{vv,z}) 或 (C_{vv,\theta})。
3. DHO 的 γ 是 finite-system、wall/thermostat/protocol-dependent effective 参数；不能以 γ₀+​(Ak^p) 的拟合截距宣称测到连续 (q=0) 横声或纵声线宽。
4. (W(k)) 在当前 10L 低-​(k) 范围内与常数连续谱兼容；Gaussian (W_0e^{-Ck^2}) 并未取得更好的 AICc/leave-one-seed-out 误差。它可作 UV regularizer，不能被当作已测得的物理 cutoff。
5. 声速 (c_s\simeq16.27\ \AA\,\mathrm{ps}^{-1}) 可以固定低-​(k) 相位标尺，却不能单独给出 γ、相位/记忆核、(F_s) 或静态投影权重。因此仅用 (L/c_s) 不足以完全重构 VACF。

---

## 10. 新 implicit-CNT 分支的建议输出合同

每个明确协议保存一份独立结果包，例如：

```text
results/collective_mode_response/<explicit_or_implicit>/<system>/<date>/
  README.md                    # 定义、协议、明确的比较边界
  metadata.json                # Lz, Reff, N, T, wall, thermostat, cadence, window
  input_audit.csv              # dump 字段/帧/ID/cadence 审计
  per_replica/
    mode_arrays.npz            # lag,k,m,CJJ,CJJ0,Phi,Fs,K,c,a,W,direct VACF
    mode_summary.csv
    vacf_msd_alpha.csv
    metadata.json, SUCCESS.txt
  aggregate/
    ensemble_arrays.npz
    matched_k_summary.csv
    reconstruction_metrics.csv
    FINISHED.txt, metadata.json
  figures/
  QA.md
```

`assets.md` 只登记具备上述最小文件证据的完成包；展示图放入 `assets/library/collective_dynamics/`，而数值、脚本、元数据和 QA 留在 `results/` 主档案。新分支应先以本文第 8.3 节的最小 pair 通过字段/协议/归一化审计，再讨论 explicit--implicit 的物理差别。

---

## 11. 一句话的当前物理基线

显式 CNT 中，低-​(k)、`m=0` 轴向模及其实测 tagged 顶角、(F_s) 和 Φ 能无自由总幅度地重构 20L 中等时间轴向 VACF；跨长度可转移的是相位载波、静态强度和 (F_s)，最不稳定且最可能暴露 explicit wall/thermostat/memory 差异的是最低模的有效耗散。对周向/螺旋模，传播方向决定 L/T；在低 (k_z) 的中期尾中交叉项很小，但在 (k_z\sim m/R) 短时粘弹性区不能忽略。
