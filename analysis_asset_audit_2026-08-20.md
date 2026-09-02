# 02_isf_collective_modes：分析资产全量审计

**审计日期：** 2026-08-20  
**范围：** 本项目下已落盘的分析结果包、展示资产和已完成的专题 archive；任务脚本、仅提交的 stage、孤立 dump 和队列记录不计作已完成分析资产。  
**文件级入口：** `assets/library/collective_dynamics/` 的每一个 `source/` 是到权威结果包的无复制目录链接。其余主题的权威结果仍在 `results/collective_mode_response/<topic>/<date>/`。

## 0. 本轮文件系统复核（2026-08-20）

- 当前有 **34 个默认结果包**与 **3 个完整隔离包**（合计原 37 个分析包）；本表逐包登记，不把同源的 PNG/PDF/SVG/TIFF、同一结果的 per-k 重绘或 block-convergence 图误计为独立模拟。隔离后的完整包位于 `_quarantine/`，不进入默认检索。
- `F:/ccfep_gcmc_archive_20260814/` 当前有 **7 个 campaign、336 个 `.dump`、约 1.09 TiB**。这是可访问的本地 archive，但除第 7 节明示 P1 的链路外，仍须以相对路径加字节数或 SHA-256 核验后，才能称为远端输入的验证副本。
- `(8,8)` 5L 全模态 CJJ 的主包 `88_5L_LA_TAr_TAtheta_dispersion/2026-08-19` 已直接复核 `metadata.json` 与 `FINISHED.txt`：3 个输入 dump、O `z,vx,vy,vz`、10 fs，分别为 347.26/371.01/349.95 ps；`N_water=665`、`Lz=504.20 Å`、`k=0.01246169..0.99693541 Å^-1`（n=1..80），相关最大 lag 100 ps。分析定义为 instantaneous O-COM-subtracted modal current；CJJ 去除每条轨迹的 temporal modal-current mean；Welch `nperseg=16384`、50% overlap、频率分辨率 0.0061035 ps^-1。它的三张大源表同时已有行数/源 SHA-256/输出 SHA-256 验证的 HDF5 紧凑工作副本，见 `derived_data/compact/`；CSV 仍为无损审计源。

## 1. 审计原则与键定义

### 1.1 不把盒长当成波矢资产的主键

对于 `CJJ(k,t)`、density/DSF `S(k,omega)`、current spectrum 与色散，资产主键为：

`(observable, physical k=2*pi*n/Lz, chirality, simulation protocol, analysis definition, usable omega/t window)`。

因此，只要相同物理 `k` 的源数据充足、协议和定义可比，来自不同盒长的曲线应进入同一 `k` bin；它们不是必须分别补齐的“盒长缺口”。盒长只用于：提供新的更低 `k_min`、提供同一 `k` 的独立复现、或改变有限时间窗/频率分辨率时的协议标签。以下所有“待补”均按这个规则判定。

### 1.2 状态术语

- **可引用**：源表、元数据/README、图件及完成/QA 证据已存在。
- **诊断/限定**：资产完整，但只能回答其写明的筛查、灵敏度或单协议问题。
- **不合并**：同一 observable 但模拟或分析协议不同；可并列比较，不可合池拟合。
- **superseded**：保留溯源，不作为结论来源。

## 2. 全量结果包清单（项目级）

下表是 37 个 `results/collective_mode_response` 根结果包的全量登记；每一行代表该包中全部图、CSV、JSON、脚本、QA 与逐 replica 文件，而非单一图片。

| 分析层 / 结果包 | observable 与模式 | chirality / L 覆盖 | 仿真与分析协议 | 状态与用途 |
|---|---|---|---|---|
| `88_5L_LA_TAr_TAtheta_dispersion/2026-08-19` | `CJJ(k,t)` n=1..80、LA/TA-r/TA-theta current spectra、色散、负瓣和 mode weight | `(8,8)`，5L，3 rep | O 10 fs、约 347--371 ps 可用段；instantaneous COM subtraction；per-replica + mean/SEM | **可引用**；5L all-mode CJJ 的主档案 |
| `88_5L_low_frequency_signed_Skw_CJJ/2026-08-19` | density/LA/TA-r/TA-theta two-sided Welch `S(k,omega)`，n=1..20 | `(8,8)`，5L，3 rep | O 10 fs；Welch 16384、50% overlap；低频 ±5 rad ps^-1 | **可引用**；heatmap 为每 channel,n 归一化，绝对强度应读 CSV |
| `_quarantine/unqualified/88_5L_per_k_semilog_Skw_LA_TAr_TAtheta/2026-08-19` | 5L per-k semilog spectra | `(8,8)`，5L | 与上述 5L source family 相同 | **unqualified，已隔离**；未见完成证据，只作待恢复候选 |
| `88_5L_semilog_Skw_LA_TAr_TAtheta/2026-08-19` | 5L semilog 汇总谱 | `(8,8)`，5L | 同一 5L source family | **可引用但非独立**；是同源重绘 |
| `88_5L_mode_weight_n1_n20/2026-08-19` | n=1..20 mode-weight 汇总 | `(8,8)`，5L | 同一 5L source family | **可引用但非独立**；谱积分/权重视图 |
| `88_5L_Wk_block_convergence/2026-08-19` | `W_k` block convergence | `(8,8)`，5L | 同一 5L source family，block estimator | **诊断**；用于稳定性，不增加新参数点 |
| `88_10L_per_k_semilog_Skw_LA_TAr_TAtheta/2026-08-19` | 10L per-k semilog `S(k,omega)` 与 phase/group velocity | `(8,8)`，10L | 10L 输入协议；频率/时间分辨率与 5L 10-fs 资产不同 | **可引用、协议分离** |
| `88_10L_per_k_semilog_Skw_pm10_LA_TAr_TAtheta/2026-08-19` | ±m/相位约定的 10L per-k 谱 | `(8,8)`，10L | 与上行同源、不同模态/相位处理 | **可引用、替代分析定义**；不得与标准版静默合并 |
| `cjj_kmin_spectrum_from_acf/2026-08-18` | n=1 `CJJ/CJJ(0)` correlation spectrum、absolute diagnostic | `(8,8)`，3L/4L/5L/10L | 300 ps Bartlett cosine transform | **可引用**；是 correlation-spectrum estimator，非原始 PSD；2L 缺本地 matching CJJ |
| `minimum_k_current_lobe_t_over_L/2026-08-10` | n=1 first complete negative lobe timing | `(8,8)`，1--5L、10L | zero-crossing bounded CJJ lobe；1--5L 10 fs/1 ns，10L 100 fs/10 ns | **可引用、协议分离** |
| `fig2_longitudinal_modes_88_rh75_330k/2026-08-11` | matched-k ISF `F/Fs/Fd`、longitudinal current modes | `(8,8)`，主要 2--5L，含限定 10L | fixed CNT、weak NH/no momentum；2--5L 10 fs/1 ns，10L 100 fs/10 ns | **可引用**；matched physical k 优先于 box label |
| `msd_numerical_alpha_cjj_phase/2026-08-11` | MSD/alpha 与 CJJ phase 数值关联 | `(8,8)` 为主 | 见包内脚本与源表 | **诊断**；不可替代 CJJ 或 MSD 主档案 |
| `absolute_cvj_first_negative_lobe/2026-08-09` | absolute/per-water/normalized CvJ、first lobe | `(7,7)/(8,8)/(9,9)/(17,0)`，1--5L、10L | 1--5L 10 fs/1 ns；10L 100 fs/10 ns；逐 replica SEM | **可引用、协议分离**；`superseded_one_decade_window/` 禁止作为主结论 |
| `stage_cjj_dispersion_crosschirality_10L_20260814/`（结果在 `results/`） | CJJ n=1..9、DHO peak、低-k 色散与 sound speed | `(7,7)/(9,9)/(17,0)`，10L，3 rep；另有 `(8,8)` 对照表 | 100 fs/10 ns，max lag 300 ps，COM-subtracted current | **可引用**；只证明 10L 同协议色散，非长度闭合 |
| `dsf_k1_all_10fs_assets/2026-08-18` | connected density `S(k1,omega)` 主矩阵、broad/low-freq semilog | `(7,7)/(9,9)`，2--5L，每格 8 velocity-seed rep | 10 fs、2 ns；Hann spectrum，按 O 数归一；SEM across velocity seeds | **可引用**；64 条 matched 轨迹；不等同独立构型 SEM |
| `_quarantine/superseded/dsf_kmin_10fs_asset_matrix/2026-08-18` | density-ISF transform matrix | 若干 kmin 资产 | interim 分析 | **superseded，已隔离**；不得引用 |
| `vdos_dsf_L2L10_weakNH_zvz/2026-08-13` | axial peculiar VDOS、minimum-k 与 matched-lambda DSF | `(8,8)`，2--5L、10L | 各 length metadata 为准；源为 z/vz | **可引用、协议分离** |
| `oxygen_vacf_spectrum_cm_L2L10/` | oxygen VACF spectrum | `(8,8)`，L2--L10 | COM 约定见源表 | **诊断**；与 VDOS 需核对窗函数/归一化后才可并表 |
| `vacf_vdos_10fs_all_available/2026-08-18` | 10-fs VACF PSD/VDOS availability matrix | 多 chirality、若干 L/rep | availability-driven；每 case 资产完整度不同 | **诊断/索引**；不是统一 parameter matrix |
| `vacf_alpha_10fs_L2L10_unified_2026-08-20` | axial peculiar Cvv、strict Cvv-ODE alpha | `(8,8)`，2--5L x3、10L x8 | 10 fs、max lag 200 ps；instantaneous oxygen-mean-vz subtraction；all-origin | **可引用**；replica 深度不等，保留 counts |
| `crosschirality_vacf_10fs_8rep_2ns_200ps_2026-08-18` | Cvv 及 complete lobe archive | `(7,7)/(8,8)/(9,9)`，2--5L x8 | 10 fs/2 ns；instantaneous oxygen mean-vz removal | **可引用**；8 个速度种子，不是独立构型起点 |
| `crosschirality_vacf_cvv_alpha_ode/2026-08-14` | Cvv/strict alpha ODE | `(7,7)/(9,9)/(17,0)`，1--5L x2、10L x3 | 1--5L 10 fs/1 ns；10L 100 fs/10 ns；lag <=200 ps | **可引用、协议分离** |
| `crosschirality_cvv_comparison/2026-08-18` | cross-chirality Cvv figures/tables | `(7,7)/(8,8)/(9,9)` 为主 | 依赖 10-fs archive | **展示/比较层**；不是独立模拟 |
| `decom_msd_cvv_alpha_L2L10_weakNH_zvz/2026-08-13` | de-COM MSD、Cvv、alpha | `(8,8)`，2--10L | weak-NH/no-momentum，z/vz 输入 | **可引用但与 2026-08-20 unified 包交叉**；优先后者的统一 10-fs定义 |
| `vacf_cjj_L3L5_weakNH_zvz/2026-08-13` | Cvv 与 n=1 CJJ 并列 | `(8,8)`，3L/5L | weak-NH/no-momentum z/vz | **诊断**；两长度关联，不是完整 k scan |
| `vacf_tail_2L_weakNH_8rep/2026-08-12` | 2L VACF tail | `(8,8)`，2L x8 | weak-NH | **可引用、单点** |
| `vacf_tail_8_8_L10_10fs_8rep_1ns_2026-08-19` | 10L 10-fs VACF tail/VDOS special | `(8,8)`，10L x8 | 10 fs/1 ns | **可引用**；与历史 100-fs 10L 分开 |
| `vacf_2L_vs_4L_oddmode_projection/2026-08-13` | odd-mode projection / counterfactual VACF | `(8,8)`，2L/4L | modal filtering | **诊断**；不能当未投影 VACF |
| `counterfactual_vacf_5L_n1_bandremoval/2026-08-13` | 5L n=1 band removal counterfactual | `(8,8)`，5L | modal filtering | **诊断**；与原 VACF 同源 |
| `alpha_z_collapse_decom_20ns_1ps/2026-08-09` | 20 ns/1 ps de-COM alpha collapse | `(8,8)`，1--5L；10L 输入曾审计 | fixed CNT NVT/momentum control；1 ps | **限定**；`superseded_unverified_10L_input/` 不可引用 |
| `alpha_z_functional_t_over_L_collapse_1ps_cadence/2026-08-10` | 统一 1 ps cadence t/L collapse | 多 chirality，1--5L/10L | 原轨迹统一 decimation 到 1 ps | **可引用、时长仍不统一** |
| `alpha_z_t_over_L_crosschirality_diagnostic/2026-08-10` | protocol-matched alpha timing | `(7,7)/(8,8)/(9,9)`，2--5L x3 | fixed CNT、20 ns、1 ps；Gaussian log-OLS | **可引用**；不含 `(17,0)` |
| `alpha_z_crosschirality_C075/2026-08-10` | C0.75 alpha sensitivity | `(7,7)/(9,9)` 2--5L x3；`(17,0)` 1--3L x1 | 20 ns/1 ps（17,0 异质） | **限定**；17,0 仅 timing evidence |
| `_quarantine/superseded/alpha_z_crosschirality_10L_100fs_10ns/2026-08-10` | 10L C0.75 sensitivity | 多 chirality，10L | 100 fs/10 ns | **已隔离**；仅保留作历史单点敏感性，不能作为 time-scaling 证据 |
| `alpha_z_trailing_decade/2026-08-09` | trailing-one-decade alpha | `(8,8)` 为主 | complete window 5--200 ps，per-replica OLS | **可引用**；是 alpha 主估计器之一 |
| `crosschirality_alpha_collapse_2L10L/2026-08-10` | cross-chirality alpha collapse figure/data | 多 chirality | 汇总层 | **展示/比较层**；溯源到对应 MSD/alpha 包 |
| `finite_size_collapse_2L10L/2026-08-10` | `(8,8)` MSD/alpha collapse、beta scan | `(8,8)`，2--5L、10L | 5--200 ps；10L cadence不同 | **可引用、协议分离** |
| `transverse_mode_88_5L_screen/2026-08-18` | transverse kz current spectrum | `(8,8)`，5L | transverse components/selected modes | **机制筛查**；不能推广到 chirality matrix |

### 同项目、但不在上述 37 包中的已完成专题资产

| 位置 | 内容 | 参数与状态 |
|---|---|---|
| `full_static_longitudinal_4chirality_20260810/archive/` | radial density、longitudinal ACF、helical density modes | `(7,7)/(8,8)/(9,9)/(17,0)`，5L rep1；**静态/筛查** |
| `full_kz_transverse_4chirality_20260809/archive/` | axial/transverse kz spectra、peak metrics | 四 chirality，5L rep1；**筛查，缺重复** |
| `exploratory_kz_transverse_localtest_20260809/` | 小规模 transverse local test | exploratory，不进入定量比较 |
| `results/flexible_fixed_collective_modes/2026-08-10/` 与 `results/flexible_fixed_4L_20260809/` | fixed vs flexible CNT 的 MSD、CJJ、DSF、cross-section/PSD | `(8,8)` 4L；flexible 为单轨机制筛查，**不可与固定壁多 replica 结论合池** |

## 3. 参数空间闭合审计（按物理 k）

### 3.0 Correction: 10L raw CJJ exists

The previous index understated the 10L CJJ coverage. `(8,8)` has three
replica-resolved 10L `CJJ_alln.csv` files with `n=1..10`, `lag_ps`, `k_inv_A`,
absolute CJJ and normalized CJJ. They are at
`assets/library/collective_dynamics/cjj_k_t/88_L2L10_allmodes_2026-08-11/source/`.
The physical-k mapping and exact matched-k bins are maintained in
`assets/library/collective_dynamics/cjj_k_t/CJJ_K_INDEX.md`. Cross-chirality
10L raw CJJ n=1..9 is separately available for 7,7/9,9/17,0.

### 3.1 CJJ(k,t)、色散与电流谱

| 目标参数空间 | 已有充分源 | 真正待补 | 不应误报为待补 |
|---|---|---|---|
| `(8,8)` all-mode CJJ | 5L n=1..80 x3，含 CJJ(t)、lobe、current spectra、LA/TA 色散 | 若需跨尺度色散，补充与 5L `k` bins 重叠的独立 source（任一 L 均可）及更低 k；当前 10L 资料须按其 cadence 单列 | 5L 以外每一个 L 的同名 n；若物理 k 已被 5L 覆盖，则不是缺口 |
| `(7,7)/(9,9)/(17,0)` low-k CJJ dispersion | 10L n=1..9 x3，100 fs/10 ns | 各 chirality 至少一个 10-fs、multi-k source，或能够覆盖相同 physical-k bins 的任意 L source；`(8,8)` 也需以同一协议做桥接 | 逐 length 的 n=1；不同 L 可对应同一 k |
| n=1 current correlation spectrum | `(8,8)` 3L/4L/5L/10L | 2L 的 matching local CJJ；若目标是跨 chirality，也缺 7,7/9,9/17,0 的同定义谱 | 在已有 3--10L 中，为“盒长整齐”重复同一 k |
| LA/TA mode splitting | `(8,8)` 5L x3 | 其它 chirality 的 >=3 replica 同定义 spectra；需要同 physical-k bins 来判断 chirality effect | 单一 `(7,7)` 5L screening 当作 matrix completion |

### 3.2 density / DSF `S(k,omega)`

| 目标参数空间 | 已有充分源 | 真正待补 | 不应误报为待补 |
|---|---|---|---|
| connected DSF k1 | `(7,7)/(9,9)`，2--5L x8，10 fs/2 ns | `(8,8)` 和 `(17,0)` 的可比较 position-resolved matrix；更低 physical-k 或 matched-k bin 的资料 | 每个盒长都做 k1；k1 随 L 改变，不能替代同-k 覆盖 |
| `(8,8)` density/current S(k,w) | 5L n=1..20 x3，10L per-k 谱 | 5L 与 10L overlap physical-k bins 的同 estimator比较；若要 chirality trend，补 7,7/9,9/17,0 | 只因为 L=2/3/4 不在表中而补；先检查 k bin 是否已有 |
| DSF--VDOS matched lambda | `(8,8)` L2--L10 包 | 各 length/source 的 metadata 一致性审计；扩展 chirality 需 position+velocity 同时存在 | 只有 velocity 的 56 条资产；它们不能“补 DSF” |

### 3.3 MSD / alpha / VACF

| 目标参数空间 | 已有覆盖 | 真正待补 | 主要限制 |
|---|---|---|---|
| `(8,8)` unified 10-fs Cvv/alpha | 2--5L x3、10L x8 | 若需 1L/8L 或高统计 lobe，按既定新计划取得独立构型起点 | 2--5L/10L replica 数不同 |
| cross-chirality high-resolution Cvv | 7,7/8,8/9,9，2--5L x8 | 17,0；1L/8L/10L 10-fs 同定义资产 | 现有 200 ps archive 不能回答 500 ps--1 ns smooth tail |
| cross-chirality alpha time scaling | 7,7/8,8/9,9，2--5L x3（20 ns/1 ps） | 17,0 的 protocol-matched多 replica 2--5L；如要 10L 关系，需要协议对齐 | 10L 100 fs/10 ns 不可直接与 1-5L 10 fs/1 ns 池化 |
| flexible vs fixed 机制 | 8,8 4L | flexible 多 independent starts 且与 fixed 同 estimator/cadence/mode | 单条 flexible trajectory 只能机制筛查 |

## 4. 重叠、冗余与协议不一致审计

| 家族 | 重叠来源 | 处置 | 原因 |
|---|---|---|---|
| 5L `(8,8)` spectral family | `88_5L_LA_TAr_TAtheta_dispersion`、`*_low_frequency_signed_Skw_CJJ`、`*_per_k_semilog`、`*_semilog`、`*_mode_weight`、`*_Wk_block_convergence` | 归为一个 source family；保留多个视图，不重复计为独立证据 | 大多来自同一 3 条 10-fs trajectory segments，只是 k/omega window、channel、normalization 或统计摘要不同 |
| 10L `(8,8)` per-k spectra | standard 与 `pm10` | 并列保存，明确 analysis-definition variant | mode/phase convention不同，不能默认等价或平均 |
| CJJ lobe assets | all-mode 5L、n=1 t/L、10L cross-chirality、absolute CvJ | 按 `CJJ` 定义、normalization、k 与 protocol 分层 | absolute/per-water/normalized 是不同 observable，不能合为一个 exponent |
| DSF | `dsf_k1_all_10fs_assets` 与 `dsf_kmin_10fs_asset_matrix` | 后者保持 superseded provenance | interim transform 已被明确拒绝；不可与主矩阵并列 |
| alpha_z | trailing-decade、C0.75、functional 1-ps、20-ns collapse、10L C0.75 | primary estimator与 sensitivity/diagnostic 分层 | OLS window、cadence、长度/时长及输入 provenance 不同 |
| VACF | 10-fs 8-rep archive、旧 weak-NH 1--5L/10L、unified 8,8、counterfactual projection/band-removal | 原 VACF、协议对照、modal counterfactual 三层分开 | COM 定义、cadence、duration、replica 深度和 filtering 均可改变曲线 |
| static/transverse | full static、full kz transverse、5L transverse screen、flex/fixed cross-section | 保留为 screening family | 多数只有 5L 或 single replica，不构成 chirality-length 定量矩阵 |

## 5. 建议的统一主索引结构

1. **波矢索引表**：未来每项 CJJ/DSF/S(k,w) 记录实际 `Lz`、`n`、physical `k`、omega resolution、window、normalization、source count；查询以 k-bin 为先。
2. **协议指纹**：固定字段为 CNT fixed/flexible、thermostat、momentum removal、rethermalization、coordinate/velocity COM 定义、cadence、duration、origin/block/window、replica 类型。
3. **主/辅版本标签**：`primary`、`protocol-distinct`、`sensitivity`、`counterfactual`、`screening`、`superseded`。展示图不再按文件名猜测其资格。
4. **补算优先级**：先补真正缺少的 physical-k、chirality 或独立 source，而不是机械补齐 L 标签；每次补算前先检查现有 k bin 与协议指纹。

## 6. 本报告的限制

- 本报告是落盘资产包的结构化清单，不把远端仅提交任务计为“已有分析”。
- 部分 2026-08-19 谱包缺少 README；本报告仅依据其 metadata、FINISHED 标志、源表和文件名给出资格，不虚构未写入的物理解释。
- `stage_*` 中的配置/脚本若没有完成输出，已从“已有分析”排除；它们可在后续单列为待执行计划。

## 7. Dump 来源与 F 盘迁移链路审计

### 7.1 证据等级

- **P1 — 显式且本地可访问：** 结果 metadata/manifest 给出 dump 路径，且该路径或已迁入 F 盘的对应文件当前可访问。
- **P2 — 显式远端、F 盘有同 campaign：** 结果 metadata 给出 CCFEP `/lustre/...` 路径，F 盘存在同 campaign/命名的 archive；本轮未逐文件 hash 比对，因此不能写为逐字节已验证副本。
- **P3 — 派生链：** 当前包以另一结果包的 CSV/JSON 为直接输入，需沿父包继续追溯至 P1/P2 dump。
- **P0 — unresolved：** 仅有图或脚本，尚未找到能证实具体 dump 的 metadata/manifest；不能作为来源闭合资产。

### 7.2 F 盘已发现的本地迁入 archive

`F:/ccfep_gcmc_archive_20260814/` 是本项目最关键的本地 archive 根。当前可访问 7 个 campaign 根、336 个 dump、约 1.09 TiB：

| F 盘 campaign 根 | dump 数 / 体量 | 可对应的分析资产 |
|---|---:|---|
| `remote_stage_88_length_scaling_20260719/` | 34 / 121.6 GiB | `(8,8)` 1--5L 原 length-scaling 的 O 10-fs/1-ns assets；5L CJJ/S(k,w) 主 family 的元数据直接引用同名 H 盘 raw-case 文件 |
| `viscfric_length_88_RH75_20260731/` | 216 / 339.3 GiB | `(8,8)` matched/rethermalized weak-NH、cross-chirality 1--5L 以及多项 thermostat/current/DSF 分支 |
| `viscfric_length_77_88_RH75_20260806/` | 24 / 36.6 GiB | `(7,7)/(8,8)` 2--5L x3 long weak-NH controls |
| `viscfric_length_all_chirality_RH75_20260806/` | 12 / 411.0 GiB | 四 chirality 10L x3 的 100-fs/10-ns long weak-NH inputs及其 restart provenance |
| `stage_vacf_tail_8_8_L2L10_8rep_weaknh_zvz_20260812/` | 40 / 59.5 GiB | `(8,8)` 2--5L/10L 8-rep oxygen `id,z,vz` VACF-tail inputs |
| `transverse_velocity_5L_10fs_weakNH_nomom_4chirality_20260808/` | 8 / 71.2 GiB | 四 chirality 5L transverse `kz` screening |
| `cnt_flexible_8_8_4L_nve_cntcsvr_20260808/` | 2 / 46.6 GiB | flexible `(8,8)` 4L mechanism-screening branch |

### 7.3 结果包到 dump 的已复原链路

| 结果家族 | 直接输入 dump / 直接父资产 | 原始远端 provenance | 当前本地位置 | 等级 / 注意事项 |
|---|---|---|---|---|
| 5L CJJ、LA/TA 色散、5L signed/per-k `S(k,w)`、`W_k` | `nvt20ns_8_8_RH75_5L_rep{1,2,3}_oxygen_10fs_1ns.dump` | `NVT20ns_5xL_8_8_RH75..._20260719` campaign | 元数据直接写 H 盘 raw-case；F 盘 `remote_stage_88_length_scaling_20260719/` 存在同名 rep1--3 | **P1**（两个本地路径可访问）；F/H 是否逐字节相同仍待 hash receipt |
| 5L CJJ mode-weight / semilog 图 | `CJJ_all_modes_per_replica.csv` 或 `low_frequency_signed_spectra_ensemble_mean_sem.csv` | 继承上一行 | 同一结果包的 `derived_data/` | **P3 -> P1**；不得作为另一套独立轨迹 |
| cross-chirality Cvv/Cvv-ODE 1--5L | `RETHERM_WEAKNH_NOMOM_10FS_<chi>_L<1..5>_rep<1..2>_oxygen_10fs_1ns.dump` | `/lustre/.../rethermalized_weakNH_nomom_cross_chirality_1L5L_2rep_20260804/` | F `viscfric_length_88_RH75_20260731/.../rethermalized_weakNH_nomom_cross_chirality_1L5L_2rep_20260804/` | **P2**；metadata 保存远端绝对路径，F archive 有同 campaign，需相对路径+字节/Hash 完成迁移验证 |
| cross-chirality 10L Cvv/CJJ 色散 | `LONG_WEAKNH_NOMOM_<chi>_L10_rep<1..3>` 的 O 100-fs/10-ns dump | `/lustre/.../long_weakNH_nomom_10L_100fs_10ns_3rep_20260806/` | F `viscfric_length_all_chirality_RH75_20260806/.../long_weakNH_nomom_10L_100fs_10ns_3rep_20260806/` | **P2**；同协议 10L，不是 1--5L 10-fs 的可合池副本 |
| `(8,8)` 10L 10-fs/8-rep Cvv tail | `VACF_8_8_L10_1ns_10fs_rep<1..8>.oxygen_id_z_vz_10fs_1ns.dump` | `/lustre/.../stage_vacf_tail_8_8_L10_8rep_weaknh_nomom_1ns_10fs_20260819/` | 结果包 metadata 已记录远端 source；本轮在 F archive 仅确认 older L2--L10 8-rep campaign，未确认该 2026-08-19 新 campaign 的完整 F 副本 | **P2/P0**；不得假定 F 中 `100fs` 文件就是本 10-fs input |
| `(8,8)` 2--5L 10-fs unified Cvv/alpha | 3 条 rethermalized O 10-fs/1-ns dump/length | 由 README/manifest 标为 locally migrated trajectories | F `viscfric_length_88_RH75_20260731/matched_10fs_rethermalized_weakNH_no_momentum_1L5L_20260804/` | **P2**；需将每 length/rep 的 metadata 补入统一包 |
| `(7,7)/(9,9)` 2--5L 10-fs/2-ns x8 Cvv/DSF | `VACF_<chi>_<L>_2ns_10fs_rep<1..8>.oxygen_id_z_vz_10fs_2ns.dump` | `/lustre/.../stage_vacf_tail_crosschirality_7_7_9_9_L2L5_8rep_weaknh_nomom_2ns_10fs_20260818/` | 当前结果 metadata 记录远端 source；本轮未在 F 盘发现同名 2026-08-18 campaign 根 | **P2（F 未定位）**；本地 F 副本状态为 unresolved |
| transverse `kz` 4 chirality | full-water 10-fs `x,y,z,ix,iy,iz,vx,vy,vz` | 2026-08-08 transverse continuation | F `transverse_velocity_5L_10fs_weakNH_nomom_4chirality_20260808/` | **P1 campaign-level**；需逐 case metadata 连接到具体 dump |
| flexible/fixed 4L | flexible NVE/CNT-CSVR 与 fixed comparator dump | flexible campaign 2026-08-08 | F `cnt_flexible_8_8_4L_nve_cntcsvr_20260808/`；fixed 输入另属其原 campaign | **P1 flexible / P0 fixed pairing**；不能把 flexible archive 当作 fixed 的来源 |
| static longitudinal / full transverse archive | 各 5L rep1 source dump | archive README/summary 指向 output | local result/archive，原 dump 需沿 `output/*/summary.json` 继续定位 | **P3** |

### 7.4 必须补齐的 provenance 工作，而非模拟工作

1. 为所有 `P2`/`P0` 包建立 `provenance.json`：`result_package`、`source_dump_remote`、`source_dump_local_F`、相对路径、字节数、SHA-256、frame count、字段列表和迁移核验日期。
2. 首先核验 5L `(8,8)` 的 F/H 三个 10-fs dump，以及 1--5L cross-chirality 与 10L 四-chirality campaign 的完整相对路径/字节清单；在 hash 或至少 byte-size 清单完成前，F archive 只能称为**候选本地副本**，不是已验证替代来源。
3. 为 `vacf_alpha_10fs_L2L10_unified_2026-08-20` 写入 2--5L 每个实际 F dump 路径和 10L 8-rep source 路径。当前 manifest 只有统计级 `Lz/k1` 与 replica 数，不能独立复原输入。
4. `source_dump` 若仍是 `/lustre/...`，不得在后续报告中改写成“F 盘输入”；应写成“远端原始来源；F 本地迁入候选/已验证副本”。

## 8. 下一步：先整合已有数据，再决定补算

### 8.1 不需要新模拟、应立即完成的统一分析

| 优先级 | 工作 | 使用已有资产 | 交付与接受条件 |
|---|---|---|---|
| P0 | CJJ physical-k master table | `(8,8)` 2--10L `CJJ_alln.csv` x3、5L n=1..20 x3、cross-chirality 10L raw CJJ x3 | 每行含 source/replica/Lz/n/k/lag/CJJ_abs/CJJ_norm/protocol fingerprint；同-k 先分 protocol 再做 mean/SEM，禁止按 L 盲目平均 |
| P0 | ISF `F/Fs/Fd` source recovery and recomputation | F archive 的 position-bearing 10-fs dumps；已有 Fig2 matched-k output | 建立 physical-k source matrix；每个 retained bin 有 total/self/distinct、replica SEM、coordinate definition、time window。先做 `(8,8)`，再扩展到已有 7,7/9,9 position assets |
| P0 | CvJ master table | `absolute_cvj_first_negative_lobe`、cross-chirality 1--10L raw outputs、F archive input manifests | 强制分 absolute/per-water/normalized CvJ；记录 current definition、COM subtraction、lobe zero crossings、units、k 或 n；不把它们压成一个“CvJ exponent” |
| P0 | F/H provenance receipt | 2026-08-14 F archive 与 H/result inputs | relative path + byte count；关键主源再 SHA-256；将 verification result 写入每个 `provenance.json` |

### 8.2 只有在 P0 结果明确显示缺口后，才值得新增模拟

- **CJJ/ISF/DSF:** 补缺的是 protocol-matched physical-k bin 或缺失 chirality，不是“每个 L 都再跑一遍”。优先 7,7/9,9/17,0 的 position+velocity 10-fs assets，使 CJJ、ISF 和 DSF 可从同一源重算。
- **CvJ:** 若现有 absolute/per-water/normalized 总表在同一协议下没有足够 replicas，补独立配置起点；速度种子只能给出条件噪声，不能替代构型独立性。
- **10L:** 只有当要比较 sub-ps morphology 或与 10-fs bins 池化时，才需要 10L 10-fs production；现有 100-fs 10L CJJ 对低频/色散仍是有效源。

任何昂贵新资产须先提出矩阵、实际 k-bins、replica/构型策略、cadence、cost 与硬 acceptance criteria，获批准后再提交。

### 8.3 明显去重/隔离候选（本轮不删除）

| 候选 | 建议 | 不能直接删除的原因 |
|---|---|---|
| `dsf_kmin_10fs_asset_matrix/2026-08-18` | 保留在 `superseded/` provenance；从默认检索和展示入口移除 | 已有 `SUPERSEDED.md`，但仍记录中间方法与输入选择 |
| 已移入 `_quarantine/in_package_superseded/` 的 unverified 10L 输入与 one-decade-window 附件 | 保留为审计附件；默认排除 | 需要能解释为何旧版本被替代 |
| `88_5L_per_k_semilog...`（无 FINISHED） | 标 `unqualified`，不纳入 primary library | 可能是有价值的可恢复渲染；不能因无完成标记而误删 |
| F 盘与 H/项目树同名 CJJ raw files | 完成 relative-path + byte/Hash receipt 后，将一个设为 canonical source、另一个标为 verified archive link | 两处目前都可访问，但未完成逐文件内容等价证明 |
| `remote_fetch/`、`heartbeat_fetch/` 与已归档 `results/.../remote_raw/` 的重复取回物 | 完成 hash 与引用搜索后，仅保留 authoritative result copy + archive receipt | 可能含不同 fetch 时点、缺失/修复版本或独有日志 |
| `recovered_corrected/` 与同名旧 projection outputs | 只让 corrected 版进入 primary；旧版归 archive | 需保留 correction provenance，不能硬删 |

**不应当去重删除：** PNG/PDF/SVG/TIFF 交付格式；absolute/per-water/normalized CvJ；standard 与 `pm10` mode convention；5L full-spectrum、low-frequency、per-k、mode-weight、block-convergence 视图。这些不是等价副本，分别承载格式、定义或诊断信息。
