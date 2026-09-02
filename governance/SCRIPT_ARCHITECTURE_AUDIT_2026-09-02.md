# 集体动力学资产与脚本架构审计（2026-09-02）

## 审计范围与结论

本审计覆盖项目根目录及 `scripts/`、`stage_*`、`remote_fetch/`、
`results/collective_mode_response/`、`assets/library/` 中可发现的脚本文本。
机器盘点记录 **11,521** 个文件、**1,579** 个脚本快照、**210** 个字节完全
相同的脚本簇和 **172** 个同 basename 变体簇；逐文件证据在
`governance/inventory/2026-09-02/`。

“完全相同”只表示同一脚本被镜像保存，最常见链路是
`stage -> results package/scripts -> assets/library/.../source/scripts`；它不是
删除历史副本的证据。本报告只对根 `scripts/` 的 107 个源码候选提出未来主线
角色。所有其他位置的副本均保留为执行/再现快照。

## 资产总图

| 层 | 当前内容与规模 | 权威用途 | 整理决策 |
|---|---|---|---|
| 原始模拟/取回数据 | `remote_fetch/` 约 1.9 GiB，另有 `results/` 内嵌大输入；CCFEP 根为 `/lustre/home/users/ewu/vb_gcmc/MD` | 原始 dump、restart、PBS 结束证据 | 不迁移、不进入 Git；只以路径、字节/哈希、fields/frame/cadence manifest 管理 |
| stage | 根 `stage_*` 与远端同名 stage | 不可变提交/实验快照 | 永久保留、Git ignore；不得以文件名决定主线 |
| 结果主档案 | `results/collective_mode_response/<topic>/<date>/`，约 54 GiB | 可引用的 README、紧凑表、QA、图、源脚本副本 | 保持按 observable/topic + 日期；以后记录 canonical script + commit |
| 展示层 | `assets/library/` 中 PNG/PDF/SVG 和 source 入口 | 直接浏览和论文使用 | 仅是主档案的展示入口，不能替代数值 authority |
| 代码候选 | 根 `scripts/` 107 个源码候选，另有根目录 2 个 legacy helper | 未来跨 case 复用的唯一代码来源 | 先版本化、再按本报告重构；目前均未自动认证为 canonical |

### 主要科学资产族

1. **轴向 density/current：** ISF-01 与 CJJ-01--05 覆盖 `(8,8)`
   matched physical-k ISF/CJJ、10L 全模态负瓣、相关谱、5L 的
   `Jz/Jr/Jtheta` 色散，以及跨手性 10L 色散。必须保留 `F/Fs/Fd`、
   raw/normalized CJJ 与实际 `k=2*pi*n/Lz`。
2. **显式 fixed-CNT 低 k 纵向链：** CJJ-06--10 是 10L DHO、beating
   筛查、长窗复核和 phase coherence；CJJ-07 是已 superseded 的短窗
   primary beating test，仍是模型识别历史证据，不能删除。
3. **圆柱/螺旋与横向链：** CJJ-11--13、CJJ-14/16/17/19--33、ISF-02
   包含 `m=0` 周向、`m>0` 螺旋、`CLT/CTL`、静态 vertex、有限盒权重和
   10L-to-20L transfer。共同核心是 `(kz,m/Rcnt)`，但结论/拟合不能混为一支。
4. **self transport：** `vacf_*`、MSD/alpha、cross-chirality Cvv/CvJ 和
   protocol-independence packages。权威量是 selected-water instantaneous
   COM-subtracted peculiar VACF 及其一致积分 MSD；seed SEM 与独立构型
   不确定度必须分开。
5. **implicit C88/C99：** C88 100 fs/6 ns matched-length current/VACF
   主链与 C99 多 cadence current/static-W/constructibility/ordered cross
   current 链。C99 的 1/10/100 fs 分层不可拼为一个未标注来源的谱或 VACF。
6. **QA 与档案：** dump field/cadence inventory、HDF5 compact conversion、
   source/display hash、NVE temperature/energy gate 和 protocol audit，是全部主线共同依赖。

## 不能被一个开关抹平的输入等级

| 输入等级 | 最低字段 | 可构造 observable | 明确禁止 |
|---|---|---|---|
| `axial_oxygen` | `id,type,z[,iz],vz` | axial `Fs/F/Fd`、axial CJJ（有 `vz` 时）、axial peculiar VACF、VACF-MSD/alpha | `Jr/Jtheta`、helical `m>0`、分子转动、CNT wall frame |
| `oxygen_3d_velocity` | `id,type,x,y,z,vx,vy,vz` | axial/cylindrical oxygen-site currents、`Jr/Jtheta`、coherence、3D O-COM subtraction | 水分子 COM/转动、CNT motion correction |
| `full_water` | `id,mol,type,x,y,z[,ix,iy,iz],vx,vy,vz` | 分子映射、O/H/水 COM、角动量、velocity/position frame audit、全部 oxygen-site量 | 无 CNT atoms/metadata 时推断 flexible-wall frame |
| `explicit_CNT_dynamic` | full-water + CNT type + CNT `x,y,z,vx,vy,vz` | 随 CNT 平移/转动的 wall frame、water-vs-wall relative velocity | 默认用 box centre 或 water COM 代替 CNT frame |
| `implicit_CNT` | water fields + `Rcnt`/axis protocol metadata | 固定 analytic cylinder 的 mode projection | 从 density peak 反推 `Rcnt`；输出 CNT momentum/force |

因此 `--implicit` / `--explicit` 不应是唯一开关。统一入口必须先读
`case.yaml`，验证 `available_fields`、`water_selection`、`wall_model`、
`axis_source`、`velocity_frame`，再拒绝物理上不可构造的命令。

## 当前脚本角色判定

### A. 优先重构为泛用主线的候选（17 个）

| 未来模块 | 当前候选 | 保留的算法价值 | 需要消除的局限 |
|---|---|---|---|
| `io.lammps_dump` + audit | `audit_dump_assets.py`, `extract_lammps_dump_frames.py`, `audit_fullwater_orbital_Lz.py`, `audit_oxygen_orbital_Lz.py` | header/field audit、有限读取、molecule/oxygen identity 检查 | 四套 reader 应成为一个 streaming parser 的 validators |
| axial ISF | `rebuild_axial_isf.py`, `rebuild_collective_isf_kseries.py`, `rebuild_isf_components_sampled.py`, `summarize_isf_2L_3L_demo.py` | all-origin `F/Fs/Fd`，ID 排序，protocol grouping，physical-k metadata | 合为 `isf`；manifest、mode list、lag/window 和输出格式 CLI 化 |
| current modes | `analyze_88_5L_full_dispersion.py`, `analyze_full_kz_static_longitudinal.py`, `analyze_kz_transverse_exploratory.py` | complex current、COM subtraction、`Jz/Jr/Jtheta`、FFT ACF/spectrum、cadence gate | 移除 oxygen type=3、fixed box centre、0.5 fs 和 expected-frame 的隐含假定 |
| self transport | `analyze_vacf_tail.py`, `analyze_vacf_integral_lockin.py`, `analyze_lowfreq_msd_loglog_decade.py`, `merge_vacf_88_L2L5_8rep.py` | axial lab/peculiar VACF、all-origin FFT、block diagnostic、lag/cadence metadata | 合为 `vacf` 与 `msd-alpha`；不静默合并 normalized/absolute/direct-MSD/protocol-distinct series |
| archive/compact | `compact_tabular_assets.py`, `inventory_collective_modes.py` | table hash/schema/row QA、源代码/资产 inventory | 保持独立工具；每个结果包必须调用 |

### B. 保留为物理专用插件或后验分析，不并入 dump 核心（25 个）

- **纵向 CJJ/DHO 拟合：** `analyze_10L_CJJ_dispersion_n1n6.py`,
  `analyze_10L_tdamp1ns_LA_gamma_n3_n8.py`,
  `analyze_implicitC88_LA_Gamma0_from_SJJ.py`,
  `analyze_implicit_matched_k_dho_damping.py`,
  `fit_88_10L_LA_linewidth_powerlaw.py`, `fit_88_10L_TAr_DHO_linewidth.py`,
  `fit_88_10L_TAtheta_zeropeak_linewidth.py`,
  `fit_explicit_LA_gamma0_plus_powerlaw.py`。可共享 `fit-current` 框架，
  但 LA/Tr/Ttheta 的 model、Nyquist 和 acceptance 必须在 profile 中明确。
- **implicit C88/C99 后验链：**
  `analyze_c99_implicit_vacf_matchedk_constructibility.py`,
  `analyze_implicitC88_transverse_SJJ.py`, `analyze_implicitCNT_TAtheta_linewidth.py`,
  `analyze_implicit_nve_thermo_drift.py`,
  `summarize_C99_static_vertex_lowk_block_audit.py`,
  `summarize_implicit_C88_static_vertex_partial.py`,
  `build_C99_CJJ53_multirate_extension.py`。C99 constructibility/static W/
  cross current 是现有 compact-table schema 的后验诊断，不能冒充 raw core。
- **explicit/flexible mechanism：** `analyze_88_5L_mode_weight_n1_n20.py`,
  `analyze_fullwater_Pz_time_acf.py`, `analyze_k0_axial_current_10L_8rep.py`,
  `analyze_multirep_tatheta_nve.py`, `make_flexible_fixed_mode_comparison_20260810.py`,
  `analyze_windowed_alpha_from_allorigin_vacf.py`。保留为 plugin/reference implementation。
- **cross-length protocol comparison：** `analyze_matched_k_5L_10L_protocolmatched.py`,
  `audit_matched_k_crosslength_LA_weights.py`, `compare_cjj_10fs_100fs_matched_k.py`,
  `compare_implicit_old_tatheta_fixed_exponents.py`。应读统一结果 schema，且拒绝 unmatched protocol pooling。

### C. 历史或一次性派生分析：保留复现，不作为新 case 入口

`aggregate_1L10L_alpha_cvj_nature.py`, `aggregate_alpha200_absolute_cvj_nature.py`,
`analyze_1L10L_lobe_timing_linearity_nature.py`,
`analyze_88_10L_LA_zvz_100fs_10ns_8rep.py`,
`analyze_88_5L_low_frequency_signed_skw.py`, `analyze_88_5L_Wk_block_convergence.py`,
`analyze_alpha_fit_sensitivity_nature.py`, `analyze_current_first_rebound_vs_k.py`,
`analyze_lowfreq_alphaz_cvj_association.py`, `analyze_lowfreq_cvj_n1.py`,
`build_vacf_alpha_88_L2L10_8rep_figures.py`, `compare_5L_10L_allmode_lobes.py`,
`compare_88_baseline_weaknh_and_10L_lobes.py`, `diagnose_88_10L_invk_lobe_scaling.py`,
`make_dsf_kmin_10fs_asset_matrix.py`, `reconcile_88_baseline_vs_weaknh_cjj_n1.py`,
`replot_88_5L_dispersion_summary.py`, `replot_88_5L_low_frequency_density_vs_LA.py`,
`test_10L_depth_log_vs_saturation.py`, `test_10L_lobe_area_lambda_collapse.py`,
`test_decom_msd_L_collapse_nature.py` 等均写死具体 case、日期、table 或假设。
它们的机制/图件 provenance 有价值，但新 case 应走 core + profile + report。

### D. 纯展示/归档：下沉为 archive snapshot（32 个）

所有 `plot_*`, `archive_*` 和带 `_nature` 的脚本主要固定一张图的 layout、
hard-coded input path 或文章叙事。保留方法是：原件留在结果包
`scripts_historical/`，新图才使用 config-driven `plot-response`；不以新画图器
追求逐像素重画旧图。

### E. 远端 shell launcher：提交快照，不进入分析主线（4 个）

`run_fullwater_Lz_audit_88L5_4rep_CCFEP.sh`,
`run_fullwater_Pz_audit_88L5_4rep_CCFEP.sh`,
`run_implicit_C88_N1600_TAtheta_CCFEP.sh`,
`run_remote_vacf_rep4_8_20260820.sh` 与根目录 legacy stage shell 均是
case-specific PBS launcher。以后应由 profile-driven PBS template 生成。

## 推荐的合并目录与命令线

```
scripts/
  collective_modes/
    io/lammps_dump.py        # single streaming parser + field/schema gate
    io/manifest.py           # case.yaml, source/path/hash/protocol records
    core/selection.py        # oxygen/water/CNT selection; id/mol stability
    core/frames.py           # box/CNT axis, unwrap, water/CNT COM frames
    core/correlation.py      # real/complex all-origin ACF, block/seed SEM
    core/modes.py            # rho, axial J, Jr/Jtheta, cylindrical L/T
    commands/{audit_dump,isf,current,vacf,msd_alpha,compact}.py
  plugins/{fit_current,static_vertex,compare_protocols}.py
  profiles/{explicit_fixed_88,explicit_flexible_88,implicit_C88,implicit_C99}.yaml
```

```text
collective-modes audit-dump --manifest cases.yaml
collective-modes isf --manifest cases.yaml --mode axial --n 1,2,3
collective-modes current --manifest cases.yaml --basis axial|cylindrical --n 1:20 --m 0:4
collective-modes vacf --manifest cases.yaml --component z|r|theta --velocity-frame lab|water-com|wall-relative
collective-modes msd-alpha --vacf-package <authority> --direct-msd <optional>
collective-modes fit-current --input <CJJ/SJJ package> --profile <yaml>
```

## 合并决策矩阵

| 候选合并 | 判断 | 理由 |
|---|---|---|
| 所有 dump readers | **必须合并底层** | 至少 ISF、full dispersion、kz transverse、VACF、orbital audit 各自手写 parser，是重复 bug 源 |
| axial CJJ + `Jr/Jtheta` | **合并为 `current` core** | 都是 selected-water velocity × phase；channel/basis 是参数 |
| `C_LT` 与 `C_TL` | **同一 engine、不同产物** | 共用 complex cross-correlation，必须分别输出，绝不以 symmetry 合并数值 |
| total/self/distinct ISF | **同一 `isf` command** | 同一 position/ID reader；`F`, `Fs`, `Fd` 必须一并输出 |
| VACF、积分 MSD、alpha | **同一 self-transport pipeline** | 一致积分是定义链；direct position-MSD 是显式第二输入/比较，不能替换 |
| full-water 与 oxygen-only | **同一 parser + feature gates** | 共有 frame/header/selection；字段不足时 full-water 专属命令必须失败 |
| explicit fixed 与 implicit | **共用 water-mode core，profile 区分** | 数学相同；`Rcnt` 来源、type mapping、wall frame/metadata 不同 |
| explicit flexible 与 fixed | **只共用 core，不共用默认 frame** | flexible 必须有 CNT group/relative frame，不能退回 fixed box centre |
| LA/Tr/Ttheta line fit | **共享 fitting framework，不共享默认模型** | radial DHO、theta Lorentzian、LA DHO 的 acceptance/Nyquist 不同 |
| C88/C99 vertex/constructibility | **暂不并入 raw core** | 是 compact-table 后验诊断，部分量不可识别 |
| 论文作图 | **不合并历史图；新增 generic style layer** | 历史图应精确复现，未来图使用 config-driven plotting |

## 实施顺序与“保留”的定义

1. 先写 `io/`、`core/` 和 `audit-dump`，以当前 compact fixture 测试；不改结果包。
2. 先迁移 ISF、axial/cylindrical current、VACF 三条 core；旧脚本留在包内并登记 canonical path + commit。
3. 再迁移 fit/compare plugins；每个 profile 有独立 acceptance 与 protocol guards。
4. 最后把一次性图/archive 脚本从根 `scripts/` 候选清单下沉为结果包历史副本；这不是删除。

`canonical` 指可扩展主线；`plugin` 指物理专用但可复用的后验分析；
`historical` 指仅复现已有包。任何未完成上述审查的文件仍是 `candidate`。
