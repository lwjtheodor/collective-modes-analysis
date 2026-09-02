# 项目执行约束：分析资产归档

适用范围：`H:/gcmc_explore/translational_anomaly/02_isf_collective_modes` 下的所有任务。

## 代码与资产治理（2026-09-02 起）

- 本项目的 Git 仓库只追踪治理文件、根目录方法/资产入口和 `scripts/` 中的可复用源代码候选；基线中的 `scripts/` 文件在通过审阅/回归前均不可称为 canonical。它不追踪 dump、restart、调度输出、`results/`、`assets/`、`remote_fetch/` 或 `stage_*` 的大规模运行产物。Git 提交是**代码历史**，不替代分析完成或科学资格。
- `stage_*`、`remote_fetch/`、`heartbeat_fetch/` 为不可变执行快照；它们可以提供 provenance，却不能因同名或被复制到结果包而自动成为 canonical script。脚本升级必须先通过 `governance/inventory/<日期>/` 的 hash/同名变体审计，再在 `scripts/` 中明确确定主线实现。
- 主线的物理边界固定为：`LONG-ISF`（`F/Fs/Fd`）、`LONG-CJJ`（去 instantaneous O axial COM 的 axial current）、`CYL-LT`（`q=(kz,m/Rcnt)`，保留 `Jr`、`Jtheta`、`CLT`、`CTL` 分支）、`SELF-VACF`（peculiar VACF 与一致积分 MSD/alpha）、`STATIC-VERTEX` 和 `QA-PROVENANCE`。具体定义、迁移顺序与禁混规则见 `governance/MAINLINE.md`。
- 新的权威结果包除原有 README/manifest/QA 外，必须记录 canonical script 相对路径、Git commit、输入 asset ID/远端绝对路径和输出 schema；历史复制脚本保留在该包内作复现快照，但应注明 `historical`、`superseded` 或 `exploratory` 状态。
- 自 2026-09-02 起，新的 raw-dump 分析主线入口为 `scripts/collective_modes_cli.py` 与 `scripts/collective_modes/`。它先执行 field/cadence/protocol capability gate，再输出可读 CSV 与 metadata：`audit`、`isf`、`current`、`vacf`、`fit-current`、`construct`、`plot`。旧脚本和结果包内 scripts 只作 provenance，除非审查后明确升级。
- `current` 以圆柱基 `(n,m)` 构造并**分别**输出 `Jz,Jr,Jtheta,L,Tinplane,Tr`、时域 CJJ、频域 periodogram 和 ordered `C_AB(tau)=Re<delta J_A(t+tau)delta J_B(t)^*>`；`C_LT` 与 `C_TL` 绝不以对称性互换。`isf` 同时输出 `F_total,F_self,F_distinct`。`construct` 只接受外部 measured `W(n,m)`，无自由总幅度地写 `sum W Fs Phi_J`，有 direct VACF 时必须输出残差表。
- 自动协议识别只能确认 dump 的字段能力；water-only dump 无法可靠判定 implicit/explicit CNT，必须在 profile 明确 wall model。implicit case 的 `Rcnt` 必须来自 CNT field/protocol metadata；flexible explicit CNT 的 wall-relative frame 必须有 CNT atom position/velocity，禁止退回 box centre 猜测。

## 完成定义

形成新的分析图、表、拟合、谱、统计汇总或派生数据后，任务**尚未完成**。在下列归档完成前，只能称为 staged / unqualified：

1. 将资产包写入 `results/collective_mode_response/<主题>/<YYYY-MM-DD>/`，包含 README、可重现脚本、紧凑源表、metadata/manifest、图件以及 QA 或完成证据。对大于约 10 MB 的数值 CSV，保留 CSV 作为可审计源，并在同一 `derived_data/compact/` 生成压缩列式 HDF5（或已有协议规定的等价格式）、manifest 与读取说明；不得只保留不可读的二进制文件。
2. 更新项目根目录 `assets.md`：登记权威路径、observable 定义、体系/length/replica、cadence/time window、输入来源、误差表示、结论边界和状态。
3. 若需要直接展示，再将 PNG/PDF/SVG/TIFF 和紧凑源表放入 `assets/library/<主题>/`；展示副本必须能追溯到主档案。
4. 对 CJJ(k,t)、ISF/DSF、S(k,omega)、色散、PSD/VDOS 或任何可复现完整分析包，额外在 `assets/library/collective_dynamics/<observable>/<版本>/source/` 建立无复制入口，并在 `assets.md` 登记版本和适用边界。
5. 对所有与波矢有关的资产，登记实际 `Lz`、mode index `n`、physical `k=2*pi*n/Lz`、频率/时间窗和 source count；补算审计以可比 physical-k bin 为主，不以 box label 的齐全度替代 k-space 覆盖。
6. 每个可引用资产包必须保存 dump provenance：远端绝对来源、实际本地输入（特别是 `F:/ccfep_gcmc_archive_20260814/`）、相对路径、字节数或 SHA-256、frame/field audit 与迁移验证日期。未经相对路径加字节/Hash 核验的 F 盘文件只能标为候选本地副本。
7. 紧凑表须记录源 CSV 的 SHA-256、行数、列/类型、压缩格式与精度；转换后至少核验行数与源哈希。默认数值工作格式为按列压缩的 HDF5，`float32` 仅可作为 CSV 仍被保留时的紧凑工作副本；需要数值无损时改用 `float64` 并登记理由。
8. 新增或重新发现 dump 时，先更新 `results/collective_mode_response/dump_asset_inventory/<日期>/` 的逐文件索引：远端/本地路径、手性、Lz/nominal L、采样间隔、时长、frame count、`ITEM: ATOMS` 字段、协议指纹、字节数/Hash 与 provenance 状态。字段不足的 dump 必须明确限制可生成的 observable；不得以 `id,z,vz` 资产补造 transverse、分子或 force 分析。

## 版本与重合处理

- 不覆盖、不静默删除历史资产；先检查相对路径、字节大小和 SHA-256。
- 同一图的多格式文件不是重复项。协议不同、输入不同、脚本修正或统计口径变化的结果也不是可互换重复项。
- 对真正被替代的产物，保留旧版本并在其 README/`assets.md` 标记 `superseded`；对未完成或探索性产物标记 `unqualified` 或 `exploratory`。
- 不可把不同 cadence、时间窗、COM/velocity 定义、replica 类型或 thermostat/momentum 协议的曲线混合为同一结论。

## 可引用门槛

只在源数据、完成日志/标志、README/manifest 和 QA 证据均已存在时，才可在 `assets.md` 标为“已归档/可引用”。队列状态、作业已提交或仅生成 dump/脚本都不是分析完成证据。

## C99 implicit-CNT 纵向动力学展示与统一重分析合同（2026-08-29 起）

- 所有新生成的 VACF、MSD、`alpha(t)`、`Phi_J(k,t)`、`F_s(k,t)` 及 VACF 构造/外推时间图，展示横轴一律使用 **log time**，统一可视范围为 **1--500 ps**。`t=0` 仅保留在原始表和积分/归一化审计中，不作为图中的对数坐标点。
- C99 的统一重分析必须覆盖 N200/N400/N800/N1600/N2400/N3200 与每个 case 的四条 velocity seeds；可比较的纵向 `m=0` 物理波数使用 `k=2*pi*n/Lz`。默认共同上限为 `k<=0.3141592654 A^-1`，对应 `nmax=5,10,20,40,60,80`，不得以相同离散 `n` 代替 physical-k 匹配。
- 每条生产轨迹须从 full-water `id mol type x y z vx vy vz` dump 中选择 oxygen type 1，并先逐帧去除 instantaneous oxygen axial COM velocity。产物至少包含 raw/normalized CJJ、`CJJ(0)`、`Phi_J`、`F_s`、`K,c,a,W`、direct peculiar VACF 及由其一致积分得到的 MSD/`alpha`；N3200 的时域 CJJ 与 vertex 不得再以 Welch-only 数据替代。
- 任何跨盒长 VACF 构造都须以 `sum_n W_n F_s(k_n,t) Phi_J(k_n,t)` 的无自由总幅度形式进行，并将 static weight、self dephasing、carrier/phase、damping 误差分别列出。只有 completed per-replica outputs、normal log endings、frame/cadence audit 与 aggregate QA 到位后，才能更新 assets 或讨论 N800/N3200 的预测性。
- 对 PBS array 作业，`jobinfo -c <jobid>` 的单-ID 查询不作为状态判据。应运行完整 `jobinfo -c` 后筛选目标 job ID，同时以预定远端 `output/<case>/rep<seed>/` 的非空成员日志、`rep_arrays.npz`、`mode_summary.csv`、`metadata.json`、`SUCCESS.txt` 和 PBS 正常末行判定真实进展与完成。
- C99 的横纵向 current 分解必须把 `T_r=J_r` 与 `T_theta=J_theta` 保持为不同支；按正时延定义的有序交叉核为 `C_AB(k,tau)=Re<delta J_A(t+tau) delta J_B(t)^*>`，因此 `C_LT` 与 `C_TL` 必须分别储存、绘制和审计，绝不得以自相关或同一时刻近似替换。不同 cadence 的 cross spectrum 也必须保留实/虚 quadrature，不能静默相加或平均。
