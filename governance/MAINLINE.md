# Scientific mainline and consolidation map

This map makes protocol boundaries explicit before source files are unified.
It is a decision map, not a claim that every listed stage is completed or
scientifically qualified.

| Mainline ID | Observable family | Canonical physical definition | Current evidence/archive entry | Do not merge with |
|---|---|---|---|---|
| `LONG-ISF` | Total/self/distinct axial ISF and DSF | `rho_k=sum_j exp(i k z_j)`, `F`, `Fs`, `Fd=F-Fs`; compare equal physical `k=2*pi*n/Lz` | `assets.md` entries `ISF-01`, current C88/C99 self-distinct stages | VACF or a same-`n` cross-box comparison |
| `LONG-CJJ` | Axial collective current, time CJJ, spectra and damping | oxygen peculiar `v_z` after instantaneous O-COM removal; retain raw and normalized CJJ, `CJJ(0)`, `k`, window and cadence | `assets.md` CJJ series; `results/collective_mode_response/` | 1 ps screening as a high-frequency/damping result, or different physical `k` |
| `CYL-LT` | Cylindrical/helical longitudinal and transverse currents | `q=(kz,m/Rcnt)`, `L=(kz Jz+kθ Jθ)/|q|`, `T_inplane=(-kθ Jz+kz Jθ)/|q|`, `Tr=Jr` | `assets.md` CJJ-11 and C88/C99 LT stage family | treating `Jtheta` and `Jr` as one transverse channel, or identifying `CLT` with `CTL` |
| `SELF-VACF` | Peculiar VACF, VACF-integrated MSD and `alpha(t)` | all-origin peculiar VACF; MSD by the consistent double integral; direct-MSD comparison retained | `assets.md` VACF/MSD series and C88/C99 multirate stages | raw/global-COM velocity, truncated-lobe metrics, or direct-MSD-only claims |
| `STATIC-VERTEX` | Static weights and collective/self coupling | retain `K,c,a,W,Fs,Phi_J` with source cadence and rank/QA | `stage_C99_static_*`, `stage_implicit_C88_static_*` | empirical reconstruction treated as identified physical vertex |
| `QA-PROVENANCE` | dump field/cadence/frame audit, result manifests and archive qualification | exact remote path + local path + byte/hash + frame/field audit + completed log | `results/collective_mode_response/dump_asset_inventory/` and `assets.md` | scheduler acceptance, directory existence or a raw dump alone |

## Priority consolidation waves

1. **Foundation (this baseline):** version the governance rules and source
   inventory; do not move original files or change result status.
2. **Longitudinal current/ISF:** review exact-copy clusters and all
   same-basename variants, then promote one CLI per `LONG-ISF` / `LONG-CJJ`
   task into `scripts/longitudinal/` with regression inputs from compact
   archived tables.
3. **Cylindrical LT:** isolate geometry/projection code in
   `scripts/cylindrical/`; retain `Jr`, `Jtheta`, `CLT`, and `CTL` separately.
4. **Self transport:** isolate COM/peculiar velocity, correlation and MSD
   integration helpers in `scripts/self_transport/`; test against a compact
   archived fixture before replacing any stage script.
5. **Archive backfill:** add `SCRIPT_PROVENANCE.yaml` and explicit status to
   each result package touched by waves 2--4.  Only then create `v0.1`.

## Remote inventory scope (read-only baseline)

The CCFEP root is `/lustre/home/users/ewu/vb_gcmc/MD`.  The current project
references many immutable campaign and `stage_*` paths below it.  Remote
scripts must continue to use `/lustre/home/users/ewu/.conda/envs/HB_analysis/bin/python`
for Python analysis.  The local inventory extracts referenced remote paths;
the next remote pass should inspect only those distinct roots and emit a
compact manifest, never download dumps as part of a code review.
