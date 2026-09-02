from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'results/collective_mode_response/isf_88_L10_100fs_10ns_k1-k2_components_1ps/2026-08-20/analysis_5ns_lag'
COL={1:'#542788',2:'#2c7fb8'}; COMP=[('self',r'$F_{\mathrm{s}}(k,t)$'),('total',r'$F(k,t)$'),('distinct',r'$F_{\mathrm{d}}(k,t)$')]
def load(n):
    with np.load(DATA/f'ISF_88_L10_k{n}_components_mean_sem.npz') as z:return {x:z[x] for x in z.files}
def style(ax):
    ax.axhline(0,color='0.45',lw=1); ax.set_xlim(0,5000); ax.tick_params(direction='out',top=False,right=False,length=3); ax.spines[['top','right']].set_visible(False); ax.set_xlabel(r'$t\ (\mathrm{ps})$')
def main():
    plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.linewidth':1,'pdf.fonttype':42})
    d={n:load(n) for n in (1,2)}; out=DATA/'figures';out.mkdir(exist_ok=True)
    fig=plt.figure(figsize=(5.5,2.55)); ax=fig.add_axes([.12,.19,.84,.69])
    for n in (1,2):
      t=d[n]['time_ps'];m=d[n]['F_self_mean'];s=d[n]['F_self_sem'];ax.fill_between(t,m-s,m+s,color=COL[n],alpha=.17,lw=0);ax.plot(t,m,color=COL[n],lw=1.15,label=fr'$k_{n}={d[n]["k_inv_A"]:.4f}\ \mathrm{{\AA}}^{{-1}}$')
    style(ax);ax.set_ylabel(r'$F_{\mathrm{s}}(k,t)$');ax.legend(frameon=False);fig.text(.5,.975,r'$(8,8)$, 10L, 100 fs / 10 ns; 1 ps analysis grid; 3 replicas',ha='center',va='top')
    fig.savefig(out/'ISF_88_L10_k1-k2_self_aggregate.png',dpi=600);fig.savefig(out/'ISF_88_L10_k1-k2_self_aggregate.pdf');plt.close(fig)
    for n in (1,2):
      fig=plt.figure(figsize=(5.5,2.55))
      for j,(key,label) in enumerate(COMP):
       ax=fig.add_axes([.10+j*.30,.20,.24,.66]);t=d[n]['time_ps'];m=d[n][f'F_{key}_mean'];s=d[n][f'F_{key}_sem'];ax.fill_between(t,m-s,m+s,color=COL[n],alpha=.17,lw=0);ax.plot(t,m,color=COL[n],lw=1.1);style(ax);ax.set_ylabel(label);fig.text(.10+j*.30-.018,.88,f"({'abc'[j]})",fontweight='bold',fontsize=9)
      fig.text(.5,.975,fr'$(8,8)$, 10L: $k_{n}={d[n]["k_inv_A"]:.6f}\ \mathrm{{\AA}}^{{-1}}$; 3 replicas',ha='center',va='top');fig.savefig(out/f'ISF_88_L10_k{n}_components.png',dpi=600);fig.savefig(out/f'ISF_88_L10_k{n}_components.pdf');plt.close(fig)
if __name__=='__main__':main()
