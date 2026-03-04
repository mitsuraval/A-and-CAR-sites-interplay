#!/usr/bin/env bash
#SBATCH --job-name=make_correlation
#SBATCH --output=out_make_correlation_%j.log
#SBATCH --error=err_make_correlation_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=exx96
#SBATCH --mem=7168M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mraval@wesleyan.edu


###
# Collects binary stacking (distance cutoff) and H-bond (presence/absence) events from all NEUTRAL_* directories nested under the running directory
# pools post-equilibration frames (>2000)
# computes φ-correlation matrices with CSV and heatmap outputs in BASE_DIR
# Expects input *_hbond.dat and *_stack.dat files with two columns (Frame, Value) located within NEUTRAL_*
###


# ─── 1. Set up Python environment ─────────────────────────────────────────
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaenv
module purge
module load python/3.12.0
echo "Using Python: $(which python3) $(python3 --version)"

# ─── 2. Define directories & expected files ───────────────────────────────
BASE_DIR="/home66/mraval/tRNAmod/unmod"
SKIP_LOG="$BASE_DIR/skipped_replicates.log"
> "$SKIP_LOG"   # clear old log

STACK_KEYS=(a b c d e f g h i j)
HBOND_KEYS=(k l m n o p q r s t u v w x y z)
EXPECTED=()
for c in "${STACK_KEYS[@]}"; do EXPECTED+=("${c}_stack.dat"); done
for c in "${HBOND_KEYS[@]}"; do EXPECTED+=("${c}_hbond.dat"); done

VALID_DIRS=()

# ─── 3. Find & validate each NEUTRAL_* directory ──────────────────────────
while IFS= read -r subdir; do
  echo "🔍 Checking: $subdir"
  missing=()
  for f in "${EXPECTED[@]}"; do
    [[ -f "$subdir/$f" ]] || missing+=("$f")
  done

  if (( ${#missing[@]} )); then
    echo "❌ Missing in $subdir:"
    for m in "${missing[@]}"; do echo "   - $m"; done
    echo "$subdir" >> "$SKIP_LOG"
    echo "⚠️ Skipping $subdir."
  else
    echo "✅ Valid: $subdir"
    VALID_DIRS+=("$subdir")
  fi
  echo "--------------------------------------------"
done < <(find "$BASE_DIR" -type d -name "NEUTRAL_*")

# ─── 4. Define stacking cutoffs ───────────────────────────────────────────
CUTOFFS=(4.5)

# ─── 5. Create Python script ───────────────────────────────────────────────
PY="$BASE_DIR/make_phi_correlations.py"
cat > "$PY" << 'PYCODE'
#!/usr/bin/env python3
"""
Compute φ‐coefficient matrices (full & cross‐site) for multiple stack cutoffs.
"""
import sys, os, glob
import pandas as pd, numpy as np, matplotlib.pyplot as plt

# read cutoff and directories
cutoff = float(sys.argv[1])
dirs = sys.argv[2:]
# base output folder
base = "/home66/mraval/tRNAmod/correlations"

# 1) Load & tag
records = []
for d in dirs:
    for fp in glob.glob(os.path.join(d, '*_stack.dat')) + glob.glob(os.path.join(d, '*_hbond.dat')):
        name = os.path.basename(fp)
        bond = name.split('_',1)[0]
        dtype = 'stack' if '_stack' in name else 'hbond'
        df = pd.read_csv(fp, sep=r'\s+', comment='#', header=None, names=['Frame','Value'])
        df = df[df.Frame > 2000]
        if dtype == 'stack':
            df.Value = (df.Value <= cutoff).astype(int)
        else:
            df.Value = (df.Value > 0).astype(int)
        df['bond'] = bond
        df['type'] = dtype
        records.append(df)

# 2) Pivot wide
data = pd.concat(records, ignore_index=True)
wide = data.pivot_table(index='Frame', columns=['bond','type'], values='Value')
wide.columns = [f"{b}_{t}" for b,t in wide.columns]

# 3) Full φ‐matrix
phi = wide.corr(method='pearson')
out_full_csv = os.path.join(base, f'all_phi_corr_cut{cutoff}.csv')
out_full_png = os.path.join(base, f'phi_heatmap_cut{cutoff}.png')
phi.to_csv(out_full_csv)
print(f"Saved full φ‐matrix: {out_full_csv}")
fig, ax = plt.subplots(figsize=(12,10))
im = ax.imshow(phi.values, cmap='bwr', vmin=-1, vmax=1, aspect='equal')
labels = [c.split('_')[0] for c in phi.columns]
ax.set_xticks(np.arange(len(labels))); ax.set_yticks(np.arange(len(labels)))
ax.set_xticklabels(labels, rotation=90, fontsize=6)
ax.set_yticklabels(labels, fontsize=6)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('φ‐coefficient', rotation=270, labelpad=15)
ax.set_title(f'Full φ‐heatmap (cutoff {cutoff} Å)', pad=10)
plt.tight_layout()
fig.savefig(out_full_png, dpi=300)
print(f"Saved full heatmap: {out_full_png}")

# 4) Cross‐site submatrix
A_site = ["a","b","f","g","k","l","m","n","o","p","q"]
CAR_site = ["d","e","i","j","t","u","v","w","x","y","z"]
cl = phi.stack().reset_index(name='rho')
cl.columns = ['full1','full2','rho']
cl['bond1'] = cl.full1.str.split('_').str[0]
cl['bond2'] = cl.full2.str.split('_').str[0]
mask = cl.bond1.isin(A_site) & cl.bond2.isin(CAR_site)
cs = cl.loc[mask]
cs_mat = cs.pivot(index='bond1', columns='bond2', values='rho')
cs_mat = cs_mat.reindex(index=A_site, columns=CAR_site)
out_cross_csv = os.path.join(base, f'cross_phi_corr_cut{cutoff}.csv')
out_cross_png = os.path.join(base, f'cross_phi_heatmap_cut{cutoff}.png')
cs_mat.to_csv(out_cross_csv)
print(f"Saved cross‐site CSV: {out_cross_csv}")
fig2, ax2 = plt.subplots(figsize=(8,6))
im2 = ax2.imshow(cs_mat.values, cmap='bwr', vmin=-1, vmax=1, aspect='equal')
ax2.set_xticks(np.arange(len(CAR_site))); ax2.set_yticks(np.arange(len(A_site)))
ax2.set_xticklabels(CAR_site, rotation=45, ha='right')
ax2.set_yticklabels(A_site)
for i,a in enumerate(A_site):
    for j,c in enumerate(CAR_site):
        v = cs_mat.loc[a,c]
        if pd.notna(v): ax2.text(j,i,f"{v:.2f}",ha='center',va='center',fontsize=8)
cbar2 = fig2.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
cbar2.set_label('φ‐coefficient', rotation=270, labelpad=15)
ax2.set_title(f'Cross φ‐heatmap (cutoff {cutoff} Å)', pad=10)
plt.tight_layout()
fig2.savefig(out_cross_png, dpi=300)
print(f"Saved cross‐heatmap: {out_cross_png}")
PYCODE

chmod +x "$PY"

# ─── 6. Run for each cutoff ────────────────────────────────────────────────
echo "📊 Running φ‐correlation for cutoffs: ${CUTOFFS[*]}"
for cut in "${CUTOFFS[@]}"; do
  python3 "$PY" "$cut" "${VALID_DIRS[@]}"
done

# ─── 7. Summary ───────────────────────────────────────────────────────────
if [[ -s "$SKIP_LOG" ]]; then
  echo "⚠️ Some dirs skipped (see $SKIP_LOG)"
else
  echo "✅ All processed with φ‐coefficient for all cutoffs!"
fi
