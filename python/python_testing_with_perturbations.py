"""
simulate_nn_plant.py
--------------------
Simuleaza planta NELINIARA (pendul inversat pe cart) cu controllerul NN.
Integrator: RK4 pur Python/NumPy (nu necesita CasADi).
Plot: 4 stari + comanda u + perturbatii w_k, pe aceeasi figura.
"""

import re
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('TkAgg')          # schimba cu 'Qt5Agg' sau 'Agg' daca nu merge
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ══════════════════════════════════════════════════════════════════
# 1. PARAMETRII PLANTEI
# ══════════════════════════════════════════════════════════════════

PLANT = dict(
    mcart = 1.0,          # [kg]  – masa carului
    m     = 0.135,        # [kg]  – masa pendulului
    b     = 0.0,          # [N·s/m] – frecare cart
    q     = 0.0,          # [N·s/rad] – frecare pendul
    g     = 9.8,          # [m/s²]
    l_tot = 0.47,         # [m]   – lungime bara
)
PLANT['l'] = PLANT['l_tot'] / 2.0                       # centru de masa
PLANT['I'] = PLANT['m'] * PLANT['l_tot']**2 / 3.0       # inertie bara

# Constrangeri comanda
U_MIN, U_MAX = -14.0, 14.0
# Pas de timp
Ts = 0.02   # [s]

# Perturbatie aditiva pe stare (identica cu w_bound din MATLAB)
W_BOUND = np.array([0.001, 0.001, 0.002, 0.002])   # [m, m/s, rad, rad/s]


# ══════════════════════════════════════════════════════════════════
# 2. ODE NELINIARA
# ══════════════════════════════════════════════════════════════════

def nonlinear_ode(x, u, p):
    mcart = p['mcart']; m = p['m']; b = p['b']
    q = p['q'];         g = p['g']; L = p['l']
    I = p['I']

    pos, dpos, theta, dtheta = x

    den = (mcart + m) * (m*L**2 + I) - (m*L*np.cos(theta))**2

    Fc = u - b*dpos + m*L*dtheta**2 * np.sin(theta)
    Fp = m*g*L*np.sin(theta) - q*dtheta

    ddpos   = ((m*L**2 + I)*Fc + m*L*np.cos(theta)*Fp) / den
    ddtheta = ((mcart + m)*Fp  - m*L*np.cos(theta)*Fc) / den

    return np.array([dpos, ddpos, dtheta, ddtheta])


# ══════════════════════════════════════════════════════════════════
# 3. INTEGRATOR RK4
# ══════════════════════════════════════════════════════════════════

def rk4_step(x, u, dt, p):
    k1 = nonlinear_ode(x,              u, p)
    k2 = nonlinear_ode(x + dt/2 * k1, u, p)
    k3 = nonlinear_ode(x + dt/2 * k2, u, p)
    k4 = nonlinear_ode(x + dt    * k3, u, p)
    return x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)


# ══════════════════════════════════════════════════════════════════
# 4. PARSER nn_weights.h
# ══════════════════════════════════════════════════════════════════

def parse_header(path):
    with open(path, 'r') as f:
        content = f.read()

    def get_define(name):
        m = re.search(rf'#define\s+{name}\s+(\d+)', content)
        return int(m.group(1)) if m else None

    def get_scalar(name):
        m = re.search(rf'const double\s+{name}\s*=\s*([+-]?\d+\.\d+(?:e[+-]?\d+)?)', content)
        return float(m.group(1)) if m else None

    def get_array_1d(name):
        m = re.search(rf'const double\s+{name}\[.*?\]\s*=\s*\{{([^}}]+)\}}', content)
        if not m:
            return None
        vals = re.findall(r'[+-]?\d+\.\d+(?:e[+-]?\d+)?', m.group(1))
        return np.array([float(v) for v in vals])

    def get_array_1d_int(name):
        m = re.search(rf'const int\s+{name}\[.*?\]\s*=\s*\{{([^}}]+)\}}', content)
        if not m:
            return None
        vals = re.findall(r'[+-]?\d+', m.group(1))
        return np.array([int(v) for v in vals])

    def get_array_2d(name, rows, cols):
        pattern = rf'const double\s+{name}\[{rows}\]\[{cols}\]\s*=\s*\{{(.*?)\}}\s*;'
        m = re.search(pattern, content, re.DOTALL)
        if not m:
            return None
        vals = re.findall(r'[+-]?\d+\.\d+(?:e[+-]?\d+)?', m.group(1))
        arr  = np.array([float(v) for v in vals])
        if len(arr) != rows * cols:
            print(f"  ATENTIE: {name} — gasit {len(arr)} valori, asteptat {rows*cols}")
        return arr.reshape(rows, cols)

    n_fc     = get_define('NN_N_FC')
    nx       = get_define('NN_NX')
    nu       = get_define('NN_NU')
    rows_arr = get_array_1d_int('NN_ROWS')
    cols_arr = get_array_1d_int('NN_COLS')

    ROWS = [int(r) for r in rows_arr]
    COLS = [int(c) for c in cols_arr]

    x_mean = get_array_1d('NN_X_MEAN')
    x_std  = get_array_1d('NN_X_STD')
    u_mean = get_scalar('NN_U_MEAN')
    u_std  = get_scalar('NN_U_STD')

    W_list, b_list = [], []
    for li in range(1, n_fc + 1):
        r = ROWS[li - 1];  c = COLS[li - 1]
        W = get_array_2d(f'NN_W{li}', r, c)
        b = get_array_1d(f'NN_B{li}')
        if W is None or b is None:
            print(f"  EROARE: NN_W{li}/NN_B{li} nu a fost gasit!")
            sys.exit(1)
        W_list.append(W)
        b_list.append(b.flatten())

    return dict(n_fc=n_fc, nx=nx, nu=nu, ROWS=ROWS, COLS=COLS,
                x_mean=x_mean, x_std=x_std, u_mean=u_mean, u_std=u_std,
                W=W_list, b=b_list)


# ══════════════════════════════════════════════════════════════════
# 5. FORWARD PASS NN
# ══════════════════════════════════════════════════════════════════

def nn_evaluate(x, nn):
    h = (x - nn['x_mean']) / nn['x_std']
    for li in range(nn['n_fc'] - 1):
        h = np.maximum(0.0, nn['W'][li] @ h + nn['b'][li])
    u_n = nn['W'][-1] @ h + nn['b'][-1]
    u   = nn['u_mean'] + nn['u_std'] * float(u_n.flatten()[0])
    return float(np.clip(u, U_MIN, U_MAX))


# ══════════════════════════════════════════════════════════════════
# 6. SIMULARE IN BUCLA INCHISA
# ══════════════════════════════════════════════════════════════════

def simulate(x0, nn, T=4.0, dt=Ts, plant=PLANT,
             with_disturbance=False, w_bound=W_BOUND, seed=42):
    rng = np.random.default_rng(seed)
    N   = int(T / dt)
    t   = np.linspace(0, T, N + 1)
    X   = np.zeros((N + 1, 4))
    U   = np.zeros(N)
    W   = np.zeros((N, 4))

    X[0] = x0
    for k in range(N):
        u      = nn_evaluate(X[k], nn)
        U[k]   = u
        x_next = rk4_step(X[k], u, dt, plant)

        if with_disturbance:
            wk   = rng.uniform(-w_bound, w_bound)
            W[k] = wk
            x_next = x_next + wk

        X[k+1] = x_next

    return t, X, U, W


# ══════════════════════════════════════════════════════════════════
# 7. PLOT
# ══════════════════════════════════════════════════════════════════

STATE_LABELS = [
    (r'Pozitie cart  $p$',          '[m]'),
    (r'Viteza cart   $\dot{p}$',    '[m/s]'),
    (r'Unghi pendul  $\theta$',     '[rad]'),
    (r'Viteza unghi  $\dot\theta$', '[rad/s]'),
]

W_COLORS  = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
W_LABELS  = [r'$w_1$ pos', r'$w_2$ dpos', r'$w_3$ $\theta$', r'$w_4$ $\dot\theta$']

def plot_results(t, X, U, W, x0, header_name, with_disturbance=False):
    n_rows = 6 if with_disturbance else 5
    fig = plt.figure(figsize=(12, 10 if with_disturbance else 9))

    dist_tag = '+ perturbatie $w_k$' if with_disturbance else 'fara perturbatie'
    fig.suptitle(
        f'Simulare NN controller ({dist_tag})\n'
    )

    gs = gridspec.GridSpec(n_rows, 1, hspace=0.60, figure=fig)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    for i in range(4):
        ax = fig.add_subplot(gs[i])
        ax.plot(t, X[:, i], color=colors[i], lw=1.5)
        ax.axhline(0, color='k', lw=0.6, ls='--', alpha=0.5)
        ax.set_ylabel(STATE_LABELS[i][1], fontsize=9)
        ax.set_title(STATE_LABELS[i][0], fontsize=9, loc='left', pad=2)
        ax.set_xlim(t[0], t[-1])
        ax.tick_params(labelsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_xticklabels([])

    t_u  = t[:-1]
    ax_u = fig.add_subplot(gs[4])
    ax_u.step(t_u, U, where='post', color='#9467bd', lw=1.5)
    ax_u.axhline( U_MAX, color='gray', lw=0.8, ls=':', alpha=0.7, label=f'+{U_MAX} N')
    ax_u.axhline(-U_MAX, color='gray', lw=0.8, ls=':', alpha=0.7, label=f'-{U_MAX} N')
    ax_u.axhline(0,      color='k',    lw=0.6, ls='--', alpha=0.5)
    ax_u.fill_between(t_u, U, 0, where=(np.abs(U) > 0), alpha=0.15, color='#9467bd')
    ax_u.set_ylabel('[N]', fontsize=9)
    ax_u.set_title(r'Comanda  $u$', fontsize=9, loc='left', pad=2)
    ax_u.set_xlim(t[0], t[-1])
    ax_u.set_ylim(U_MIN*1.15, U_MAX*1.15)
    ax_u.tick_params(labelsize=8)
    ax_u.legend(fontsize=8, loc='upper right', framealpha=0.6)
    ax_u.grid(True, alpha=0.3)

    if with_disturbance:
        ax_u.set_xticklabels([])
        ax_w = fig.add_subplot(gs[5])
        for i in range(4):
            ax_w.plot(t_u, W[:, i], lw=0.9, alpha=0.8,
                      color=W_COLORS[i], label=W_LABELS[i])
        for i, wb in enumerate(W_BOUND):
            ax_w.axhline( wb, color=W_COLORS[i], lw=0.6, ls='--', alpha=0.4)
            ax_w.axhline(-wb, color=W_COLORS[i], lw=0.6, ls='--', alpha=0.4)
        ax_w.axhline(0, color='k', lw=0.5, ls='--', alpha=0.4)
        ax_w.set_ylabel('[wk]', fontsize=9)
        ax_w.set_xlabel('Timp [s]', fontsize=9)
        ax_w.set_xlim(t[0], t[-1])
        ax_w.tick_params(labelsize=8)
        ax_w.legend(fontsize=7, loc='upper right', ncol=4, framealpha=0.6)
        ax_w.grid(True, alpha=0.3)
    else:
        ax_u.set_xlabel('Timp [s]', fontsize=9)

    fname = 'simulare_nn_plant_dist.png' if with_disturbance else 'simulare_nn_plant.png'
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    print(f"\n  Plot salvat: {fname}")
    plt.show()


# ══════════════════════════════════════════════════════════════════
# 8. MAIN
# ══════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Simulare NN controller pe planta neliniara (pendul pe cart)'
    )
    parser.add_argument('--header',
                        default='nn_weights.h',
                        help='Cale catre nn_weights.h')
    parser.add_argument('--x0', nargs=4, type=float,
                        default=[0.5, 0.0, 0.15, 0.0],
                        metavar=('POS', 'DPOS', 'THETA', 'DTHETA'),
                        help='Stare initiala: pos dpos theta dtheta')
    parser.add_argument('--T',
                        type=float, default=4.0,
                        help='Durata simulare [s]')
    parser.add_argument('--mcart',
                        type=float, default=1.0,
                        help='Masa carului [kg]')
    args = parser.parse_args()

    PLANT['mcart'] = args.mcart
    x0 = np.array(args.x0)

    print(f"\n  Citire header : {args.header}")
    nn = parse_header(args.header)

    print(f"  Arhitectura   : {nn['n_fc']} straturi FC, "
          f"intrare {nn['nx']}, iesire {nn['nu']}")
    print(f"  Stare initiala: {x0}")
    print(f"  Durata        : {args.T} s  ({int(args.T/Ts)} pasi Ts={Ts}s)")
    print(f"  Masa cart     : {PLANT['mcart']} kg")
    print(f"  Masa pendul   : {PLANT['m']} kg,  l_tot = {PLANT['l_tot']} m")

    # ─── ACTIVARE ZGOMOT ───────────────────────────────────────────
    print("\n  Simulare in curs (CU perturbatii) ...", end='', flush=True)
    t, X, U, W = simulate(x0, nn, T=args.T, with_disturbance=True)
    print(" gata.")
    # ───────────────────────────────────────────────────────────────

    theta_max = np.max(np.abs(X[:, 2]))
    t_stab    = t[np.argmax(np.abs(X[:, 2]) < 0.01)] if np.any(np.abs(X[:, 2]) < 0.01) else None
    print(f"\n  |theta|_max   : {np.degrees(theta_max):.2f} deg")
    if t_stab is not None:
        print(f"  t_stab (~1deg): {t_stab:.3f} s")
    else:
        print(f"  t_stab        : pendul nu a ajuns sub 1 deg in {args.T}s")
    print(f"  |u|_max       : {np.max(np.abs(U)):.2f} N")

    # Generare plot cu zgomot inclus
    plot_results(t, X, U, W, x0, args.header, with_disturbance=True)