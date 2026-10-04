"""
simulate_nn_plant.py
--------------------
Simuleaza planta NELINIARA (pendul inversat pe cart) cu controllerul NN.
Integrator: RK4 pur Python/NumPy
Plot: 4 stari + comanda u, pe aceeasi figura.

Rulare:
    python simulate_nn_plant.py --header nn_weights.h
    python simulate_nn_plant.py --header nn_weights.h --x0 0.0 0.0 0.15 0.0
    python simulate_nn_plant.py --header nn_weights.h --x0 0.05 0.1 0.1 0.3 --T 5.0
"""

import re
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('TkAgg')         
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ══════════════════════════════════════════════════════════════════
# 1. PARAMETRII PLANTEI
# ══════════════════════════════════════════════════════════════════

PLANT = dict(
    mcart = 1.0,          # [kg]  – masa carului  (nu e in snippet-ul tau, valoare tipica)
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


# ══════════════════════════════════════════════════════════════════
# 2. ODE NELINIARA
# ══════════════════════════════════════════════════════════════════

def nonlinear_ode(x, u, p):
    """
    Starea: x = [pos, dpos, theta, dtheta]
    Iesire: dx/dt
    """
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

def simulate(x0, nn, T=4.0, dt=Ts, plant=PLANT):
    N   = int(T / dt)
    t   = np.linspace(0, T, N + 1)
    X   = np.zeros((N + 1, 4))
    U   = np.zeros(N)

    X[0] = x0
    for k in range(N):
        u      = nn_evaluate(X[k], nn)
        U[k]   = u
        X[k+1] = rk4_step(X[k], u, dt, plant)

    return t, X, U


# ══════════════════════════════════════════════════════════════════
# 7. PLOT
# ══════════════════════════════════════════════════════════════════

STATE_LABELS = [
    (r'Pozitie cart  $p$',          '[m]'),
    (r'Viteza cart   $\dot{p}$',    '[m/s]'),
    (r'Unghi pendul  $\theta$',     '[rad]'),
    (r'Viteza unghi  $\dot\theta$', '[rad/s]'),
]

def plot_results(t, X, U, x0, header_name):
    fig = plt.figure(figsize=(12, 9))
    fig.suptitle(
        f'Simulare NN controller — planta neliniara\n'
        f'$x_0 = [{", ".join(f"{v:.3f}" for v in x0)}]$   '
        f'({header_name})',
        fontsize=12, fontweight='bold'
    )

    gs = gridspec.GridSpec(5, 1, hspace=0.55, figure=fig)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    # ── State plots ──────────────────────────────────────────────
    for i in range(4):
        ax = fig.add_subplot(gs[i])
        ax.plot(t, X[:, i], color=colors[i], lw=1.5)
        ax.axhline(0, color='k', lw=0.6, ls='--', alpha=0.5)
        ax.set_ylabel(STATE_LABELS[i][1], fontsize=9)
        ax.set_title(STATE_LABELS[i][0], fontsize=9, loc='left', pad=2)
        ax.set_xlim(t[0], t[-1])
        ax.tick_params(labelsize=8)
        ax.grid(True, alpha=0.3)
        if i < 3:
            ax.set_xticklabels([])

    # ── Control plot ─────────────────────────────────────────────
    t_u = t[:-1]   # u definit pe intervalele [k, k+1)
    ax_u = fig.add_subplot(gs[4])
    ax_u.step(t_u, U, where='post', color='#9467bd', lw=1.5)
    ax_u.axhline( U_MAX, color='gray', lw=0.8, ls=':', alpha=0.7, label=f'+{U_MAX} N')
    ax_u.axhline(-U_MAX, color='gray', lw=0.8, ls=':', alpha=0.7, label=f'-{U_MAX} N')
    ax_u.axhline(0,      color='k',    lw=0.6, ls='--', alpha=0.5)
    ax_u.fill_between(t_u, U, 0,
                      where=(np.abs(U) > 0),
                      alpha=0.15, color='#9467bd')
    ax_u.set_ylabel('[N]', fontsize=9)
    ax_u.set_title(r'Comanda  $u$', fontsize=9, loc='left', pad=2)
    ax_u.set_xlabel('Timp [s]', fontsize=9)
    ax_u.set_xlim(t[0], t[-1])
    ax_u.set_ylim(U_MIN*1.15, U_MAX*1.15)
    ax_u.tick_params(labelsize=8)
    ax_u.legend(fontsize=8, loc='upper right', framealpha=0.6)
    ax_u.grid(True, alpha=0.3)

    plt.savefig('simulare_nn_plant.png', dpi=150, bbox_inches='tight')
    print("\n  Plot salvat: simulare_nn_plant.png")
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
                        help='Stare initiala: pos dpos theta dtheta  (default: 0 0 0.15 0)')
    parser.add_argument('--T',
                        type=float, default=4.0,
                        help='Durata simulare [s]  (default: 4.0)')
    parser.add_argument('--mcart',
                        type=float, default=1.0,
                        help='Masa carului [kg]  (default: 1.0)')
    args = parser.parse_args()

    # Actualizeaza masa carul daca e specificata
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

    print("\n  Simulare in curs ...", end='', flush=True)
    t, X, U = simulate(x0, nn, T=args.T)
    print(" gata.")

    # Statistici rapide
    theta_max = np.max(np.abs(X[:, 2]))
    t_stab    = t[np.argmax(np.abs(X[:, 2]) < 0.01)] if np.any(np.abs(X[:, 2]) < 0.01) else None
    print(f"\n  |theta|_max   : {np.degrees(theta_max):.2f} deg")
    if t_stab is not None:
        print(f"  t_stab (~1deg): {t_stab:.3f} s")
    else:
        print(f"  t_stab        : pendul nu a ajuns sub 1 deg in {args.T}s")
    print(f"  |u|_max       : {np.max(np.abs(U)):.2f} N")

    plot_results(t, X, U, x0, args.header)