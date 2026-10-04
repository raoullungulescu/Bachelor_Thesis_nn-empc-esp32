import matplotlib.pyplot as plt

# Inițializăm structurile de date
stari = {0: [], 1: [], 2: [], 3: []}
comanda = []

file_name = 'raoul_hai.txt'

# Citirea și parsarea fișierului complet
with open(file_name, 'r', encoding='utf-8', errors='ignore') as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) == 5:
            try:
                vals = [float(p) for p in parts]
                for i in range(4):
                    stari[i].append(vals[i])
                comanda.append(vals[4])
            except ValueError:
                continue

# Definim intervalul de eșantioane dorit
start_idx = 8100
end_idx = 9236

# Extragem doar porțiunea de date solicitată
stari_filtrate = {i: stari[i][start_idx:end_idx + 1] for i in range(4)}
comanda_filtrata = comanda[start_idx:end_idx + 1]

# Resetăm axa X: începe de la 0 și numără eșantioanele din interval
axa_x = list(range(len(comanda_filtrata)))


# =========================================================================
# FEREASTRA 1: GRAFICUL PENTRU STĂRI (4 Subploturi separate, toate albastre)
# =========================================================================
fig1, axs = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
etichete_stari = ['Poziție', 'Viteză', 'Unghi', 'Viteză unghiulară']

for i in range(4):
    # Plotăm fiecare stare cu linie albastră continuuă
    axs[i].plot(axa_x, stari_filtrate[i], color='blue', linewidth=1.5, label=etichete_stari[i])
    axs[i].set_ylabel(etichete_stari[i], fontweight='bold', fontsize=10)
    axs[i].grid(True, linestyle='--', alpha=0.5)
    axs[i].legend(loc='upper right')

# Punem eticheta de timp/eșantioane doar pe ultimul subplot de jos
axs[3].set_xlabel(f'Număr eșantioane', fontweight='bold', fontsize=11)
fig1.suptitle('Vizualizarea celor 4 stări', fontsize=14, fontweight='bold', y=0.99)
fig1.tight_layout()


# =========================================================================
# FEREASTRA 2: GRAFICUL PENTRU COMANDĂ (Alt grafic separat, tot albastru)
# =========================================================================
fig2, ax_com = plt.subplots(figsize=(11, 4))

# Plotăm comanda separat cu linie albastră
ax_com.plot(axa_x, comanda_filtrata, color='blue', linewidth=1.5, label='Comandă')
ax_com.set_title('Vizualizarea comenzii', fontweight='bold', fontsize=12)
ax_com.set_ylabel('Valoare comandă', fontweight='bold')
ax_com.set_xlabel(f'Număr eșantioane', fontweight='bold', fontsize=11)
ax_com.grid(True, linestyle='--', alpha=0.5)
ax_com.legend(loc='upper right')
fig2.tight_layout()


# Afișăm ambele ferestre în același timp
plt.show()