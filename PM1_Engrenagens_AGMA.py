import math
import matplotlib.pyplot as plt

# =============================================================================
# PROJETO DE MÁQUINAS 1 - DIMENSIONAMENTO DE ENGRENAGENS CILÍNDRICAS (AGMA)
# Dentes RETOS e HELICOIDAIS em um único código
# Aluno: João Vicente Rosal Giovannetti Daros
# Referência: Shigley - Projeto de Engenharia Mecânica (Caps. 13 e 14, SI)
# Unidades: mm, kW, rpm, MPa
#
# O programa pergunta apenas: tipo de dente, número de dentes do pinhão,
# material e largura de face. Os demais parâmetros são os PADRÕES abaixo
# (edite esta seção para mudar o problema).
# Saídas: resumo no terminal + tabela gráfica dos fatores (com a localização
# de cada um no Shigley) + mapa 2D "número de dentes x largura de face"
# mostrando as regiões PASSA (verde) / NÃO PASSA (vermelho).
# =============================================================================

# ---------------------- PARÂMETROS PADRÃO DO PROBLEMA -----------------------
PADROES = {
    "reto": dict(
        H_motor=0.339,     # Potência do motor (kW)
        n_motor=11.0,      # Rotação do motor (rpm)
        n_acionado=2.2,    # Rotação do eixo acionado (rpm)
        K_o=1.0,           # Fator de sobrecarga (Seção 14-8)
        FS=1.0,            # Fator de segurança mínimo do projeto
        Q_v=7,             # Número de qualidade (Seção 14-7)
        phi_n=20.0,        # Ângulo de pressão normal (graus)
        psi=0.0,           # Ângulo de hélice (graus) - 0 para dentes retos
        m_n=3.0,           # Módulo normal (mm)
        ciclos_p=1e9,      # Vida útil do pinhão (ciclos)
        R=0.95,            # Confiabilidade
        Yj_p=0.33,         # Fator geométrico J do pinhão (Fig. 14-6)
        Yj_g=0.42,         # Fator geométrico J da engrenagem (Fig. 14-6)
        N_p_padrao=20,
    ),
    "helicoidal": dict(
        H_motor=9.0,
        n_motor=1000.0,
        n_acionado=500.0,
        K_o=1.25,
        FS=2.0,
        Q_v=7,
        phi_n=20.0,
        psi=30.0,
        m_n=2.5,
        ciclos_p=1e8,
        R=0.95,
        Yj_p=0.41,         # (Figs. 14-7 e 14-8)
        Yj_g=0.444,
        N_p_padrao=14,
    ),
}

# Fatores fixos (hipóteses usuais do curso)
K_b = 1     # Fator de espessura de borda (Seção 14-16, Fig. 14-16), mB >= 1,2
K_t = 1     # Fator de temperatura (Seção 14-15), T_oleo < 120 C
C_f = 1     # Fator de condição de superfície (Seção 14-9)
C_h = 1     # Fator de razão de dureza (Seção 14-12, Eq. 14-36)
C_p = 191   # Coeficiente elástico aço-aço, raiz(MPa) (Eq. 14-13, Tabela 14-8)

# Materiais disponíveis: (nome, grau AGMA, dureza HB)
# St pela Fig. 14-2 e Sc pela Fig. 14-5 (aços endurecidos por completo, MPa)
MATERIAIS = {
    1: ("Aço grau 1 - HB 300", 1, 300.0),
    2: ("Aço grau 2 - HB 350", 2, 350.0),
    3: ("Aço grau 2 - HB 397", 2, 397.0),
}

# Cores de status do mapa passa/não passa
COR_PASSA = "#2E7D32"
COR_FALHA = "#C62828"
COR_LINHA = "#37474F"


def pede_int(msg, padrao):
    txt = input(f"{msg} [{padrao}]: ").strip()
    return int(txt) if txt else padrao


def pede_float(msg, padrao):
    txt = input(f"{msg} [{padrao}]: ").strip().replace(",", ".")
    return float(txt) if txt else padrao


# Tabela 14-2 (Shigley) - Fator de forma de Lewis Y (interpolação linear)
TABELA_LEWIS = [(12, 0.245), (13, 0.261), (14, 0.277), (15, 0.290), (16, 0.296),
                (17, 0.303), (18, 0.309), (19, 0.314), (20, 0.322), (21, 0.328),
                (22, 0.331), (24, 0.337), (26, 0.346), (28, 0.353), (30, 0.359),
                (34, 0.371), (38, 0.384), (43, 0.397), (50, 0.409), (60, 0.422),
                (75, 0.435), (100, 0.447), (150, 0.460), (300, 0.472), (400, 0.480)]


def lewis_Y(N):
    """Fator de Lewis Y interpolado da Tabela 14-2."""
    if N <= TABELA_LEWIS[0][0]:
        return TABELA_LEWIS[0][1]
    if N >= TABELA_LEWIS[-1][0]:
        return TABELA_LEWIS[-1][1]
    for (N1, Y1), (N2, Y2) in zip(TABELA_LEWIS, TABELA_LEWIS[1:]):
        if N1 <= N <= N2:
            return Y1 + (Y2 - Y1) * (N - N1) / (N2 - N1)


def fator_tamanho(b, m_t, Y):
    """Ks - Fator de tamanho (Seção 14-10)."""
    return 0.8433 * (b * m_t * math.sqrt(Y)) ** 0.0535


def fator_distribuicao(b, d_p):
    """Kh - Fator de distribuição de carga (Seção 14-11, Eqs. 14-30 a 14-35).
    Hipóteses: dentes não coroados, montagem entre mancais,
    unidade fechada comercial, condições normais de montagem."""
    C_mc = 1
    razao = b / (10 * d_p)
    if razao < 0.05:
        razao = 0.05
    if b <= 25:
        C_pf = razao - 0.025
    else:
        C_pf = razao - 0.0375 + 4.92e-4 * b
    C_pm = 1
    C_ma = 0.127 + 0.622e-3 * b - 1.69e-7 * b**2
    C_e = 1
    return 1 + C_mc * (C_pf * C_pm + C_ma * C_e)


# =============================================================================
# ENTRADAS DO USUÁRIO (apenas 4 perguntas)
# =============================================================================
print("=" * 70)
print(" DIMENSIONAMENTO DE ENGRENAGENS - MÉTODO AGMA (Shigley)")
print("=" * 70)
tipo = ""
while tipo not in ("1", "2"):
    tipo = input("Tipo de engrenagem (1 = dentes retos, 2 = helicoidal): ").strip()
helicoidal = (tipo == "2")
cfg = PADROES["helicoidal" if helicoidal else "reto"]

H_motor = cfg["H_motor"]
n_motor = cfg["n_motor"]
n_acionado = cfg["n_acionado"]
K_o = cfg["K_o"]
FS = cfg["FS"]
Q_v = cfg["Q_v"]
phi_n = cfg["phi_n"]
psi = cfg["psi"]
m_n = cfg["m_n"]
ciclos_p = cfg["ciclos_p"]
R = cfg["R"]
Yj_p = cfg["Yj_p"]
Yj_g = cfg["Yj_g"]

N_p = pede_int("Número de dentes do pinhão", cfg["N_p_padrao"])

print("Materiais disponíveis:")
for i, (nome, _, _) in MATERIAIS.items():
    print(f"  {i} = {nome}")
print("  4 = Personalizado (digitar St e Sc)")
escolha_mat = pede_int("Material", 2)
if escolha_mat == 4:
    nome_mat = "Personalizado"
    S_t = pede_float("Tensão de flexão admissível St (MPa)", 359.0)
    Sc = pede_float("Resistência ao contato admissível Sc (MPa)", 1080.5)
else:
    nome_mat, grau, HB = MATERIAIS[escolha_mat]
    if grau == 2:
        S_t = 0.703 * HB + 113        # Fig. 14-2, grau 2
        Sc = 2.41 * HB + 237          # Fig. 14-5, grau 2
    else:
        S_t = 0.533 * HB + 88.3       # Fig. 14-2, grau 1
        Sc = 2.22 * HB + 200          # Fig. 14-5, grau 1

# =============================================================================
# FATORES INDEPENDENTES DA GEOMETRIA (calculados uma única vez)
# =============================================================================
reducao = n_motor / n_acionado
ciclos_g = ciclos_p / reducao

K_r = 0.658 - 0.0759 * math.log(1 - R)          # Eq. 14-38 (Seção 14-14)
Yn_p = 1.3558 * ciclos_p**-0.0178               # Fig. 14-14 (Seção 14-13)
Yn_g = 1.3558 * ciclos_g**-0.0178
Zn_p = 1.4488 * ciclos_p**-0.023                # Fig. 14-15 (Seção 14-13)
Zn_g = 1.4488 * ciclos_g**-0.023


def prepara(N):
    """Geometria, carregamento e fatores que dependem do nº de dentes N."""
    g = {"N_p": N}
    if helicoidal:
        g["phi_t"] = math.degrees(math.atan(math.tan(math.radians(phi_n)) / math.cos(math.radians(psi))))
        g["m_t"] = m_n / math.cos(math.radians(psi))
        g["p_x"] = math.pi * m_n / math.sin(math.radians(psi))
    else:
        g["phi_t"] = phi_n
        g["m_t"] = m_n
        g["p_x"] = None
    g["N_g"] = max(round(N * reducao), N + 1) if reducao > 1 else round(N * reducao)
    g["m_g"] = g["N_g"] / N
    g["d_p"] = g["m_t"] * N
    g["d_g"] = g["m_t"] * g["N_g"]
    g["Y_p"] = lewis_Y(N)
    g["Y_g"] = lewis_Y(g["N_g"])

    # Nº mínimo de dentes do pinhão contra interferência (Eq. 13-11 / Eq. 13-22)
    k = 1
    sin2 = math.sin(math.radians(g["phi_t"])) ** 2
    fh = math.cos(math.radians(psi)) if helicoidal else 1.0
    mg = g["m_g"]
    g["N_p_min"] = 2 * k * fh * (mg + math.sqrt(mg**2 + (1 + 2 * mg) * sin2)) / ((1 + 2 * mg) * sin2)

    # Carregamento e fator dinâmico (Seção 14-7)
    g["V"] = math.pi * g["d_p"] * n_motor / 60000
    g["W_t"] = 1000 * H_motor / g["V"]
    B = 0.25 * (12 - Q_v) ** (2 / 3)
    A = 50 + 56 * (1 - B)
    g["K_v"] = ((A + math.sqrt(200 * g["V"])) / A) ** B      # Eq. 14-27
    g["V_max"] = (A + (Q_v - 3)) ** 2 / 200                  # Eq. 14-29

    # Fator geométrico de contato I (Seção 14-5, Eqs. 14-21 a 14-25)
    if helicoidal:
        rb_p = (g["d_p"] / 2) * math.cos(math.radians(g["phi_t"]))
        rb_g = (g["d_g"] / 2) * math.cos(math.radians(g["phi_t"]))
        a = m_n
        Z = (math.sqrt((g["d_p"] / 2 + a)**2 - rb_p**2) + math.sqrt((g["d_g"] / 2 + a)**2 - rb_g**2)
             - (g["d_p"] / 2 + g["d_g"] / 2) * math.sin(math.radians(g["phi_t"])))
        P_N = math.pi * m_n * math.cos(math.radians(phi_n))
        m_N = P_N / (0.95 * Z)
    else:
        m_N = 1
    g["m_N"] = m_N
    g["Z_I"] = (math.cos(math.radians(g["phi_t"])) * math.sin(math.radians(g["phi_t"])) / (2 * m_N)) * mg / (mg + 1)
    return g


def b_requerido(g):
    """Largura de face requerida por flexão e desgaste (iterativa: Ks e Kh dependem de b)."""
    b = 30.0
    for _ in range(100):
        Ks_p = fator_tamanho(b, g["m_t"], g["Y_p"])
        K_h = fator_distribuicao(b, g["d_p"])
        b_flex = FS * g["W_t"] * K_o * g["K_v"] * Ks_p * (1 / g["m_t"]) * (K_h * K_b / Yj_p) * (K_t * K_r) / (S_t * Yn_p)
        b_desg = ((C_p * K_t * K_r / (Sc * Zn_p * C_h)) ** 2) * FS * g["W_t"] * K_o * g["K_v"] * Ks_p * (K_h / g["d_p"]) * (C_f / g["Z_I"])
        b_novo = max(b_flex, b_desg)
        if abs(b_novo - b) < 1e-4:
            break
        b = b_novo
    return b_flex, b_desg, max(b_flex, b_desg)


def tensoes(g, b):
    """Tensões e fatores de segurança AGMA com a largura de face b."""
    r = {}
    r["Ks_p"] = fator_tamanho(b, g["m_t"], g["Y_p"])
    r["Ks_g"] = fator_tamanho(b, g["m_t"], g["Y_g"])
    r["K_h"] = fator_distribuicao(b, g["d_p"])
    comum = g["W_t"] * K_o * g["K_v"]
    r["Tflex_p"] = comum * r["Ks_p"] * (1 / (b * g["m_t"])) * (r["K_h"] * K_b / Yj_p)     # Eq. 14-15
    r["Tflex_g"] = comum * r["Ks_g"] * (1 / (b * g["m_t"])) * (r["K_h"] * K_b / Yj_g)
    r["Tcont_p"] = C_p * math.sqrt(comum * r["Ks_p"] * (r["K_h"] / (g["d_p"] * b)) * (C_f / g["Z_I"]))  # Eq. 14-16
    r["Tcont_g"] = C_p * math.sqrt(comum * r["Ks_g"] * (r["K_h"] / (g["d_p"] * b)) * (C_f / g["Z_I"]))
    r["SF_p"] = (S_t * Yn_p) / (K_t * K_r * r["Tflex_p"])       # Eq. 14-41
    r["SF_g"] = (S_t * Yn_g) / (K_t * K_r * r["Tflex_g"])
    r["SH_p"] = (Sc * Zn_p * C_h) / (K_t * K_r * r["Tcont_p"])  # Eq. 14-42
    r["SH_g"] = (Sc * Zn_g * C_h) / (K_t * K_r * r["Tcont_g"])
    return r


# =============================================================================
# CÁLCULO DO CASO DO USUÁRIO
# =============================================================================
g = prepara(N_p)
b_flex, b_desg, b_req = b_requerido(g)

if helicoidal:
    b_min_rec = 2 * g["p_x"]
    b_max_rec = None
    intervalo_txt = f"b >= {b_min_rec:.1f} mm (2 passos axiais)"
else:
    b_min_rec = 3 * math.pi * g["m_t"]
    b_max_rec = 5 * math.pi * g["m_t"]
    intervalo_txt = f"{b_min_rec:.1f} a {b_max_rec:.1f} mm (3*pi*m a 5*pi*m)"

print(f"\nLargura de face requerida: {b_req:.2f} mm "
      f"(flexão: {b_flex:.2f} mm | desgaste: {b_desg:.2f} mm)")
print(f"Intervalo recomendado: {intervalo_txt}")
b_sugerido = float(math.ceil(max(b_req, b_min_rec)))
b = pede_float("Largura de face adotada (mm)", b_sugerido)

r = tensoes(g, b)
passa = (r["SF_p"] >= FS and r["SF_g"] >= FS
         and r["SH_p"]**2 >= FS and r["SH_g"]**2 >= FS
         and N_p >= g["N_p_min"])

# ------------------------------ RESUMO NO TERMINAL --------------------------
print("\n" + "=" * 70)
print(f" RESULTADOS - DENTES {'HELICOIDAIS' if helicoidal else 'RETOS'} | "
      f"N_p = {N_p} | b = {b:.1f} mm | {nome_mat}")
print("=" * 70)
print(f"Redução: {reducao:.3f} | N_g = {g['N_g']} | d_p = {g['d_p']:.2f} mm | d_g = {g['d_g']:.2f} mm")
if helicoidal:
    print(f"Ângulo de pressão transversal: {g['phi_t']:.2f} graus | Passo axial: {g['p_x']:.2f} mm")
print(f"N mínimo de dentes (interferência): {math.ceil(g['N_p_min'])}"
      + ("  <-- VIOLADO!" if N_p < g["N_p_min"] else ""))
print(f"V = {g['V']:.3f} m/s (máx. p/ Qv={Q_v}: {g['V_max']:.1f} m/s) | Wt = {g['W_t']:.1f} N")
print(f"St = {S_t:.1f} MPa | Sc = {Sc:.1f} MPa")
print("-" * 70)
print(f"Tensão de flexão: pinhão {r['Tflex_p']:.1f} MPa | coroa {r['Tflex_g']:.1f} MPa")
print(f"Tensão de contato: pinhão {r['Tcont_p']:.1f} MPa | coroa {r['Tcont_g']:.1f} MPa")
print(f"SF (flexão): pinhão {r['SF_p']:.2f} | coroa {r['SF_g']:.2f}")
print(f"SH (contato): pinhão {r['SH_p']:.2f} (SH²={r['SH_p']**2:.2f}) | "
      f"coroa {r['SH_g']:.2f} (SH²={r['SH_g']**2:.2f})")
print("-" * 70)
print(f" VEREDITO: {'PASSA (SF e SH² >= FS = %.1f)' % FS if passa else 'NÃO PASSA (algum fator < FS = %.1f)' % FS}")
print("=" * 70)

# =============================================================================
# FIGURA 1 - TABELA DOS FATORES + ONDE ENCONTRAR NO SHIGLEY
# =============================================================================
linhas = [
    ("Ko",   f"{K_o:.2f}",           "Fator de sobrecarga",              "Seção 14-8"),
    ("Kv",   f"{g['K_v']:.4f}",      "Fator dinâmico",                   "Seção 14-7, Eqs. 14-27/28, Fig. 14-9"),
    ("Ks,p", f"{r['Ks_p']:.4f}",     "Fator de tamanho (pinhão)",        "Seção 14-10 (Y: Tabela 14-2)"),
    ("Ks,g", f"{r['Ks_g']:.4f}",     "Fator de tamanho (coroa)",         "Seção 14-10 (Y: Tabela 14-2)"),
    ("Kh",   f"{r['K_h']:.4f}",      "Distribuição de carga",            "Seção 14-11, Eqs. 14-30 a 14-35"),
    ("Kb",   f"{K_b:.2f}",           "Espessura de borda",               "Seção 14-16, Fig. 14-16"),
    ("Kt",   f"{K_t:.2f}",           "Temperatura",                      "Seção 14-15"),
    ("Kr",   f"{K_r:.4f}",           "Confiabilidade",                   "Seção 14-14, Eq. 14-38, Tab. 14-10"),
    ("Cf",   f"{C_f:.2f}",           "Condição de superfície",           "Seção 14-9"),
    ("Ch",   f"{C_h:.2f}",           "Razão de dureza",                  "Seção 14-12, Eq. 14-36"),
    ("Cp",   f"{C_p:.0f}",           "Coef. elástico (raiz de MPa)",     "Eq. 14-13, Tabela 14-8"),
    ("Yn,p", f"{Yn_p:.4f}",          "Ciclagem de flexão (pinhão)",      "Seção 14-13, Fig. 14-14"),
    ("Yn,g", f"{Yn_g:.4f}",          "Ciclagem de flexão (coroa)",       "Seção 14-13, Fig. 14-14"),
    ("Zn,p", f"{Zn_p:.4f}",          "Ciclagem de contato (pinhão)",     "Seção 14-13, Fig. 14-15"),
    ("Zn,g", f"{Zn_g:.4f}",          "Ciclagem de contato (coroa)",      "Seção 14-13, Fig. 14-15"),
    ("J,p",  f"{Yj_p:.3f}",          "Geométrico de flexão (pinhão)",    "Seção 14-5, " + ("Figs. 14-7/14-8" if helicoidal else "Fig. 14-6")),
    ("J,g",  f"{Yj_g:.3f}",          "Geométrico de flexão (coroa)",     "Seção 14-5, " + ("Figs. 14-7/14-8" if helicoidal else "Fig. 14-6")),
    ("I",    f"{g['Z_I']:.4f}",      "Geométrico de contato",            "Seção 14-5, Eqs. 14-21 a 14-25"),
    ("St",   f"{S_t:.1f} MPa",       "Flexão admissível",                "Fig. 14-2, Tabs. 14-3/14-4"),
    ("Sc",   f"{Sc:.1f} MPa",        "Contato admissível",               "Fig. 14-5, Tabela 14-6"),
    ("SF",   f"{r['SF_p']:.2f} / {r['SF_g']:.2f}", "Segurança de flexão (pinhão/coroa)",  "Seção 14-17, Eq. 14-41"),
    ("SH",   f"{r['SH_p']:.2f} / {r['SH_g']:.2f}", "Segurança de contato (pinhão/coroa)", "Seção 14-17, Eq. 14-42"),
]

fig1, ax1 = plt.subplots(figsize=(11, 0.30 * len(linhas) + 1.2))
fig1.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.03)
ax1.axis("off")
ax1.set_title(f"Fatores AGMA - dentes {'helicoidais' if helicoidal else 'retos'} "
              f"(N_p = {N_p}, b = {b:.1f} mm, {nome_mat})", fontsize=12, pad=14)
tabela = ax1.table(cellText=[list(l) for l in linhas],
                   colLabels=["Fator", "Valor", "O que é", "Onde está no Shigley"],
                   colWidths=[0.09, 0.15, 0.35, 0.41],
                   cellLoc="left", loc="center")
tabela.auto_set_font_size(False)
tabela.set_fontsize(9.5)
tabela.scale(1, 1.35)
for (lin, col), cel in tabela.get_celld().items():
    cel.set_edgecolor("#CFD8DC")
    if lin == 0:
        cel.set_facecolor(COR_LINHA)
        cel.set_text_props(color="white", weight="bold")
    elif lin % 2 == 0:
        cel.set_facecolor("#F5F5F5")
fig1.tight_layout()

# =============================================================================
# FIGURA 2 - MAPA "NÚMERO DE DENTES x LARGURA DE FACE" (PASSA / NÃO PASSA)
# =============================================================================
N_ini, N_fim = 12, 60
Ns = list(range(N_ini, N_fim + 1))
fronteira = [b_requerido(prepara(N))[2] for N in Ns]

y_max = max(max(fronteira) * 1.15, b * 1.25, (b_max_rec or 0) * 1.15)

fig2, ax2 = plt.subplots(figsize=(10, 6))
ax2.fill_between(Ns, fronteira, y_max, color=COR_PASSA, alpha=0.15)
ax2.fill_between(Ns, 0, fronteira, color=COR_FALHA, alpha=0.15)
ax2.plot(Ns, fronteira, color=COR_LINHA, lw=2.2,
         label="b mínimo requerido (flexão e desgaste, SF = SH² = FS)")

# Rótulos das regiões (não depender só da cor)
N_meio = Ns[len(Ns) // 2]
b_meio = fronteira[len(Ns) // 2]
ax2.text(N_meio, (b_meio + y_max) / 2, "PASSA", color=COR_PASSA,
         fontsize=15, weight="bold", ha="center")
ax2.text(N_meio, b_meio / 2, "NÃO PASSA", color=COR_FALHA,
         fontsize=15, weight="bold", ha="center")

# Intervalo recomendado de largura de face
ax2.axhline(b_min_rec, color=COR_LINHA, ls="--", lw=1.2, alpha=0.7,
            label=("b mínimo recomendado (2 passos axiais)" if helicoidal
                   else "Intervalo recomendado (3·pi·m a 5·pi·m)"))
if b_max_rec is not None:
    ax2.axhline(b_max_rec, color=COR_LINHA, ls="--", lw=1.2, alpha=0.7)

# Limite de interferência (nº mínimo de dentes)
N_min_interf = math.ceil(g["N_p_min"])
if N_min_interf > N_ini:
    ax2.axvline(N_min_interf, color=COR_FALHA, ls=":", lw=1.5, alpha=0.8,
                label=f"N mínimo contra interferência ({N_min_interf} dentes)")

# Ponto do projeto do usuário
cor_ponto = COR_PASSA if passa else COR_FALHA
ax2.scatter([N_p], [b], s=130, color=cor_ponto, edgecolor="black", zorder=5,
            label=f"Seu projeto (N={N_p}, b={b:.0f} mm) - {'PASSA' if passa else 'NÃO PASSA'}")

ax2.set_xlim(N_ini, N_fim)
ax2.set_ylim(0, y_max)
ax2.set_xlabel("Número de dentes do pinhão")
ax2.set_ylabel("Largura de face b (mm)")
ax2.set_title(f"Mapa de projeto - dentes {'helicoidais' if helicoidal else 'retos'} "
              f"(módulo {m_n} mm, {nome_mat}, FS = {FS:.1f})")
ax2.grid(alpha=0.3)
ax2.legend(loc="upper right", fontsize=9)
fig2.text(0.01, 0.01,
          "Obs.: fatores J fixados nos valores padrão; ao variar muito o nº de dentes, "
          "confira J nas Figs. 14-6/14-7/14-8.",
          fontsize=8, color="#607D8B")
fig2.tight_layout()

plt.show()
