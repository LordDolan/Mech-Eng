import math

# =============================================================================
# PROJETO DE MÁQUINAS 1 - DIMENSIONAMENTO DE ENGRENAGENS CILÍNDRICAS (AGMA)
# Dentes RETOS e HELICOIDAIS em um único código
# Aluno: João Vicente Rosal Giovannetti Daros
# Referência: Shigley - Projeto de Engenharia Mecânica (Caps. 13 e 14, SI)
# Unidades: mm, kW, rpm, MPa
# =============================================================================
# Pressione Enter em qualquer pergunta para aceitar o valor padrão [entre colchetes].


def pede_float(msg, padrao):
    txt = input(f"{msg} [{padrao}]: ").strip().replace(",", ".")
    return float(txt) if txt else padrao


def pede_int(msg, padrao):
    txt = input(f"{msg} [{padrao}]: ").strip()
    return int(txt) if txt else padrao


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
    """Kh - Fator de distribuição de carga (Eqs. 14-30 a 14-35).
    Hipóteses: dentes não coroados, montagem entre mancais,
    unidade fechada comercial, condições normais de montagem."""
    C_mc = 1                                    # Eq. 14-31: dentes não coroados
    razao = b / (10 * d_p)
    if razao < 0.05:
        razao = 0.05
    if b <= 25:
        C_pf = razao - 0.025                    # Eq. 14-32 (b <= 25 mm)
    else:
        C_pf = razao - 0.0375 + 4.92e-4 * b     # Eq. 14-32 (25 < b <= 425 mm)
    C_pm = 1                                    # Eq. 14-33: montagem entre mancais
    C_ma = 0.127 + 0.622e-3 * b - 1.69e-7 * b**2  # Eq. 14-34: unidade fechada comercial
    C_e = 1                                     # Eq. 14-35: montagem normal
    return 1 + C_mc * (C_pf * C_pm + C_ma * C_e)


# =============================================================================
# SELEÇÃO DO TIPO DE ENGRENAGEM
# =============================================================================
print("=" * 70)
print(" DIMENSIONAMENTO DE ENGRENAGENS - MÉTODO AGMA (Shigley)")
print("=" * 70)
tipo = ""
while tipo not in ("1", "2"):
    tipo = input("Tipo de engrenagem (1 = dentes retos, 2 = helicoidal): ").strip()
helicoidal = (tipo == "2")

# =============================================================================
# DADOS DE ENTRADA (padrões = exemplos usados em aula para cada tipo)
# =============================================================================
print("\n--- Dados de entrada ---")
H_motor = pede_float("Potência do motor (kW)", 9.0 if helicoidal else 0.339)
n_motor = pede_float("Rotação do motor (rpm)", 1000.0 if helicoidal else 11.0)
n_acionado = pede_float("Rotação do eixo acionado (rpm)", 500.0 if helicoidal else 2.2)
reducao = n_motor / n_acionado

R = pede_float("Confiabilidade", 0.95)
ciclos_p = pede_float("Vida útil do pinhão (ciclos)", 1e8 if helicoidal else 1e9)
ciclos_g = ciclos_p / reducao

K_o = pede_float("Fator de sobrecarga Ko (Tabela em aula / Shigley)", 1.25 if helicoidal else 1.0)
FS = pede_float("Fator de segurança mínimo do projeto", 2.0 if helicoidal else 1.0)
Q_v = pede_int("Número de qualidade Qv", 7)

print("\n--- Decisões de projeto ---")
phi_n = pede_float("Ângulo de pressão normal (graus)", 20.0)
psi = pede_float("Ângulo de hélice (graus)", 30.0) if helicoidal else 0.0
N_p = pede_int("Número de dentes do pinhão", 14 if helicoidal else 20)
m_n = pede_float("Módulo normal (mm)", 2.5 if helicoidal else 3.0)

print("\n--- Material ---")
H_brinell = pede_float("Dureza do núcleo (HB)", 397.0 if helicoidal else 350.0)
grau = pede_int("Grau do material AGMA (1 ou 2)", 2)
# Figs. 14-2 e 14-5 (aços endurecidos por completo, em MPa)
if grau == 2:
    S_t_padrao = 0.703 * H_brinell + 113
    S_c_padrao = 2.41 * H_brinell + 237
else:
    S_t_padrao = 0.533 * H_brinell + 88.3
    S_c_padrao = 2.22 * H_brinell + 200
S_t = pede_float("Tensão de flexão admissível St (MPa) (Fig. 14-2 ou Tabela 14-3/14-4)",
                 round(S_t_padrao, 1))
Sc = pede_float("Resistência ao contato admissível Sc (MPa) (Fig. 14-5 ou Tabela 14-6)",
                round(S_c_padrao, 1))
C_p = pede_float("Coeficiente elástico Cp (raiz(MPa)) (Tabela 14-8, aço-aço = 191)", 191.0)

# =============================================================================
# GEOMETRIA (Cap. 13)
# =============================================================================
if helicoidal:
    phi_t = math.degrees(math.atan(math.tan(math.radians(phi_n)) / math.cos(math.radians(psi))))
    m_t = m_n / math.cos(math.radians(psi))     # módulo transversal
    p_x = math.pi * m_n / math.sin(math.radians(psi))  # passo axial
else:
    phi_t = phi_n
    m_t = m_n
    p_x = None

N_g_exato = N_p * reducao
N_g = round(N_g_exato)
if abs(N_g - N_g_exato) > 1e-6:
    print(f"\nAVISO: N_g = {N_g_exato:.2f} não é inteiro; arredondado para {N_g} "
          f"(relação de transmissão real = {N_g / N_p:.4f}).")
m_g = N_g / N_p

d_p = m_t * N_p
d_g = m_t * N_g

# Número mínimo de dentes do pinhão para evitar interferência
# (Eq. 13-11 dentes retos; Eq. 13-22 helicoidais), k = 1 (profundidade completa)
k = 1
sin2 = math.sin(math.radians(phi_t)) ** 2
fator_helice = math.cos(math.radians(psi)) if helicoidal else 1.0
N_p_min = 2 * k * fator_helice * (m_g + math.sqrt(m_g**2 + (1 + 2 * m_g) * sin2)) / ((1 + 2 * m_g) * sin2)

print("\n--- Geometria ---")
print(f"Redução de velocidade: {reducao:.3f}")
print(f"Ângulo de pressão transversal: {phi_t:.3f} graus")
print(f"Módulo transversal: {m_t:.3f} mm")
print(f"Diâmetro primitivo do pinhão: {d_p:.2f} mm")
print(f"Diâmetro primitivo da engrenagem: {d_g:.2f} mm")
if helicoidal:
    print(f"Passo axial: {p_x:.2f} mm")
print(f"Número mínimo de dentes do pinhão (interferência): {math.ceil(N_p_min)}")
if N_p < N_p_min:
    print(f"AVISO: N_p = {N_p} está ABAIXO do mínimo de {math.ceil(N_p_min)} dentes -> risco de interferência!")

# =============================================================================
# CARREGAMENTO E FATOR DINÂMICO
# =============================================================================
V = math.pi * d_p * n_motor / 60000            # velocidade da linha primitiva (m/s)
W_t = 1000 * H_motor / V                       # carregamento tangencial (N)

B = 0.25 * (12 - Q_v) ** (2 / 3)               # Eq. 14-28
A = 50 + 56 * (1 - B)
K_v = ((A + math.sqrt(200 * V)) / A) ** B      # Eq. 14-27 (SI)
V_max = (A + (Q_v - 3)) ** 2 / 200             # Eq. 14-29: velocidade máxima p/ este Qv

print("\n--- Carregamento ---")
print(f"Velocidade da linha primitiva: {V:.3f} m/s (máxima para Qv={Q_v}: {V_max:.2f} m/s)")
if V > V_max:
    print(f"AVISO: V > V_max! Aumente o número de qualidade Qv.")
print(f"Carregamento tangencial Wt: {W_t:.2f} N")
print(f"Fator dinâmico Kv: {K_v:.4f}")

# =============================================================================
# FATORES AGMA QUE NÃO DEPENDEM DA LARGURA DE FACE
# =============================================================================
K_b = 1     # Fator de borda (Fig. 14-16), mB >= 1,2
K_t = 1     # Fator de temperatura (Seção 14-15), T_oleo < 120 C
C_f = 1     # Fator de condição de superfície (Seção 14-2)
C_h = 1     # Fator de razão de dureza (Seção 14-12)

Y_p = lewis_Y(N_p)                              # Tabela 14-2
Y_g = lewis_Y(N_g)
print("\n--- Fatores ---")
print(f"Fator de Lewis: Y_p = {Y_p:.3f} | Y_g = {Y_g:.3f} (Tabela 14-2)")

# Fatores geométricos de flexão J (gráficos: Fig. 14-6 p/ retos;
# Figs. 14-7/14-8 p/ helicoidais) -> entrar com o valor lido
Yj_p = pede_float("Fator geométrico de flexão J do pinhão (Fig. 14-6 ou 14-7/14-8)",
                  0.41 if helicoidal else 0.33)
Yj_g = pede_float("Fator geométrico de flexão J da engrenagem", 0.444 if helicoidal else 0.42)

K_r = 0.658 - 0.0759 * math.log(1 - R)          # Eq. 14-38 (0,5 < R < 0,99)
print(f"Fator de confiabilidade Kr: {K_r:.4f}")

# Fatores de ciclagem (Figs. 14-14 e 14-15) - curvas usuais p/ N > 10^7
curva = pede_int("Curva de ciclagem (1 = nominal, 2 = conservadora)", 1)
if curva == 2:
    Yn_p, Yn_g = 1.6831 * ciclos_p**-0.0323, 1.6831 * ciclos_g**-0.0323
    Zn_p, Zn_g = 2.466 * ciclos_p**-0.056, 2.466 * ciclos_g**-0.056
else:
    Yn_p, Yn_g = 1.3558 * ciclos_p**-0.0178, 1.3558 * ciclos_g**-0.0178
    Zn_p, Zn_g = 1.4488 * ciclos_p**-0.023, 1.4488 * ciclos_g**-0.023
print(f"Fator de ciclagem de flexão: Yn_p = {Yn_p:.4f} | Yn_g = {Yn_g:.4f}")
print(f"Fator de ciclagem de contato: Zn_p = {Zn_p:.4f} | Zn_g = {Zn_g:.4f}")

# =============================================================================
# FATOR GEOMÉTRICO DE CONTATO I (Eqs. 14-21 a 14-25)
# =============================================================================
if helicoidal:
    rb_p = (d_p / 2) * math.cos(math.radians(phi_t))
    rb_g = (d_g / 2) * math.cos(math.radians(phi_t))
    a = m_n                                     # adendo (Tabela 13-4)
    Z = (math.sqrt((d_p / 2 + a)**2 - rb_p**2) + math.sqrt((d_g / 2 + a)**2 - rb_g**2)
         - (d_p / 2 + d_g / 2) * math.sin(math.radians(phi_t)))   # Eq. 14-25
    P_N = math.pi * m_n * math.cos(math.radians(phi_n))           # Eq. 14-24
    m_N = P_N / (0.95 * Z)                      # Eq. 14-21: razão de compartilhamento de carga
else:
    m_N = 1                                     # dentes retos
Z_I = (math.cos(math.radians(phi_t)) * math.sin(math.radians(phi_t)) / (2 * m_N)) * m_g / (m_g + 1)  # Eq. 14-23
print(f"Fator geométrico de contato I: {Z_I:.4f}")

# =============================================================================
# LARGURA DE FACE REQUERIDA (iterativa, pois Ks e Kh dependem de b)
# =============================================================================
def largura_requerida(b):
    Ks_p = fator_tamanho(b, m_t, Y_p)
    K_h = fator_distribuicao(b, d_p)
    # Flexão: Eq. 14-15 com sigma = St*Yn/(FS*Kt*Kr), isolando b
    b_flex = FS * W_t * K_o * K_v * Ks_p * (1 / m_t) * (K_h * K_b / Yj_p) * (K_t * K_r) / (S_t * Yn_p)
    # Contato: Eq. 14-16 com sigma_c = Sc*Zn*Ch/(Kt*Kr) e SH^2 = FS, isolando b
    b_desg = ((C_p * K_t * K_r / (Sc * Zn_p * C_h)) ** 2) * FS * W_t * K_o * K_v * Ks_p * (K_h / d_p) * (C_f / Z_I)
    return b_flex, b_desg

b_iter = 30.0
for _ in range(100):
    b_flex, b_desg = largura_requerida(b_iter)
    b_req = max(b_flex, b_desg)
    if abs(b_req - b_iter) < 1e-4:
        break
    b_iter = b_req

# Intervalo recomendado de largura de face
if helicoidal:
    b_min = 2 * p_x        # recomendação: no mínimo 2 passos axiais (Seção 13-11)
    b_max = None
    intervalo_txt = f"b >= {b_min:.2f} mm (2 passos axiais)"
else:
    b_min = 3 * math.pi * m_t
    b_max = 5 * math.pi * m_t
    intervalo_txt = f"{b_min:.2f} a {b_max:.2f} mm (3*pi*m a 5*pi*m)"

print("\n" + "=" * 70)
print(" LARGURA DE FACE")
print("=" * 70)
print(f"Requerida pela flexão: {b_flex:.2f} mm")
print(f"Requerida pelo desgaste (contato): {b_desg:.2f} mm")
print(f"Requerida (a maior das duas): {b_req:.2f} mm")
print(f"Intervalo recomendado: {intervalo_txt}")

b_sugerido = math.ceil(max(b_req, b_min))
if b_max is not None and b_req > b_max:
    print(f"AVISO: a largura requerida ({b_req:.2f} mm) excede o máximo recomendado "
          f"({b_max:.2f} mm) -> aumente o módulo ou o número de dentes e recalcule.")
print(f"Sugestão: b = {b_sugerido} mm")

b = pede_float("Largura de face adotada (mm)", float(b_sugerido))
if b < b_req:
    print(f"AVISO: b adotado é menor que o requerido ({b_req:.2f} mm) -> fatores de segurança abaixo do especificado.")
if b_max is not None and not (b_min <= b <= b_max):
    print(f"AVISO: b adotado está fora do intervalo recomendado ({intervalo_txt}).")
elif b_max is None and b < b_min:
    print(f"AVISO: b adotado está abaixo do mínimo recomendado ({b_min:.2f} mm).")

# =============================================================================
# TENSÕES E FATORES DE SEGURANÇA COM O b ADOTADO
# =============================================================================
Ks_p = fator_tamanho(b, m_t, Y_p)
Ks_g = fator_tamanho(b, m_t, Y_g)
K_h = fator_distribuicao(b, d_p)

# Tensão de flexão (Eq. 14-15)
Tflex_p = W_t * K_o * K_v * Ks_p * (1 / (b * m_t)) * (K_h * K_b / Yj_p)
Tflex_g = W_t * K_o * K_v * Ks_g * (1 / (b * m_t)) * (K_h * K_b / Yj_g)

# Tensão de contato (Eq. 14-16)
Tcont_p = C_p * math.sqrt(W_t * K_o * K_v * Ks_p * (K_h / (d_p * b)) * (C_f / Z_I))
Tcont_g = C_p * math.sqrt(W_t * K_o * K_v * Ks_g * (K_h / (d_p * b)) * (C_f / Z_I))

# Fatores de segurança AGMA (Eqs. 14-41 e 14-42)
SF_p = (S_t * Yn_p) / (K_t * K_r * Tflex_p)
SF_g = (S_t * Yn_g) / (K_t * K_r * Tflex_g)
SH_p = (Sc * Zn_p * C_h) / (K_t * K_r * Tcont_p)
SH_g = (Sc * Zn_g * C_h) / (K_t * K_r * Tcont_g)

print("\n" + "=" * 70)
print(f" RESULTADOS - DENTES {'HELICOIDAIS' if helicoidal else 'RETOS'} | b = {b:.1f} mm")
print("=" * 70)
print(f"Ks pinhão = {Ks_p:.4f} | Ks engrenagem = {Ks_g:.4f} | Kh = {K_h:.4f}")
print("-" * 70)
print(f"Tensão de flexão no pinhão: {Tflex_p:.2f} MPa")
print(f"Tensão de flexão na coroa: {Tflex_g:.2f} MPa")
print(f"Fator de segurança AGMA de flexão (SF) - pinhão: {SF_p:.2f}")
print(f"Fator de segurança AGMA de flexão (SF) - coroa: {SF_g:.2f}")
print("-" * 70)
print(f"Tensão de contato no pinhão: {Tcont_p:.2f} MPa")
print(f"Tensão de contato na coroa: {Tcont_g:.2f} MPa")
print(f"Fator de segurança AGMA de contato (SH) - pinhão: {SH_p:.2f} (SH^2 = {SH_p**2:.2f})")
print(f"Fator de segurança AGMA de contato (SH) - coroa: {SH_g:.2f} (SH^2 = {SH_g**2:.2f})")
print("-" * 70)
print("Obs.: para comparar flexão e contato na mesma base, compare SF com SH^2 (Seção 14-19).")
print("=" * 70)
