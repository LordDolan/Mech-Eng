# Mech-Eng

Repositório com códigos desenvolvidos ao longo da graduação em Engenharia Mecânica :)

---

## 📁 Projeto_Transcal.cpp — Transferência de Calor em Tubo

Simulação numérica, em C++, do resfriamento de água quente escoando por um tubo exposto a um ambiente frio, considerando **convecção interna**, **convecção externa** e **radiação** na superfície.

### O problema físico

Água a **90 °C** escoa a **1,2 m/s** por um tubo de **40 mm** de diâmetro e **15 m** de comprimento, com parede de **1 mm** de espessura, em um ambiente a **10 °C**. O programa calcula, ao longo do comprimento do tubo:

- **T_w(x)** — temperatura da superfície do tubo;
- **T_m(x)** — temperatura média do fluido;
- **Q_dto(x)** — taxa de transferência de calor.

### Etapas do cálculo

1. **Número de Reynolds** (`Re_D = U·D/ν`) para classificar o regime de escoamento:
   - Laminar (`Re < 2300`): `Nu_D = 3,66` (temperatura de superfície constante);
   - Turbulento: correlação de **Dittus-Boelter** para resfriamento, `Nu_D = 0,023·Re^0,8·Pr^0,3`.
2. **Coeficiente de convecção interno** `h_i = Nu_D·k/D` e **número de Biot** da parede `Bi = h_i·e/k_parede`.
3. **Radiação linearizada**: o coeficiente `h_rad = ε·σ·(T_w² + T_c²)·(T_w + T_c)` transforma o termo de radiação em algo com a mesma forma da convecção.
4. **Método de Newton-Raphson** para resolver o balanço de energia na superfície e obter a temperatura inicial da parede `T_w(0)`:

   ```
   f(T_w) = T_mi − T_w − (1 + Bi)·[(h_e/h_i)·(T_w − T_∞) + (h_rad/h_i)·(T_w − T_c)] = 0
   ```

5. **Método de Runge-Kutta de 4ª ordem (RK4)** para integrar a EDO do balanço de energia por unidade de comprimento ao longo do tubo:

   ```
   dT_w/dx = p·[h_i·(T_mi − T_w) − h_e·(T_w − T_∞) − h_rad·(T_w − T_c)] / (ṁ·c_p)
   ```

6. **Gráficos** de `T_w(x)`, `T_m(x)` e `Q_dto(x)` gerados via [matplotlib-cpp](https://github.com/lava/matplotlib-cpp).

### Como compilar e executar

O código depende do header [`matplotlibcpp.h`](https://github.com/lava/matplotlib-cpp) (arquivo único), do **Python 3 com dev headers** e do pacote **matplotlib** instalado no Python:

```bash
# Dependências (exemplo em Ubuntu/Debian)
sudo apt install g++ python3-dev python3-matplotlib
wget https://raw.githubusercontent.com/lava/matplotlib-cpp/master/matplotlibcpp.h

# Compilação (ajuste a versão do Python se necessário)
g++ -std=c++17 -DWITHOUT_NUMPY -I. -I/usr/include/python3.12 Projeto_Transcal.cpp -o transcal -lpython3.12

# Execução
./transcal
```

> Se o NumPy estiver instalado com os headers disponíveis, pode-se omitir `-DWITHOUT_NUMPY` e adicionar o include do NumPy.

### Saída esperada

```
Re_D: 152075
h_i: 6641.92
Bi: 0.127729
A temperatura inicial da superfície Tw_i é: 362.764 K
x: 0 m, T_w: 362.764 K, T_m: 363 K, Q_dto: 2954.56 W
...
```

Como `h_i ≫ h_e`, a resistência dominante é a externa: a parede fica praticamente na temperatura do fluido (≈ 362,8 K) e varia muito pouco ao longo dos 15 m.

### 🔧 Correções aplicadas em relação à versão original

1. **Chamada de plotagem que não compilava** — `plt::plot(x, y, "r-", {{"label", ...}})` não corresponde a nenhuma sobrecarga do matplotlib-cpp (erro de compilação verificado). Substituída por `plt::named_plot("T_w(x)", x, y, "r-")`, que aplica o rótulo da legenda corretamente.
2. **Viscosidade cinemática errada** — `nu` estava definida como `0.001 m²/s`, que na verdade é o valor da viscosidade *dinâmica* da água a 20 °C em Pa·s. Isso resultava em `Re ≈ 48` (laminar), fisicamente implausível para água a 1,2 m/s em tubo de 40 mm. Corrigida para `nu = mu/rho ≈ 3,16·10⁻⁷ m²/s`, o que leva a `Re ≈ 1,5·10⁵` (turbulento) e muda completamente o valor de `h_i`.
3. **Derivada errada no Newton-Raphson** — `function_derivative` não era a derivada de `function`: continha termos da EDO (`h_i·p/(ṁ·c_p)`) e um erro de digitação na derivada da radiação (`T_c³` no lugar de `T_c²`). Substituída pela derivada analítica correta:
   `f'(T_w) = −1 − (1+Bi)/h_i · [h_e + h_rad + (dh_rad/dT_w)·(T_w − T_c)]`, com `dh_rad/dT_w = ε·σ·(3T_w² + 2T_c·T_w + T_c²)`.
4. **Perímetro faltando na EDO** — em `dTwdx`, o fluxo de calor (W/m²) precisa ser multiplicado pelo perímetro `p` (m) para virar taxa por unidade de comprimento (W/m) antes de dividir por `ṁ·c_p` (W/K); sem isso a equação era dimensionalmente inconsistente (a variável `p` era declarada mas nunca usada).
5. **Comentários das propriedades** — os valores de `mu`, `c_p`, `Pr` e `k` correspondem à água a ~90 °C (temperatura de operação), não a 25 °C como indicavam os comentários; comentários ajustados. A expressão de `h_rad` foi extraída para a função `h_radiation()` para eliminar código repetido em 4 lugares.

### ⚠️ Pontos de atenção (não alterados)

- `rho = 997 kg/m³` é o valor a 25 °C; a ~90 °C seria ≈ 965 kg/m³. Isso afeta `ṁ` e `Re` em ~3 %.
- `k_iso = 52 W/m·K` está comentado como "isolante", mas esse valor é típico de **aço-carbono** — coerente com o uso no Biot da parede do tubo, apenas o nome/comentário que pode confundir.
- `Q_dto(x)` usa a área **total** de troca (`A_troca`, os 15 m inteiros) com a temperatura local `T_w(x)`, ou seja, representa "a taxa total caso todo o tubo estivesse à temperatura daquele ponto", e não a taxa acumulada até `x`.
