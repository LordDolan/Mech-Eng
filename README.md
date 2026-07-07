# Mech-Eng

Repositório com códigos desenvolvidos ao longo da graduação em Engenharia Mecânica :)

---

## Projeto_Transcal.cpp

Simula o escoamento de água quente em um tubo, representando o comportamento da água em um ambiente espacial: o tubo perde calor para o ambiente frio por convecção e radiação, e o programa calcula como a temperatura evolui ao longo do comprimento.

**Métodos aplicados:** correlação de Dittus-Boelter (convecção interna), Newton-Raphson (temperatura inicial da parede) e Runge-Kutta de 4ª ordem (integração da temperatura ao longo do tubo).

**Entrada:** constantes definidas no início do código — água a 90 °C e 1,2 m/s, tubo de 40 mm de diâmetro e 15 m de comprimento, ambiente a 10 °C.

**Saída:** Reynolds, coeficiente de convecção `h_i`, Biot, tabela de `T_w(x)`, `T_m(x)` e `Q_dto(x)` no terminal, e gráficos das temperaturas e da taxa de calor ao longo do tubo.

**Compilação** (requer [matplotlibcpp.h](https://github.com/lava/matplotlib-cpp) e Python 3 com matplotlib):

```bash
g++ -std=c++17 -DWITHOUT_NUMPY -I. -I/usr/include/python3.12 Projeto_Transcal.cpp -o transcal -lpython3.12
```

---

## PM1_Engrenagens_AGMA.py

Dimensiona engrenagens cilíndricas de dentes retos ou helicoidais pelo método AGMA (Shigley, cap. 14), verificando falha por flexão e por desgaste (contato).

**Métodos aplicados:** correlações AGMA de tensão e resistência, fator de Lewis interpolado (Tabela 14-2), cálculo iterativo da largura de face mínima e verificação de interferência (cap. 13).

**Entrada:** 4 perguntas — tipo de dente (reto/helicoidal), número de dentes do pinhão, material (menu com aços graus 1 e 2 ou St/Sc personalizados) e largura de face. Os demais parâmetros (potência, rotações, módulo, FS...) são padrões editáveis no topo do código.

**Saída:** resumo no terminal com tensões e fatores de segurança AGMA (veredito passa/não passa), tabela gráfica dos fatores com a localização de cada um no Shigley, e mapa 2D "nº de dentes × largura de face" com as regiões passa (verde) e não passa (vermelho).

**Execução** (requer matplotlib):

```bash
pip install matplotlib
python3 PM1_Engrenagens_AGMA.py
```
