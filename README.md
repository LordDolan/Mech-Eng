# Mech-Eng

Repositório com códigos desenvolvidos ao longo da graduação em Engenharia Mecânica :)

---

## Projeto_Transcal.cpp

Simulação em C++ do resfriamento de água quente escoando por um tubo exposto a um ambiente frio, considerando convecção interna, convecção externa e radiação. Usa Newton-Raphson para achar a temperatura inicial da parede e Runge-Kutta de 4ª ordem para integrar a temperatura ao longo do tubo.

**Entrada:** as propriedades definidas como constantes no início do código — água a 90 °C escoando a 1,2 m/s num tubo de 40 mm de diâmetro, 15 m de comprimento e parede de 1 mm, com ambiente a 10 °C.

**Saída:** no terminal, o número de Reynolds, o coeficiente de convecção interno `h_i`, o número de Biot e uma tabela com `T_w(x)` (temperatura da parede), `T_m(x)` (temperatura do fluido) e `Q_dto(x)` (taxa de calor) ao longo do tubo; além de dois gráficos (temperaturas e taxa de calor vs. comprimento).

**Como compilar** (requer o header [matplotlibcpp.h](https://github.com/lava/matplotlib-cpp), Python 3 com dev headers e matplotlib):

```bash
g++ -std=c++17 -DWITHOUT_NUMPY -I. -I/usr/include/python3.12 Projeto_Transcal.cpp -o transcal -lpython3.12
./transcal
```
