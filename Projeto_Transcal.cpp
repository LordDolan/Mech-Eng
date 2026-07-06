#include <iostream>
#include <cmath>
#include <vector>
#include "matplotlibcpp.h" // Inclua a biblioteca matplotlib-cpp

namespace plt = matplotlibcpp; // Namespace para matplotlib-cpp

// Propriedades Termofluidodinâmicas
const float rho = 997.0; // Massa específica da água (kg/m^3), 25C
const float mu = 3.146878e-4; // Viscosidade dinâmica da água (Pa*s), 90C
const float nu = mu / rho; // Viscosidade cinemática da água (m^2/s) = mu/rho ~ 3.16e-7
const float c_p = 4206; // Calor específico da água (J/kg*K), 90C
const float Pr = 1.96; // Número de Prandtl da água, 90C
const float k = 0.675; // Condutividade térmica da água (W/m*K), 90C
const float U = 1.2; // Velocidade do fluido (m/s)
const float epsilon = 0.7; // Emissividade do isolante
const float D = 0.04; // Diâmetro do tubo (m)
const float L = 15; // Comprimento do tubo (m)
const float k_iso = 52; // Condutividade térmica da parede do tubo (W/m*K)
const float e = 0.001; // Espessura do tubo (m)
const float h_e = 12; // Coeficiente de transferência de calor externo (W/m^2*K)
const float T_mi = 273 + 90; // Temperatura média interna do fluído (K)
const float T_c = 273 + 10; // Temperatura dos ambientes (K)
const float T_inf = 273 + 10; // Temperatura ambiente (K)
const float sigma = 5.67e-8; // Constante de Stefan-Boltzmann
float A_troca = 2 * M_PI * (D/2) * L; // Área de troca térmica (m^2)
float A = M_PI * pow(D, 2) / 4; // Área da seção transversal do tubo (m^2)
float m_dot = rho * U * A; // Vazão mássica do fluido (kg/s)
float p = M_PI * D; // Perímetro do tubo (m)

float Nu_D; // Número de Nusselt

// Coeficiente de transferência de calor por radiação linearizado:
// h_rad = epsilon*sigma*(T_w^2 + T_c^2)*(T_w + T_c) = epsilon*sigma*(T_w^3 + T_c*T_w^2 + T_c^2*T_w + T_c^3)
float h_radiation(float T_w) {
    return epsilon * sigma * (pow(T_w,3) + T_c * pow(T_w,2) + pow(T_c,2) * T_w + pow(T_c,3));
}

// Função de Twi: f(T_w) = 0 no balanço de energia da superfície
float function(float T_w, float h_i, float Bi) {
    float h_rad = h_radiation(T_w);
    float conv = (h_e / h_i) * (T_w - T_inf);
    float rad = (h_rad / h_i) * (T_w - T_c);

    return T_mi - T_w - (1 + Bi) * (conv + rad);
}

// Derivada da função f(T_w) em relação a T_w (usada no Newton-Raphson)
float function_derivative(float T_w, float h_i, float Bi) {
    float h_rad = h_radiation(T_w);
    float dhrad_dTw = epsilon * sigma * (3 * pow(T_w,2) + 2 * T_c * T_w + pow(T_c,2));
    return -1 - ((1 + Bi) / h_i) * (h_e + h_rad + dhrad_dTw * (T_w - T_c));
}

// Método de Newton-Raphson para encontrar a raiz
float newton_raphson(float initial_guess, float h_i, float Bi, float tolerance = 1e-6, int max_iterations = 100) {
    float T_w = initial_guess;
    for (int i = 0; i < max_iterations; ++i) {
        float f = function(T_w, h_i, Bi);
        float f_prime = function_derivative(T_w, h_i, Bi);

        if (fabs(f_prime) < 1e-12) {
            std::cerr << "Derivada muito pequena, pode não convergir." << std::endl;
            break;
        }

        float T_w_next = T_w - f / f_prime;

        if (fabs(T_w_next - T_w) < tolerance) {
            return T_w_next;
        }

        T_w = T_w_next;
    }

    std::cerr << "O método de Newton-Raphson não convergiu." << std::endl;
    return T_w;
}

// Função que define a ODE dT_w/dx = f(x, T_w)
// Balanço de energia por unidade de comprimento: multiplica-se pelo perímetro p
float dTwdx(float x, float T_w, float h_i, float Bi) {
    float h_rad = h_radiation(T_w);
    return p * (h_i * (T_mi - T_w) - h_e * (T_w - T_inf) - h_rad * (T_w - T_c)) / (m_dot * c_p);
}

// Implementação do método de Runge-Kutta de 4a ordem
void runge_kutta(float h_i, float Bi, float dx, int n_steps, std::vector<float>& x_values, std::vector<float>& T_w_values) {
    float x = 0.0;
    float T_w = T_w_values[0]; // Temperatura inicial (valor de chute inicial)

    for (int i = 0; i < n_steps; ++i) {
        float k1 = dx * dTwdx(x, T_w, h_i, Bi);
        float k2 = dx * dTwdx(x + dx / 2.0, T_w + k1 / 2.0, h_i, Bi);
        float k3 = dx * dTwdx(x + dx / 2.0, T_w + k2 / 2.0, h_i, Bi);
        float k4 = dx * dTwdx(x + dx, T_w + k3, h_i, Bi);

        T_w += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        x += dx;

        x_values.push_back(x);
        T_w_values.push_back(T_w);
    }
}

// Função para calcular T_m(x)
float calculate_Tm(float T_w, float h_i, float Bi) {
    float h_rad = h_radiation(T_w);
    return T_w + (1 + Bi) * ( (h_e / h_i) * (T_w - T_inf) + (h_rad / h_i) * (T_w - T_c) );
}

// Função para calcular Q_dto(x)
float calculate_Qdto(float T_w, float h_i) {
    return h_i * A_troca * (T_mi - T_w);
}

int main() {
    // Cálculo do Reynolds
    float Re_D = U * D / nu;
    std::cout << "Re_D: " << Re_D << std::endl;

    // Cálculo do Nusselt
    if (Re_D < 2300) { // Regime laminar
        Nu_D = 3.66; // Assumindo (T_s constante)
    } else { // Regime turbulento
        Nu_D = 0.023 * pow(Re_D, 0.8) * pow(Pr, 0.3);
    }

    // Cálculo do Coeficiente de transferência de calor por convecção interno (h_i)
    float h_i = (Nu_D * k) / D;
    std::cout << "h_i: " << h_i << std::endl;

    // Cálculo do Biot (Bi)
    float Bi = (h_i * e) / k_iso;
    std::cout << "Bi: " << Bi << std::endl;

    // Determinação da Temperatura da superfície interna do tubo (Tw_i) - Newton Raphson
    float initial_guess = (T_mi + T_c) * 0.5; // Chute inicial de temperatura em Kelvin
    float T_w_initial = newton_raphson(initial_guess, h_i, Bi);
    std::cout << "A temperatura inicial da superfície Tw_i é: " << T_w_initial << " K" << std::endl;

    // Parâmetros para o método de Runge-Kutta
    float dx = 0.1; // Passo no comprimento do tubo
    int n_steps = static_cast<int>(L / dx); // Número de passos
    std::vector<float> x_values = {0.0}; // Valores de x
    std::vector<float> T_w_values = {T_w_initial}; // Valores de Tw(x)
    std::vector<float> T_m_values; // Valores de Tm(x)
    std::vector<float> Q_dto_values; // Valores de Q_dto(x)

    // Executa o método de Runge-Kutta
    runge_kutta(h_i, Bi, dx, n_steps, x_values, T_w_values);

    // Calcula T_m(x) e Q_dto(x) para cada valor de x
    for (size_t i = 0; i < x_values.size(); ++i) {
        T_m_values.push_back(calculate_Tm(T_w_values[i], h_i, Bi));
        Q_dto_values.push_back(calculate_Qdto(T_w_values[i], h_i));
    }

    // Exibe os resultados na mesma linha
    for (size_t i = 0; i < x_values.size(); ++i) {
        std::cout << "x: " << x_values[i] << " m, T_w: " << T_w_values[i] << " K, ";
        std::cout << "T_m: " << T_m_values[i] << " K, ";
        std::cout << "Q_dto: " << Q_dto_values[i] << " W" << std::endl;
    }

    // Plota os gráficos
    plt::figure();
    plt::named_plot("T_w(x)", x_values, T_w_values, "r-");
    plt::named_plot("T_m(x)", x_values, T_m_values, "b-");
    plt::xlabel("Comprimento do tubo (x) [m]");
    plt::ylabel("Temperatura [K]");
    plt::legend();
    plt::title("Temperaturas ao longo do comprimento do tubo");
    plt::grid(true);

    plt::figure();
    plt::named_plot("Q_dto(x)", x_values, Q_dto_values, "g-");
    plt::xlabel("Comprimento do tubo (x) [m]");
    plt::ylabel("Taxa de transferência de calor (Q_dto) [W]");
    plt::legend();
    plt::title("Taxa de transferência de calor ao longo do comprimento do tubo");
    plt::grid(true);

    plt::show();

    return 0;
}
