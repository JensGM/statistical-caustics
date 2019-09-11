#include <array>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <random>

#include <glm/glm.hpp>

using namespace std::literals::complex_literals;
using
    std::complex,
    std::conj,
    std::exp,
    std::imag,
    std::real,
    std::sqrt,

    glm::vec2,
    glm::dot,
    glm::length,
    glm::normalize;

constexpr auto g = 9.81f; // gravity
constexpr auto pi = 3.14159265359f;
constexpr auto recp_root_2 = 0.70710678118f; // 1 / sqrt(2)
constexpr auto A = 1.0f; // arbitrary amplitude
constexpr auto L = 70.0f; // field dimension (i think in meters)
constexpr auto N = 256; // grid resolution
constexpr auto w = vec2{6.0f, 5.0f}; // horizontal wind vector

using Field = std::array<std::array<complex<float>, N>, N>;

template<typename T>
auto remap_range(T v, T lo0, T hi0, T lo1, T hi1) {
    return lo1 + (v - lo0) * (hi1 - lo1) / (hi0 - lo0);
}

Field fft(const Field& field) {
    auto* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    auto* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    auto plan = fftw_plan_dft_2d(N, N, in, out, FFTW_BACKWARD, FFTW_MEASURE);

    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            auto idx = m * N + n;
            in[idx][0] = real(field[n][m]);
            in[idx][1] = imag(field[n][m]);
        }
    }

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    Field result;
    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            auto idx = m * N + n;
            float r = out[idx][0],
                  i = out[idx][1];
            result[n][m] = {r, i};

            // HACK, something something grid translation.
            // https://www.keithlantz.net/2011/11/ocean-simulation-part-two-using-the-fast-fourier-transform/
            auto sign = (n + m) % 2 == 0 ? 1.0f : -1.0f;
            result[n][m] *= sign;
        }
    }

    fftw_free(in);
    fftw_free(out);

    return result;
}

auto idx2k(int i, int k) {
    auto n = remap_range<float>(i, 0, N - 1, -N / 2, N / 2);
    auto m = remap_range<float>(k, 0, N - 1, -N / 2, N / 2);
    return 2.0f * pi * vec2{n, m} / L;
}

auto gaussian() {
    static std::random_device rd{};
    static std::mt19937 gen{rd()};
    static std::normal_distribution dist{0.0f, 1.0f};
    return dist(gen);
}

auto 🌊(vec2 k) {
    // max_wave is the second L in the paper. There are two. In the same font.
    auto windspeed = length(w);
    auto max_wave = windspeed * windspeed / g;

    auto k_unit = normalize(k);
    auto w_unit = normalize(w);
    auto kw = dot(k_unit, w_unit);
    auto kw2 = kw * kw;
    auto kw4 = kw2 * kw2;
    auto kw6 = kw4 * kw2; // Modification in the paper

    auto k_length = length(k);
    auto k2 = k_length * k_length;
    auto k4 = k2 * k2;

    // k_max_ws2 is (kL)^2 in the paper.
    auto k_max_ws = k_length * max_wave;
    auto k_max_ws2 = k_max_ws * k_max_ws;

    // We return a complex number to allow sqrt of negative numbers
    return complex{A * kw6 * exp(-1.0f / k_max_ws2) / k4, 0.0f};
}

auto hₒ(vec2 k) {
    complex ξ{gaussian(), gaussian()};
    return recp_root_2 * ξ * sqrt(🌊(k)); // √ is not in the allowed range...
}

auto w₀(vec2 k) {
    return sqrt(g * length(k));
}

int main() {
    Field spectrum;
    Field specconj;

    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            auto k = idx2k(n, m);
            spectrum[n][m] =      hₒ( k);
            specconj[n][m] = conj(hₒ(-k));
        }
    }

    float t = 1.0f;

    for (int f = 0; f < 240; ++f) {
        t += 60.f / 1000.f;
        Field frequencies;

        for (int n = 0; n < N; ++n) {
            for (int m = 0; m < N; ++m) {
                auto k = idx2k(n, m);
                auto wt = w₀(k) * t;
                auto c₀ = cos( wt) + 1.0if * sin( wt);
                auto c₁ = cos(-wt) + 1.0if * sin(-wt);
                frequencies[n][m] = spectrum[n][m] * c₀ + specconj[n][m] * c₁;
            }
        }

        Field heights = fft(frequencies);

        for (int n = 0; n < N; ++n) {
            for (int m = 0; m < N; ++m) {
                std::cout << real(heights[n][m]) << "\t";
            }
            std::cout << '\n';
        }
    }
}

static_assert(std::is_same_v<complex<float>, decltype(hₒ(vec2{}))>,
              "The final return type is not of type complex<float>");
