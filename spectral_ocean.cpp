#include <array>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <random>

#include <glm/glm.hpp>

using namespace std::literals::complex_literals;
using
    std::abs,
    std::complex,
    std::conj,
    std::exp,
    std::imag,
    std::pow,
    std::real,
    std::sqrt,

    glm::vec2,
    glm::dot,
    glm::length,
    glm::normalize;

constexpr auto 𝑔 = 9.81f; // gravity
constexpr auto π = 3.14159265359f;
constexpr auto recp_root_2 = 0.70710678118f; // 1/sqrt(2),sqrt is not constexpr
constexpr auto 𝐴 = 1.0f; // arbitrary amplitude
constexpr auto L = 70.0f; // field dimension (i think in meters)
constexpr auto N = 256; // grid resolution
constexpr auto 𝑤 = vec2{6.0f, 5.0f}; // horizontal wind vector

using Field = std::array<std::array<complex<float>, N>, N>;

template<typename T>
auto remap_range(T v, T lo0, T hi0, T lo1, T hi1) {
    return lo1 + (v - lo0) * (hi1 - lo1) / (hi0 - lo0);
}

auto fft(const Field& field) {
    auto* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    auto* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);

    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            auto idx = m * N + n;
            in[idx][0] = real(field[n][m]);
            in[idx][1] = imag(field[n][m]);
        }
    }

    auto plan = fftw_plan_dft_2d(N, N, in, out, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_execute(plan);

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

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return result;
}

auto idx2k(int i, int k) {
    auto n = remap_range<float>(i, 0, N - 1, -N / 2, N / 2);
    auto m = remap_range<float>(k, 0, N - 1, -N / 2, N / 2);
    return 2.0f * π * vec2{n, m} / L;
}

auto gaussian() {
    static std::random_device rd{};
    static std::mt19937 gen{rd()};
    static std::normal_distribution dist{0.0f, 1.0f};
    return dist(gen);
}

/**
 * Phillips frequency for 𝐤
 */
auto Ph(vec2 𝐤) {
    auto 𝑉 = length(𝑤);
    auto 𝐿 = 𝑉 * 𝑉 / 𝑔; // Largets possible wave

    auto ḱ = normalize(𝐤);
    auto ẃ = normalize(𝑤);

    // Modification in the paper ḱ·ẃ² -> ḱ·ẃ⁶
    auto ḱ·ẃ⁶ = pow(abs(dot(ḱ, ẃ)), 6.0f); // Note, dot product to power 6
    auto 𝑘 = length(𝐤);
    auto 𝑘⁴ = pow(𝑘, 4.0f);
    auto 𝑘𝐿² = pow(𝑘 * 𝐿, 2.0f); // Note, product to power of 2

    return 𝐴 * exp(-1.0f / 𝑘𝐿²) / 𝑘⁴ * ḱ·ẃ⁶;
}

/**
 * Fourier amplitude, ℎₒ with a tilde in the paper
 */
auto ℎₒ(vec2 𝐤) {
    complex ξ{gaussian(), gaussian()};
    return recp_root_2 * ξ * sqrt(Ph(𝐤));
}

/**
 * Dispersion relation
 */
auto ω(vec2 𝐤) {
    return sqrt(𝑔 * length(𝐤));
}

int main() {
    Field 𝐇₀₊;
    Field 𝐇₀₋; // Complex conjugate

    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            auto 𝐤 = idx2k(n, m);
            𝐇₀₊[n][m] =      ℎₒ( 𝐤);
            𝐇₀₋[n][m] = conj(ℎₒ(-𝐤));
        }
    }

    float t = 0.0f;

    for (int f = 0; f < 240; ++f) {
        t += 60.f / 1000.f;
        Field 𝐇;

        for (int n = 0; n < N; ++n) {
            for (int m = 0; m < N; ++m) {
                auto 𝐤 = idx2k(n, m);
                auto ωt = ω(𝐤) * t;
                auto c₀ = cos( ωt) + 1.0if * sin( ωt);
                auto c₁ = cos(-ωt) + 1.0if * sin(-ωt);
                𝐇[n][m] = 𝐇₀₊[n][m] * c₀ + 𝐇₀₋[n][m] * c₁;
            }
        }

        Field heights = fft(𝐇);

        for (int n = 0; n < N; ++n) {
            for (int m = 0; m < N; ++m) {
                std::cout << real(heights[n][m]) << "\t";
            }
            std::cout << '\n';
        }
    }
}

static_assert(std::is_same_v<complex<float>, decltype(ℎₒ(vec2{}))>,
              "The final return type is not of type complex<float>");
