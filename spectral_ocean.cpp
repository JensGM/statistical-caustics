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
constexpr auto w = vec2{10.0f, 0.0f}; // horizontal wind vector

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

            // HACK
            result[n][m] *= n % 2 == 0 ? 1.0f : -1.0f;
            result[n][m] *= m % 2 == 0 ? 1.0f : -1.0f;
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

auto phillips(vec2 k) {
    // max_wave is the second L in the paper. There are two. In the same font.
    auto windspeed = length(w);
    auto max_wave = windspeed * windspeed / g;

    auto k_unit = normalize(k);
    auto w_unit = normalize(w);
    auto kw = dot(k_unit, w_unit);

    auto k_length = length(k);
    auto k_2 = k_length * k_length;
    auto k_4 = k_2 * k_2;

    // k_max_ws2 is (kL)^2 in the paper.
    auto k_max_ws = k_length * max_wave;
    auto k_max_ws2 = k_max_ws * k_max_ws;
    auto k_max_ws4 = k_max_ws2 * k_max_ws2;
    auto k_max_ws8 = k_max_ws4 * k_max_ws4;
    auto k_max_ws16 = k_max_ws8 * k_max_ws8;

    // We return a complex number to allow sqrt of negative numbers
    return complex{A * kw * exp(-1.0f / k_max_ws16) / k_4, 0.0f};
}

auto h0(vec2 k) {
    complex E{gaussian(), gaussian()};
    return recp_root_2 * E * sqrt(phillips(k));
}

auto w0(vec2 k) {
    return sqrt(g * length(k));
}

int main() {
    Field spectrum;
    Field specconj;

    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            auto k = idx2k(n, m);
            spectrum[n][m] =      h0( k);
            specconj[n][m] = conj(h0(-k));
        }
    }

    Field frequencies;

    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            auto k = idx2k(n, m);
            frequencies[n][m] = spectrum[n][m] * exp(-1.0if * w0(k) * 1.0f)
                              + specconj[n][m] * exp( 1.0if * w0(k) * 1.0f);
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

static_assert(std::is_same_v<complex<float>, decltype(h0(vec2{}))>,
              "The final return type is not of type complex<float>");
