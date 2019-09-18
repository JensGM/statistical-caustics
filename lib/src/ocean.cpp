#include <cmath>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <random>
#include <utility>

#include <ocean/ocean.hpp>

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

namespace {

constexpr auto 𝑔 = 9.81f; // gravity
constexpr auto π = 3.14159265359f;
constexpr auto recp_root_2 = 0.70710678118f; // 1/sqrt(2),sqrt is not constexpr

template<typename T>
auto remap_range(T v, T lo0, T hi0, T lo1, T hi1) {
    return lo1 + (v - lo0) * (hi1 - lo1) / (hi0 - lo0);
}

auto fft(const ocean::field& field, int N, int M) {
    auto* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * M);
    auto* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * M);

    for (auto k = 0; k < N * M; ++k) {
        in[k][0] = real(field[k]);
        in[k][1] = imag(field[k]);
    }

    auto plan = fftw_plan_dft_2d(N, M, in, out, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_execute(plan);

    auto result = ocean::field(N * M);
    for (auto k = 0; k < N * M; ++k) {
        float r = out[k][0],
              i = out[k][1];
        result[k] = {r, i};

        // HACK, something something grid translation.
        // https://www.keithlantz.net/2011/11/ocean-simulation-part-two-using-the-fast-fourier-transform/
        auto sign = (k + k / N) % 2 == 0 ? 1.0f : -1.0f;
        result[k] *= sign;
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return result;
}

auto idx2k(int k, int N, int M, int Lx, int Ly) {
    auto row = k / N;
    auto col = k % N;
    auto n = remap_range<float>(row, 0, N - 1, -N / 2, N / 2);
    auto m = remap_range<float>(col, 0, M - 1, -M / 2, M / 2);
    return 2.0f * π * vec2{n / Lx, m / Ly};
}

auto gaussian() {
    static std::random_device rd{};
    static std::mt19937 gen{rd()};
    static std::normal_distribution dist{0.0f, 1.0f};
    return dist(gen);
}

} // ^ Unnamed namespace

float ocean::Ph(vec2 𝐤) const {
    auto 𝑘 = length(𝐤);
    if (𝑘 < 0.000001)
        return 0.0f;

    auto 𝑉 = length(𝑤);
    auto 𝐿 = 𝑉 * 𝑉 / 𝑔; // Largets possible wave

    auto ḱ = normalize(𝐤);
    auto ẃ = normalize(𝑤);

    // Modification in the paper ḱ·ẃ² -> ḱ·ẃ⁶
    auto ḱ·ẃ⁶ = pow(abs(dot(ḱ, ẃ)), 6.0f); // Note, dot product to power 6
    auto 𝑘⁴ = pow(𝑘, 4.0f);
    auto 𝑘𝐿² = pow(𝑘 * 𝐿, 2.0f); // Note, product to power of 2

    return 𝐴 * exp(-1.0f / 𝑘𝐿²) / 𝑘⁴ * ḱ·ẃ⁶;
}

complex<float> ocean::ℎₒ(vec2 𝐤) const {
    complex ξ{gaussian(), gaussian()};
    return recp_root_2 * ξ * sqrt(Ph(𝐤));
}

float ocean::ω(vec2 𝐤) const {
    return sqrt(𝑔 * length(𝐤));
}

ocean::field ocean::height_field(float t) const noexcept (true) {
    ocean::field 𝐇(K);
    for (auto k = 0; k < K; ++k) {
        auto 𝐤 = idx2k(k, N, M, Lx, Ly);
        auto ωt = ω(𝐤) * t;
        auto c₀ = cos( ωt) + 1.0if * sin( ωt);
        auto c₁ = cos(-ωt) + 1.0if * sin(-ωt);
        𝐇[k] = 𝐇₀₊[k] * c₀ + 𝐇₀₋[k] * c₁;
    }
    return fft(𝐇, N, M);
}

ocean::ocean(int N, int M, float Lx, float Ly, float 𝐴, vec2 𝑤)
    : N(N)
    , M(M)
    , K(N * M)
    , Lx(Lx)
    , Ly(Ly)
    , 𝐴(𝐴)
    , 𝑤(𝑤)
    , 𝐇₀₊(K)
    , 𝐇₀₋(K)
{
    for (auto k = 0; k < K; ++k) {
        auto 𝐤 = idx2k(k, N, M, Lx, Ly);
        𝐇₀₊[k] =      ℎₒ( 𝐤);
        𝐇₀₋[k] = conj(ℎₒ(-𝐤));
    }
}

using final_type = decltype(std::declval<ocean>().ℎₒ(vec2{}));
static_assert(std::is_same_v<complex<float>, final_type>,
              "The final return type is not of type complex<float>");
