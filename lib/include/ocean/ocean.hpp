/* Implementation of the ocean water simulation described by Jerry Tessendorf in
 * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.161.9102&rep=rep1&type=pdf
 *
 * The idea of this implementation is to keep the code as close to the paper as
 * possible.
 *
 * Warning, unicode in source code 🌊
 */

#include <complex>
#include <vector>

#include <glm/glm.hpp>

class ocean {
public:
    using field = std::vector<std::complex<float>>;
    ocean(int N, int M, float Lx, float Ly, float 𝐴, glm::vec2 𝑤);
    field height_field(float t) const noexcept (true);

private:
    /* Ocean properties */

    int N, M; // grid resolution
    int K;
    float Lx, Ly;  // field dimension (i think in meters)
    float 𝐴; // arbitrary amplitude
    glm::vec2 𝑤; // horizontal wind vector

    /* Spectra */

    field 𝐇₀₊; // fourier amplitudes
    field 𝐇₀₋; // complex conjugate

public:
    /**
     * Phillips frequency for 𝐤
     */
    float Ph(glm::vec2 𝐤) const;

    /**
     * Fourier amplitude, ℎₒ with a tilde in the paper
     */
    std::complex<float> ℎₒ(glm::vec2 𝐤) const;

    /**
     * Dispersion relation
     */
    float ω(glm::vec2 𝐤) const;
};
