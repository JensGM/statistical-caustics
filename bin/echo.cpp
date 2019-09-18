#include <iostream>

#include <ocean/ocean.hpp>


int main() {
    ocean 🌊{
        256,         // N
        256,         // M
        70.0f,       // Lx
        70.0f,       // Ly
        1.0f,        // 𝐴
        {6.0f, 5.0f} // ω
    };
    auto heights = 🌊.height_field(0.0);
    for (const auto& height : heights) {
        std::cout << real(height) << "\t";
    }
}
