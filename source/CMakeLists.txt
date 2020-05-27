#
# core
#
add_library(core INTERFACE)
target_compile_features(core INTERFACE cxx_std_17)
target_include_directories(core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}) # not cool, but okay for toy projects!
target_compile_options(core INTERFACE -Wall -Wextra -Wpedantic -fno-omit-frame-pointer -Wnarrowing)
target_sources(core INTERFACE
        core/types.hpp
        core/random.hpp
        core/matrix.hpp
        )

#
# analyse
#
add_library(analyse INTERFACE)
target_sources(analyse INTERFACE
        analyse/fourier.hpp
        )
target_link_libraries(analyse INTERFACE core)

#
# sampling
#
add_library(sampling INTERFACE)
target_sources(sampling INTERFACE
        sampling/sampler.hpp
        sampling/random.hpp
        sampling/stratified.hpp
        sampling/elementary.hpp

        sampling/imageio.hpp
        analyse/fourier.hpp

        sampling/integrator.hpp
        )

target_link_libraries(sampling INTERFACE core)