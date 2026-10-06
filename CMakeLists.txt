cmake_minimum_required(VERSION 3.20)
project(DeepNAV LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
enable_testing()

# ------------------------------
# 1. Configuration & Standards
# ------------------------------
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(APPLE)
    add_compile_options("-gdwarf-4")
endif()

# Generic radio targets share this file; numerical use requires only Eigen.
function(add_radiometric_targets)
    find_package(Eigen3 CONFIG REQUIRED)
    add_library(deepnav_radiometric STATIC
        src/observations/radiometric/TwoWayLink.cpp
        src/observations/radiometric/PropagationCorrections.cpp
        src/observations/radiometric/TwoWayMeasurement.cpp
        src/opnav/core/LightTimeSolver.cpp)
    target_include_directories(deepnav_radiometric PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(deepnav_radiometric PUBLIC Eigen3::Eigen)
    target_compile_features(deepnav_radiometric PUBLIC cxx_std_23)

    if(DEEPNAV_RADIOMETRIC_WITH_CSPICE)
        find_path(CSPICE_INCLUDE_DIR NAMES SpiceUsr.h
            HINTS "$ENV{CSPICE_HOME}/include" "$ENV{CSPICE_HOME}/cspice/include")
        find_library(CSPICE_LIB NAMES cspice cspice.a
            HINTS "$ENV{CSPICE_HOME}/lib" "$ENV{CSPICE_HOME}/cspice/lib")
        if(NOT CSPICE_INCLUDE_DIR OR NOT CSPICE_LIB)
            message(FATAL_ERROR "CSPICE adapter requested: set CSPICE_HOME or CSPICE_INCLUDE_DIR/CSPICE_LIB.")
        endif()
        add_library(deepnav_radiometric_spice STATIC
            src/dynamics/SpiceStateProvider.cpp
            src/dynamics/SpiceKernelSet.cpp
            src/observations/radiometric/SpiceLightTime.cpp)
        target_include_directories(deepnav_radiometric_spice PUBLIC ${CSPICE_INCLUDE_DIR})
        target_link_libraries(deepnav_radiometric_spice PUBLIC deepnav_radiometric ${CSPICE_LIB})
    endif()
endfunction()

# Numerical dynamics/filter targets are also available without CSPICE or graphics.
find_package(Eigen3 CONFIG REQUIRED)
set(CLOCK_SOURCES
    src/Clocks/Clocks.cpp
    src/Clocks/ClockHistory.cpp
    src/Clocks/ClockModel.cpp
    src/Clocks/LocalOscillator.cpp
    src/Clocks/DSAC.cpp
    src/Clocks/ClockTruthSimulator.cpp
    src/Clocks/Allan.cpp
    src/Clocks/Calibration.cpp
)
set(CLOCK_CSV_SOURCES src/utils/ClockCsvWriter.cpp)

add_library(deepnav_clocks STATIC ${CLOCK_SOURCES})
target_include_directories(deepnav_clocks PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(deepnav_clocks PUBLIC Eigen3::Eigen)

set(EKF_SOURCES
    src/dynamics/EphemerisInterpolator.cpp
    src/dynamics/IntegratedStateProvider.cpp
    src/RKF45Integrator.cpp
    src/filters/filter.cpp
    src/filters/EKF.cpp
    src/filters/CartesianPropagator.cpp
    src/observations/GeometricRadiometricModel.cpp
    src/perturbations/ForceModel.cpp
    src/perturbations/GravitationalModel.cpp
    src/perturbations/SRPModel.cpp
    src/perturbations/Eclipse.cpp
    src/perturbations/J2.cpp
)

# Reusable numerical filtering target without CSPICE or graphics linkage.
add_library(deepnav_ekf STATIC ${EKF_SOURCES})
target_include_directories(deepnav_ekf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(deepnav_ekf PUBLIC Eigen3::Eigen deepnav_clocks)


# Build the reusable radio model without configuring the full application's
# graphics/SPICE dependencies. Voyager's demo is not part of these targets.
option(DEEPNAV_RADIOMETRIC_ONLY "Build only reusable radio targets and their tests" OFF)
option(DEEPNAV_RADIOMETRIC_WITH_CSPICE "Build the optional CSPICE radio adapter" ON)
if(DEEPNAV_RADIOMETRIC_ONLY)
    add_radiometric_targets()
    add_executable(test_two_way_link tests/test_two_way_link.cpp src/dynamics/LinearStateProvider.cpp)
    target_link_libraries(test_two_way_link PRIVATE deepnav_radiometric)
    add_test(NAME test_two_way_link COMMAND test_two_way_link)
    add_executable(test_radiometric_navigation tests/test_radiometric_navigation.cpp
        src/dynamics/LinearStateProvider.cpp)
    target_link_libraries(test_radiometric_navigation PRIVATE deepnav_radiometric deepnav_ekf)
    add_test(NAME test_radiometric_navigation COMMAND test_radiometric_navigation)
    add_executable(test_ekf tests/test_ekf.cpp)
    target_link_libraries(test_ekf PRIVATE deepnav_ekf)
    add_test(NAME test_ekf COMMAND test_ekf)
    if(DEEPNAV_RADIOMETRIC_WITH_CSPICE)
        add_executable(cassini_two_way_demo cassini_two_way_demo.cpp)
        target_link_libraries(cassini_two_way_demo PRIVATE deepnav_radiometric_spice)
        target_compile_definitions(cassini_two_way_demo PRIVATE DEEPNAV_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
        add_test(NAME test_cassini_two_way_benchmark COMMAND cassini_two_way_demo --samples 3)
    endif()
    return()
endif()

# ------------------------------
# 2. Dependencies (Eigen, GLFW, GLEW, OpenGL)
# ------------------------------
find_package(Eigen3 CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(GLEW REQUIRED)
find_package(OpenGL REQUIRED)

message(STATUS "Graphics Libraries Found:")
message(STATUS "  GLFW: ${glfw3_DIR}")
message(STATUS "  GLEW: ${GLEW_INCLUDE_DIRS}")
message(STATUS "  OpenGL: ${OPENGL_LIBRARIES}")

# ------------------------------
# 3. Dependencies (CSPICE)
# ------------------------------
if(NOT DEFINED ENV{CSPICE_HOME})
    message(FATAL_ERROR "Errore: CSPICE_HOME non definita!")
endif()

set(CSPICE_DIR "$ENV{CSPICE_HOME}")

find_path(CSPICE_INCLUDE_DIR
    NAMES SpiceUsr.h
    PATHS "${CSPICE_DIR}/include" "${CSPICE_DIR}/cspice/include"
    NO_DEFAULT_PATH
)

find_library(CSPICE_LIB 
    NAMES cspice cspice.a
    PATHS "${CSPICE_DIR}/lib" "${CSPICE_DIR}/cspice/lib"
    NO_DEFAULT_PATH
)

# ------------------------------
# 4. ★ NEW: ImGui + ImPlot
# ------------------------------
set(IMGUI_DIR  ${CMAKE_SOURCE_DIR}/third_party/imgui)
set(IMPLOT_DIR ${CMAKE_SOURCE_DIR}/third_party/implot)

add_library(imgui_lib STATIC
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
    ${IMPLOT_DIR}/implot.cpp
    ${IMPLOT_DIR}/implot_items.cpp
)

target_include_directories(imgui_lib SYSTEM PUBLIC
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
    ${IMPLOT_DIR}
)

# ImGui deliberately uses raw memory operations for some internal types.
# Newer Clang versions diagnose this vendor implementation pattern, so keep
# that warning scoped to the third-party target rather than disabling it for
# DeepNAV sources.
target_compile_options(imgui_lib PRIVATE
    $<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang>:-Wno-nontrivial-memcall>
)

# ImGui needs GLFW + OpenGL — reuse what you already found
target_link_libraries(imgui_lib PUBLIC
    glfw
    ${OPENGL_LIBRARIES}
)

if(APPLE)
    target_link_libraries(imgui_lib PUBLIC
        "-framework Cocoa"
        "-framework IOKit"
        "-framework CoreVideo"
    )
endif()

# ------------------------------
# 5. Target: project_dependencies (Interface)
# ------------------------------
add_library(project_dependencies INTERFACE)

target_include_directories(project_dependencies INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    "${CSPICE_INCLUDE_DIR}"
    ${GLEW_INCLUDE_DIRS}
    "${CMAKE_SOURCE_DIR}/third_party/glm"
)

target_link_libraries(project_dependencies INTERFACE 
    Eigen3::Eigen 
    "${CSPICE_LIB}"
    glfw
    GLEW::GLEW
    ${OPENGL_LIBRARIES}
    imgui_lib          
    m
)

if(APPLE)
    target_link_libraries(project_dependencies INTERFACE 
        "-framework Cocoa"
        "-framework IOKit"
        "-framework CoreVideo"
    )
endif()

# ------------------------------
# 6. Static Library (Astro functions)
# ------------------------------
add_radiometric_targets()

set(BIG_FUNCTIONS_SOURCE_CANDIDATES
    ${CLOCK_SOURCES}
    ${CLOCK_CSV_SOURCES}
    src/stations/StationCatalog.cpp
    src/stations/ElevationMask.cpp
    src/RKF45Integrator.cpp
    src/dynamics/EphemerisInterpolator.cpp
    src/dynamics/SpiceInterpolator.cpp
    src/dynamics/ConstantStateProvider.cpp
    src/dynamics/SpiceEarthEnvironment.cpp
    src/simulation/EarthOrbitScenario.cpp
    src/dynamics/LinearStateProvider.cpp
    src/filters/filter.cpp
    src/filters/WLS.cpp
    src/filters/BatchLeastSquaresDriver.cpp
    src/filters/EKF.cpp
    src/filters/CartesianPropagator.cpp
    src/observations/GeometricRadiometricModel.cpp
    src/observations/synth/obs_synth.cpp
    src/observations/synth/TargetStateProvider.cpp
    src/observations/synth/RangeSynth.cpp
    src/observations/synth/RangeRateSynth.cpp
    src/observations/synth/VLBISynth.cpp
    src/perturbations/Gravitational.cpp
    src/perturbations/GravitationalModel.cpp
    src/perturbations/SRP.cpp
    src/perturbations/SRPModel.cpp
    src/perturbations/ForceModel.cpp

    src/opnav/CameraModel.cpp
    src/opnav/CameraAttitude.cpp
    src/opnav/DistortionModel.cpp
    src/opnav/core/GeometricLineOfSight.cpp
    src/opnav/core/LightTimeSolver.cpp
    src/opnav/core/ApparentDirection.cpp
    src/opnav/image/CircularGaussian.cpp
    src/opnav/image/CentroidEstimator.cpp
    src/opnav/image/PhotocenterCorrection.cpp
    src/opnav/image/CassiniImage.cpp
    src/opnav/image/EllipticalGaussian.cpp
    src/opnav/image/EllipticalCentroidEstimator.cpp
)

set(BIG_FUNCTIONS_SOURCES)
foreach(candidate_src IN LISTS BIG_FUNCTIONS_SOURCE_CANDIDATES)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${candidate_src}")
        list(APPEND BIG_FUNCTIONS_SOURCES "${candidate_src}")
    endif()
endforeach()

list(REMOVE_ITEM BIG_FUNCTIONS_SOURCES ${EKF_SOURCES} ${CLOCK_SOURCES}
    src/opnav/core/LightTimeSolver.cpp)
add_library(big_functions STATIC ${BIG_FUNCTIONS_SOURCES})
target_link_libraries(big_functions PUBLIC project_dependencies deepnav_ekf deepnav_radiometric)

# ------------------------------
# 7. Executables
# ------------------------------
set(APP_SOURCE_CANDIDATES
    Try3D.cpp
    main.cpp
    station_catalog_demo.cpp
    test_stations.cpp
    test_rkf45.cpp
    tests/test_ephemeris_interpolator.cpp
    tests/test_synth_observations.cpp
    tests/test_synth_source_comparison.cpp
    tests/test_voyager_position_od.cpp
    tests/test_optical.cpp
    tests/test_apparent_direction.cpp
    tests/test_centroid_estimator.cpp
    tests/test_cassini_image.cpp
    tests/test_elliptical_centroid.cpp
    cassini_centroid_demo.cpp
    tests/test_geometric_line_of_sight.cpp
    tests/test_clocks.cpp
    tests/test_clock_history.cpp
    tests/test_ekf.cpp
    tests/test_perturbation_models.cpp
    tests/test_spice_environment.cpp
    tests/test_earth_orbit_scenario.cpp
    ekf_demo.cpp
    cassini_two_way_demo.cpp
    cassini_two_way_ekf_demo.cpp
    tests/test_two_way_link.cpp
    tests/benchmark_two_way.cpp
    tests/test_radiometric_navigation.cpp
    tests/test_spice_radio_clock.cpp
    montecarlo/ekf_montecarlo.cpp
    tests/test_montecarlo.cpp

    tests/Optical/test_opnav_geometry.cpp
)

set(APP_SOURCES)
foreach(candidate_src IN LISTS APP_SOURCE_CANDIDATES)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${candidate_src}")
        list(APPEND APP_SOURCES "${candidate_src}")
    endif()
endforeach()

foreach(app_src IN LISTS APP_SOURCES)
    get_filename_component(app_name "${app_src}" NAME_WE)
    add_executable(${app_name} "${app_src}")

    if(app_name STREQUAL "cassini_two_way_demo" OR app_name STREQUAL "cassini_two_way_ekf_demo")
        if(NOT DEEPNAV_RADIOMETRIC_WITH_CSPICE)
            message(FATAL_ERROR "The full demo build requires DEEPNAV_RADIOMETRIC_WITH_CSPICE=ON.")
        endif()
        target_link_libraries(${app_name} PRIVATE deepnav_radiometric_spice deepnav_ekf)
        target_compile_definitions(${app_name} PRIVATE
            DEEPNAV_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    elseif(app_name STREQUAL "ekf_demo" OR app_name STREQUAL "ekf_montecarlo"
       OR app_name STREQUAL "test_spice_environment"
       OR app_name STREQUAL "test_earth_orbit_scenario")
        target_link_libraries(${app_name} PRIVATE big_functions)
        target_compile_definitions(${app_name} PRIVATE
            DEEPNAV_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    elseif(app_name STREQUAL "benchmark_two_way")
        target_sources(${app_name} PRIVATE src/dynamics/LinearStateProvider.cpp)
        target_link_libraries(${app_name} PRIVATE deepnav_radiometric_spice)
        target_compile_definitions(${app_name} PRIVATE DEEPNAV_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    elseif(app_name STREQUAL "test_spice_radio_clock")
        target_link_libraries(${app_name} PRIVATE deepnav_radiometric_spice)
        target_compile_definitions(${app_name} PRIVATE DEEPNAV_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    elseif(app_name STREQUAL "test_radiometric_navigation")
        target_sources(${app_name} PRIVATE src/dynamics/LinearStateProvider.cpp)
        target_link_libraries(${app_name} PRIVATE deepnav_radiometric deepnav_ekf)
    elseif(app_name STREQUAL "test_two_way_link")
        target_sources(${app_name} PRIVATE src/dynamics/LinearStateProvider.cpp)
        target_link_libraries(${app_name} PRIVATE deepnav_radiometric)
    elseif(app_name STREQUAL "test_montecarlo")
        target_include_directories(${app_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/montecarlo)
        target_link_libraries(${app_name} PRIVATE Eigen3::Eigen)
    elseif(app_name STREQUAL "test_ekf"
       OR app_name STREQUAL "test_perturbation_models")
        target_link_libraries(${app_name} PRIVATE deepnav_ekf)
    elseif(app_name STREQUAL "test_clocks" OR app_name STREQUAL "test_clock_history")
        target_sources(${app_name} PRIVATE ${CLOCK_CSV_SOURCES})
        target_include_directories(${app_name} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
        )
        target_link_libraries(${app_name} PRIVATE deepnav_clocks)
        if(app_name STREQUAL "test_clock_history")
            target_compile_definitions(${app_name} PRIVATE
                DEEPNAV_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
            )
        endif()
    elseif(app_name STREQUAL "test_apparent_direction")
        target_sources(${app_name} PRIVATE
            src/opnav/core/LightTimeSolver.cpp
            src/opnav/core/ApparentDirection.cpp
            src/opnav/CameraAttitude.cpp
            src/opnav/CameraModel.cpp
            src/opnav/DistortionModel.cpp
            src/dynamics/LinearStateProvider.cpp
        )
        target_include_directories(${app_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
        target_link_libraries(${app_name} PRIVATE Eigen3::Eigen)
    elseif(app_name STREQUAL "test_centroid_estimator")
        target_compile_definitions(${app_name} PRIVATE
            DEEPNAV_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
        )
        target_sources(${app_name} PRIVATE
            src/opnav/image/CircularGaussian.cpp
            src/opnav/image/CentroidEstimator.cpp
            src/opnav/image/PhotocenterCorrection.cpp
            src/opnav/core/LightTimeSolver.cpp
            src/opnav/core/ApparentDirection.cpp
            src/opnav/CameraAttitude.cpp
            src/opnav/CameraModel.cpp
            src/opnav/DistortionModel.cpp
            src/dynamics/LinearStateProvider.cpp
        )
        target_include_directories(${app_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
        target_link_libraries(${app_name} PRIVATE Eigen3::Eigen)
    elseif(app_name STREQUAL "cassini_centroid_demo" OR app_name STREQUAL "test_cassini_image"
           OR app_name STREQUAL "test_elliptical_centroid")
        target_compile_definitions(${app_name} PRIVATE DEEPNAV_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
        target_sources(${app_name} PRIVATE src/opnav/image/CassiniImage.cpp)
        if(app_name STREQUAL "cassini_centroid_demo" OR app_name STREQUAL "test_elliptical_centroid")
            target_sources(${app_name} PRIVATE
                src/opnav/image/CircularGaussian.cpp
                src/opnav/image/CentroidEstimator.cpp
                src/opnav/image/EllipticalGaussian.cpp
                src/opnav/image/EllipticalCentroidEstimator.cpp)
        endif()
        target_include_directories(${app_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
        target_link_libraries(${app_name} PRIVATE Eigen3::Eigen)
    elseif(app_name STREQUAL "test_geometric_line_of_sight")
        # Keep this geometry-only test independent of CSPICE and the graphics
        # stack. The implementation remains part of big_functions above.
        target_sources(${app_name} PRIVATE
            src/opnav/core/GeometricLineOfSight.cpp
        )
        target_include_directories(${app_name} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
        )
        target_link_libraries(${app_name} PRIVATE Eigen3::Eigen)
    else()
        target_link_libraries(${app_name} PRIVATE big_functions)
    endif()

    if(app_name STREQUAL "test_synth_observations"
       OR app_name STREQUAL "test_synth_source_comparison"
       OR app_name STREQUAL "test_voyager_position_od")
        add_test(NAME ${app_name} COMMAND ${app_name} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
        if(app_name STREQUAL "test_voyager_position_od")
            set_tests_properties(${app_name} PROPERTIES DEPENDS test_synth_observations)
        endif()
    elseif(app_name MATCHES "^test_")
        add_test(NAME ${app_name} COMMAND ${app_name} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
endforeach()

add_test(NAME test_montecarlo_runner COMMAND ekf_montecarlo --runs 2
    --output "${CMAKE_BINARY_DIR}/Output_montecarlo_smoke")

add_test(NAME test_cassini_two_way_benchmark COMMAND cassini_two_way_demo --samples 3)

add_test(NAME test_cassini_two_way_ekf_smoke COMMAND cassini_two_way_ekf_demo --hours 0.5 --truth matched
    --acceleration-noise-density 0
    --output "${CMAKE_BINARY_DIR}/Output_cassini_ekf_smoke")
