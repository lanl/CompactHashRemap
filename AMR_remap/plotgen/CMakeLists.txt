set (PLOT_SRC distribution_plot_gen.cc)
set (PLOT_HDR distribution_plot_gen.h)

add_executable(plot_gen ${PLOT_SRC} ${PLOT_HDR})
target_link_libraries(plot_gen meshgen)
