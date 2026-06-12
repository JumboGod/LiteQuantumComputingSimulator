foreach(example bell_state ghz_state qft grover shor algorithms_demo vqe_h2 noisy_bell trajectory_noise)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE lqcs)
endforeach()
