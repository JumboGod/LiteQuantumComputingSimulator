foreach(example bell_state ghz_state qft grover shor algorithms_demo)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE lqcs)
endforeach()
