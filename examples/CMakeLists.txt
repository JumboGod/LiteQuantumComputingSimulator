foreach(example bell_state ghz_state)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE lqcs)
endforeach()
