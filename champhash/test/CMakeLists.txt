cmake_minimum_required(VERSION 3.10)

add_executable(test_champhash test_champhash.cpp)
target_link_libraries(test_champhash PRIVATE champhash::champhash)
add_test(NAME champhash COMMAND test_champhash)
