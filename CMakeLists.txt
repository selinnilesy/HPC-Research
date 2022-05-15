project(Yunan_Implementation)

set(CMAKE_CXX_STANDARD 17)

set(HEADERS header.h)
set(SOURCE
        SSSconflictFree.cpp
        readSSS.cpp
        CSR_encoder.cpp
        )
add_executable(all SSSconflictFree.cpp readSSS.cpp CSR_encoder.cpp)