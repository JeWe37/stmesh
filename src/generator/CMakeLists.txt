add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/src/generator/table.cpp
                   COMMAND python marching_hypercubes_table.py ${CMAKE_BINARY_DIR}/src/generator/table.cpp
                   DEPENDS marching_hypercubes_table.py
                   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/generator
                   COMMENT "Generating marching hypercubes table")

add_custom_target(RunGenerator DEPENDS ${CMAKE_BINARY_DIR}/src/generator/table.cpp
                  COMMENT "Checking if re-generation is required")
