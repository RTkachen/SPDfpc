cmake_minimum_required(VERSION 3.10)
project(SPD)

# 🔧 Nazwa programu (łatwa do zmiany w jednym miejscu)
set(EXECUTABLE_NAME fpcmax)

# 📂 Katalog wyjściowy dla binarki
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# 📌 Ustawienie standardu C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 📥 Dodanie katalogu z nagłówkami
include_directories(${CMAKE_SOURCE_DIR}/inc)

# 🔍 Znalezienie plików źródłowych i nagłówkowych
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB HEADERS "${CMAKE_SOURCE_DIR}/inc/*.hpp")

# 🛠️ Tworzenie pliku wykonywalnego
add_executable(${EXECUTABLE_NAME} ${SOURCES} ${HEADERS})

# 🔗 Tworzenie symlinka w katalogu głównym po buildzie
add_custom_command(
  TARGET ${EXECUTABLE_NAME} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXECUTABLE_NAME}
          ${CMAKE_SOURCE_DIR}/${EXECUTABLE_NAME}
)
