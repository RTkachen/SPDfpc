#!/bin/bash
set -e  # Przerwij skrypt, jeśli jakiś krok się nie powiedzie

# Nazwa programu – musi być zgodna z EXECUTABLE_NAME w CMakeLists.txt
EXECUTABLE_NAME="fpcmax"

# Budowanie projektu
mkdir -p build
cd build
cmake ..
make

# Wróć do katalogu głównego po kompilacji
cd ..

if [ "$1" == "run" ]; then
    if [ -z "$2" ]; then
        echo "🔧 Użycie: ./build.sh run nazwa_pliku.txt"
        exit 1
    fi

    INPUT_FILE="$2"
    echo "🚀 Uruchamianie ./${EXECUTABLE_NAME} z plikiem '$INPUT_FILE'..."
    ./${EXECUTABLE_NAME} "$INPUT_FILE"
fi
