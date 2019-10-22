SRC="$(dirname $(realpath $0))/../src"

clang-format -i $SRC/*.cpp
clang-format -i $SRC/*.h
