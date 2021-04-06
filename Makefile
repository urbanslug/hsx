CXX=g++

CXXFLAGS=-Wall -g -D_GLIBCXX_DEBUG
ASAN_FLAGS=-fsanitize=address -fsanitize=undefined -fsanitize=leak -fno-omit-frame-pointer -Wno-format-security
SOURCE_FILES=src/cmd.cpp src/hsx.cpp

all:  hsx-cpp

hsx-cpp: $(SOURCR_FILES)
	$(CXX) $(CXXFLAGS) $(ASAN_FLAGS) -o hsx-cpp src/cmd.cpp

clean:
	rm -f *.o hsx-cpp
