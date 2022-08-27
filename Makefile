CPPFLAGS += -std=c++17 -W -Wall -g -Wno-unused-parameter
CPPFLAGS += -I include -I src

HPPFILES := $(shell find include/ -type f -name "*.h")
CPPFILES := $(shell find src/ -type f -name "*.cpp")
OBJS = $(patsubst %.cpp,%.o,$(CPPFILES))

src/%.o: src/%.cpp $(HPPFILES)
	g++ $(CPPFLAGS) -c -o $@ $<


bin/main : $(OBJS)
	mkdir -p bin
	g++ $(CPPFLAGS) $^ -o $@

.PHONY: clean

clean :
	rm -rf bin/*
	find src/ -type f -name '*.o' -delete
	find src/ -type f -name '*.o' -delete
	


