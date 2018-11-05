src = src/main.cpp

obj = $(src:.cpp=.o)

LDFLAGS = -lm
CXXFLAGS = -std=c++11 -O3

analogsat : $(obj)
	$(CXX) $(CXXFLAGS) -o analogsat $(obj) $(LDFLAGS)

.PHONY : clean
clean :
	rm -f analogsat $(obj) core *~
