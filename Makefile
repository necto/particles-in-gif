
POTENTIAL=U3.cpp
INTEGRATOR=verle.cpp

all : generator reactor

simulation.gif: start.cfg reactor
	./reactor start.cfg

start.cfg : generator example.cfg
	./generator example.cfg

generator : generator.cpp common.h Makefile
	g++ generator.cpp -lconfig++ -g -std=c++0x -o generator

reactor: reactor.cpp *.h $(POTENTIAL) $(INTEGRATOR) Makefile
	g++ reactor.cpp $(POTENTIAL) $(INTEGRATOR) \
		  -std=c++0x -g -L/usr/lib \
          `Magick++-config --cxxflags --cppflags` \
		  `Magick++-config --ldflags --libs` -lconfig++ -o reactor

clean :
	rm generator reactor State/*.txt
