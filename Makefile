generator : generator.cpp
	g++ generator.cpp -lconfig++ -g -std=c++0x -o generator

reactor: reactor.cpp
	g++ reactor.cpp -std=c++0x -g -L/usr/lib -I /usr/include/ImageMagick/ -l Magick++ -l MagickCore -lconfig++ -o reactor

clean :
	rm generator reactor