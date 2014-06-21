generator : generator.cpp
	g++ generator.cpp -o generator

reactor: reactor.cpp
	g++ reactor.cpp -L/usr/lib -I /usr/include/ImageMagick/ -l Magick++ -l MagickCore -o reactor