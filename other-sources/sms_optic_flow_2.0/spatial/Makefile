
CFLAGS= -O4
LFLAGS=iio.o -lpng -ljpeg -ltiff -lm

main: main.cpp brox_optic_flow.h mask.h bicubic_interpolation.h  zoom.h gaussian.h iio.o
	$(CXX) -std=c++98 $(CFLAGS) -o main main.cpp  $(LFLAGS)


iio.o: iio.c
	$(CC) -std=c99 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE -c iio.c

clean:
	rm -f iio.o main.o main
