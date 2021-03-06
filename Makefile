# base library
TARGET = TiltFilterBase.swc

# wrapper
SRC = main.cpp tiltfilter.cpp
OBJ = $(foreach b,$(basename $(SRC)), $s.o)


# build options
CFLAGS = $(BASE_CFLAGS) -O3 -DNDEBUG -DWEBP_HAVE_PNG -DWEBP_HAVE_JPEG -DWEBP_HAVE_TIFF -DWEBP_USE_THREAD \
         -Wextra -Wshadow -Wall -flto-api=exports.txt

.PHONY: all
all: $(TARGET)

.PHONY: clean
clean:
	rm -f *.swf *.swc *.as3 *.abc *.o

$(TARGET): $(SRC)
	"$(FLASCC)/usr/bin/g++" $(CFLAGS) -O4 $(SRC) -emit-swc=info.smoche.TiltFilterBase -lFlash++ -lAS3++ -o $@

include Makefile.common

