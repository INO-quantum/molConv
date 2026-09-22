APP = molConv

# Add any other object files to this list below
APP_OBJS = MolecularConversion.o
APP_OBJS += threads.o

#CPPFLAGS += -I/usr/local/include        include directory
CPPFLAGS += -O3 -Wall
#LDFLAGS += -L/usr/local/lib             library directory
LDFLAGS += -pthread
#LDLIBS += -lusb-1.0                     #libary without 'lib'

all: clean build

build: $(APP)

.PHONY: clean

clean:
	rm -f *.o

debug: CXXFLAGS += -D_DEBUG -g
debug: clean
debug: $(APP)

$(APP): $(APP_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(APP_OBJS) $(LDLIBS)
