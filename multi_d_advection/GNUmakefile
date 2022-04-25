
# Program name

# find all the source files and header files
SOURCES := $(wildcard *.cpp)
HEADERS := $(wildcard *.H)

# list of object files by replacing .cpp with .o
OBJECTS := $(SOURCES:.cpp=.o)

# Multiple main functions
MAINS := test.o

# Filter out mains from objects
OBJECTS := $(filter-out $(MAINS), $(OBJECTS))

#Allow user to choose debug mode, default is false, turn off debug
DEBUG ?= FALSE

#CFLAGS
# -w: inhibit all warning messages
# -Werror: make all warnings into errors
# -Wpedantic:
# -Wall: enable all warnings about constructions that some users consider questionable.
# -Wextra: enables some extra warning flags not enabled by -Wall
# -Wshadow: warn whenever a local variable or declaration shows another variabe, parameter, type, class member, or instance variable.

CXXFLAGS := -Wall -Wextra -Wpedantic -Wshadow 

ifeq (${DEBUG}, TRUE)
	CXXFLAGS += -g
else 
# -O0 disbale all optimization: faster compilation but slow runtime
# -O0 -O1 -O2 -O3 -Os -Oz -Ofast -Omax 
	CXXFLAGS += -DNDEBUG -O3 #Turn on debug

endif

# Create object file with same name as .cpp file with all header files
%.o : %.cpp ${HEADERS}
	g++ ${CXXFLAGS} -c $<

# Needs to include test.o and timing.o since they're excluded
test: test.o ${OBJECTS}
	g++ -o $@ test.o ${OBJECTS}


#all: 

clean:
	rm -f *.o test *~ *.dat
