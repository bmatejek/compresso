CXX = g++-4.8
CPPFLAGS = -O4 -std=c++11
LDFLAGS = -O4
OBJS = compresso.o
INC_DIR = 
LIB_DIR = 
EXTERNAL_LIBS = 

compresso: $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) $(LIB_DIR) $(EXTERNAL_LIBS) -o compresso

compresso.o: compresso.cpp
	$(CXX) $(CPPFLAGS) $(INC_DIR) -c compresso.cpp

clean:
	rm -f *~ *.o \#* compresso