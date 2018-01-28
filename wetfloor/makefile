CXX	=	g++
CCFLAGS	=   -Wall -O3 -I../engine
EXTRA	=	#-std=c++0x
OBJS	=	main.o
TARGET	=	wet

$(TARGET):	$(OBJS)
		$(CXX) $(CCFLAGS) $(EXTRA) -o $(TARGET) $(OBJS)

clean:
		rm -f $(OBJS) $(TARGET) *~ core

.cc.o:
		$(CXX) $(CCFLAGS) $(EXTRA) -c $<

