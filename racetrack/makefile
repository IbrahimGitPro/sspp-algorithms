CXX	=	g++-6
CCFLAGS	=	-Wall -O3 -I../engine

OBJS	=	main.o parsing.o
TARGET	=	race

$(TARGET):	$(OBJS)
		$(CXX) $(CCFLAGS) $(EXTRA) -o $(TARGET) $(OBJS)

clean:
		rm -f $(OBJS) $(TARGET) *~ core

.cc.o:
		$(CXX) $(CCFLAGS) $(EXTRA) -c $<

