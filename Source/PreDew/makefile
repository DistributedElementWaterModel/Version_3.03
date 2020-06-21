OBJSD = DateTime.o Date.o CTime.o
OBJSP = sourcePreDew.o classPreDew.o
OBJSR = sourceRiverFlow.o classRiverFlow.o
C++C=g++

predew:	$(OBJSP)
	$(C++C) $(C++FLAGS) $(OBJSP) -o predew 

riverFlow:	$(OBJSR)
		$(C++C) $(C++FLAGS) $(OBJSR) -o riverFlow

stationMask:	sourceStationMask.cpp 
		$(C++C) $(C++FLAGS) sourceStationMask.cpp -o stationMask

clean:
	rm -f $(OBJS) *.o core
