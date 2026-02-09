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

stationMaskAll:	sourceStationMaskAll.cpp 
		$(C++C) $(C++FLAGS) sourceStationMaskAll.cpp -o stationMaskAll

stationMaskElevation:	sourceStationMaskElevation.cpp 
			$(C++C) $(C++FLAGS) sourceStationMaskElevation.cpp -o stationMaskElevation

mask:	sourceMask.cpp 
		$(C++C) $(C++FLAGS) sourceMask.cpp -o mask

gridMaskNorge:	$(OBJSD) sourceGridMaskNorge.cpp
		$(C++C) $(C++FLAGS) $(OBJSD) sourceGridMaskNorge.cpp -o gridMaskNorge 

gridMaskVariable:	$(OBJSD) sourceGridMaskVariable.cpp
			$(C++C) $(C++FLAGS) $(OBJSD) sourceGridMaskVariable.cpp -o gridMaskVariable 

gridRectangleVariable:	$(OBJSD) sourceGridRectangleVariable.cpp
			$(C++C) $(C++FLAGS) $(OBJSD) sourceGridRectangleVariable.cpp -o gridRectangleVariable 

clean:
	rm -f $(OBJS) *.o core
