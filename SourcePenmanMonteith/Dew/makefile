# Define the compiler
CXX = g++

# Define the source files, object files, and executable name
SOURCES = CTime.cpp \
	Date.cpp \
	DateTime.cpp \
	Dew.cpp \
	DistributedElement.cpp \
	Glacier.cpp \
	GlacierElements.cpp \
	GlacierIce.cpp \
	HbvAquifer.cpp \
	HBV.cpp \
	InputElement.cpp \
	InputTimeSeries.cpp \
	KinematicWave.cpp \
	KwaAquifer.cpp \
	Lake.cpp \
	LakeWaterBalance.cpp \
	MeteorologicalStations.cpp \
	ModelControl.cpp \
	Parameters.cpp \
	ParametersGeneral.cpp \
	ParametersGlacierRetreat.cpp \
	ParametersKiWa.cpp \
	ParametersLandSurface.cpp \
	ParametersSubSurfaceHbv.cpp \
	Routing.cpp \
	SelectedHbvTimeSeriesElements.cpp \
	SelectedKiWaTimeSeriesElements.cpp \
	SelectedSubCatchmentTimeSeriesElements.cpp \
	Snow.cpp \
	stdafx.cpp \
	SubCatchment.cpp \
	TotalReservoirStorage.cpp \
	Util.cpp \
	Vegetation.cpp

OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = dew

CXXFLAGS = $(shell pkg-config --cflags netcdf)
LDFLAGS = $(shell pkg-config --libs netcdf)

# Default target to build the project
all: $(EXECUTABLE)

# Rule to link the object files into the executable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

# Generic rule to compile a .cpp file into a .o file (using implicit rules)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to clean up generated files
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
