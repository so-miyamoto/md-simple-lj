

CXX=clang++
CXXSTD=-std=c++14
OPT=-O3 -flto -Wall

BUILD=./build



all:
	mkdir -p $(BUILD)
	$(CXX) -c $(OPT) $(CXXSTD) simbox.cpp -o $(BUILD)/simbox.o
	$(CXX) -c $(OPT) $(CXXSTD) observer.cpp -o $(BUILD)/observer.o
	$(CXX) -c $(OPT) $(CXXSTD) systemIO.cpp -o $(BUILD)/systemIO.o
	$(CXX) -c $(OPT) $(CXXSTD) -DNVT_ENSEMBLE md.cpp -o $(BUILD)/md_nvt.o 
	$(CXX) $(BUILD)/simbox.o $(BUILD)/observer.o $(BUILD)/systemIO.o $(BUILD)/md_nvt.o -o $(BUILD)/nvt.out
	$(CXX) -c $(OPT) $(CXXSTD) -DNVE_ENSEMBLE md.cpp -o $(BUILD)/md_nve.o 
	$(CXX) $(BUILD)/simbox.o $(BUILD)/observer.o $(BUILD)/systemIO.o $(BUILD)/md_nve.o -o $(BUILD)/nve.out

clean:
	rm -rf $(BUILD)
