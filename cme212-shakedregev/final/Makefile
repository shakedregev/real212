#EXEC = main
EXEC = test_legendre
####################
## Makefile Setup ##
####################

# Define the C++ compiler to use
CXX := g++ 
# Dependency directory and flags
DEPSDIR := $(shell mkdir -p .deps; echo .deps)
# MD: Dependency as side-effect of compilation
# MF: File for output
# MP: Include phony targets
DEPSFILE = $(DEPSDIR)/$(notdir $*.d)
DEPSFLAGS = -MD -MF $(DEPSFILE) #-MP


# Define cxx compile flags
CXXFLAGS  := -Wall -Wconversion  -Wextra  -pedantic  -std=c++11

####################
## Makefile Rules ##
####################

# Suffix replacement rules
#   $^: the name of the prereqs of the rule
#   $<: the name of the first prereq of the rule
#   $@: the name of the target of the rule

# 'make' - default rule
all: $(EXEC)

# Default rule for creating an exec of $(EXEC) from a .o file
$(EXEC): main.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LDLIBS)

# Default rule for creating a .o file from a .cpp file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEPSFLAGS) -o $@ -c $<

# 'make clean' - deletes all .o and temp files, exec, and dependency file
clean:
	-$(RM) *.o $(EXEC)
	$(RM) -r $(DEPSDIR)

# Define rules that do not actually generate the corresponding file
.PHONY: clean all

# Include the dependency files
-include $(wildcard $(DEPSDIR)/*.d)