IDIR = /home/myname/ellipsoid_functions
CC = gcc #the C compiler
CFLAGS = -I. -Wall -O3 -lm
DEPS = ellipsoid_functions.h $(IDIR)

#Common source files
COMMON_SRC = zeros.c symmetric.c multiply.c \
             cholesky.c digitc.c max_abs_column.c \
	     display.c alpha_sort.c initial_values.c \
	     direct_calculation.c

COMMON_OBJ = $(patsubst %.c, $(IDIR)/%.o, $(COMMON_SRC))

#Main program 1 (Ellipsoid fitting using the separation in groups technique)
MAIN1_SRC = separation_in_groups.c matrix_summary.c
MAIN1_OBJ = $(patsubst %.c, $(IDIR)/%.o, $(MAIN1_SRC))
EXEC1 = separation_in_groups

#Main program 2 (Ellipsoid fitting using the sequential adjustments technique)
MAIN2_SRC = sequential_adjustments.c sequential.c
MAIN2_OBJ = $(patsubst %.c, $(IDIR)/%.o, $(MAIN2_SRC))
EXEC2 = sequential_adjustments

all : $(EXEC1) $(EXEC2) #all the executables in one target

#Rule to compile object files
$(IDIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXEC1): $(COMMON_OBJ) $(MAIN1_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

$(EXEC2): $(COMMON_OBJ) $(MAIN2_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(IDIR)/*.o $(EXEC1) $(EXEC2)

