EXE_DIR := bin/
SRC_DIR := src/
OBJ_DIR := obj/
PARSER_DIR:=../parser/src/

TARGET1=$(EXE_DIR)traj_hdf5.out

FC := gfortran
FXX := h5fc
FFLAGS := -ffree-line-length-none -fopenmp # gfortran
#FFLAGS = -O3 -convert big_endian #-ffpe-trap=invalid,zero,overflow -fbacktrace
#FFLAGS = -O3 -convert big_endian #-ffpe-trap=invalid,zero,overflow -fbacktrace
#FFLAGS = -h byteswapio,list=a -eZ -hpic -dynamic
FFLAGS += -J$(OBJ_DIR) 

#FFLAGS := -mcmodel=large -shared-intel -xHost -qopenmp
#FFLAGS += $(FFLAGS) #-O0 -CB -traceback -g -fpe0

#FFLAGS0= -fconvert=big-endian -ffree-line-length-none -fopenmp
#FFLAGS = $(FFLAGS0) -Wall -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace -g
# FFLAGS += -module $(OBJ_DIR) 

.SUFFIXES: .f90 .F90 .o .mod
%.o: %.mod

MOD_PARSER := input_parser.f90
SRC_PARSER :=

MOD:= input_parser.f90 module_tdata.f90 module_trajectory.f90 # ejecta_ye_time.f90

SRC := main.f90 #index_rank.f90

MOD := $(addprefix $(PARSER_DIR), $(MOD_PARSER)) $(addprefix $(SRC_DIR), $(MOD))
SRC := $(addprefix $(PARSER_DIR), $(SRC_PARSER)) $(addprefix $(SRC_DIR), $(SRC))

SRC := $(MOD) $(SRC)

OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC:.f90=.o))) 

.PHONY: all dirs test

all : dirs $(TARGET1)
objs : dirs $(OBJ_FILES)
dirs : $(EXE_DIR) $(OBJ_DIR)
$(EXE_DIR):
	mkdir -p $(EXE_DIR)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR) 


# 1st rule
$(TARGET1): $(OBJ_FILES)
	$(FXX) $(FFLAGS) -o $@ $(OMPFLAGS) $(LIBS) $^

$(OBJ_DIR)%.o: $(SRC_DIR)%.F90
	$(FXX) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $< -o $@
$(OBJ_DIR)%.o: $(PARSER_DIR)%.F90
	$(FXX) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $< -o $@
$(OBJ_DIR)%.o: $(SRC_DIR)%.f90
	$(FXX) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $< -o $@
$(OBJ_DIR)%.o: $(PARSER_DIR)%.f90
	$(FXX) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $< -o $@
$(OBJ_DIR)%.mod: $(SRC_DIR)%.F90
	$(FXX) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $< -o $@
$(OBJ_DIR)%.mod: $(PARSER_DIR)%.F90
	$(FXX) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $< -o $@
$(OBJ_DIR)%.mod: $(SRC_DIR)%.f90
	$(FXX) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $< -o $@
$(OBJ_DIR)%.mod: $(PARSER_DIR)%.f90
	$(FXX) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $< -o $@



# clean rule
.PHONY: clean
clean:
	$(RM) $(OBJ_FILES) fort* *.mod
