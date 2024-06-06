CC      = clang++
CFLAGS  = -O3 -std=c++11

SRC     = S_optim.cpp ../adjoint.cpp ../tracer.cpp ../inputs.cpp
OBJ     = $(SRC:.cpp=.o)
LIBS    = -lc 
OUT     = S_optim


$(OUT): $(OBJ) #  : this a rule that Out requires Obj
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

.PHONY: clean
clean:
	rm -rf ${OUT} $(OBJ)