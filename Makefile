CC      = clang++
CFLAGS  = -O3 -std=c++11

SRC     = src/optim.cpp src/adjoint.cpp src/tracer.cpp src/inputs.cpp
OBJ     = $(SRC:.cpp=.o)
LIBS    = -lc 
OUT     = optim


$(OUT): $(OBJ) #  : this a rule that Out requires Obj
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

.PHONY: clean
clean:
	rm -rf ${OUT} $(OBJ)