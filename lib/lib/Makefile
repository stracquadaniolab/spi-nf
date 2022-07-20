all:spi_matrix_utils.o
	gcc -lm -shared spi_matrix_utils.o -o spi_matrix_utils.so

spi_matrix_utils.o:spi_matrix_utils.h spi_matrix_utils.c
	gcc -c -Wall -Werror -fpic spi_matrix_utils.c -o spi_matrix_utils.o

clean:
	rm *.o *.so
