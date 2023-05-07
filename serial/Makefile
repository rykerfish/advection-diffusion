default:
	gcc -g -o main main.c serial.c -lm 

run:
	./main 128

clean:
	rm main

test:
	gcc -g -o test test.c serial.c -lm
	./test
	rm ./test