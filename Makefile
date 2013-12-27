all: clean
	g++ -I lib source/*.cpp

clean:
	rm -f ./a.out
