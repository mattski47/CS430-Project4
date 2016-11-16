all: raytrace.c
	gcc raytrace.c -o raytrace

clean:
	rm -rf raytrace *~