all:
	gcc -g -Wall -O2 -I${SAMDIR} calcRoiCovg.c -o calcRoiCovg -L${SAMDIR} -lbam -lm -lz
clean:
	rm -f calcRoiCovg
