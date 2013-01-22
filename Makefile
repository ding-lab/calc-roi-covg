all:
    ifndef SAMDIR
	@echo "Please define environment variable SAMDIR to point to your samtools libraries"
    else
	gcc -g -Wall -O2 -I${SAMDIR} calcRoiCovg.c -o calcRoiCovg -L${SAMDIR} -lbam -lm -lz
    endif
clean:
	rm -f calcRoiCovg
