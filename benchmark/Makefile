
CASE := test-tiny.inp
CASE := test-tiny-long.inp
FRAME := 4

# Use single cpu for omp to avoid round-off errors
RUN_TEST = rm -f *.0 sys.msg output.asc; OMP_NUM_THREADS=1 ../src/flac ${CASE}

ORIG = orig-${CASE}

all: cmp

set:
	# run test, store test output to old/
	${RUN_TEST}
	rm -rf ${ORIG}
	mkdir -p ${ORIG}
	mv *.0 *.rs sys.msg output.asc _contents.save vbc.s ${ORIG}
	git log -1 HEAD > ${ORIG}/version

cmp:
	${RUN_TEST}
	python3 compare.py ${FRAME} ${ORIG}

.PHONY: all set cmp
