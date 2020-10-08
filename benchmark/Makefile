
#CASE := test-tiny.inp
CASE := test-tiny-long.inp
#CASE := test-big.inp
FRAME := 4
EXE := ../src/flac

# Use single cpu for omp to avoid round-off errors
RUN_TEST = rm -f *.0 sys.msg output.asc; OMP_NUM_THREADS=1 ${EXE} ${CASE}

ORIG = orig-${CASE}

all: cmp

set:
	# run test, store test output to old/
	${RUN_TEST}
	${MAKE} store

cmp:
	${RUN_TEST}
	python3 compare.py ${ORIG} ${FRAME}

store:
	rm -rf ${ORIG}
	mkdir -p ${ORIG}
	mv *.0 *.rs sys.msg output.asc _contents.save vbc.s ${ORIG}
	cp ../src/snapshot.diff ${ORIG}

.PHONY: all set cmp store
