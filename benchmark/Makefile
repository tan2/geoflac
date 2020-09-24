all: cmp

set:
	# run test, store test output to old/
	rm -f *.0 sys.msg output.asc
	OMP_NUM_THREADS=1 ../src/flac test-tiny.inp
	rm -rf orig/
	mkdir -p orig/
	mv *.0 *.rs sys.msg output.asc _contents.save vbc.s orig/

cmp:
	rm -f *.0 sys.msg output.asc
	OMP_NUM_THREADS=1 ../src/flac test-tiny.inp
	python3 compare.py 2

.PHONY: all set cmp
