all:
	gcc sds.c hashmap.c newpde2d_paper.c -g -Wall -o newpde -lm -fopenmp
	time ./newpde newentryfile1.txt
	python plot_graphs2d.py teste6
