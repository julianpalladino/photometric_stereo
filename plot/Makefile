all: plots

plots:
	mkdir -p plot
	gnuplot < plots.gpi
	cp -a plot/ ../doc/

clean:
	rm -rf plot/
	rm -f *.dat

.PHONY: all plots clean
