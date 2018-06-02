.PHONY: doc clean src

all: doc

doc: plot
	@$(MAKE) -C doc

plot: src test
	@$(MAKE) -C plot

src:
	@cmake src/
	@$(MAKE) -C src

test:
	@./run_tests

clean:
	@$(MAKE) clean -C doc
	@$(MAKE) clean -C src
