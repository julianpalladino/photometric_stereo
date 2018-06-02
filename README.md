# Photometric Stereo

## Documentation
* All documentation can be found in **Documentation - Informe.pdf**. This includes development and experimentation


## Dependencies

* pdflatex
* latexmk
* gcc

## Installation

```
sudo apt-get install -y \
  biber \
  build-essential \
  latexmk \
  texlive-bibtex-extra \
  texlive-lang-spanish \
  texlive-latex-recommended
```

## Compilation

* `make clean` deletes temporal files
* `make test` runs unit tests
