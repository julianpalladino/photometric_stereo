# Laboratorio de Métodos Numéricos

## Dependencias

* pdflatex
* latexmk
* gcc

## Instalación

```
sudo apt-get install -y \
  biber \
  build-essential \
  latexmk \
  texlive-bibtex-extra \
  texlive-lang-spanish \
  texlive-latex-recommended
```

## Uso

* `make doc` genera el informe
* `make clean` borra archivos temporales
* `make test` corre los tests de unidad
