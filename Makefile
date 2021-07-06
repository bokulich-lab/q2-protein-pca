.PHONY: all lint test test-cov install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	py.test --cov=q2_protein_pca

install: all
	maturin build --release -m ranking/Cargo.toml
	$(PYTHON) setup.py install
	$(PYTHON) -m pip install --no-deps aln_ranking --find-links ranking/target/wheels/

dev: all
	pip install -e .

clean: distclean

distclean: ;
