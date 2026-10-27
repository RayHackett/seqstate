[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-green.svg?maxAge=2678400&style=flat)](https://github.com/rayHackett/seqstate/blob/main/CHANGELOG.md)
[![Docs](https://img.shields.io/readthedocs/seqstate/latest?style=flat&maxAge=600)](https://seqstate.readthedocs.io)
[![Issues](https://img.shields.io/github/issues/RayHackett/seqstate.svg?style=flat&maxAge=600)](https://github.com/RayHackett/enseqstatezymm/issues)
<!-- [![version](https://img.shields.io/github/v/tag/rayhackett/seqstate?label=version&sort=semver)](https://github.com/rayhackett/seqstate/tags) -->
<!-- [![Actions](https://img.shields.io/github/actions/workflow/status/RayHackett/enzymm/test.yml?branch=main&style=flat&maxAge=300)](https://github.com/RayHackett/Enzymm/actions/workflows/test.yml)
[![Coverage](https://img.shields.io/codecov/c/gh/rayhackett/enzymm?logo=codecov&style=flat&maxAge=3600)](https://codecov.io/gh/rayhackett/enzymm/) -->
<!-- [![Python Versions](https://img.shields.io/pypi/pyversions/enzymm.svg?style=flat&maxAge=600&logo=python)](https://pypi.org/project/enzymm/#files)
[![PyPI](https://img.shields.io/pypi/v/enzymm.svg?style=flat&maxAge=3600)](https://pypi.python.org/pypi/enzymm)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/enzymm?logo=anaconda&style=flat&maxAge=3600)](https://anaconda.org/bioconda/enzymm)
[![Wheel](https://img.shields.io/pypi/wheel/enzymm?style=flat&maxAge=3600)](https://pypi.org/project/enzymm/#files) -->
<!-- [![PyPI Downloads](https://static.pepy.tech/personalized-badge/enzymm?period=total&units=INTERNATIONAL_SYSTEM&left_color=grey&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/enzymm) -->
<!-- [![AUR](https://img.shields.io/aur/version/python-pyhmmer?logo=archlinux&style=flat&maxAge=3600)](https://aur.archlinux.org/packages/python-pyhmmer) -->
<!-- [![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat)](https://github.com/RayHackett/enzymm/) -->

# seqstate [![Star me](https://img.shields.io/github/stars/rayhackett/seqstate.svg?style=social&label=Star&maxAge=3600)](https://github.com/rayhackett/seqstate/stargazers)

starting...

# Repo TODOs
- check over docstrings
- Actually write a README.md
- generate docs with sphinx
- add github workflows for tests
- publish to pypi


# Feature implementaton TODOs
- can I improve typing for sequencecollection.load ? not sure...
- Handle the sequence accession names better
- embl and ncbi input
- Decide which other information to optionally get from a genbank or even gff.
- CDS finder
    - non translated cds genes?
    - other RNA genes
    - pseudogenes?
    - locus tag?
    - gene name?
- Which other bonus info can we collect from pyrodigal?
    - RBS
- Gerneralize labels which dont have i_evalues necessarily
- Implement a general feature type (Protein would be a sub-part of that)
- Consider renaming the domain type to something generic which is the result of an HMM
- Consider adding Rfams
- Fix the preprocessing cli to handle various hmm sources via keyword better
- Link protein type with handling from peptides.py
- write tests for HMMCollection class
- Test speed maybe with lz4 compression of hmms. nominally, decompression rate should be better
- can we save space in genetable and feature table by index protein_id - joining would always have to be on seq_id + protein_id?

## Features for Gecco

## Features for Chamois

## Features for Sismis