starting...

# Repo TODOs
- Fix the cli call in the pyproject.toml
- check over docstrings
- Actually write a README.md
- generate docs with sphinx
- add github workflows for tests
- publish to pypi


# Feature implementaton TODOs
- Handle the sequence accession names better
- Decide which other information to optionally get from a genbank.
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

## Features for Gecco

## Features for Chamois

## Features for Sismis