starting...

# Repo TODOs
- Fix the cli call in the pyproject.toml
- check over docstrings
- Actually write a README.md
- generate docs with sphinx
- add github workflows for tests
- publish to pypi


# Feature implementaton TODOs
- can I improve typing for sequencecollection.load ? not sure...
- test faa only
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