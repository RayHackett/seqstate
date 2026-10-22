CREATE TABLE Gene (
    id     INTEGER PRIMARY KEY ASC,
    contig TEXT NOT NULL,
    name   TEXT NOT NULL,
    start  INTEGER,
    end    INTEGER,
    strand INTEGER
);

CREATE TABLE Vocab (
    id   INTEGER PRIMARY KEY ASC,
    name TEXT NOT NULL
);

CREATE TABLE Domain  (
    gene    INTEGER NOT NULL,
    domain  INTEGER NOT NULL,
    pvalue  REAL,
    start   INTEGER,
    FOREIGN KEY (gene)   REFERENCES Gene(id)
    FOREIGN KEY (domain) REFERENCES Vocab(id)
);

CREATE TABLE Label (
    gene    INTEGER,
    label   INTEGER,
    FOREIGN KEY (gene)  REFERENCES Gene(id),
    FOREIGN KEY (label) REFERENCES LabelVocab(id),
    PRIMARY KEY (gene, label)
);

CREATE TABLE LabelVocab (
    id INTEGER PRIMARY KEY ASC,
    name TEXT NOT NULL
);