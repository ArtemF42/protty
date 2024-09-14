# Protty: A pretty simple tool for protease prediction
Protty utilizes profile Hidden Markov Models (HMMs) constructed from
[MEROPS](https://www.ebi.ac.uk/merops/) database to predict putative
proteases in query protein sequences

# Installation
```bash
git clone https://github.com/ArtemF42/protty.git
cd protty
pip install -e .
```

```bash
# Not available at the moment
pip install protty
```

# Building profile HMMs
## Running `protty-build`
**Note:** `protty-build` uses [Clustal Omega](http://www.clustal.org/omega/)
to perform multiple sequence alignment and requires it to be installed.
By default, Protty assumes `clustalo` is in your `PATH`. If this is not the
case, you should specify the `--clustalo` parameter

```bash
protty-build [options] $DATABASE
```

## Output
Once the process is complete, the profile HMMs will be available in
`$DATABASE/profiles`. The easiest way to merge them into the database is
to use `cat`

```bash
cat $DATABASE/profiles/*.hmm > $DATABASE/merops.hmm
```
## Advanced usage
The `protty-build` pipeline consists of 3 major steps:

1. Downloading MEROPS data
2. Filtering raw FASTA files
3. Building profile HMMs

Use `--skip` option if you want to manually control the pipeline. For example,
the command below will only download data from the MEROPS server

```bash
protty-build --skip 2,3 $DATABASE
```

# Scanning query sequences against the database
## Running `protty-scan`
```bash
protty-scan [options] $DATABASE/merops.hmm proteins.faa
```

## Output
By default, `protty-scan` generates two files named `predicted_proteases.tsv`
and `predicted_proteases.faa`, located in the working directory. Use `--tsv`
and `--faa` options to change the default behavior