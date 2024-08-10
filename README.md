# Protty: A pretty simple tool for protease prediction
Protty utilizes profile Hidden Markov Models (HMMs) constructed from
[MEROPS](https://www.ebi.ac.uk/merops/) database to predict putative
proteases in query protein sequences

# Installation
```bash
pip install protty
```

# Building profile HMMs
## Running `protty-build`
**Note:** `protty-build` uses [Clustal Omega](http://www.clustal.org/omega/)
to perform multiple sequence alignment and requires it to be installed.
By default, Protty assumes `clustalo` is in your `PATH`. If this is not the
case, you should specify the `--clustalo` parameter

```bash
protty-build path/to/database
```

## Output

