# Protty: A pretty tool for protease prediction
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
pip install protty
```

# ...
## Running `protty-build`
**Note:** `protty-build` uses [Clustal Omega](http://www.clustal.org/omega/) to perform multiple sequence alignment and requires it to be installed. By default, Protty assumes `clustalo` is in your `PATH`. If this is not the case, you should specify the `--clustalo` parameter

```bash
protty-build
```

## Output

# ...
## Running `protty-scan`
```bash
protty-scan
```

## Output