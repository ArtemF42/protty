# Protty: A pretty simple tool for protease prediction
Protty utilizes profile Hidden Markov Models (HMMs) constructed from [MEROPS](https://www.ebi.ac.uk/merops/) database to predict putative proteases in query protein sequences

## Installation
```bash
git clone https://github.com/ArtemF42/protty.git
pip install ./protty
```

## Scanning query sequences against the database
### Running `protty-scan`
```bash
protty-scan [options] /PATH/TO/DATABASE/merops.hmm proteins.faa
```

### Output
By default, `protty-scan` generates two files named `predicted_proteases.tsv` and `predicted_proteases.faa`, located in the working directory. Use `--tsv` and `--faa` options to change the default behavior
