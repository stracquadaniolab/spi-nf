# spi-nf: The SCRaMbLE Polymer Interaction model (SPI) workflow

![](https://img.shields.io/badge/current_version-0.1.0-blue)
![](https://github.com/stracquadaniolab/spi-nf/workflows/build/badge.svg)
## Overview

SPI is a statical mechanics model to conduct \acrshort{scramble} experiments in
silico. SPI uses polymer physics to model recombinations between loxpSym sites,
and implements efficient rejection sampling and histogram reweighting algorithms
to study synthetic genome evolution properties.

## Configuration

- `name`: mnemonic name of the run (default: spi)
- `outdir`: basedir for storing results (default: results/)
- `genome.structure`: genome structure file (default: syn9r_structure.txt, bundled)
- `genome.reconstruction`: genome reconstruction file (default: syn9r_reconstructions.txt, bundled)
- `model.lambda`: list of lambda parameters for the model (default: [4,5])
- `model.nu`: list of nu parameters for the model (default: [0.3,0.4])
- `model.configs`: maximum number of genome configurations (samples) to generate (default: 10)
## Running the workflow

### Install or update the workflow

```bash
nextflow pull stracquadaniolab/spi-nf
```

### Run the analysis

```bash
nextflow run stracquadaniolab/spi-nf
```

### Run the analysis with singularity on a slurm cluster

```bash
nextflow run stracquadaniolab/spi-nf -profile slurm,singularity
```

### Run the analysis with docker on a slurm cluster

```bash
nextflow run stracquadaniolab/spi-nf -profile docker,singularity
```

### Run the analysis with docker locally

```bash
nextflow run stracquadaniolab/spi-nf -profile docker
```
## Results

- `loglik.txt`: file with the loglik estimates for each parameter setting.
- `profiles.gz`: compressed file with the probability of deletion of each segment.
- `trajectories.gz`: compressed file with the recombination events for sampled genome.
- `run.info.txt`: log of the workflow parameters.

## Authors

- Giovanni Stracquadanio
