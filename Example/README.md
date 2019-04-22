# Demonstration

Users may follow this tutorial to see the performance and output of the function *findPhysLink* performs.

## Data set

The test data were derived from an *Escherichia coli* data set published by D. Ingle et al. and can be downloaded (DOI: [10.26188/5cbd1f18b2602](https://melbourne.figshare.com/ndownloader/files/14939444)) from a FigShare server maintained by the University of Melbourne.

Citation: Ingle, D. J. et al. [Evolution of atypical enteropathogenic *E. coli* by repeated acquisition of LEE pathogenicity island variants](https://www.nature.com/articles/nmicrobiol201510). *Nat. Microbiol*. 1, 15010 (2016).

## Procedure

When [SLURM Workload Manager](https://slurm.schedmd.com/documentation.html) is used, please specify arguments indicated by square brackets in the following SLURM scripts before running the scripts.

```bash
sbatch findPhysLink.slurm  # to produce a linkage network
sbatch lmm.slurm  # to produce an association network
```

Users can also launch local runs with R commands:

```bash
Rscript --vanilla findPhysLink.R  # to produce a linkage network
Rscript --vanilla lmm.R  # to produce an association network
```

Script process_results.R shows basic commands for processing GeneMates outputs.

## Results

Outputs from *findPhysLink* for two runs (with and without allelic physical distances) are archived and compressed into file Output.tar.gz (DOI: [10.26188/5cbd22cd00780](https://melbourne.figshare.com/ndownloader/files/14939450)), which can be downloaded from the FigShare server as well. Since the LMM-only run (without physical distances) shares most of its outputs with the full run (with physical distances), only the result file lmm.rds is included in Output.tar.gz.

Output: the directory for GeneMates outputs for the test data set

- lmm.rds: produced with lmm.slurm.
- findPhysLink: a directory containing complete results of findPhysLink.slurm.
  - physlink.rds: evidence of physical linkage.
  - snps.rds (processed SNP tables), which is the same as that produced by lmm.slurm.

