# Expe_SPSP
Every experiments made for [SuperSampler's paper](link to biorxive)

The performance comparison experiment can be reproduced by using the snakefile in the folder 'Performance comparison'.
In the meantime, here are the basic informations needed to reproduce our experiments:

## Tools used

- [Simka](https://github.com/GATB/simka)
- [Sourmash](https://github.com/sourmash-bio/sourmash)
- [SuperSampler](https://github.com/TimRouze/SuperSampler) Commit number for performance experiments is: [33548bc](https://github.com/TimRouze/supersampler/commit/33548bc332cc19e997f5e99a0f377f751d494d36). Commit number for abundance experiments is [187fabc](https://github.com/TimRouze/supersampler/commit/187fabca555934d16a99a301da49f9be597f3c7c).

## Data

- [Refseq](fof_refseq.txt)
- [Salmonellas](fof_salmonellas.txt)
- C. elegans base genome: GCA_000002985.3 WBcel235 bristol n2 haploid 07/02/2013.

Every genome used for experiments where taken from these sets. Always in the order of appearance in the files. e.g. When 100 salmonellas are used, the firs 100 genomes in the file of file for salmonellas are selected.

## Command lines
### Simka
```sh
./simka -in {input file of file} -out {folder for output} -out-tmp {folder for temporary files} -abundance-min 1 -kmer-size {k-mer size}
```
/!\ Simka requires a special formating for input files of file, see [Simka's repository](https://github.com/GATB/simka) for details /!\

### Sourmash
```sh
conda activate sourmash_env
sourmash sketch dna -p scaled={subsampling rate},k={k-mer size} {abund} --from-file {input file of file} -o {output name for sketch}
Sourmash compare {input sketch} {--containment} --csv {output filename} --ksize {k-mer size}
```
Sourmash results were sorted to match the input file of file order as Simka and SuperSampler keep this order.
SortCSV is present on [SuperSampler's repository](https://github.com/TimRouze/supersampler)
```sh
./sortCSV {input comparison matrix} {output name} {input file of file (to get original order)}
```

### SuperSampler
```sh
./sub_sampler -f {input file of file} -s {subsampling rate} -p {prefix for output sketches}_ -k {k-mer size} -m {minimizer size}
./comparator -f {input file of file} -o {prefix for output}
```
other options for ./sub_sampler, -j to avoid saving k-mer abundances. -l to save log(abundances) -g to save one abundance per super-k-mer. -gl for log and grouped at the same time.


## Values tested:

### Performance comparison

- K-mer size = [31, 63]
- Subsampling rate = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
- Minimizer size (specific to SuperSampler) = [11, 13, 15]
- 1024 RefSeq genomes and 1024 Salmonellas.


### Abundance experiment
As tools were compared to their own result with no subsampling applied, Simka was not run for this experiment.
- K-mer size = [31, 63]
- Subsampling rates = [2^{0} - 2^{15}]
- Minimizer sizes = [13, 15]
- 20 refseq and salmonellas genomes respectively.

### Scalability experiment
As only computational time and ram were monitored, we did not launch Simka on these experiments.
- K-mer size = 63
- Subsampling rate = 1000
- Minimizer size = 15
- From 100 to 128,000 RefSeq genomes.
