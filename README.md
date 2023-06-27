# supersampler

## Fractionnal hitting set implementation for lightweight genomic data sketching

## Bird-eye view
SuperSampler (SPSP) is an implementation for a novel k-mer selection scheme we called Fractional Hitting Sets (FHS) which is a generalisation of Universal Hitting Sets (UHS).
It allows to quickly create sketches of genomes/ metagenomes and to compare such sketches to obtain Containment or Jaccard indices of the input data.

SuperSampler uses super-k-mers instead of k-mers which allows for lighter sketches, less RAM usage and less computational time when performing comparison than traditional subsampling methods. Thanks to a clever sketch organisation allowed by the super-k-mers structure.

Sketch creation is an application of FracMinHash on the selection of minimizers (a m-mer of a k-mer which hash value is minimal). When a minimizer is selected, every k-mer around it which shares the same minimizer is selected and will form a super-k-mer.

[Preprint](https://www.biorxiv.org/content/10.1101/2023.06.21.545875v1)

## Key properties
SuperSampler:
- allows to output k-mers from the sketches it created in plain text.
- Produces a majority of Maximal super-k-mers (super-k-mers which size is 2k - m).
- Gives equivalent results as State-of-the-art methods in less time and less memory consumption.


## Compilation
After the simple installation using:

```sh
git clone --depth 1  https://github.com/TimRouze/SuperSampler
make -j 8
```

It is possible to create sketches from fasta files and to compare such sketches with each other to perform analysis.

## Sketch creation
### Basic example
```sh
./sub_sampler -f my_fasta.fa.gz
```

will output a sketch file named sub_my_fasta.gz in the current directory. 
The sketch will be composed of super-31-mers with minimizers of size 11.
The number of k-mers will be divided by 1000 compared to the number present in the input fasta file.

### Main parameters
We can tune the sketch creation with several parameters impacting the size of output sketches.

#### Input file -i
The input fasta file.

#### Input File of File -f
If creating sketches for several fasta files at once, use this option and give the file-of-file for the inputs (txt format).
When using this option, SuperSampler will also output a file-of-file for the sketches it just created so they can be compared easily afterwards.

#### K-mer size -k
length of the k-mers to index.
The default value is 31.

#### Minimizer size -m
The minimizer of a kmer is its minimal word of length m  (according to a hash function).
SuperSampler will select k-mers according to the hash value of minimizers it encounters.
If the hash value is below the selction threshold, the k-mer is selected and a super-k-mer is constructer around this minimizer.
A smaller minimizer size improves the maximal lenght of super-k-mers which implies a gain in memory and computational time.
However, minimizer size shouldn't be too small (below 11 in general) as it makes comparison perform badly.
THe default value is 11.

#### Subsampling rate -s
The Subsampling rate will dictate how many k-mers will be kept in the created sketch.
a value of 10 means roughly one k-mer out of 10 will be selected in the output.
Higher values improve performances but too high value will make comparisons difficult because of the lack of k-mers.
The default value is 1000.

#### Threads used -t
This parameter defines the number of threads used by SuperSampler to construct its sketch.
The default value is 8 threads.

#### Prefix for output -p
The prefix to give to the output sketche(s). 
SuperSampler names sketches this way: "sub_"+input_filename+".gz".
In the case of using -f, the output file of file will also have this prefix.
The default value is "subsampled_"

## Sketch comparison
### Basic example
```sh
./comparator -f subsampled_file_of_file.txt
```
This command will output two csv files named "results_containment.csv.gz" and "results_jaccard.csv.gz" in the current directory. 
Basic behavior is an All vs All comparison of the sketches.
The csv files contain respectively the containment index and jaccard index matrices.

### Main parameters
We can tune the sketch comparison with several parameters depending on the user's needs.

#### Input file -f
The input sketch file of files.
If only using this parameter without -q, an All vs All comparison will be performed.

#### Input query file -q
The sketches given with this parameter will be compared to all the sketches given in the -f parameter.
This parameter allows to perform N vs All comparisons, N being the number of files given with this parameter.

#### Required precision to be output -p
The number of digits in output matrices.
Default value is 6.

#### Minimum value to be output -m
The minimum containment or jaccard value to output in the matrices.
THe default value is 0.0.

#### Prefix for output -o
The prefix to give to the output matrices. 
SuperSampler names matrices this way: "results"+input_filename+".csv.gz"
The default value is "results"

## Experiments:
Every commands, tools and data used to test this work are available on [this repo](https://github.com/TimRouze/Expe_SPSP)
