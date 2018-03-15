# Filterreads

This module is built to convert a fastq file into occupancy matrix of kmers per reads. The result of running the script will be matrix file with the following format_type

| read name | AA | AC | AG | AT | CA | CC | CG | CT | GA | GC | GG | GT | TA | TC | TG | TT |
|:---------:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| read1     |  0 |  1 | 1  |  1 |  0 | 1  | 0  | 0  |  0 | 0  | 1  | 0  | 1  |  0 |  0 | 0  |
| read2     | 1  | 1  | 0  | 0  | 0  | 0  | 0  | 1  |  1 | 1  | 1  | 0  | 0  |  1 |  0 | 1  |

The script can be run with the following command:

```{sh}
python3 filterreads.py -i <input fastq path> -o <output kmatrix file path> -k <klen>
```
