# YP_Hotspots
This script is compatible with both Python 2 and Python 3, having been tested on Python 2.7+ and Python 3.6 & 3.7. 

## OS Requirements

This package is supported for macOS and Linux. The package has been tested on the following systems:

macOS: Monterey (12.5.1)
Linux: CentOS (7.5.1804)

## Python Dependencies
```
scipy  
numpy  
statsmodels  
multiprocessing  
```

## Installation Guide
```
pip install scipy  
pip install numpy  
pip install statsmodels  
pip install multiprocessing
```

## Usage
```
python SlidingWindow_PermutationTest_HotSpots_Identification.py [-h] -s SNP_MATRIX -l CHR_LEN -o OUTFILE
                                                                [-w WINDOW_LEN]
                                                                [-d MERGE_DISTANCE]  
                                                                [-r REPEAT_NUM] [-t THREADS]  
                                                                [-p PVALUE]
```

`
python3 SlidingWindow_PermutationTest_HotSpots_Identification.py dataset/3318S_SNP-INDEL.matrix
-l 4653728 -o d1k_r10k_w500.hotregion -d 1000 -r 10000 -w 500 -t 20
`

Running the complete test dataset will take 10 hours or more a "normal" desktop computer. If you just want to test whether it can run successfully, you can either use a quarter of the test dataset (while proportionally reducing the reference genome length of 4,653,728 bp) or reduce the "-r" parameter.

## License
This project is covered under the GUN General Public License, version 3.0 (GPL-3.0).
