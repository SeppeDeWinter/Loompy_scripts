# [Loom_scripts]

**Command line interface (CLI) for directly manipulating loom files**

## Requirements

* Python 3

* Loompy

* Numpy

* Pandas

## Installation

```bash
git clone https://github.com/SeppeDeWinter/Loompy_scripts.git
```

## Usage

**Calculate average gene expression over defined set(s) of column attributes.**

```bash
python Loom_scripts.py average_expression \
                       -i <INPUT> \
                       -o <OUTPUT> \
                       -ca <COLUMN_ATTRIBUTES \
                       -r <RA_USE>
                       [COLUMN_ATTRIBUTES ...]> \
                       [-ra <ROW_ATTRIBUTE>] \
                       [-g <GENES>] \
                       [-b <BATCH_SIZE>] \
                       [-ly ,LAYERS [LAYERS ...]>]
                       [-h]

Arguments:
    -i <INPUT>, --input <INPUT>
                        <Required> Path to loom file.
    -o <OUTPUT>, --output <OUTPUT>
                        <Required> Path to write output csv file containing average gene expression.
    -ca <COLUMN_ATTRIBUTES [COLUMN_ATTRIBUTES ...]>, --column_attributes <COLUMN_ATTRIBUTES [COLUMN_ATTRIBUTES ...]>
                        <Required> Column attribute or combination of column attributes to calculate average gene expression on. 
                        For example -ca ClusterID to calculate average gene expression across all clusters.
                        For example -ca ClusterID TimePoint to calculate average gene expression across all combinations of ClusterIDs and timepoints (e.g. cluster_1 timepoint_1, cluster_1 timepoint_2, cluster_2 timepoint_1 and cluster_2 timepoint_2 when the loom files has two ClusterIDs and two TimePoints).
    -r <RA_USE>, --ra_use <RA_USE>
                        <Required> which row attribute is used to calculate average exoression on (e.g. gene).
    -ra <ROW_ATTRIBUTE>, --row_attribute <ROW_ATTRIBUTE>
                        <Optional> Boolean row attribute (in loom file) specifying on which genes average expression should be calculated (e.g> selected), leave empty for all genes.
    -g <GENES>, --genes <GENES>
                        <Optional> new line seperated file containing gene names from which the average expression has to be calculated, leave empty for all genes
    -b <BATCH_SIZE>, --batch_size <BATCH_SIZE>
                        <Optional> batch size for scanning through loom file, default is 512
    -ly <LAYERS [LAYERS ...]>, --layers <LAYERS [LAYERS ...]>
                        <Optional> Layers to calculate average expression on (e.g. spliced)
    -h, --help            show help message.
                                 
```
**subsample loom file down to a random sample of cells (without replacement)**

```bash
python Loom_scripts.py subsample_cells 
                       -i <INPUT> \
                       -o <OUTPUT> \
                       -n <NUM_CELLS>

Arguments:
    -i <INPUT>, --input <INPUT>
                        <Required> Path to loom file, make sure loom file is writable
    -o <OUTPUT>, --output <OUTPUT>
                        <Required> Path to write output loom file
    -n <NUM_CELLS>, --num_cells <NUM_CELLS>
                        <Required> Number of cells to sample
    -h, --help            show this help message
```
**Map row attribute to different name (e.g. map gene names to orthologous gene names).**

```bash

python Loom_scripts.py map_ra 
                        -i <INPUT>
                        -ra <ROW_ATTR>
                        -m <MAPPING>
                        -n <NEW_RA_NAME>
                        [-h]   

Arguments:
    -i INPUT, --input INPUT
                        <Required> Path to loom file (Loom file shoud be
                        writable)
    -ra ROW_ATTR, --row_attr ROW_ATTR
                        <Required> Row attribute to map from (e.g. Gene)
    -m MAPPING, --mapping MAPPING
                        <Required> Path to csv file containing mapping (A 2-column data frame defining variable name mapping. First column is source variable name and second column is target variable name
    -n NEW_RA_NAME, --new_ra_name NEW_RA_NAME
                        <Required> New name of row attribute for the new mapping.

  -h, --help            show this help message.
```



