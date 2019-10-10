# Using Barcode_split

`Barcode_split` is a simple tool to split summary sequencing files according to barcodes.

## User interface

`Barcode_split` was designed to be used either through a python Application programming interface (API) for Jupyter notebook or a command line interface (CLI).

### Jupyter API

* [Barcode_split API usage notebook](https://a-slide.github.io/pycoQC/demo/Barcode_split_API_demo/)

### Shell CLI

* [Barcode_split CLI usage notebook](https://a-slide.github.io/pycoQC/demo/Barcode_split_CLI_demo/)


## IO and options

`Barcode_split` requires at least one sequencing_summary.txt file generated with either Albacore or Guppy, but it can also take multiple files similar to pycoQC

If barcodes information are not contained in the summary file it is also possible to pass externally provided barcodes from either earlier versions Guppy and Deepbinner.

Output summary files split by barcodes are created in a user defined directory (default = current dir) and prefixed according to the barcode name.

By default `Barcode_split` doesn't generate files for unclassified and low frequency barcodes, which are likely to be false positive. Users can control this behaviour be requestion unclassified reads (`output_unclassified`),
and changing the low frequency threshold (`min_barcode_percent`).

At the end of execution `Barcode_split` print a summary of the barcodes found.

```
                  Counts  Write
    barcode02          2  False
    barcode07          1  False
    barcode08         30  False
    barcode09       9945   True
    barcode10      12644   True
    barcode11      13594   True
    barcode12       9813   True
    unclassified    3971  False  
```
