# Fast Cohen's d calculator
![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/mrhedmad/fast-cohen/rust.yml)

A fast implementation of a Cohen's d calculator.

## Installation
You will need [cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) installed.
```bash
cargo install --git https://github.com/MrHedmad/fast-cohen.git
```

## Usage
You will need two `.csv` files, each with the two sample groups.
The first row of each input file should be the column names.
The first column of each file should contain gene names (or generically any item name).
The rest of the columns should be the data, and the order of these does not matter.
The data should be numeric.

Then, simply call:
```bash
fast-cohen "control_samples" "case_samples" "output_csv"
```

You can use `fast-cohen --help` to see the help message:
```
Calculate cohen's d of expression values.

Usage: fast-cohen [OPTIONS] <CASE_EXPRESSION_MATRIX> <CONTROL_EXPRESSION_MATRIX> <OUTPUT_PATH>

Arguments:
  <CASE_EXPRESSION_MATRIX>     Path to the expression matrix with the 'case' expression matrix
  <CONTROL_EXPRESSION_MATRIX>  Path to the expression matrix with the 'control' expression matrix
  <OUTPUT_PATH>                Path and filename of the output file

Options:
  -d, --delimiter <DELIMITER>  Delimiter of the input files [default: "\t"]
  -h, --help  
```

The output csv will have a `row_names` column with the row names, and a `cohen_d` column with the Cohen's d values.

NOTE: The order of the samples in the two input files MUST be the same between the two input files.
If you need to sort the rows, and the first column is the sample names, you can use [`xsv`](https://github.com/BurntSushi/xsv) to sort the columns with the following command:
```bash
xsv sort "some_file" > "sorted_file"
```
