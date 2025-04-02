# TaxdumpMaker

TaxdumpMaker is a Python-based tool designed to generate NCBI-like taxdump files (`names.dmp` and `nodes.dmp`) from a user-provided list of bacterial species names. It utilizes the `ete3` library to fetch taxonomy data from the NCBI database and processes it into the required format.

## Features

- Fetches taxonomy data for a list of species names.
- Generates `names.dmp` and `nodes.dmp` files in the NCBI taxdump format.
- Supports optional imputation of missing taxonomic ranks.
- Allows updating the local NCBI taxonomy database.

## Requirements

- Python 3.7 or higher
- Required Python packages:
  - `pandas`
  - `ete3`

## Installation

1. Clone the repository or download the script files.
2. Install the required Python packages:
   
   ```bash
   pip install pandas ete3
   ```

## Usage

Run the `taxdumpMaker.py` script with the following arguments:

```bash
python3 taxdumpMaker.py --species-list <species_list.txt> --output-dir <output_file> [--update-db] [--impute]
```

### Arguments:
- `<species_list.txt>`: Path to a text file containing a list of bacterial species (one per line).
- `<output_file>`: Name of the output directory where the generated `names.dmp` and `nodes.dmp` files will be saved.
- `--update-db` (optional): If specified, the script will update the local NCBI taxonomy database before processing.
- `--impute` (optional): If specified, the script will attempt to impute missing taxonomic ranks if there are any missing.

### Example:

```bash
python3 taxdumpMaker.py --species-list species_list_example.txt --output-dir test_output --update-db --impute
```

## Output

After successful execution, the output directory (`<output_file>`) will contain:
- `names.dmp`
- `nodes.dmp`

For more details, refer to the documentation or check the script's help message:

```bash
python3 taxdumpMaker.py --help
```

