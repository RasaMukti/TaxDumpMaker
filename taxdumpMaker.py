import argparse
import pandas as pd
from ete3 import NCBITaxa
import os
import sys
from typing import List, Dict, Optional

# Assuming dump_functions.py is in the same directory or accessible via PYTHONPATH
from dump_functions import impute_taxons, create_names_dump, create_nodes_dump

# Define standard rank order
RANK_ORDER = ["no rank", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

def read_species_list(filepath: str) -> List[str]:
    """Reads a list of species names from a file."""
    try:
        with open(filepath, "r") as f:
            speciesList = [line.strip() for line in f if line.strip()] # Ignore empty lines
        if not speciesList:
            print(f"Error: Species list file '{filepath}' is empty.")
            sys.exit(1)
        return speciesList
    except FileNotFoundError:
        print(f"Error: Species list file not found at '{filepath}'")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading species list file '{filepath}': {e}")
        sys.exit(1)

def get_taxonomy_data(speciesList: List[str], ncbi: NCBITaxa) -> Optional[pd.DataFrame]:
    """Fetches and processes taxonomy data for the given species list."""
    print("Fetching TaxIDs...")
    name_translator = ncbi.get_name_translator(speciesList)
    tax_ids_dict = dict(name_translator.items()) # Convert to dict for easier processing

    # Identify species not found
    found_species = set(tax_ids_dict.keys())
    not_found_species = [sp for sp in speciesList if sp not in found_species]
    if not_found_species:
        print("\nWarning: The following species were not found in the NCBI database:")
        for sp in not_found_species:
            print(f"- {sp}")
        print("-" * 20) # Separator

    if not tax_ids_dict:
        print("Error: No TaxIDs found for any species in the list.")
        return None

    tax_ids = [item[0] for item in tax_ids_dict.values()] # Extract the first ID if multiple exist

    print("Fetching lineages and ranks...")
    dflist = []
    for tax_id in tax_ids:
        try:
            lineage = ncbi.get_lineage(tax_id)
            ranks = ncbi.get_rank(lineage)
            # Translate taxids in ranks dict to names
            translated_ranks = {ncbi.get_taxid_translator([k]).get(k, 'unknown_rank_name'): v for k, v in ranks.items()}
            # Invert dict: {rank: name}
            rank_to_name = {v: k for k, v in translated_ranks.items()}
            dflist.append(rank_to_name)
        except Exception as e:
            species_name = ncbi.get_taxid_translator([tax_id]).get(tax_id, f"TaxID {tax_id}")
            print(f"Warning: Could not process lineage/ranks for {species_name}: {e}")

    if not dflist:
        print("Error: Failed to retrieve rank information for all found species.")
        return None

    df = pd.DataFrame(dflist)

    # Reorder columns based on standard rank order
    ordered_columns = [col for col in RANK_ORDER if col in df.columns]
    df = df[ordered_columns]

    return df

def sort_taxonomy_data(df: pd.DataFrame) -> pd.DataFrame:
    """Sorts the taxonomy dataframe based on the defined rank order."""
    sort_columns = [col for col in RANK_ORDER if col in df.columns]
    if not sort_columns:
        print("Warning: Cannot sort DataFrame as no standard rank columns are present.")
        return df
    print("Sorting taxonomy data...")
    # Use Categorical for robust sorting even with NaNs or mixed types if necessary
    # For now, standard sort_values should work if data is consistent
    try:
        sorted_df = df.sort_values(by=sort_columns)
    except Exception as e:
        print(f"Warning: Could not sort DataFrame by ranks: {e}. Returning unsorted.")
        return df
    return sorted_df

def main():
    parser = argparse.ArgumentParser(description="Generate NCBI-like taxdump files (names.dmp, nodes.dmp) from a species list.")
    parser.add_argument("--species-list", required=True, help="Path to the input file containing species names (one per line).")
    parser.add_argument("--output-dir", required=True, help="Path to the directory where output files (names.dmp, nodes.dmp) will be saved.")
    parser.add_argument("--update-db", action="store_true", help="Update the local NCBI taxonomy database before processing.")
    parser.add_argument("--impute", action="store_true", help="Impute missing taxonomic ranks with 'unknownX' placeholders.")

    args = parser.parse_args()

    # Validate output directory
    if not os.path.isdir(args.output_dir):
        try:
            os.makedirs(args.output_dir)
            print(f"Created output directory: {args.output_dir}")
        except OSError as e:
            print(f"Error: Could not create output directory '{args.output_dir}': {e}")
            sys.exit(1)

    # Initialize NCBITaxa
    print("Initializing NCBI Taxonomy database...")
    try:
        ncbi = NCBITaxa()
        if args.update_db:
            print("Updating NCBI Taxonomy database (this may take a while)...")
            ncbi.update_taxonomy_database()
    except Exception as e:
        print(f"Error initializing or updating NCBI Taxonomy database: {e}")
        sys.exit(1)

    species_list = read_species_list(args.species_list)
    taxonomy_df = get_taxonomy_data(species_list, ncbi)
    if taxonomy_df is None or taxonomy_df.empty:
        print("Exiting due to errors in fetching taxonomy data.")
        sys.exit(1)

    # Optional imputation
    if args.impute:
        print("Imputing missing taxonomic ranks...")
        taxonomy_df = impute_taxons(taxonomy_df) # Use the function from dump_functions

    # Sort the dataframe
    sorted_taxonomy_df = sort_taxonomy_data(taxonomy_df)

    # Define output file paths
    names_dump_path = os.path.join(args.output_dir, "names.dmp")
    nodes_dump_path = os.path.join(args.output_dir, "nodes.dmp")

    # Create taxdump files
    print(f"Creating names dump file: {names_dump_path}")
    create_names_dump(sorted_taxonomy_df, ncbi, names_dump_path)

    print(f"Creating nodes dump file: {nodes_dump_path}")
    create_nodes_dump(sorted_taxonomy_df, ncbi, nodes_dump_path)

    print("\nTaxdump files generated successfully.")
    print(f"Output directory: {args.output_dir}")

if __name__ == "__main__":
    main()