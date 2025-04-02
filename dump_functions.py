"""Functions for creating nodes.dmp and names.dmp files in NCBI format."""

import numpy as np
import pandas as pd
from ete3 import NCBITaxa


def impute_taxons(taxon_df: pd.DataFrame) -> pd.DataFrame:
    """
    Replace nan values in taxon dataframe with
    unknown1, unknown2, ...
    """
    # Ensure we are working with a copy to avoid SettingWithCopyWarning
    taxon_df = taxon_df.copy()
    taxon_df = taxon_df.replace('unknown_rank_name', None)
    row, col = np.where(pd.isnull(taxon_df) == True)
    i = 1
    for r, c in zip(row, col):
        # Use .iloc for setting values
        taxon_df.iloc[r, c] = "unknown" + str(i)
        i += 1
    return taxon_df


def getTaxIDFromNCBI(taxa: str, ncbi_instance: NCBITaxa) -> str:
    """
    Return NCBI tax ID given the name of the taxon using a pre-initialized NCBITaxa instance.
    If the NCBI tax ID does not exist for that taxon
    then return the name of the taxon itself.
    """
    taxid = ncbi_instance.get_name_translator([taxa])
    if taxid:
        return str(list(taxid.values())[0][0])
    else:
        print(f"Warning: NCBI: Unable to match tax ID for '{taxa}'")
        return taxa


def create_names_dump(speciesSorted: pd.DataFrame, ncbi_instance: NCBITaxa, output_filepath: str) -> None:
    """
    Create a names.dmp file.
    """

    def write_line(level, name):
        tax_id = getTaxIDFromNCBI(name, ncbi_instance)
        scientific_name_entry = f"{level[0]}__{name}" if name != "root" else "root"
        return f"{tax_id}\t|\t{scientific_name_entry}\t|\t\t|\tscientific name\n"
    
    rank_order = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    lines = []

    # Add root entry
    lines.append("1\t|\troot\t|\t\t|\tscientific name\n")

    first_entry = speciesSorted.iloc[0]
    for level in rank_order:
        lines.append(write_line(level[0], first_entry[level]))

    previous = first_entry
    for _, row in speciesSorted.iloc[1:].iterrows():
        for level in rank_order:
            if row[level] != previous[level]:
                lines.append(write_line(level[0], row[level]))
        previous = row

    with open(output_filepath, "w") as outfile:
        outfile.writelines(lines)


def create_nodes_dump(speciesSorted: pd.DataFrame, ncbi_instance: NCBITaxa, output_filepath: str) -> None:
    """
    Create a nodes.dmp file.
    """

    def write_line(child, parent, rank):
         child_id = getTaxIDFromNCBI(child, ncbi_instance)
         parent_id = getTaxIDFromNCBI(parent, ncbi_instance) if parent != "root" else "1"
         # Ensure rank is lowercase as per convention seen in NCBI files
         formatted_rank = rank.lower() if rank else "no rank"
         return f"{child_id}\t|\t{parent_id}\t|\t{formatted_rank}\t|\n"

    # Define the standard NCBI rank order, ensure columns exist
    rank_order = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    available_ranks = [rank for rank in rank_order if rank in speciesSorted.columns]

    lines = []
    # Add root entry
    lines.append("1\t|\t1\t|\tno rank\t|\n")

    # Use a set to keep track of written parent-child relationships to avoid duplicates
    written_nodes = {("1", "1")}

    for _, row in speciesSorted.iterrows():
        parent = "root"
        for rank in available_ranks:
            child = row[rank]
            if pd.notna(child): # Ensure child is not NaN
                child_id = getTaxIDFromNCBI(child, ncbi_instance)
                parent_id = getTaxIDFromNCBI(parent, ncbi_instance) if parent != "root" else "1"

                node_pair = (child_id, parent_id)
                if node_pair not in written_nodes:
                    lines.append(write_line(child, parent, rank))
                    written_nodes.add(node_pair)
                parent = child # Current child becomes parent for the next level
            else:
                # If a rank is missing, stop processing lower ranks for this row
                break
        
    with open(output_filepath, "w") as outfile:
        outfile.writelines(lines)
