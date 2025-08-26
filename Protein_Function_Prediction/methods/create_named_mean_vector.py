# Verison 1.0.0
# Author: Aaron Schroeder
import pandas as pd
import numpy as np

def process_data(csv_path, pkl_path, output_path):
    """
    Reads protein IDs and mean vector saved from centering the SVD matrix and combines them into one file

    Args:
        csv_path (str): path to input csv containing proteins IDs
        pkl_path (str): Path to pickle file with mean vector
        output_path (str): path to output csv
    """
    try:
        # read file
        df_ids = pd.read_csv(csv_path, usecols=[0])

        # read protein_ids
        protein_ids = df_ids.iloc[:, 0].tolist()
        print(f"Successfully read {len(protein_ids)} Protein-IDs from {csv_path} .")

        # read mean vector
        mean_vector = np.load(pkl_path, allow_pickle=True)
        print(f"Successfully read {len(mean_vector)} mean vector entries from {pkl_path}.")

        # check lengths
        if len(protein_ids) != len(mean_vector):
            print("Warning: ProteinIDs number does not match mean vector size; This should not happen")
            print(f"IDs: {len(protein_ids)}, Vector entries: {len(mean_vector)}")

        # create dataframe
        result_df = pd.DataFrame({
            'Protein_ID': protein_ids,
            'Mean_Vector_Value': mean_vector
        })
        print("Successfully combined dataframe.")

        # Schritt 4: Ergebnis als neue CSV-Datei speichern
        result_df.to_csv(output_path, index=False)
        print(f"result saved in {output_path}.")

    except FileNotFoundError as e:
        print(f"Error: At least one file not found {e}")
    except Exception as e:
        print(f"An unexpected error occured: {e}")


if __name__ == "__main__":
    # Use SVD result to ensure Protein_ID order is the same as in the mean vector
    input_csv_file = '../results/centered_svd_reduced_matrix.csv'
    input_pkl_file = '../results/mean_vector.pkl'
    output_csv_file = '../results/protein_mean_vector.csv'

    process_data(input_csv_file, input_pkl_file, output_csv_file)
