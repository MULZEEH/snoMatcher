import os
import sys


def rename_csv_files(directory):
    """
    Renames CSV files in the format 'a@b.csv' to 'a.csv'

    Args:
        directory: Path to the directory containing CSV files
    """
    # Check if directory exists
    if not os.path.exists(directory):
        print(f"Error: Directory '{directory}' does not exist")
        return

    if not os.path.isdir(directory):
        print(f"Error: '{directory}' is not a directory")
        return

    # Get all CSV files in the directory
    files = [f for f in os.listdir(directory) if f.endswith('.csv')]

    if not files:
        print(f"No CSV files found in '{directory}'")
        return

    renamed_count = 0

    for filename in files:
        # Check if filename contains '@'
        if '@' in filename:
            # Split by '@' and take the part before it
            new_name = filename.split('@')[0] + '.csv'

            old_path = os.path.join(directory, filename)
            new_path = os.path.join(directory, new_name)

            # Check if target file already exists
            if os.path.exists(new_path):
                print(f"Warning: '{new_name}' already exists. Skipping '{filename}'")
                continue

            try:
                os.rename(old_path, new_path)
                print(f"Renamed: '{filename}' -> '{new_name}'")
                renamed_count += 1
            except Exception as e:
                print(f"Error renaming '{filename}': {e}")
        else:
            print(f"Skipping '{filename}' (no '@' found)")

    print(f"\nTotal files renamed: {renamed_count}")


if __name__ == "__main__":
    # Get directory from command line argument or prompt user
    if len(sys.argv) > 1:
        directory = sys.argv[1]
    else:
        directory = input("Enter the directory path: ").strip()

    rename_csv_files(directory)