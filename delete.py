import os

def delete_small_files(folder_path):
    try:
        # List all files in the specified folder
        files = os.listdir(folder_path)

        for file in files:
            file_path = os.path.join(folder_path, file)

            # Check if the file is a regular file and its size is less than 30KB
            if os.path.isfile(file_path) and os.path.getsize(file_path) < 30 * 1024:
                print(f"Deleting {file_path}...")
                os.remove(file_path)
                print(f"{file_path} deleted successfully.")

        print("Deletion of small files completed.")

    except Exception as e:
        print(f"An error occurred: {str(e)}")

# Specify the folder path
folder_path = "D:\Jatkasaari\original_data\preprocessed_lidar\T\T52"

# Call the function to delete small files
delete_small_files(folder_path)
