import os

def safe_fetch(file_path):
    """ Function to check if a file_path exists and return the contents
        if so.
    """
    if os.path.exists(file_path):
        with open(file_path, 'r', encoding='UTF-8') as file_contents:
            return file_contents.read()
    return ''