#!/usr/bin/env python
import sys, os, time, gzip, bz2, pickle, json, logging, functools, hashlib
from datetime import datetime
from memory_profiler import memory_usage
from pandas.errors import EmptyDataError

# Read/Write
# ==========
# Get file object
def open_file_reader(filepath: str, compression="auto", binary=False):
    """
    Opens a file for reading with optional compression.

    Args:
        filepath (str): Path to the file.
        compression (str, optional): Type of compression {None, 'gzip', 'bz2'}. Defaults to "auto".
        binary (bool, optional): Whether to open the file in binary mode. Defaults to False.

    Returns:
        file object: A file-like object.
    """
    # Determine compression type based on the file extension if 'auto' is specified
    if compression == "auto":
        ext = filepath.split(".")[-1].lower()
        if ext == "gz":
            compression = "gzip"
        elif ext == "bz2":
            compression = "bz2"
        else:
            compression = None

    # Determine the mode based on the 'binary' flag
    mode = "rb" if binary else "rt"

    # Open the file with or without compression
    if not compression:
        return open(filepath, mode)
    elif compression == "gzip":
        return gzip.open(filepath, mode)
    elif compression == "bz2":
        return bz2.open(filepath, mode)
    else:
        raise ValueError(f"Unsupported compression type: {compression}")
            
# Get file object
def open_file_writer(filepath: str, compression="auto", binary=False):
    """
    Args:
        filepath (str): path/to/file
        compression (str, optional): {None, gzip, bz2}. Defaults to "auto".
        binary (bool, optional): Whether to open the file in binary mode. Defaults to False.
    
    Returns:
        file object
    """
    if compression == "auto":
        ext = filepath.split(".")[-1].lower()
        if ext == "gz":
            compression = "gzip"
        elif ext == "bz2":
            compression = "bz2"
        else:
            compression = None

    if binary:
        mode = "wb"
    else:
        mode = "wt"

    if not compression:
        return open(filepath, mode)
    elif compression == "gzip":
        return gzip.open(filepath, mode)
    elif compression == "bz2":
        return bz2.open(filepath, mode)
    else:
        raise ValueError(f"Unsupported compression type: {compression}")

# Pickle I/O
def read_pickle(filepath, compression="auto"):
    with open_file_reader(filepath, compression=compression, binary=True) as f:
        return pickle.load(f)
    
def write_pickle(obj, filepath, compression="auto"):
    with open_file_writer(filepath, compression=compression, binary=True) as f:
        pickle.dump(obj, f)
        
# Json I/O
def read_json(filepath):
    with open_file_reader(filepath, compression=None, binary=False) as f:
        return json.load(f)
    
def write_json(obj, filepath, indent=4):
    with open_file_writer(filepath, compression=None, binary=False) as f:
        return json.dump(obj, f)
    
# Formatting
# ==========
# Get duration
def format_duration(duration):
    """
    Format the elapsed time since `t0` in hours, minutes, and seconds.
    
    Adapted from @john-fouhy:
    https://stackoverflow.com/questions/538666/python-format-timedelta-to-string
    """
    hours, remainder = divmod(int(duration), 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02}:{minutes:02}:{seconds:02}"

# Format header for printing
def format_header(text, line_character="=", n=None):
    if n is None:
        n = len(text)
    line = n*line_character
    return "{}\n{}\n{}".format(line, text, line)

# Format memory
def format_bytes(B, unit="auto", return_units=True):
    """
    Return the given bytes as a human-readable string in KB, MB, GB, or TB.
    1 KB = 1024 Bytes

    Adapted from the following source (@whereisalext):
    https://stackoverflow.com/questions/12523586/python-format-size-application-converting-b-to-kb-mb-gb-tb/52379087
    """
    KB = 1024
    MB = KB ** 2  # 1,048,576
    GB = KB ** 3  # 1,073,741,824
    TB = KB ** 4  # 1,099,511,627,776

    def format_with_unit(size, unit_name):
        return f"{size:.2f} {unit_name}" if return_units else size

    unit = unit.lower()
    if unit != "auto":
        unit = unit.lower()
        if unit == "b":
            return format_with_unit(B, "B")
        elif unit == "kb":
            return format_with_unit(B / KB, "KB")
        elif unit == "mb":
            return format_with_unit(B / MB, "MB")
        elif unit == "gb":
            return format_with_unit(B / GB, "GB")
        elif unit == "tb":
            return format_with_unit(B / TB, "TB")
        else:
            raise ValueError(f"Unknown unit: {unit}")
    else:
        if B < KB:
            return format_with_unit(B, "B")
        elif KB <= B < MB:
            return format_with_unit(B / KB, "KB")
        elif MB <= B < GB:
            return format_with_unit(B / MB, "MB")
        elif GB <= B < TB:
            return format_with_unit(B / GB, "GB")
        else:
            return format_with_unit(B / TB, "TB")
        
# Logging
# =======
def build_logger(logger_name=__name__, stream=sys.stdout):
    # Create a logger object
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)  # Set the logging level
    
    # Create a stream handler to output logs to stdout
    stream_handler = logging.StreamHandler(stream)
    stream_handler.setLevel(logging.DEBUG)  # Set the level for the handler
    
    # Create a formatter and set it to the handler
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stream_handler.setFormatter(formatter)
    
    # Add the handler to the logger
    logger.addHandler(stream_handler)

    return logger
    
def reset_logger(logger):
    # Remove all existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()
    
    # Set a new handler (for example, to output to stdout)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    
    # Optionally set a new level
    logger.setLevel(logging.DEBUG)
    
# Timestamp
def get_timestamp(format_string:str="%Y-%m-%d %H:%M:%S"):
    # Get the current date and time
    now =  datetime.now()
    # Create a timestamp string
    return now.strftime(format_string)

# Profiling
# =========
def profile_peak_memory(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Measure memory usage
        mem_usage = memory_usage((func, args, kwargs), max_usage=True, retval=True, max_iterations=1)
        peak_memory, result = mem_usage[0], mem_usage[1]
        print(f"Peak memory usage for {func.__name__}: {format_bytes(peak_memory)}")
        return result
    return wrapper

# Directory
# =========
def get_file_size(filepath:str, format=False):
    size_in_bytes = os.stat(filepath).st_size
    if format:
        return format_bytes(size_in_bytes)
    else:
        return size_in_bytes
    
def check_file(filepath:str, empty_ok=False, minimum_filesize=1): # Doesn't handle empty gzipped files
    if not os.path.exists(filepath):
        raise FileNotFoundError(filepath)
    if not empty_ok:
        if get_file_size(filepath) < minimum_filesize:
            raise EmptyDataError(filepath)

# md5 hash from file
def get_md5hash_from_file(filepath:str, block_size=65536):
    """
    Calculate the MD5 hash of a file.

    Parameters:
    - file_path: The path to the file.
    - block_size: The size of each block read from the file (default is 64KB).

    Returns:
    - A string containing the MD5 hash.
    """
    md5 = hashlib.md5()
    with open(filepath, 'rb') as f:
        for block in iter(lambda: f.read(block_size), b''):
            md5.update(block)
    return md5.hexdigest()


# Check argument choices
def check_argument_choice(query, choices:set):
    """_summary_

    Args:
        query (_type_): Query option
        choices (set): Acceptable options

    Raises:
        ValueError: _description_
    """
    choices = set(choices)
    if query not in choices:
        raise ValueError(f"Invalid option '{query}'. Allowed choices are: {choices}")
