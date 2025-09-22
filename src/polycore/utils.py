import logging, psutil, numpy as np
from typing import List, Tuple, Dict, Optional

IUPAC_BITS = {
    'A':1, 'C':2, 'G':4, 'T':8,
    'R':1|4, 'Y':2|8, 'S':4|2, 'W':1|8, 'K':4|8, 'M':1|2,
    'B':2|4|8, 'D':1|4|8, 'H':1|2|8, 'V':1|2|4
}
ALLELES = [1,2,4,8]
POPCOUNT16 = np.array([bin(i).count("1") for i in range(16)], dtype=np.uint8)

def ambiguity_size(char: str) -> int:
    return bin(IUPAC_BITS.get(char.upper(), 0)).count('1')

def set_up_logging():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(message)s',
                        datefmt='%H:%M:%S')

def auto_chunk_size(n: int, safety_fraction: float = 0.8) -> int:
    avail_bytes = psutil.virtual_memory().available
    max_bytes = avail_bytes * safety_fraction
    bytes_per_site = n * n if n > 0 else max_bytes
    return max(1000, int(max_bytes / bytes_per_site))
