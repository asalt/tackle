# experimental. not in use yet
import os
import pandas as pd
import threading
from typing import Dict, Any, Type, Optional

class ThreadSafeLazyLoader:
    def __init__(self, file_path: str, read_kwargs: Dict[str, Any] = None):
        self._file_path = file_path
        self._read_kwargs = read_kwargs or {}
        self._df = None
        self._lock = threading.Lock()

    def _read_file(self) -> Optional[pd.DataFrame]:
        """
        Read file with proper error handling and logging
        """
        if not os.path.exists(self._file_path):
            print(f"Error: Could not find {self._file_path}")
            return None

        try:
            # Determine reader based on file extension
            if self._file_path.endswith(".tsv"):
                reader = pd.read_table
            elif self._file_path.endswith(".xlsx"):
                reader = pd.read_excel
            else:
                print(f"Unsupported file format: {self._file_path}")
                return None

            # Read file with provided kwargs
            return reader(self._file_path, **self._read_kwargs)
        
        except Exception as e:
            print(f"Error reading {self._file_path}: {e}")
            return None

    @property
    def df(self) -> Optional[pd.DataFrame]:
        """
        Thread-safe lazy loading of dataframe
        """
        if self._df is None:
            with self._lock:
                # Double-checked locking to prevent multiple reads
                if self._df is None:
                    self._df = self._read_file()
        return self._df


class SingletonMeta(type):
    """
    Thread-safe Singleton Metaclass
    """
    _instances = {}
    _lock = threading.Lock()

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            with cls._lock:
                # Double-checked locking
                if cls not in cls._instances:
                    cls._instances[cls] = super().__call__(*args, **kwargs)
        return cls._instances[cls]


class GeneMapper(ThreadSafeLazyLoader, metaclass=SingletonMeta):
    def __init__(self):
        file_path = os.path.join(
            "data", 
            "genetable20201208_median_isoform_mass.tsv"
        )
        read_kwargs = {
            "dtype": {"GeneID": str}, 
            "index_col": "GeneID"
        }
        super().__init__(file_path, read_kwargs)

        # Memoized derived dictionaries
        self._symbol = None
        self._funcat = None
        self._description = None
        self._taxon = None

    def _derive_mappings(self):
        """
        Derive mappings only when first accessed and dataframe is loaded
        """
        if self.df is not None:
            self._symbol = self.df["GeneSymbol"].to_dict()
            self._funcat = self.df["FunCats"].fillna("").to_dict()
            self._description = self.df["GeneDescription"].to_dict()
            self._taxon = self.df["TaxonID"].to_dict()

    @property
    def symbol(self):
        if self._symbol is None:
            self._derive_mappings()
        return self._symbol

    @property
    def funcat(self):
        if self._funcat is None:
            self._derive_mappings()
        return self._funcat

    # Similar methods for description and taxon
