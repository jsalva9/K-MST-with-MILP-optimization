import os
import yaml
from pathlib import Path


class Config:
    def __init__(self):
        # Define root directory, one level up from current directory
        self.root = Path(__file__).parent.parent
        self.config_path = f'{self.root}/config.yaml'

        self.config = self.load_config()
        self.data_path = f'{self.root}/{self.config["data_path"]}'

    def load_config(self):
        # Read yaml file in self.config_path
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config
