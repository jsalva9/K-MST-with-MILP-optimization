import yaml
from pathlib import Path


class Config:
    def __init__(self):
        """
        Initialize Config class. Load config.yaml file and define paths.
        """
        # Define root directory, one level up from current directory
        self.root = Path(__file__).parent.parent
        self.config_path = f'{self.root}/config.yaml'

        self.config = self.load_config()
        self.data_path = f'{self.root}/{self.config["basic"]["data_path"]}'
        self.output_path = f'{self.root}/{self.config["basic"]["output_path"]}'
        self.config["basic"]['data_path'] = self.data_path
        self.config['output_path'] = self.output_path

    def load_config(self):
        """
        Load config.yaml file
        Returns:
            config (dict): Dictionary with config.yaml file
        """

        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config
