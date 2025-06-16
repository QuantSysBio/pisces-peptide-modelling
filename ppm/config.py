""" Definition of Config class.
"""
import os
import glob
import shutil
from pathlib import Path

import yaml


class Config:
    """ Holder for configuration of the SPI-ART pipeline.
    """
    def __init__(self, config_file):
        with open(config_file, 'r', encoding='UTF-8') as stream:
            config_dict = yaml.safe_load(stream)

        self.path = config_file
        self._load_data(config_dict)
        self._clean_file_paths()

    def _clean_file_paths(self):
        """ Function to clean the file paths given to inspire.
        """
        home = str(Path.home())

        self.output_folder = self.output_folder.replace('~', home).replace('%USERPROFILE%', home)
        if self.output_folder.endswith('/'):
            self.output_folder = self.output_folder[:-1]

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        if not os.path.exists(f'{self.output_folder}/trainingDatasets'):
            os.makedirs(f'{self.output_folder}/trainingDatasets')
        if not os.path.exists(f'{self.output_folder}/mmDatasets'):
            os.makedirs(f'{self.output_folder}/mmDatasets')
        if not os.path.exists(f'{self.output_folder}/models'):
            os.makedirs(f'{self.output_folder}/models')
        if not os.path.exists(f'{self.output_folder}/imgs'):
            os.makedirs(f'{self.output_folder}/imgs')
        if not os.path.exists(f'{self.background_folder}'):
            os.makedirs(f'{self.background_folder}')
        if not os.path.exists(f'{self.background_folder}/frequency_dfs'):
            os.makedirs(f'{self.background_folder}/frequency_dfs')
        if not os.path.exists(f'{self.background_folder}/sample_ratios'):
            os.makedirs(f'{self.background_folder}/sample_ratios')
        if not os.path.exists(f'{self.background_folder}/random_dfs'):
            os.makedirs(f'{self.background_folder}/random_dfs')
        if not os.path.exists(f'{self.background_folder}/remapped'):
            os.mkdir(f'{self.background_folder}/remapped')
    
    def _load_data(self, config_dict):
        self.title = config_dict['title']
        self.output_folder = config_dict['outputFolder']
        self.meta_df = config_dict.get('metaDf')
        self.model = config_dict.get('model')
        self.cell_line = config_dict.get('cellLine')
        self.peptides_pq = config_dict.get('peptidesParquet')
        self.antigen_folder = config_dict.get('antigenFolder')
        self.background_folder = config_dict.get('backgroundFolder')
        self.n_cores = config_dict.get('nCores', 60)
        self.proteome = config_dict.get('proteome')
        self.cryptic_folder = config_dict.get('crypticFolder')
        self.canonical_results = config_dict.get('canonicalResults')

