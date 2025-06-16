""" Main Script from which the whole program runs.
"""
from argparse import ArgumentParser

from ppm.config import Config
from ppm.create_background import create_bg
from ppm.create_pisces_db import create_pisces_db
from ppm.analysis import analyse_model
from ppm.analysis_spliced import analyse_spliced_model
from ppm.preprocess import preprocess_canonical, preprocess_cryptic
from ppm.preprocess_spliced import process_spliced
from ppm.score import score_multi_mappers
from ppm.train import train_models
from ppm.train_spliced import train_spliced_models

PPM_PIPELINES = [
    'createPiscesDB',
    'bg',
    'all',
    'preprocess',
    'train', 'train+',
    'analysis',
    'score',
]

def get_arguments():
    """ Function to collect command line arguments.

    Returns
    -------
    args : argparse.Namespace
        The parsed command line arguments.
    """
    parser = ArgumentParser(description='inSPIRE Pipeline for MS Search Results.')

    parser.add_argument(
        '--config_file',
        default=None,
        help='Config file to be read from.',
        required=False,
    )
    parser.add_argument(
        '--pipeline',
        choices=PPM_PIPELINES,
        default='all',
        help='Pipeline to run.',
        required=False,
    )
    parser.add_argument(
        '--pep_length',
        default=0,
        help='Pipeline to run.',
        required=False,
        type=int,
    )

    return parser.parse_args()

def run_ppm(pipeline=None, config_file=None):
    """ Main function.
    """
    if pipeline is None:
        args = get_arguments()
        pipeline = args.pipeline
        config_file = args.config_file

    config = Config(config_file)


    if pipeline == 'createPiscesDB':
        create_pisces_db(config)

    if pipeline == 'bg':
        create_bg(config, args.pep_length)

    if pipeline in ('preprocess', 'all'):
        if config.model == 'cryptic':
            preprocess_cryptic(config)
        elif config.model == 'canonical':
            preprocess_canonical(config)
        elif config.model == 'spliced':
            process_spliced(config)

    if pipeline in ('train', 'train+', 'all'):
        if config.model in ('cryptic', 'canonical'):
            train_models(config)
        if config.model == 'spliced':
            train_spliced_models(config)

    if pipeline in ('analysis', 'train+', 'all'):
        print(f'Running analysis for model {config.model}')
        if config.model in ('canonical', 'cryptic'):
            analyse_model(config)
        elif config.model == 'spliced':
            analyse_spliced_model(config)

    if pipeline == 'score':
        score_multi_mappers(config)

if __name__ == '__main__':
    run_ppm()



