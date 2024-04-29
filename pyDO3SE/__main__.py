#!/usr/bin/env python3
from pyDO3SE.tools.cli.root import cli
from multiprocessing import freeze_support


if __name__ == "__main__":
    ''' Get the config and input data locations from the cli
    Then loads these into the config and runs the model

    '''
    freeze_support()
    cli()
