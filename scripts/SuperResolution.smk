import sys
sys.path.append("/home/rcarlson/mountb/Channel/OpticalPooledScreens")
from functools import partial
import ops.annotate
import ops.firesnake
from ops.firesnake import Snake
import ops.io
import pandas as pd

WELLS = config['WELLS']

# display options for saved .tif files (view in ImageJ)
channels = ('SEP', 'mRuby3', 'STING', 'Ratio')
LUTS = [getattr(ops.io, config['LUTS'][x]) for x in channels]
DISPLAY_RANGES = [config['DISPLAY_RANGES'][x] for x in channels]

# naming convention for input and processed files
input_files = partial(ops.firesnake.input_files_superres,
                      directory=config['INPUT_DIRECTORY'])

processed_input = partial(ops.firesnake.processed_file,
                         directory=config['PROCESS_DIRECTORY'],magnification=config['MAGNIFICATION']
                         )

processed_output = partial(ops.firesnake.processed_file,
                         directory=config['PROCESS_DIRECTORY'],
                         temp_tags=config['TEMP_TAGS'],magnification=config['MAGNIFICATION']
                         )

rule all:
    input:
        # request individual files or list of files
        [expand(processed_input(x), well=WELLS, tile = 0)
            for x in config['REQUESTED_TAGS']],
           [config['PROCESS_DIRECTORY'] + '/' + x for x in config['REQUESTED_FILES']],

rule segment_puncta:
    input:
        input_files(config['INPUT_TAG'])
    output:
        processed_output('puncta.tif')
    run:
        Snake.segment_puncta_stack(output=output, data=input,
            channel_axis=config['GOLGI_AXIS'],
            channel=config['GOLGI_CHANNEL'], threshold=config['GOLGI_THRESHOLD'],
            radius=config['GOLGI_RADIUS'], area_min=config['GOLGI_AREA'][0], area_max=config['GOLGI_AREA'][1])


rule call_puncta:
     input:
        input_files(config['INPUT_TAG']),
        processed_input('puncta_masked.tif'),
     output:
        processed_output('pheno_golgi.csv')
     run:
        Snake.extract_phenotype_golgi_puncta_live(output=output, data_phenotype=input[0], 
              puncta=input[1], wildcards=wildcards)


rule mask_puncta:
     input:
        processed_input('puncta.tif'),
        input_files('mask.tif')
     output:
        processed_output('puncta_masked.tif')
     run:
        Snake.mask_data(output=output, data=input[0],
              mask=input[1], wildcards=wildcards)


