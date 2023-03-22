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

rule adjust_luts:
     input:
        input_files(config['INPUT_TAG']),
     output:
        processed_output('processed.tif')
     run:
        Snake.nd2_to_tif(output=output, data=input, axis = 0, display_ranges=DISPLAY_RANGES, luts=LUTS,
              wildcards=wildcards)


rule ratio:
    input:
        input_files(config['INPUT_TAG'])
    output:
        processed_output('ratio.tif')
    run:
        Snake.sep_mruby3_ratio_sr(output=output, data=input, display_ranges=DISPLAY_RANGES, luts=LUTS)

