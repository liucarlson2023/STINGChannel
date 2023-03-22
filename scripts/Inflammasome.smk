import sys
sys.path.append("/home/rcarlson/mountb/Channel/OpticalPooledScreens")
from functools import partial
import ops.annotate
import ops.firesnake
from ops.firesnake import Snake
import ops.io
import pandas as pd

WELLS = config['WELLS']
TILES = config['TILES']
TILES = [str(item).zfill(4) for item in TILES]

# display options for saved .tif files (view in ImageJ)
channels = ('DAPI', 'NLRP3', 'PSTING', 'STING')
LUTS = [getattr(ops.io, config['LUTS'][x]) for x in channels]
DISPLAY_RANGES = [config['DISPLAY_RANGES'][x] for x in channels]

# set paramspaces if a paramsearch mode is selected
if config['MODE'] == 'paramsearch_segmentation':
    (config,
        nuclei_segmentation_paramspace,
        cell_segmentation_paramspace) = ops.firesnake.initialize_paramsearch(config)
elif config['MODE']!='process':
    raise ValueError(f'MODE="{config["MODE"]}" not recognized, use either "process" or "paramsearch"')
else:
    if isinstance(config['NUCLEUS_AREA'][0],list):
        raise ValueError('NUCLEUS_AREA cannot be a list of lists for MODE="process"')

# naming convention for input and processed files
input_files = partial(ops.firesnake.input_files_live,
                      directory=config['INPUT_DIRECTORY'], channels=config['CHANNELS'])

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
        [expand(processed_input(x), well=WELLS, tile=TILES)
            for x in config['REQUESTED_TAGS']],
           [config['PROCESS_DIRECTORY'] + '/' + x for x in config['REQUESTED_FILES']],

rule nd2_to_tif:
    input:
        input_files(config['INPUT_TAG'])
    output:
        processed_output('inflammasome.tif')
    run:
        Snake.nd2_to_tif(output=output, data=input, display_ranges=DISPLAY_RANGES, luts=LUTS)

rule segment:
    input:
        processed_output('inflammasome.tif')
    output:
        processed_output('nuclei.tif'),
        processed_output('cells.tif'),
    run:
        if config['SEGMENT_METHOD'] == 'cell_2019':
            Snake.segment_cell_2019(
                output=output, 
                data=input[0],
                nuclei_threshold=config['THRESHOLD_DAPI'],
                nuclei_area_min=config['NUCLEUS_AREA'][0],
                nuclei_area_max=config['NUCLEUS_AREA'][1],
                cell_threshold=config['THRESHOLD_CELL'],
                smooth=config['SMOOTH'],
                radius=config['RADIUS']
            )
        elif config['SEGMENT_METHOD'] == 'cell_2019_select_channels':
            Snake.segment_cell_2019_select_channels(
                output=output,
                data=input[0],
                nuclei_threshold=config['THRESHOLD_DAPI'],
                nuclei_area_min=config['NUCLEUS_AREA'][0],
                nuclei_area_max=config['NUCLEUS_AREA'][1],
                cell_threshold=config['THRESHOLD_CELL'],
                smooth=config['SMOOTH'],
                radius=config['RADIUS'],
                chstart=config['SEGMENT_CHSTART'], 
                chend=config['SEGMENT_CHEND'])
        elif config['SEGMENT_METHOD'] == 'cellpose':
            # last cycle
            cycle = config['CELLPOSE']['CYTO_CYCLE']
            data = ops.io.read_stack(input[0])[cycle]
            Snake.segment_cellpose(
                output=output, 
                data=data, 
                dapi_index=0, 
                cyto_index=config['CELLPOSE']['CYTO_CHANNEL'],
                diameter=config['CELLPOSE']['DIAMETER'],
                )
        else:
            error = ('config entry SEGMENT_METHOD must be "cell_2019" or "cellpose", '
                     f'not {config["SEGMENT_METHOD"]}')
            raise ValueError(error)

rule find_psting_mask:
    input:
        processed_input('inflammasome.tif'),
        processed_input('cells.tif'),
    output:
        processed_output('cells_psting.tif')
    run:
        Snake.find_mask(output=output, data = input[0], input_mask=input[1], channel=config['PUNCTA_CHANNEL'], threshold=config['PUNCTA_THRESHOLD'])


rule call_pheno_cells:
     input:
        processed_input('inflammasome.tif'),
        processed_input('cells_psting.tif'),
        processed_input('cells.tif'),
     output:
        processed_output('pheno_cells.csv')
     run:
        Snake.extract_phenotype_psting(output=output, data_phenotype=input[0], puncta=input[1],
              cells=input[2], wildcards=wildcards)


