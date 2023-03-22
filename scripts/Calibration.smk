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
CYCLES = [config['CYCLE_FORMAT'].format(cycle=x) for x in config['CYCLES']]

# display options for saved .tif files (view in ImageJ)
channels = ('HOECHST', 'SEP', 'mRuby3', 'Ratio')
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
input_files = partial(ops.firesnake.input_files_calibration,
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

rule align:
    priority: -1
    input:
        input_files(config['INPUT_TAG'], CYCLES)
    output:
        processed_output('aligned.tif')
    run:
        Snake.align_SBS(output=output, data=input,
            display_ranges=DISPLAY_RANGES, luts=LUTS, wildcards = wildcards)

rule ratio:
    input:
        processed_output('aligned.tif')
    output:
        processed_output('ratio.tif')
    run:
        Snake.sep_mruby3_ratio(output=output, data=input, display_ranges=DISPLAY_RANGES, luts=LUTS)

rule segment:
    input:
        input_files(config['INPUT_TAG'], CYCLES[0]),
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

rule find_cytoplasm:
    input:
        processed_input('nuclei.tif'),
        processed_input('cells.tif'),
    output:
        processed_output('cytoplasm.tif')
    run:
        Snake.find_cytoplasm(output=output, nuclei=input[0],
            cells=input[1])

rule call_pheno:
     input:
        processed_input('aligned.tif'),
        processed_input('nuclei.tif'),
        processed_input('cells.tif'),
        processed_input('cytoplasm.tif'),
     output:
        processed_output('pheno.csv')
     run:
        Snake.extract_phenotype_hoechst_sep_mruby3_live(output=output, data=input[0], nuclei=input[1], 
              cells=input[2], cytoplasm=input[3], wildcards=wildcards)


