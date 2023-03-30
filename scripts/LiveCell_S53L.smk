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
        processed_output('live.tif')
    run:
        Snake.nd2_to_tif(output=output, data=input, display_ranges=DISPLAY_RANGES, luts=LUTS)

rule ratio:
    input:
        input_files(config['INPUT_TAG'])
    output:
        processed_output('ratio.tif')
    run:
        Snake.sep_mruby3_ratio(output=output, data=input, display_ranges=DISPLAY_RANGES, luts=LUTS)

rule segment_live:
    input:
        processed_output('live.tif')
    output:
        processed_output('nuclei.tif'),
        processed_output('cells.tif'),
    run:
        Snake.segment_cell_2019_select_channels_stack(
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

rule find_cytoplasm:
    input:
        processed_input('nuclei.tif'),
        processed_input('cells.tif'),
    output:
        processed_output('cytoplasm.tif')
    run:
        Snake.find_cytoplasm(output=output, nuclei=input[0],
            cells=input[1])

rule track_nuclei:
    input:
        processed_input('nuclei.tif'),
    output:
        processed_output('nuclei_tracked.tif')
    run:
        Snake.track_live_nuclei(output=output, nuclei=input[0], tolerance_per_frame=150, cutoff=200)

rule track_cells:
    input:
        processed_input('cells.tif'),
    output:
        processed_output('cells_tracked.tif')
    run:
        Snake.track_live_nuclei(output=output, nuclei=input[0], tolerance_per_frame=150, cutoff=200)

rule track_cytoplasm:
    input:
        processed_input('cytoplasm.tif'),
    output:
        processed_output('cytoplasm_tracked.tif')
    run:
        Snake.track_live_nuclei(output=output, nuclei=input[0], tolerance_per_frame=150, cutoff=200)

rule call_pheno_tracked:
     input:
        processed_input('live.tif'),
        processed_input('nuclei_tracked.tif'),
        processed_input('cells_tracked.tif'),
        processed_input('cytoplasm_tracked.tif'),
     output:
        processed_output('pheno_tracked.csv')
     run:
        Snake.extract_phenotype_hoechst_sep_mruby3_live(output=output, data=input[0], nuclei=input[1], 
              cells=input[2], cytoplasm=input[3], wildcards=wildcards)

rule call_pheno:
     input:
        processed_input('live.tif'),
        processed_input('nuclei.tif'),
        processed_input('cells.tif'),
        processed_input('cytoplasm.tif'),
     output:
        processed_output('pheno.csv')
     run:
        Snake.extract_phenotype_hoechst_sep_mruby3_live(output=output, data=input[0], nuclei=input[1],
              cells=input[2], cytoplasm=input[3], wildcards=wildcards)
