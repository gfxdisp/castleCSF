# castleCSF - A Contrast Sensitivity Function of Color, Area, Spatio-Temporal frequency, Luminance and Eccentricity

This repository contains both the code and the data for a contrast sensitivity model, as the function of
* spatial frequency in cycles per degree (for gratings)
* temporal frequency in Hz
* eccentricity in visual degrees
* lumunance in cd/m^2 (or nit)
* area in squared visual degrees
* Chromaticity of the background
* Chromatic modulation direction of the stimulus
* Shape of stimulus (gratings or disks)

The details about the model and the dataset can be found on the project [website](https://www.cl.cam.ac.uk/research/rainbow/projects/castleCSF/). 

## Code

Currently, the code is provided as a Matlab class in the directory `matlab`. 

### Example

Please check the directory `examples` to use castleCSF for predicting and plotting contrast sensitivity.

## Data

Each datapoint represents a Gabor patch at the detection threshold, either for individual observers, or averaged across all observers. The sensitivity is averaged in the log-contrast space. 

The data is stored in the CSV files:

* `data/data_individual.csv` - measurements for individual participants. 

If no individual data is available in a dataset/paper, it is excluded from `data_individual.csv`. 

* `data/data_aggregated.csv` - the sensitivities averaged across the participants.

The columns are identical as in `data_individual.csv` but without the columns `observer` and `age`. 

* `data/backgrounds.csv` - the LMS coordinates of the background/adaptation colour

The background IDs are unique across all datasets so that the tables `data*` and `backgrounds` can be merged using `bkg_id` as the key. 

bkg_id - the unique ID of the background colour 
L, M, S - LMS colour coordinates
R - rod response or scotopic luminance (CIE 1951 scotopic luminous function)
bkg_label - string label of the background (e.g. 'red', 'white', 'd65')
dataset - the ID of the dataset

* `data/color_direction.csv` - the LMS colour vector representing the colour axis along which the stimulus (Gabor) was modulated

The colour directions IDs are unique across all datasets so that the tables `data*` and `color_direction` can be merged using `color_direction` as the key. 

col_dir_id - the unique ID of the colour direction
L_delta , M_delta, S_delta - the vector defining the direction in the LMS colour space
dataset - the ID of the dataset

* `data/data_individual_merged.csv` - the same as `data_individual.csv` but merged with both `backgrounds.csv` and `color_direction.csv`. 

* `data/data_aggregated_merged.csv` - the same as `data_aggregated.csv` but merged with both `backgrounds.csv` and `color_direction.csv`. 

* `data/datasets.json` - metadata of each dataset

Columns:
- observer - The anonymized unique ID of an observer. Some observers can be common across the datasets. 
- age - age of the observer in years. NaN if it is unknown.
- luminance - luminance in cd/m^2, using standard CIE luminous efficiency function (or Y of CIE 1931 XYZ)
- s_frequency - spatial frequency in cycles per degree
- t_frequency - temporal frequency in cycles per degree
- orientation - spatial orientation of the stimulus in degrees
- col_dir_id - the ID of the colour direction. The LMS vectors of the colour directions are stored in `color_direction.csv`. 
- bkg_id - the ID of the background colour. The LMS coordinates of the background colour can be found in `backgrounds.csv`.
- ge_sigma - the standard deviation of the Gaussian envelope that limits the size of the Gabor patch
- ge_lambda - the number of cycles within the 1 standard deviation radius, computed as 2 * ge_sigma * frequency
- stimulus - shape of the stimulus
- area - the area of the stimulus in deg^2. For regular Gabor patches, it is computed as pi*ge_sigma^2
- log_cone_contrast - the log10 of cone contrast at the detection threshold
- var_log_cone_contrast - the variance of the log_cone_contrast. NaN if the variance for individual measurements is not available. The variance is typically estimated when fitting a psychometric function to nAFC measurements. 
- dataset - the ID of the dataset
- eccentricity - distance from the central foveal point in visual degrees
- vis_field - the angle that defines the position in the visual field at a certain eccentricity. The values:
  0 - temporal (horizontal, away from the nose);
  90 - superior (vertical, top)
  180 - nasa (horizontal, toward the nose)
  240 - inferior (vertical, bottom)
  The values are equivalent to the polar coordinates for the right eye. 
