########## ICESAT-2 ATL08 PRE-PROCESSING ###########
## take list of ATL08 h5 files and create data tables with necessary variables
## restructure 20m sub-segments to single columns
## filter by a range of quality assurance flags
## create a canopy present/absent variable and filter by land cover type

## github repository for icesat2R package: https://github.com/mqueinnec/icesat2R
library(icesat2R); library(dplyr); library(data.table); library(terra)

## directory containing .h5 files to process
dir = "data/"

# list of all ATL08 files in directory
atl08_flist = list.files(path = dir,
                         pattern = "ATL08", 
                         full.names = TRUE)

# read .h5 as data frames save as .csv
lapply(atl08_flist, 
       function(x)read_ATL08(atl08_h5 = x,
                             beam_strength = 'strong', # only process strong beams
                             odir = "data/", 
                             atl03_pc = FALSE)) # option to include ATL03 photon classification

## combine all csv and filter by quality assurance flags
#get filepaths of all ATL08 csv
in_paths = dir(path = "data/",
                pattern = ".csv",
                full.names = TRUE)

# function that reads csv to data table and appends source filename 
read_csv_filename = function(filename){
  ret = fread(filename)
  ret$Source = filename
  ret
}

# combine all csv to one large data.table
df = lapply(in_paths, function(x) read_csv_filename(filename = x)) %>%
  rbindlist()

# create variables for date and time from filename
df$date = substr(sub(".*ATL08_", "", df$Source), 1, 8)

# filter by QA flags....
df_QA = df %>% filter(., cloud_flag_atm == 0) %>% #cloud flag - 0-10, 0 = no cloud/aerosol present
  filter(., layer_flag == 0) %>% #layer flag: 0-1, 0 = absence of clouds or blowing snow
  filter(., segment_watermask == 0) %>% #filter water 
  filter(., terrain_flg == 0) %>% #terrain flag: quality check on deviation from reference DEM
  filter(., msw_flag == 0) %>% #multiple scattering flag: 0 = no layers
  filter(., night_flag == 1 ) #remove daytime

# list of columns to subset  
# X1:X5 = 20m sub-segments within 100m segments
cols = c('rgt', 'cycle','beam','UTM_zone', 'date' , 'segment_landcover', 'segment_id_beg', 'segment_id_end',
          'latitude_20m_X1', 'latitude_20m_X2','latitude_20m_X3','latitude_20m_X4', 'latitude_20m_X5', 
          'longitude_20m_X1','longitude_20m_X2','longitude_20m_X3','longitude_20m_X4', 'longitude_20m_X5', 
          'h_canopy_20m_X1', 'h_canopy_20m_X2', 'h_canopy_20m_X3','h_canopy_20m_X4', 'h_canopy_20m_X5', 'Source')

# subset columns
df_QA = select(df_QA, cols)

# pivot 20m sub-segments into a single column per variable
# create sub-segment ids using using the first id of each 100m segment
df_QA$segment_id_X2 = df_QA$segment_id_beg + 1
df_QA$segment_id_X3 = df_QA$segment_id_beg + 2
df_QA$segment_id_X4 = df_QA$segment_id_beg + 3

# for each 20m variable select the relevent columns (5) and stack into a single column 
# the 100m variables are replicated across each 20m sub-segment 5 times
h20 = cbind(stack(select(df_QA, h_canopy_20m_X1: h_canopy_20m_X5)),
             stack(select(df_QA, segment_id_beg, segment_id_X2: segment_id_X4, segment_id_end)),
             stack(select(df_QA, latitude_20m_X1: latitude_20m_X5)),
             stack(select(df_QA, longitude_20m_X1: longitude_20m_X5)),
             rep(df_QA$rgt, 5),
             rep(df_QA$cycle, 5),
             rep(df_QA$beam, 5),
             rep(df_QA$UTM_zone, 5),
             rep(df_QA$segment_landcover, 5),
             rep(df_QA$date, 5),
             rep(df_QA$Source, 5))

# rename columns
# the stack function also adds an index indicating which vector the observation originated from
names(h20) = c('h_canopy_20m', 'h_canopy_20m_ind',
                'ph_segment_id', 'ph_segment_id_ind',
                'latitude_20m', 'latitude_20m_ind',
                'longitude_20m', 'longitude_20m_ind',
                'rgt', 'cycle', 'beam', 'UTM_zone', 'segment_landcover', 'date', 'Source')

# drop index columns 
h20 = select(h20, !(contains('ind')))

# create variable for canopy presence or absence
# if 20m segment has fewer than 4 photons classified as vegetation then it is flagged with large integer 
# If height > large integer then assign 0 for canopy absence
h20$Canopy  =  ifelse(h20$h_canopy_20m < 1000, 1, 0)

# Assign height of 0 segments with absent canopy
h20[h20$Canopy == 0, ]$h_canopy_20m = 0

# filter canopy absent segments by landcover type using segment_landcover flag
# only keep canopy absent segments if they also fall within herbaceous/graminoid landcover

# land cover codes for low stature vegetation
# 30 = Herbaceous
# 60 = Bare_sparse_veg
# 90 = Herbaceous wetland
# 100 = moss_lichen

# all land cover codes
values = c(0, 111, 113, 112, 114, 115, 116, 121, 123, 122, 124, 125, 126,
            20, 30, 90, 100, 60, 40, 50, 70, 80, 200, 255)
keep_values = c(30, 60, 90, 100)
values = values[!values %in% keep_values]

# Remove canopy absent segments if they are not in land cover type with low stature veg
h20 = h20[!(h20$Canopy == 0 & h20$segment_landcover %in% values), ]

#filter segments with invalid lat/long
h20 = h20[h20$latitude_20m < 90,]
