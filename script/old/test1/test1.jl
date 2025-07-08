# Import required Julia packages for oceanographic data analysis
using NCDatasets      # For reading and writing NetCDF files
using PhysOcean       # Physical oceanography utilities
using DataStructures  # For ordered dictionaries and other data structures
using DIVAnd          # Data-Interpolating Variational Analysis in n-dimensions
using PyPlot          # Plotting library (matplotlib wrapper)
using Dates           # Date and time handling
using Statistics      # Statistical functions (mean, etc.)
using Random          # Random number generation
using Printf          # String formatting with printf-style syntax

# Define spatial grid parameters for the Adriatic Sea analysis
dx, dy = 0.125, 0.125  # Grid resolution in degrees (longitude, latitude)
lonr = 11.5:dx:20      # Longitude range from 11.5° to 20° E with 0.125° spacing
latr = 39:dy:46        # Latitude range from 39° to 46° N with 0.125° spacing
timerange = [Date(1950,1,1),Date(2017,12,31)];  # Time period for analysis

# Define depth levels for 3D analysis (in meters)
# Full depth range commented out, using simplified 3-level version
depthr = [0.,5., 10., 15., 20., 25., 30., 40., 50., 66, 
    75, 85, 100, 112, 125, 135, 150, 175, 200, 225, 250, 
    275, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 
    800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 
    1300, 1350, 1400, 1450, 1500, 1600, 1750, 1850, 2000];
depthr = [0.,10.,20.];  # Simplified to 3 depth levels: surface, 10m, 20m

# Define analysis parameters
varname = "Salinity"    # Variable being analyzed
yearlist = [1900:2017]; # Years to include in analysis
monthlist = [[1,2,3],[4,5,6],[7,8,9],[10,11,12]]; # Seasonal groupings (quarters)

# Create time selector for seasonal analysis
TS = DIVAnd.TimeSelectorYearListMonthList(yearlist,monthlist);
@show TS;

# ========================================================================
# DATA DOWNLOAD AND PREPARATION SECTION
# ========================================================================

# Set up directory paths for data files
datadir = "./"
smalldatafile = joinpath(datadir, "AdriaticSea_SDC_1000.nc")  # Test dataset with 1000 observations
datafile = joinpath(datadir, "AdriaticSea_SDC.nc")             # Full dataset

# Create data directory if it doesn't exist
if !isdir(datadir)
    @info("Creating data directory")
    mkdir(datadir) 
end

# Download small test dataset (1000 observations) if not already present
if !isfile(smalldatafile)
    @info("Downloading test data file (1000 lines)")
    download("https://dox.ulg.ac.be/index.php/s/1CevuhrnDW18fJT/download", smalldatafile)
else
    @info("Data file already downloaded")
end

# Download full dataset if not already present
if !isfile(datafile)
    @info("Downloading data file")
    download("https://dox.ulg.ac.be/index.php/s/IRYJyNZ5KuKVoQL/download", datafile)
else
    @info("Data file already downloaded")
end

# Load observational data from NetCDF files
# First load from small dataset (for testing)
@time obsval,obslon,obslat,obsdepth,obstime,obsid = NCODV.load(Float64, smalldatafile, 
    "Water body salinity");

# Then load from full dataset (overwrites the small dataset variables)
@time obsval,obslon,obslat,obsdepth,obstime,obsid = NCODV.load(Float64, datafile, 
    "Water body salinity");

# ========================================================================
# PLOTTING OBSERVATIONAL DATA DISTRIBUTION
# ========================================================================

# Create a figure showing the geographic distribution of observation points
figure("Adriatic-Data")
ax = subplot(1,1,1)
plot(obslon, obslat, "ko", markersize=.1)  # Plot observation locations as small black dots
aspectratio = 1/cos(mean(latr) * pi/180)   # Calculate proper aspect ratio for latitude
ax.tick_params("both",labelsize=6)
gca().set_aspect(aspectratio)

# Check quality and consistency of observations
checkobs((obslon,obslat,obsdepth,obstime),obsval,obsid)

# ========================================================================
# BATHYMETRY DATA DOWNLOAD AND PROCESSING
# ========================================================================

# Download bathymetry data (seafloor depth) for the region
bathname = "gebco_30sec_8.nc"
if !isfile(bathname)
    download("https://dox.ulg.ac.be/index.php/s/U0pqyXhcQrXjEUX/download",bathname)
else
    @info("Bathymetry file already downloaded")
end

# Load bathymetry data and interpolate to our grid
@time bx,by,b = load_bath(bathname,true,lonr,latr);

# Plot the bathymetry data
figure("Adriatic-Bathymetry")
ax = subplot(1,1,1)
pcolor(bx, by, permutedims(b, [2,1]));  # Create colored map of bathymetry
colorbar(orientation="vertical", shrink=0.8).ax.tick_params(labelsize=8)
contour(bx, by, permutedims(b, [2,1]), [0, 0.1], colors="k", linewidths=.5)  # Add coastline contour
gca().set_aspect(aspectratio)
ax.tick_params("both",labelsize=6)

# ========================================================================
# MASK CREATION AND EDITING FOR ANALYSIS DOMAIN
# ========================================================================

# Create a 3D mask for the analysis domain
# This mask determines which grid points are valid for analysis (water vs land)
mask = falses(size(b,1),size(b,2),length(depthr))
for k = 1:length(depthr)
    for j = 1:size(b,2)
        for i = 1:size(b,1)
            mask[i,j,k] = b[i,j] >= depthr[k]  # True where water depth >= analysis depth
        end
    end
end
@show size(mask)

# Plot the initial mask (surface level)
figure("Adriatic-Mask")
ax = subplot(1,1,1)
gca().set_aspect(aspectratio)
ax.tick_params("both",labelsize=6)
pcolor(bx,by, transpose(mask[:,:,1])); 

# Create coordinate grids for mask editing
grid_bx = [i for i in bx, j in by];
grid_by = [j for i in bx, j in by];

# Edit the mask to remove specific regions
mask_edit = copy(mask);
sel_mask1 = (grid_by .<= 42.6) .& (grid_bx .<= 14.);  # Remove northern shallow areas
sel_mask2 = (grid_by .<= 41.2) .& (grid_bx .<= 16.2); # Remove central/southern regions
mask_edit = mask_edit .* .!sel_mask1 .* .!sel_mask2;   # Apply logical NOT and multiply
@show size(mask_edit)

# Plot the edited mask
figure("Adriatic-Mask-Edited")
ax = subplot(1,1,1)
ax.tick_params("both",labelsize=6)
pcolor(bx, by, transpose(mask[:,:,1])); 
gca().set_aspect(aspectratio)


# ========================================================================
# DATA FILTERING AND QUALITY CONTROL
# ========================================================================

# Filter observational data to keep only realistic salinity values
sel = (obsval .<= 40) .& (obsval .>= 25);  # Typical Adriatic Sea salinity range

# Apply the filter to all observation arrays
obsval = obsval[sel]
obslon = obslon[sel]
obslat = obslat[sel]
obsdepth = obsdepth[sel]
obstime = obstime[sel]
obsid = obsid[sel];

# ========================================================================
# DIVAND ANALYSIS PARAMETERS SETUP
# ========================================================================

# Optional: Calculate observation weights based on data density
# Uncommented code would create spatially varying error estimates
#@time rdiag=1.0./DIVAnd.weight_RtimesOne((obslon,obslat),(0.03,0.03));
#@show maximum(rdiag),mean(rdiag)

# Define grid dimensions for parameter arrays
sz = (length(lonr),length(latr),length(depthr));

# Set correlation lengths (influence radius) for each dimension
lenx = fill(100_000.,sz)   # 100 km correlation length in longitude direction
leny = fill(100_000.,sz)   # 100 km correlation length in latitude direction  
lenz = fill(25.,sz);       # 25 m correlation length in depth direction
len = (lenx, leny, lenz);  # Combine into tuple for DIVAnd

# Set noise-to-signal ratio (regularization parameter)
epsilon2 = 0.1;            # Controls smoothness vs data fidelity tradeoff
#epsilon2 = epsilon2 * rdiag;  # Optional: spatially varying epsilon

# ========================================================================
# OUTPUT FILE SETUP AND METADATA CONFIGURATION
# ========================================================================

# Set up output directory and filename
outputdir = "./"
if !isdir(outputdir)
    mkpath(outputdir)
end
filename = joinpath(outputdir, "Water_body_$(replace(varname," "=>"_"))_Adriatic.4Danl.nc")

# Define comprehensive metadata for NetCDF file following SeaDataNet standards
metadata = OrderedDict(
    # Name of the project (SeaDataCloud, SeaDataNet, EMODNET-chemistry, ...)
    "project" => "SeaDataCloud",

    # URN code for the institution EDMO registry,
    # e.g. SDN:EDMO::1579
    "institution_urn" => "SDN:EDMO::1579",

    # Production group
    #"production" => "Diva group",

    # Name and emails from authors
    "Author_e-mail" => ["Your Name1 <name1@example.com>", "Other Name <name2@example.com>"],

    # Source of the observation
    "source" => "observational data from SeaDataNet and World Ocean Atlas",

    # Additional comment
    "comment" => "Duplicate removal applied to the merged dataset",

    # SeaDataNet Vocabulary P35 URN
    # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p35
    # example: SDN:P35::WATERTEMP
    "parameter_keyword_urn" => "SDN:P35::EPC00001",

    # List of SeaDataNet Parameter Discovery Vocabulary P02 URNs
    # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p02
    # example: ["SDN:P02::TEMP"]
    "search_keywords_urn" => ["SDN:P02::PSAL"],

    # List of SeaDataNet Vocabulary C19 area URNs
    # SeaVoX salt and fresh water body gazetteer (C19)
    # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=C19
    # example: ["SDN:C19::3_1"]
    "area_keywords_urn" => ["SDN:C19::3_3"],

    "product_version" => "1.0",
    
    "product_code" => "something-to-decide",
    
    # bathymetry source acknowledgement
    # see, e.g.
    # * EMODnet Bathymetry Consortium (2016): EMODnet Digital Bathymetry (DTM).
    # https://doi.org/10.12770/c7b53704-999d-4721-b1a3-04ec60c87238
    # 
    # taken from
    # http://www.emodnet-bathymetry.eu/data-products/acknowledgement-in-publications
    #
    # * The GEBCO Digital Atlas published by the British Oceanographic Data Centre on behalf of IOC and IHO, 2003
    #
    # taken from
    # https://www.bodc.ac.uk/projects/data_management/international/gebco/gebco_digital_atlas/copyright_and_attribution/
        
    "bathymetry_source" => "The GEBCO Digital Atlas published by the British Oceanographic Data Centre on behalf of IOC and IHO, 2003",

    # NetCDF CF standard name
    # http://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
    # example "standard_name" = "sea_water_temperature",
    "netcdf_standard_name" => "sea_water_salinity",

    "netcdf_long_name" => "sea water salinity",

    "netcdf_units" => "1e-3",

    # Abstract for the product
    "abstract" => "...",

    # This option provides a place to acknowledge various types of support for the
    # project that produced the data
    "acknowledgement" => "...",

    "documentation" => "https://doi.org/doi_of_doc",

    # Digital Object Identifier of the data product
    "doi" => "...");

# Convert metadata to NetCDF-compatible attributes
ncglobalattrib, ncvarattrib = SDNMetadata(metadata, filename, varname, lonr, latr)

# Remove any existing analysis file to start fresh
if isfile(filename)
    rm(filename) # delete the previous analysis
    @info "Removing file $filename"
end

# ========================================================================
# PLOTTING FUNCTION DEFINITION
# ========================================================================

# Set up figure output directory
figdir = "./"

# Define a function to plot interpolation results for each time step
function plotres(timeindex,sel,fit,erri)
    tmp = copy(fit)                            # Copy the fitted data to avoid modifying original
    nx,ny,nz = size(tmp)                       # Get dimensions of the fitted data array
    
    for i in 1:nz                             # Loop through each depth level
        figure("Adriatic-Additional-Data")     # Create or select figure window
        ax = subplot(1,1,1)                   # Create subplot
        ax.tick_params("both",labelsize=6)    # Set tick parameters
        ylim(39.0, 46.0);                     # Set latitude limits
        xlim(11.5, 20.0);                     # Set longitude limits
        title("Depth: $(depthr[i]) \n Time index: $(timeindex)", fontsize=6)  # Add title with depth and time info
        
        # Create colored plot of the interpolated salinity field
        pcolor(lonr.-dx/2.,latr.-dy/2, permutedims(tmp[:,:,i], [2,1]);
               vmin = 33, vmax = 40)           # Set color scale limits for salinity
        colorbar(extend="both", orientation="vertical", shrink=0.8).ax.tick_params(labelsize=8)

        # Add land mask as gray contour 
        contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0],colors = [[.5,.5,.5]])
        aspectratio = 1/cos(mean(latr) * pi/180)  # Calculate proper aspect ratio
        gca().set_aspect(aspectratio)
        
        # Save the figure with formatted filename
        figname = varname * @sprintf("_%02d",i) * @sprintf("_%03d.png",timeindex)
        PyPlot.savefig(joinpath(figdir, figname), dpi=600, bbox_inches="tight");
        PyPlot.close_figs()                   # Close figure to free memory
    end
end

# ========================================================================
# MAIN DIVAND ANALYSIS EXECUTION
# ========================================================================

# Execute the main DIVAnd 3D analysis
@time dbinfo = diva3d((lonr,latr,depthr,TS),        # Grid coordinates and time selector
    (obslon,obslat,obsdepth,obstime), obsval,        # Observation coordinates and values
    len, epsilon2,                                    # Correlation lengths and regularization
    filename,varname,                                 # Output file and variable name
    bathname=bathname,                               # Bathymetry file for land/sea mask
    plotres = plotres,                               # Plotting function to call during analysis
    mask = mask_edit,                                # Edited mask for analysis domain
    fitcorrlen = false,                              # Don't fit correlation lengths automatically
    niter_e = 2,                                     # Number of iterations for error estimation
    ncvarattrib = ncvarattrib,                       # NetCDF variable attributes
    ncglobalattrib = ncglobalattrib,                 # NetCDF global attributes
    surfextend = true                                # Extend surface values to deeper levels if needed
    );

# Save observation metadata to the output file
DIVAnd.saveobs(filename,(obslon,obslat,obsdepth,obstime),obsid);