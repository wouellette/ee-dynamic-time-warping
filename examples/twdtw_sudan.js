// ****************************************************************************************************************** //
// *********************************** TW-DTW Example for Sennar - Sudan ******************************************** //
// ****************************************************************************************************************** //

// Import external dependencies
var palettes = require('users/gena/packages:palettes');
var wrapper = require('users/adugnagirma/gee_s1_ard:wrapper');
var S2Masks = require('users/soilwatch/soilErosionApp:s2_masks.js');
var composites = require('users/soilwatch/soilErosionApp:composites.js');

// Import the Dynamic Time Warping script
var DTW = require('users/soilwatch/functions:dtw.js');

// Input data parameters
var CLASS_NAME = 'lc_class'; // Property name of the feature collection containing the crop type class attribute
var AGG_INTERVAL = 30; // Number of days to use to create the temporal composite for 2020
var TIMESERIES_LEN = 6; // Number of timestamps in the time series
var PATTERNS_LEN = 6; // Number of timestamps for the reference data points
var CLASS_NO = 7; // Number of classes to map. Those are: builtup, water, trees, baresoil, cropland, rangeland, wetland
var S2_BAND_LIST = ['B2', 'B3', 'B11', 'B12', 'ndvi']; // S2 Bands to use as DTW input
var S1_BAND_LIST = ['VV', 'VH']; // S1 Bands to use as DTW input
var BAND_NO = S1_BAND_LIST.concat(S2_BAND_LIST).length; // Number of bands to use for the DTW. Currently,
                 // The default 7 bands are: S2 NDVI, S2 B2, S2 B3, S2 B11, S2 B12, S1 VV, S1 VH
var DOY_BAND = 'doy'; // Name of the Day of Year band for the time-aware DTW implementation

// DTW time parameters
var BETA = 50; // Beta parameter for the Time-Weighted DTW, controlling the tolerance (in days) of the weighting.
var ALPHA = 0.1; // ALPHA parameter for the Time-Weighted DTW, controlling steepness of the logistic distribution

// Import external water mask dataset
var not_water = ee.Image("JRC/GSW1_2/GlobalSurfaceWater").select('max_extent').eq(0); // JRC Global Surface Water mask

// Import ALOS AW3D30 latest DEM version v3.2
var dem = ee.ImageCollection("JAXA/ALOS/AW3D30/V3_2").select("DSM");
dem = dem.mosaic().setDefaultProjection(dem.first().select(0).projection());

//Remove mountain areas that are not suitable for crop growth
var slope = ee.Terrain.slope(dem); // Calculate slope from the DEM data
var dem_mask = dem.lt(3600); // Mask elevation above 3600m, where no crops grow.
var slope_mask = slope.lt(30); // Mask slopes steeper than 30°, where no crops grow.
var crop_mask = dem_mask.and(slope_mask); // Combine the two conditions

// Function to calculate the NDVI for planet mosaics
var addNDVI = function(img){
  return img.addBands(img.normalizedDifference(['B8','B4']).multiply(10000).toInt16().rename('ndvi'));
};

// Define the are of interest to use
var adm0_name = 'Sudan';
var adm1_name = 'Sennar';
var adm2_name = 'Sennar';

// A dictionary that will be iterated over for multi-year land cover mapping.
// Comment out the years you do not wish to produce.
var year_dict = {//'2020': 'COPERNICUS/S2_SR',
                 '2019': 'COPERNICUS/S2',
                 '2018': 'COPERNICUS/S2',
                 '2017': 'COPERNICUS/S2'
                };

// Ensure that the retrieved county geometry is unique
var county = ee.Feature(
ee.FeatureCollection(ee.FeatureCollection("FAO/GAUL/2015/level2")
                    .filter(ee.Filter.and(ee.Filter.equals('ADM0_NAME', adm0_name),
                                          ee.Filter.equals('ADM2_NAME', adm2_name)))
                    ).first());

// Center map on the county and plot it on the map
Map.centerObject(county.geometry());
Map.layers().reset([ui.Map.Layer(county, {}, adm2_name)]);

// Function that performs DTW land classification for a given year.
var DTWClassification = function(year, collection_type){

  var date_range = ee.Dictionary({'start': year + '-07-01', 'end': year + '-12-30'}); // Second half of year used only.
  // Load the Sentinel-2 collection for the time period and area requested
  var s2_cl = S2Masks.loadImageCollection(collection_type, date_range, county.geometry());

  // Perform cloud masking using the S2 cloud probabilities assets from s2cloudless,
  // courtesy of Sentinelhub/EU/Copernicus/ESA
  var masked_collection = s2_cl
                          .filterDate(date_range.get('start'), date_range.get('end'))
                          .map(S2Masks.addCloudShadowMask(not_water, 1e4))
                          .map(S2Masks.applyCloudShadowMask)
                          .map(addNDVI); // Add NDVI to band list

  // Generate a list of time intervals for which to generate a harmonized time series
  var time_intervals = composites.extractTimeRanges(date_range.get('start'), date_range.get('end'), 30);

  // Generate harmonized monthly time series of FCover as input to the vegetation factor V
  var s2_ts = composites.harmonizedTS(masked_collection, S2_BAND_LIST, time_intervals, {agg_type: 'geomedian'});

  // Replace masked pixels by the mean of the previous and next timestamps
  // And add a Day of Year (DOY) band
  // This is the linear interpolation approach,
  // which is more lightweight than other gap-filling approaches like Savitzky-Golay or harmonic regression
  var s2_stack = ee.Image(s2_ts
                          .map(function(image){
                            var currentDate = ee.Date(image.get('system:time_start'));
                            var meanImage = s2_ts.filterDate(currentDate.advance(-AGG_INTERVAL-1, 'day'),
                                                                 currentDate.advance(AGG_INTERVAL+1, 'day')).mean();
                            // replace all masked values
                            var ddiff = currentDate.difference(ee.Date(ee.String(date_range.get('start'))), 'day');
                            return meanImage.where(image, image).addBands(ee.Image(ddiff).rename(DOY_BAND).toInt16());
                          })
                          .iterate(function(image, previous){return ee.Image(previous).addBands(image)}, ee.Image([])));

  // Define S1 preprocessing parameters, as per:
  // Version: v1.2
  // Date: 2021-03-10
  // Authors: Mullissa A., Vollrath A., Braun, C., Slagter B., Balling J., Gou Y., Gorelick N.,  Reiche J.
  // Sentinel-1 SAR Backscatter Analysis Ready Data Preparation in Google Earth Engine. Remote Sensing 13.10 (2021): 1954.
  // Description: This script creates an analysis ready S1 image collection.
  // License: This code is distributed under the MIT License.
  var parameter = {//1. Data Selection
                   START_DATE: date_range.get('start'),
                   STOP_DATE: date_range.get('end'),
                   POLARIZATION:'VVVH', // The polarization available may differ depending on where you are on the globe
                   ORBIT : 'DESCENDING', // The orbit availability may differ depending on where you are on the globe
                   // Check out this page to find out what parameters suit your area:
                   // https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-1/observation-scenario
                   GEOMETRY: county.geometry(),
                   //2. Additional Border noise correction
                   APPLY_ADDITIONAL_BORDER_NOISE_CORRECTION: true,
                   //3.Speckle filter
                   APPLY_SPECKLE_FILTERING: true,
                   SPECKLE_FILTER_FRAMEWORK: 'MULTI',
                   SPECKLE_FILTER: 'LEE',
                   SPECKLE_FILTER_KERNEL_SIZE: 9,
                   SPECKLE_FILTER_NR_OF_IMAGES: 10,
                   //4. Radiometric terrain normalization
                   APPLY_TERRAIN_FLATTENING: true,
                   DEM: dem,
                   TERRAIN_FLATTENING_MODEL: 'VOLUME', // More desirable for vegetation monitoring.
                                                       //Use "SURFACE" if working on urban or bare soil applications
                   TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER: 0,
                   //5. Output
                   FORMAT : 'DB',
                   CLIP_TO_ROI: false,
                   SAVE_ASSETS: false
  }

  //Preprocess the S1 collection
  var s1_ts = wrapper.s1_preproc(parameter)[1]
              .map(function(image){return image.multiply(1e4).toInt16() // Convert to Int16 using 10000 scaling factor
                                          .set({'system:time_start': image.get('system:time_start')})});

  // Create equally-spaced temporal composites covering the date range and convert to multi-band image
  var s1_stack = ee.Image(composites.harmonizedTS(s1_ts, S1_BAND_LIST, time_intervals, {agg_type: 'geomedian'})
                          .iterate(function(image, previous){return ee.Image(previous).addBands(image)}, ee.Image([])));

  // Re-order the stack order before converting it to a DTW-ready input array
  var s1s2_stack = s1_stack.addBands(s2_stack)
                   .select(ee.List(S1_BAND_LIST.concat(S2_BAND_LIST)).add(DOY_BAND) // Make sure DOY band goes last
                   .map(function(band){return ee.String(band).cat('.*')})) // Add regex for band selection
                   .unmask(0) // DTW does not tolerate null values, so gap fill to 0 if gaps remain
                   .clip(county.geometry()); // Clip again to remask unmasked values outside of area of interest

  var band_names = s1s2_stack.bandNames();

  // Sample the band values for each training data points
  // If reference signatures are already defined, it uses those signatures rather than sampling them again.
  var reference_signatures = reference_signatures || s1s2_stack
                                                     .sampleRegions({
                                                       collection: signatures,
                                                       properties: [CLASS_NAME],
                                                       scale : 10,
                                                       geometries: true
                                                     });

  // Wrapper function for the DTW implementation, intended to iterate over each land cover/crop class
  // provided in the reference signatures
  var dtw_min_dist = function(key, val){
    key = ee.Number.parse(key);
    // Function to format the signatures to a DTW-ready EE array
    var training_data_list = DTW.prepareSignatures(reference_signatures,
                                                   CLASS_NAME,
                                                   key,
                                                   BAND_NO,
                                                   PATTERNS_LEN,
                                                   band_names);

    // Compute the class-wise DTW distance
    return ee.ImageCollection(DTW.DTWDist(training_data_list,
                                          ndvi_image_list,
                                          {patterns_no: val,
                                          band_no: BAND_NO,
                                          timeseries_len: TIMESERIES_LEN,
                                          patterns_len: PATTERNS_LEN,
                                          constraint_type: 'time-weighted',
                                          beta: BETA,
                                          alpha: ALPHA
                                          })
           ).min()
           .rename('dtw')
           // Add class band corresponding to the land cover/crop class computed.
           // This is useful/necessary to generate the hard classification map from the dissimilarity values
           .addBands(ee.Image(key).toByte().rename('band'));
  };

  // Prepare the input time series of images into a DTW-ready EE array
  var ndvi_image_list = ee.Image(DTW.prepareBands(s1s2_stack,
                                                  BAND_NO,
                                                  TIMESERIES_LEN,
                                                  band_names));

  // Map the DTW distance function over each land cover class
  var dtw_image_list = reference_signatures_agg.map(dtw_min_dist);

  // Turn image collection into an array
  var array = ee.ImageCollection(dtw_image_list.values()).toArray();

  // Sort array by the first band, keeping other bands
  var axes = {image:0, band:1};
  var sort = array.arraySlice(axes.band, 0, 1);  // select bands from index 0 to 1 (DTW dissimilarity score)
  var sorted = array.arraySort(sort);

  // Take the first image only
  var values = sorted.arraySlice(axes.image, 0, 1);

  // Convert back to an image
  var min = values.arrayProject([axes.band]).arrayFlatten([['dtw', 'band']]);

  // Extract the DTW dissimilarity score
  var dtw_score = min.select(0).rename('score_' + year);
  // Extract the DTW hard classification
  var dtw_class = min.select(1).rename('classification_' + year);

  // 1. DTW outputs (dissimilarity score + classification map), 2. reference signatures, 3. stack of bands used as input to DTW
  return [dtw_class.addBands(dtw_score), reference_signatures, s1s2_stack];
};

// Manually captured signatures from photo-interpretation of 6 land cover classes in Sennar, Sudan
var signatures = ee.FeatureCollection('users/soilwatch/Sudan/SennarSignatures')
                 .filterBounds(county.geometry());

// Create a dictionary mapping land cover class to number of reference signatures per sample
var reference_signatures_agg = signatures.aggregate_histogram(CLASS_NAME);
print('Number of reference signatures per land cover class:')
print(reference_signatures_agg);

// Class list
var classification_names = ['built-up',
                            'water',
                            'trees',
                            'bare soil',
                            'active cropland',
                            'rangelands',
                            'wetland',
                            'abandoned/long-term fallow cropland (> 3 years)',
                            'short-term fallow cropland (<= 3 years)'
                            ];

// Corresponding color palette for the mapped land cover/crop classes
var classification_palette = ['#d63000',
                              'blue',
                              '#4d8b22',
                              'grey',
                              'ebff65',
                              '#98ff00',
                              'purple',
                              '#00ce6f',
                              '#FFA33F'
                              ];

// DTW Dissimilarity score palette
var score_palette = palettes.colorbrewer.RdYlGn[9].reverse();

// Compute the DTW classification for the year 2020.
var dtw_outputs = DTWClassification('2020', 'COPERNICUS/S2');
var dtw = dtw_outputs[0]; // Extract the DTW dissimilarity score and hard classification
var reference_signatures = dtw_outputs[1]; // Extract the reference signatures for 2020.
var s1s2_stack = dtw_outputs[2]; // Extract the bands stack used as input for DTW.
var s1s2_list = ee.List([s1s2_stack]); // Convert input bands to list to enable appending data from other years

// Iterate over each land cover/crop class to compute the DTW for each
Object.keys(year_dict).forEach(function(i) {
  var dtw_outputs = DTWClassification(i, year_dict[i])
  dtw = dtw.addBands(dtw_outputs[0]);
  s1s2_list = s1s2_list.add(dtw_outputs[2]);
});

// Image Visualization Parameters for the multi-temporal ndvi composite
var imageVisParam = {bands: ["ndvi_5", "ndvi_3", "ndvi_1"],
                     gamma: 1,
                     max: 7000,
                     min: 1000,
                     opacity: 1
};

// The NDVI/EVI images are masked with the 2020 crop mask generated as part of the TCP project.
Map.addLayer(s1s2_stack, imageVisParam, 'ndvi stack 2020');

imageVisParam['bands'] = ["VV_5", "VV_3", "VV_1"];
Map.addLayer(s1s2_stack, imageVisParam, 'VV stack 2020');

// Load the pre-generated layers from assets as producing the data live takes a while to load and display, and may result in memory errors
var dtw = ee.Image('users/soilwatch/Sudan/dtw_sennar_s1s2_2017_2020'); // Comment out this line to produce live data
var dtw_score = dtw.select('score_2020');
var dtw_class = dtw.select('classification_2020');

// VITO Land Cover dataset is used to generate a cropland mask.
// The dataset uses a 3-year epoch of time series data to generate the dominant land cover at any given location.
// This gives us a good estimation of cropland extent for the period 2016-2019, able to disregard single year fallowing events.
var vito_lulc_crop = ee.ImageCollection("ESA/WorldCover/v100").filterBounds(county).mosaic().eq(40);

Map.addLayer(dtw_score,
             {palette: score_palette, min: 0, max: 20000},
             'DTW dissimilarity score (30 days signature, 30 days images)');

Map.addLayer(dtw_class.updateMask(crop_mask),
             {palette: classification_palette.slice(0, classification_palette.length-2), min: 1, max: CLASS_NO},
             'DTW classification (30 days signature, 30 days images)');

// Generate the additional abandoned/long-term fallowed classes and short-term fallowed classes
var dtw_class_strat = dtw_class.where(dtw_class.eq(6) // 3 previous years identified as rangelands to be considered abandoned
                                      .and(dtw.select('classification_2019').eq(6))
                                      .and(dtw.select('classification_2018').eq(6))
                                      .and(dtw.select('classification_2017').eq(6))
                                      .and(vito_lulc_crop),
                                      ee.Image(8))
                               .where(dtw_class.eq(6) // Current year labelled as rangeland, any of the 3 previous years labeled as active cropland
                                      .and(dtw.select('classification_2019').eq(5)
                                      .or(dtw.select('classification_2018').eq(5))
                                      .or(dtw.select('classification_2017').eq(5))
                                      ).and(vito_lulc_crop),
                                      ee.Image(9));

Map.addLayer(dtw_class_strat.updateMask(crop_mask),
             {palette: classification_palette, min: 1, max: CLASS_NO+2},
             'DTW classification: active/abandoned cropland distinction');

// Generate are histogram (in hectares) of the respective land cover classes for 2020
var area_histogram = ee.Dictionary(dtw_class_strat.updateMask(crop_mask).reduceRegion(
  {reducer: ee.Reducer.frequencyHistogram(),
  geometry:county.geometry(),
  scale: 100,
  maxPixels:1e13,
  tileScale: 4
  }).get('classification_2020')).values();

// Assign class names with each area value.
var area_list = ee.List([]);
for (var i = 0; i <= CLASS_NO+1; i++) {
  area_list = area_list.add(ee.Feature(null, {
        'area': area_histogram.get(i),
        'class': classification_names[i]
      }));
}

var options = {title: 'Land Cover Classes Distribution',
               colors: classification_palette,
               sliceVisibilityThreshold: 0 // Don't group small slices.
              };

// Create the summary pie chart of the land cover class distribution
var area_chart = ui.Chart.feature.byFeature({
    features: ee.FeatureCollection(area_list),
    xProperty: 'class',
    yProperties: ['area']
  })
  .setChartType('PieChart')
  .setOptions(options);

print(area_chart);

// Define arguments for animation function parameters.
var video_args = {
  'dimensions': 2500,
  'region': county.geometry(),
  'framesPerSecond': 1,
  'crs': 'EPSG:4326',
  'min': 1,
  'max': 9,
  'palette': classification_palette
}

// Function to export GIF
var generateGIF = function(img, band_name){
  print(band_name+' GIF')
  print(img.select(band_name).getVideoThumbURL(video_args));
}

// Export GIF of the DTW classification + DTW classification distinguishing between abandoned/active cropland
generateGIF(ee.ImageCollection(ee.List([dtw.select('classification_2020').rename('classification')])
               .add(dtw_class_strat.rename('classification'))), 'classification');

// Export Video of the temporal NDVI composites for each year produced
Export.video.toDrive({
  collection: ee.ImageCollection(s1s2_list.map(function(img){return ee.Image(img).divide(100).toByte()}))
              .select(["ndvi_5", "ndvi_3", "ndvi_1"]),
  description:'temporal ndvi rgb',
  region: county.geometry(),
  scale: 100,
  crs:'EPSG:4326',
  maxPixels: 1e13
});

// Plot signatures with corresponding color code to see their locations
var setPointProperties = function(f){
  var class_val = f.get("lc_class");
  var mapDisplayColors = ee.List(classification_palette);

  // use the class as index to lookup the corresponding display color
  return f.set({style: {color: mapDisplayColors.get(ee.Number(class_val).subtract(1))}})
}

// apply the function and view the results on map
var styled_td = signatures.map(setPointProperties);

// Add reference sample points location to the map with the appropriate color legend
Map.addLayer(styled_td.style({styleProperty: "style"}), {}, 'crop type sample points');

// Export classified data. This is recommended to get the full extent of the data generated and saved,
// so it can be explored and consulted seamlessly.
Export.image.toAsset({
  image: dtw,
  description:'dtw_sennar_s1s2_2017_2020',
  assetId: 'users/soilwatch/Sudan/dtw_sennar_s1s2_2017_2020',
  region: county.geometry(),
  crs: 'EPSG:4326',
  scale: 10,
  maxPixels:1e13
});

// Create a legend for the different crop types
// set position of panel
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

// Create legend title
var legendTitle = ui.Label({
  value: 'Legend',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});

// Add the title to the panel
legend.add(legendTitle);

// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {

      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor: color,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });

      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });

      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};

// Add color and and names
for (var i = 0; i <= CLASS_NO+1; i++) {
  legend.add(makeRow(classification_palette[i], classification_names[i]));
  }

// add legend to map (alternatively you can also print the legend to the console)
Map.add(legend);
