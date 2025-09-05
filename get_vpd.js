/////////////////// add info below //////////////////////

//The vpd values will differ from the report because I used a different method to mask the clouds. Now, both the solar radiation and VPD scripts use the same method. Also, NDVI in the old VPD script was calculated for only one year


// source: https://gis.stackexchange.com/questions/436605/vapour-pressure-deficit-vpd-calculation-from-era5-google-earth-engine

// Load your polygon FeatureCollection
var polygons = ee.FeatureCollection("projects/crops-461709/assets/ind_cotton_only_usda_kharif_may_oct"); // path to your crop area

// Define year range
var startYear = 2019; // start year for the ET retrieval 
var endYear = 2023; // end year for the ET retrieval
var years = ee.List.sequence(startYear, endYear);

// Crop season period
var cropStartDay = 20; // start day of the crop season (ESA WorldCereal https://esa-worldcereal.org/en/products/crop-calendars or USDA Crop Explorer for non cereals)
var cropStartMonth = 4; // start month of the crop season

var cropEndDay = 6; // end day of the crop season
var cropEndMonth = 10; // start month of the crop season

// mid-season period (used to retrieve areas with NDVI > 0.25), adjust according to each crop type
var ndviMonth = 7;  // mid-season month for the selected crop
var ndviStartDay = 1; // start day of the mid-season month
var ndviEndDay = 30; // end day of the mid-season month

var ndvi_value = 0.25; // minimum NDVI value used as a threshold to select areas with active crops within cropland regions
var max_cloud_cover = 70;  // This sets the cloud cover threshold for retrieving Sentinel-2 imagery used in NDVI calculation during the mid-season.

/////////////////////////////////


var mergedGeometry = polygons.union().geometry();

// Load cropland mask (ESA WorldCover)
var croplandMask = ee.ImageCollection('ESA/WorldCover/v200')
  .first()
  .select('Map')
  .eq(40);

// compute vpd for each year
var yearlyStats = years.map(function(year) {
  year = ee.Number(year);

  var start = ee.Date.fromYMD(year, cropStartMonth, cropStartDay);

  // Adjust end year if crop season crosses calendar years
  var endYearAdjusted = ee.Algorithms.If(
    ee.Number(cropEndMonth).lt(cropStartMonth)
      .or(ee.Number(cropEndMonth).eq(cropStartMonth)
      .and(ee.Number(cropEndDay).lt(cropStartDay))),
    year.add(1),
    year
  );
  var end = ee.Date.fromYMD(endYearAdjusted, cropEndMonth, cropEndDay);

  // Adjust NDVI year if NDVI month is before crop start month
  var ndviYear = year.add(ee.Number(ndviMonth).lt(cropStartMonth));
  var ndviStart = ee.Date.fromYMD(ndviYear, ndviMonth, ndviStartDay);
  var ndviEnd = ee.Date.fromYMD(ndviYear, ndviMonth, ndviEndDay);

// Cloud masking using SCL band for Sentinel-2 SR
function maskS2clouds(image) {
  var scl = image.select('SCL');
  var cloudFree = scl.eq(4) // Vegetation
    .or(scl.eq(5))           // Not Vegetated
    .or(scl.eq(6))           // Water
    .or(scl.eq(7))           // Unclassified
    .or(scl.eq(1));          // Saturated/Defective (optional)

  return image.updateMask(cloudFree)
              .normalizedDifference(['B8', 'B4']).rename('NDVI');
}

  // Load Sentinel-2 imagery, apply cloud masking, and select areas with NDVI > 0.25 (or other informed value)
  var s2ndvi = ee.ImageCollection('COPERNICUS/S2_SR')
  .filterDate(ndviStart, ndviEnd)
  .filterBounds(mergedGeometry)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', max_cloud_cover))
  .map(maskS2clouds);
  

  var meanNDVI = s2ndvi.mean();
  var ndviMask = meanNDVI.gte(ndvi_value);


  var combinedMask = croplandMask.and(ndviMask);

// Load ERA5-Land data
// source: https://gis.stackexchange.com/questions/436605/vapour-pressure-deficit-vpd-calculation-from-era5-google-earth-engine
// https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics?tab=overview
// https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR
  
  var era5 = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
    .select(['temperature_2m', 'dewpoint_temperature_2m'])
    .filterDate(start, end);

  var vpdCollection = era5.map(function(img) {
    var Ta = img.select('temperature_2m');
    var Td = img.select('dewpoint_temperature_2m');

    var rh100 = img.expression(
      'exp((17.269 * (Td - 273.15)) / (273.3 + (Td - 273.15)) - (17.269 * (Ta - 273.15)) / (237.3 + (Ta - 273.15)))',
      {Td: Td, Ta: Ta}
    ).rename('rh100');

    var TaC = Ta.subtract(273.15);
    var es = TaC.multiply(17.27)
                .divide(TaC.add(237.3))
                .exp()
                .multiply(0.6108)
                .rename('es');

    var vpd = es.multiply(ee.Image(1).subtract(rh100)).rename('vpd');
    return vpd.set('system:time_start', img.get('system:time_start'));
  });

  var meanVPD = vpdCollection.mean()
    .updateMask(combinedMask)
    .clip(polygons);

  //var mergedGeometry = polygons.union().geometry();

  var meanVPDvalue = meanVPD.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: mergedGeometry,
    scale: 10000,  // appropriate for ERA5
    maxPixels: 1e13
  }).get('vpd');


  // Return a Feature with year and VPD value
  return ee.Feature(null, {
    'year': year,
    'meanVPD': meanVPDvalue
  });
});

// Convert to FeatureCollection and print
var vpdResults = ee.FeatureCollection(yearlyStats);

print('VPD Results by Year', vpdResults);


/* old version - used in the project 

// source: https://gis.stackexchange.com/questions/436605/vapour-pressure-deficit-vpd-calculation-from-era5-google-earth-engine

// date range
var start = ee.Date('2022-01-13');
var end = ee.Date('2022-08-03');

var cloudcover = 5;

// NDVI filtering period
var ndviMonth = 4;  
var ndviStartDay = 1;
var ndviEndDay = 30;
var ndviYear = 2022; 

var polygons = ee.FeatureCollection("projects/crops-461709/assets/uzb_wheat_only_usda_17166_wint_01-13_08-03");

var croplandMask = ee.ImageCollection('ESA/WorldCover/v200')
  .first()
  .select('Map')
  .eq(40);  // 40 = cropland

var ndviStartDate = ee.Date.fromYMD(ndviYear, ndviMonth, ndviStartDay);
var ndviEndDate = ee.Date.fromYMD(ndviYear, ndviMonth, ndviEndDay);

// Load Sentinel-2 imagery, apply cloud masking, and select areas with NDVI > 0.25 (or other informed value)
var s2ndvi = ee.ImageCollection('COPERNICUS/S2_SR')
  .filterDate(ndviStartDate, ndviEndDate)
  .filterBounds(polygons)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloudcover))
  .map(function(img) {
    var qa = img.select('QA60');
    var cloudFree = qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0));
    return img.updateMask(cloudFree);
  })
  .map(function(img) {
    return img.normalizedDifference(['B8','B4']).rename('NDVI')
              .copyProperties(img, ['system:time_start']);
  });

var meanNDVI = s2ndvi.mean();
var growingMask = meanNDVI.gte(0.25);

var combinedMask = croplandMask.and(growingMask);


// Load ERA5-Land daily data

// https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics?tab=overview
// https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR

var era5 = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")    //('ECMWF/ERA5/DAILY')
  .select(['temperature_2m', 'dewpoint_temperature_2m'])    //(['mean_2m_air_temperature', 'dewpoint_2m_temperature'])
  .filterDate(start, end);

// Compute RH/100 and VPD for each image
// tem in Kelvin https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics?tab=overview
var vpdCollection = era5.map(function(img) {
  var Ta = img.select('temperature_2m'); // in Kelvin 
  var Td = img.select('dewpoint_temperature_2m'); // in Kelvin

  // Compute RH/100
  // T is in Celsius in the formula https://bmcnoldy.earth.miami.edu/Humidity.html
  var rh100 = img.expression(
    'exp((17.269 * (Td - 273.15)) / (273.3 + (Td - 273.15)) - (17.269 * (Ta - 273.15)) / (237.3 + (Ta - 273.15)))',
    {Td: Td, Ta: Ta} 
  ).rename('rh100');

  // saturation vapor pressure (es) in kPa 
  // chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf unit:  mb/hPa 
  var TaC = Ta.subtract(273.15);  
  
  // saturation vapor pressure (es) in kPa
  var es = TaC.multiply(17.27)
              .divide(TaC.add(237.3))
              .exp()
              .multiply(0.6108)
              .rename('es');
  
  //  VPD = es * (1 - RH)
  var vpd = es.multiply(ee.Image(1).subtract(rh100)).rename('vpd');

  return vpd.set('system:time_start', img.get('system:time_start'));
});

//  mean VPD over the season

var meanVPD = vpdCollection.mean()
  .updateMask(combinedMask)
  .clip(polygons);


Map.addLayer(combinedMask.updateMask(combinedMask), {palette: ['white', 'green']}, 'Cropland + NDVI ≥ 0.25');


// Visualization
var vpdVis = {min: 0, max: 3, palette: ['blue', 'green', 'yellow', 'red']};
Map.centerObject(combinedMask, 6);
Map.addLayer(meanVPD, vpdVis, 'Mean VPD (kPa)');

print('Mean VPD image:', meanVPD);

var mergedGeometry = polygons.union().geometry();

var meanVPDvalue = meanVPD.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: mergedGeometry,
  scale: 1000,
  maxPixels: 1e13
}).get('vpd');

print('Mean VPD:', meanVPDvalue);


*/
