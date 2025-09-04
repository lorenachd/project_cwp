/////////////////// add info below //////////////////////


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

// compute solar radiation for each year
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
  

  /* not working, but works in the ET code
  var s2ndvi = ee.ImageCollection('COPERNICUS/S2_SR')
    .filterDate(ndviStart, ndviEnd)
    .filterBounds(mergedGeometry)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', max_cloud_cover))
    .map(function(img) {
      var qa = img.select('QA60'); // source: https://gis.stackexchange.com/questions/445556/different-methods-for-masking-clouds-of-sentinel-2-images-in-gee
      var cloudFree = qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0));
      return img.updateMask(cloudFree)
                .normalizedDifference(['B8', 'B4']).rename('NDVI');
    }); 
    
  */  

  var meanNDVI = s2ndvi.mean();
  var ndviMask = meanNDVI.gte(ndvi_value);


  var combinedMask = croplandMask.and(ndviMask);

  // Load ERA5 solar radiation and apply mask # https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR
  var era5 = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
                .filterDate(start, end);

  var solar = era5.select('surface_net_solar_radiation_sum')
                  .sum()
                  .rename('solar_rad')
                  .updateMask(combinedMask)
                  .clip(mergedGeometry);

  var meanSolar = solar.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: mergedGeometry,
    scale: 1000,
    maxPixels: 1e13
  }).get('solar_rad');

  var stats = ee.Feature(null, {
    'year': year,
    'start_date': start.format('YYYY-MM-dd'),
    'end_date': end.format('YYYY-MM-dd'),
    'mean_solar_radiation': meanSolar
  });

  return ee.FeatureCollection([stats]);
});

// print mean surface net short-wave solar radiation for each year for the selected area
var allStats = ee.FeatureCollection(yearlyStats).flatten();
print('yearly mean solar radiation for crop area:', allStats);

