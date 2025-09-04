

/// Retrieve ET from a filtered area during the crop growing season


/////////// add info below

// Load crop area feature collection
var crop_country = ee.FeatureCollection("projects/crops-461709/assets/kyrg_wheat_only_usda"); // path to your crop area

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

var ndvi = 0.25; // minimum NDVI value used as a threshold to select areas with active crops within cropland regions
var max_cloud_cover = 5;  // This sets the cloud cover threshold for retrieving Sentinel-2 imagery used in NDVI calculation during the mid-season.
// A value of 5% is used here, but sometimes no images are available for the selected period.// In that case, you may need to increase this value.

//////////////////////////////////////////////////////////


// Load cropland mask (ESA WorldCover)
var croplandMask = ee.ImageCollection('ESA/WorldCover/v200')
  .first()
  .select('Map')
  .eq(40);

// Adjust start and end dates to match the correct dekadal periods
function getDekadStart(date) {
  var d = ee.Date(date);
  var day = ee.Number(d.get('day'));
  var month = ee.Number(d.get('month'));
  var year = ee.Number(d.get('year'));

  return ee.Algorithms.If(day.lte(5),
    d.update({'day': 1}),
    ee.Algorithms.If(day.lte(15),
      d.update({'day': 11}),
      ee.Algorithms.If(day.lte(25),
        d.update({'day': 21}),
        ee.Algorithms.If(month.eq(12),
          ee.Date.fromYMD(year.add(1), 1, 1),
          ee.Date.fromYMD(year, month.add(1), 1)
        )
      )
    )
  );
}

function getDekadEnd(date) {
  var d = ee.Date(date);
  var day = ee.Number(d.get('day'));
  var month = ee.Number(d.get('month'));
  var year = ee.Number(d.get('year'));

  return ee.Algorithms.If(day.lte(5),
    ee.Date.fromYMD(year, month, 1).advance(-1, 'month').update({'day': 21}),
    ee.Algorithms.If(day.lte(14),
      d.update({'day': 1}),
      ee.Algorithms.If(day.lte(24),
        d.update({'day': 11}),
        d.update({'day': 21})
      )
    )
  );
}

// Loop over each year
var results = years.map(function(year) {
  year = ee.Number(year);

  // Define crop season dates
  
  var rawStart = ee.Date.fromYMD(year, cropStartMonth, cropStartDay);

  // If the end month is before the start month, assume the crop season ends in the next year
  var endYearAdjusted = ee.Algorithms.If(
    cropEndMonth < cropStartMonth,
    year.add(1),
    year
  );
  var rawEnd = ee.Date.fromYMD(endYearAdjusted, cropEndMonth, cropEndDay);
  
  var startDate = ee.Date(getDekadStart(rawStart)); // First day of the season adjusted to match ET data from WaPOR dekadal periods
  var endDate = ee.Date(getDekadEnd(rawEnd)); // First day of the season adjusted to match ET data from WaPOR dekadal periods

  var ndviStartDate = ee.Date.fromYMD(year, ndviMonth, ndviStartDay);
  var ndviEndDate = ee.Date.fromYMD(year, ndviMonth, ndviEndDay);

  // Load Sentinel-2 imagery, apply cloud masking, and select areas with NDVI > 0.25 (or other informed value)
  var s2ndvi = ee.ImageCollection('COPERNICUS/S2_SR')
    .filterDate(ndviStartDate, ndviEndDate)
    .filterBounds(crop_country)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', max_cloud_cover))
    .map(function(img) {
      var qa = img.select('QA60'); // source: https://gis.stackexchange.com/questions/445556/different-methods-for-masking-clouds-of-sentinel-2-images-in-gee
      var cloudFree = qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0));
      return img.updateMask(cloudFree);
    })
    .map(function(img) {
      return img.normalizedDifference(['B8','B4']).rename('NDVI')
                .copyProperties(img, ['system:time_start']);
    })
    .map(function(img) {
      return img.updateMask(croplandMask).copyProperties(img, ['system:time_start']);
    });

  var meanNDVI = s2ndvi.mean();
  var growingMask = meanNDVI.gte(ndvi);

  // Load WaPOR ET data and filter it for the crop area using the adjusted start and end dates the crop area for the adjusted start date and end dates end date
  var waporET = ee.ImageCollection("projects/UNFAO/wapor/v3/L1-AETI-D")
    .filterDate(startDate, endDate)
    .filterBounds(crop_country);

  var totalET = waporET
    .map(function(img) {
      return img.updateMask(growingMask);
    })
    .sum()
    .clip(crop_country);

  var bandName = totalET.bandNames().get(0);
  
  // calculates the mean ET value
  var meanET = totalET.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: crop_country.geometry(),
    scale: 300, // same as ET WaPOR data (300m resolution)
    maxPixels: 1e13
  }).get(bandName);

  return ee.Dictionary({
    'year': year,
    'start': rawStart.format('yyyy-MM-dd'),
    'end': rawEnd.format('yyyy-MM-dd'),
    'mean_ET': meanET
  });
});

// Print mean ET for the selected crop area per year
print('Average ET per year:', results);




