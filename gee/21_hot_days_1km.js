// Daytime LST on hot pre-monsoon days (ERA5-Land Tmax >= local 90th/95th percentile), MYD11A1, 2019-2023.
var era = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').filterDate('2019-01-01', '2024-01-01')
  .filter(ee.Filter.calendarRange(3, 5, 'month')).select('temperature_2m_max');
var pct = era.reduce(ee.Reducer.percentile([90, 95])).rename(['p90', 'p95']);
var lst = ee.ImageCollection('MODIS/061/MYD11A1').filterDate('2019-01-01', '2024-01-01')
  .filter(ee.Filter.calendarRange(3, 5, 'month'));
var joined = ee.Join.saveFirst('era').apply(lst, era, ee.Filter.maxDifference({difference: 12 * 3600 * 1000,
  leftField: 'system:time_start', rightField: 'system:time_start'}));
var daily = ee.ImageCollection(joined).map(function(i){
  i = ee.Image(i);
  var qc = i.select('QC_Day');
  var q = qc.bitwiseAnd(3).lte(1).and(qc.rightShift(6).bitwiseAnd(3).lte(1));
  var t = i.select('LST_Day_1km').multiply(0.02).subtract(273.15).updateMask(q);
  var ta = ee.Image(i.get('era'));
  var h90 = ta.gte(pct.select('p90')), h95 = ta.gte(pct.select('p95'));
  return ee.Image.cat([t.updateMask(h90).rename('LSTd_hot90'), t.updateMask(h95).rename('LSTd_hot95'), t.rename('LSTd_allMAM')]);
});
var m = daily.mean();
var n = daily.count().rename(['n_hot90', 'n_hot95', 'n_all']);
var region = ee.Geometry.Rectangle([68.0, 6.5, 97.5, 37.5], null, false);
Export.image.toAsset({image: m.addBands(n).toFloat(), description: 'india_1km_hotdays',
  assetId: 'projects/YOUR_PROJECT/assets/india_1km_hotdays', crs: 'SR-ORG:6974', scale: 926.625433055833,
  region: region, maxPixels: 1e10, pyramidingPolicy: {'.default': 'mean'}});
