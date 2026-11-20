// Global Forest Change v1.13 loss and gain (%) on the 1/720 deg grid.
var india = ee.FeatureCollection('projects/YOUR_PROJECT/assets/india_pov_boundary').geometry();
var nf = ee.ImageCollection('projects/nature-trace/assets/forest_typology/natural_forest_2020_v1_0_collection').mosaic();
var nfU = nf.unmask(0);                                                // natural forest probability 0-250; 125 = 0.5
var h = ee.Image('UMD/hansen/global_forest_change_2025_v1_13');
var f2000 = h.select('treecover2000').gte(30);
var ly = h.select('lossyear');
var hans = ee.Image.cat([
  f2000,                                             // forest (>=30% canopy) in 2000
  f2000.and(ly.gte(1)),                              // loss 2001-2025
  f2000.and(ly.gte(1)).and(ly.lte(10)),              // loss 2001-2010
  f2000.and(ly.gte(11)).and(ly.lte(20)),             // loss 2011-2020
  f2000.and(ly.gte(21)).and(nfU.gte(125)),           // loss 2021-2025 in natural forest (2020 map)
  f2000.and(ly.gte(21)).and(nfU.lt(125)),            // loss 2021-2025 in other tree cover
  h.select('gain').eq(1),                            // gain 2000-2012
  h.select('gain').eq(1).and(nfU.lt(125))            // gain classed non natural in 2020
]).rename(['hF2000','hLoss0125','hLoss0110','hLoss1120','hLossNat2125','hLossOth2125','hGain','hGainNonNat'])
  .unmask(0).toFloat()
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024});
Export.image.toAsset({image: hans.multiply(100).round().toUint8().clip(india), description: 'india_hansen_154m_v113',
  assetId: 'projects/YOUR_PROJECT/assets/india_hansen_154m_v113',
  crs: 'EPSG:4326', crsTransform: [1/720, 0, 68.0, 0, -1/720, 37.5], region: india.bounds(), maxPixels: 1e10,
  pyramidingPolicy: {'.default': 'mean'}});
