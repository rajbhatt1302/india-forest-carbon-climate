// Land cover with a natural forest threshold of 0.52 instead of 0.50 (sensitivity analysis).
var india = ee.FeatureCollection('projects/YOUR_PROJECT/assets/india_pov_boundary').geometry();
var nfU = ee.ImageCollection('projects/nature-trace/assets/forest_typology/natural_forest_2020_v1_0_collection').mosaic().unmask(0);
var wc = ee.ImageCollection('ESA/WorldCover/v200').first();
var tree = wc.eq(10);
var lc10 = ee.Image.cat([tree.and(nfU.gte(130)), tree.and(nfU.lt(130)), wc.eq(40), wc.eq(50),
                         wc.eq(80).or(wc.eq(90)).or(wc.eq(95))])
  .rename(['fNat','fPlant','fCrop','fBuilt','fWater']).toFloat()
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024});
Export.image.toAsset({image: lc10.multiply(100).round().toUint8().clip(india), description: 'india_lc_154m_nf052',
  assetId: 'projects/YOUR_PROJECT/assets/india_lc_154m_nf052',
  crs: 'EPSG:4326', crsTransform: [1/720, 0, 68.0, 0, -1/720, 37.5], region: india.bounds(), maxPixels: 1e10,
  pyramidingPolicy: {'.default': 'mean'}});
