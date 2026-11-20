// Land cover fractions (%) on a 1/720 deg grid. WorldCover 2021 trees are split into natural forest and
// plantation/other trees with the Natural Forests of the World 2020 map (p >= 0.5).
var india = ee.FeatureCollection('projects/YOUR_PROJECT/assets/india_pov_boundary').geometry();
var nf = ee.ImageCollection('projects/nature-trace/assets/forest_typology/natural_forest_2020_v1_0_collection').mosaic();
var wc = ee.ImageCollection('ESA/WorldCover/v200').first();          // 2021, 10 m
var tree = wc.eq(10);
var nfU = nf.unmask(0);                                                // natural forest probability 0-250; 125 = 0.5
var lc10 = ee.Image.cat([tree.and(nfU.gte(125)), tree.and(nfU.lt(125)), wc.eq(40), wc.eq(50),
                         wc.eq(80).or(wc.eq(90)).or(wc.eq(95))])
  .rename(['fNat','fPlant','fCrop','fBuilt','fWater']).toFloat()
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024});
var h = ee.Image('UMD/hansen/global_forest_change_2024_v1_12');
var f2000 = h.select('treecover2000').gte(30);
var ly = h.select('lossyear');
var hans = ee.Image.cat([
  f2000,                                             // forest (>=30% canopy) in 2000
  f2000.and(ly.gte(1)),                              // loss 2001-2024
  f2000.and(ly.gte(1)).and(ly.lte(10)),              // loss 2001-2010
  f2000.and(ly.gte(11)).and(ly.lte(20)),             // loss 2011-2020
  f2000.and(ly.gte(21)).and(nfU.gte(125)),           // loss 2021-2024 in natural forest (2020 map)
  f2000.and(ly.gte(21)).and(nfU.lt(125)),            // loss 2021-2024 in other tree cover
  h.select('gain').eq(1),                            // gain 2000-2020
  h.select('gain').eq(1).and(nfU.lt(125))            // gain classed non natural in 2020
]).rename(['hF2000','hLoss0124','hLoss0110','hLoss1120','hLossNat2124','hLossOth2124','hGain','hGainNonNat'])
  .unmask(0).toFloat()
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024});
var stack = lc10.addBands(hans).multiply(100).round().toUint8().clip(india);
Export.image.toAsset({image: stack, description: 'india_lc_154m', assetId: 'projects/YOUR_PROJECT/assets/india_lc_154m',
  crs: 'EPSG:4326', crsTransform: [1/720, 0, 68.0, 0, -1/720, 37.5], region: india.bounds(), maxPixels: 1e10,
  pyramidingPolicy: {'.default': 'mean'}});
