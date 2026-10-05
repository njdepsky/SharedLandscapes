/**************************************************
 * 03_extract_counts.js - 10 m land-cover pixel counts per DGG cell (Earth Engine Code Editor)
 *
 * For every DGG cell, counts 10 m pixels (EPSG:4326, scale 10 m) in these classes:
 *   n_total       valid non-water pixels: WorldCover v200 (2021) != 80 AND >= 1 Dynamic World obs
 *   n_land        n_total minus persistent water (DW water >= 95% of obs, clusters >= 10 px)
 *                 and persistent snow/ice (DW snow >= 99% of obs, clusters >= 10 px)
 *   n_freshwater, n_ice       the excluded persistent water / snow-ice pixels
 *   n_crop, n_built, n_barren DW class frequency > ANTHRO_THRESHOLD % of observations (land only)
 *   n_pasture                 Global Pasture Watch cultivated grassland (land only)
 *   n_other_intensive_treecrop  Spatial Database of Planted Trees (SDPT) Version 2.0, class 2 (land only)
 *   n_nonhabitat  land pixels that are crop OR built OR pasture OR tree crop
 *   n_habitat     n_land - n_nonhabitat
 * Year: Dynamic World composite for YEAR.
 *
 * Inputs: one Earth Engine TABLE ASSET per DGG layer (zipped shapefiles from R/01_build_dgg.R),
 * all in ASSET_FOLDER. Only layers in the project domain (tile lat0 >= MIN_LAT) are processed.
 * The 'seqnum' property is a STRING (9-digit ids overflow shapefile numbers): parsed here.
 *
 * Outputs: CSVs in Google Drive folder EXPORT_FOLDER, one per layer and batch, columns
 *   seqnum16, n_total, n_land, n_freshwater, n_ice, n_barren, n_built, n_crop, n_pasture,
 *   n_other_intensive_treecrop, n_nonhabitat, n_habitat
 * Download them all into one folder; python/05_postprocess.py reads every file matching
 * GEE_INPUT_PATTERN in that folder (see python/config.py).
 *
 * Running: the Code Editor creates one task per layer x batch; start them from the Tasks tab
 * (process the layer list in slices with FIRST_LAYER / N_LAYERS to keep task lists manageable).
 * Cost: the band stack is computed per export, so cost scales with the number of exports and
 * the area each covers; keep N_BATCHES = 1 unless a layer fails on memory or time limits.
 * Use your own Cloud project so compute bills to your quota.
 **************************************************/

/*****************************
 * 1. Settings
 *****************************/
var YEAR = 2024;
var ASSET_FOLDER = 'projects/YOUR_PROJECT/assets/nri_dgg';   // one table asset per DGG layer
var LAYER_PREFIX = 'dgg_isea3h16_land_';
var MIN_LAT = -60;                // project domain: layers with tile lat0 >= MIN_LAT
var FIRST_LAYER = 0;              // process sorted domain layers [FIRST_LAYER, FIRST_LAYER + N_LAYERS)
var N_LAYERS = 50;
var N_BATCHES = 1;                // exports per layer (seqnum mod N_BATCHES); raise only if needed

var EXPORT_FOLDER = 'nri_dgg_counts';
var EXPORT_PREFIX = 'dgg_counts_' + YEAR;   // file: <prefix>_<layer>_b<batch>.csv

var THRESHOLD_PCT_WATER = 95;
var THRESHOLD_PCT_SNOW = 99;
var MIN_CLUSTER = 10;
var ANTHRO_THRESHOLD = 40;        // DW frequency threshold (%) for crop, built and barren
var EXTRACT_SCALE = 10;           // metres
var SEQNUM_FIELD = 'seqnum';      // STRING property in the layer assets

/*****************************
 * 2. Layers in the project domain
 *****************************/
function layerLat0(name) {
  // <prefix><n|s>LL_<e|w>LLL ; polar caps: <prefix>north_cap / south_cap
  if (name.indexOf('south_cap') >= 0) return -90;
  if (name.indexOf('north_cap') >= 0) return 80;
  var m = name.match(/_([ns])(\d{2})_[ew]\d{3}$/);
  if (!m) return null;
  return (m[1] === 'n' ? 1 : -1) * parseInt(m[2], 10);
}

var listing = ee.data.listAssets(ASSET_FOLDER, {pageSize: 1000});
var layers = (listing.assets || [])
  .filter(function(a) { return a.type === 'TABLE'; })
  .map(function(a) { return a.id || a.name; })
  .filter(function(id) {
    var name = id.split('/').pop();
    var lat0 = layerLat0(name);
    return name.indexOf(LAYER_PREFIX) === 0 && lat0 !== null && lat0 >= MIN_LAT;
  })
  .sort();
var selected = layers.slice(FIRST_LAYER, FIRST_LAYER + N_LAYERS);
print('Domain layers in ' + ASSET_FOLDER + ': ' + layers.length +
      ' | this run: ' + FIRST_LAYER + '..' + (FIRST_LAYER + selected.length - 1), selected);

/*****************************
 * 3. Datasets
 *****************************/
var startDate = YEAR + '-01-01';
var endDate = YEAR + '-12-31';

var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
  .filterDate(startDate, endDate)
  .select('label');

var worldCover = ee.Image('ESA/WorldCover/v200/2021').select('Map');

var gpw = ee.ImageCollection('projects/global-pasture-watch/assets/ggc-30m/v1/grassland_c')
  .filterDate('2022-01-01', '2023-01-01').first();
var cultivatedGrass = gpw.eq(1).unmask(0).toUint8()
  .reproject({crs: 'EPSG:4326', scale: EXTRACT_SCALE})
  .rename('pasture');

var treeCrop = ee.Image('projects/ee-duzhenrong02/assets/SDPTV3')
  .eq(2).unmask(0).toUint8()
  .reproject({crs: 'EPSG:4326', scale: EXTRACT_SCALE})
  .rename('other_intensive_treecrop');

/*****************************
 * 4. Dynamic World frequencies
 *****************************/
var totalObs = dw.count().rename('total_obs');
var hasObs = totalObs.gt(0);

var dwCounts = dw.map(function(img) {
  return ee.Image.cat([
    img.eq(0).rename('water'),
    img.eq(8).rename('snow'),
    img.eq(4).rename('crop'),
    img.eq(6).rename('built'),
    img.eq(7).rename('bare')
  ]).toUint8();
}).sum();

function pct(band) {
  return dwCounts.select(band).divide(totalObs).multiply(100).updateMask(hasObs);
}
var pctWater = pct('water');
var pctSnow = pct('snow');
var pctCrop = pct('crop');
var pctBuilt = pct('built');
var pctBare = pct('bare');

/*****************************
 * 5. Masks
 *****************************/
var oceanMask = worldCover.neq(80);   // WorldCover permanent water (incl. ocean) excluded

var waterMaskRaw = pctWater.gte(THRESHOLD_PCT_WATER).unmask(0);
var waterMask = waterMaskRaw.and(waterMaskRaw.connectedPixelCount(100, true).gte(MIN_CLUSTER)).not();

var snowMaskRaw = pctSnow.gte(THRESHOLD_PCT_SNOW).unmask(0);
var snowMask = snowMaskRaw.and(snowMaskRaw.connectedPixelCount(100, true).gte(MIN_CLUSTER)).not();

var validMask = oceanMask.and(hasObs);
var landMask = waterMask.and(snowMask).and(oceanMask);

/*****************************
 * 6. Binary layers
 *****************************/
var isTotal = ee.Image(1).rename('n_total').toUint8().updateMask(validMask);
var isLand = ee.Image(1).rename('n_land').toUint8().updateMask(validMask.and(landMask));
var isFreshwater = waterMask.not().and(snowMask).rename('n_freshwater').toUint8().unmask(0).updateMask(validMask);
var isIce = snowMask.not().rename('n_ice').toUint8().unmask(0).updateMask(validMask);
var landOnly = validMask.and(landMask);
var isCrop = pctCrop.gt(ANTHRO_THRESHOLD).rename('n_crop').toUint8().unmask(0).updateMask(landOnly);
var isBuilt = pctBuilt.gt(ANTHRO_THRESHOLD).rename('n_built').toUint8().unmask(0).updateMask(landOnly);
var isBarren = pctBare.gt(ANTHRO_THRESHOLD).rename('n_barren').toUint8().unmask(0).updateMask(landOnly);
var isPasture = cultivatedGrass.rename('n_pasture').toUint8().unmask(0).updateMask(landOnly);
var isTreeCrop = treeCrop.rename('n_other_intensive_treecrop').toUint8().unmask(0).updateMask(landOnly);
var isNonhabitat = isCrop.or(isBuilt).or(isPasture).or(isTreeCrop)
  .rename('n_nonhabitat').toUint8().unmask(0).updateMask(landOnly);

var countStack = ee.Image.cat([
  isTotal, isLand, isFreshwater, isIce, isBarren, isBuilt, isCrop, isPasture, isTreeCrop, isNonhabitat
]);

var SELECTORS = [
  'seqnum16', 'n_total', 'n_land', 'n_freshwater', 'n_ice', 'n_barren', 'n_built', 'n_crop',
  'n_pasture', 'n_other_intensive_treecrop', 'n_nonhabitat', 'n_habitat'
];

/*****************************
 * 7. Extraction
 *****************************/
function extract(hexes) {
  return countStack.reduceRegions({
    collection: hexes,
    reducer: ee.Reducer.sum().unweighted(),
    scale: EXTRACT_SCALE,
    crs: 'EPSG:4326',
    tileScale: 8
  }).map(function(f) {
    function num(key) { return ee.Number(ee.Algorithms.If(f.get(key), f.get(key), 0)); }
    var nLand = num('n_land');
    var nNonhabitat = num('n_nonhabitat');
    return ee.Feature(null, {
      'seqnum16': ee.Number.parse(f.get(SEQNUM_FIELD)),
      'n_total': num('n_total'),
      'n_land': nLand,
      'n_freshwater': num('n_freshwater'),
      'n_ice': num('n_ice'),
      'n_barren': num('n_barren'),
      'n_built': num('n_built'),
      'n_crop': num('n_crop'),
      'n_pasture': num('n_pasture'),
      'n_other_intensive_treecrop': num('n_other_intensive_treecrop'),
      'n_nonhabitat': nNonhabitat,
      'n_habitat': nLand.subtract(nNonhabitat).max(0)
    });
  });
}

/*****************************
 * 8. Exports (one per layer x batch)
 *****************************/
selected.forEach(function(assetId) {
  var layer = assetId.split('/').pop();
  var hexes = ee.FeatureCollection(assetId).map(function(f) {
    return f.set('batch_id', ee.Number.parse(f.get(SEQNUM_FIELD)).mod(N_BATCHES));
  });
  for (var b = 0; b < N_BATCHES; b++) {
    var name = EXPORT_PREFIX + '_' + layer + '_b' + b;
    Export.table.toDrive({
      collection: extract(hexes.filter(ee.Filter.eq('batch_id', b))),
      description: name.slice(0, 100),
      folder: EXPORT_FOLDER,
      fileNamePrefix: name,
      fileFormat: 'CSV',
      selectors: SELECTORS
    });
  }
});

/*****************************
 * 9. Sanity check before starting the tasks
 *****************************/
if (selected.length) {
  var sample = extract(ee.FeatureCollection(selected[0]).limit(20));
  print('Sample rows (first selected layer)', sample.limit(5));
}
