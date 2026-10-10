/*
 * export_tma_grid.groovy -- QuPath script: TMA grid -> CSV for 02_tma_dearray.R
 *
 * 1. Open the Xenium MORPHOLOGY image of one slide in QuPath
 *    (<xenium outs>/morphology_focus/ch0000_dapi.ome.tif, or morphology_focus.ome.tif).
 *    Only this image has the same coordinates as the Xenium cells. An H&E or a
 *    CellScape image would give coordinates that do not match.
 * 2. TMA > TMA dearrayer: set the core diameter (e.g. 0.6 mm), the column and row
 *    labels so the names match tma_map.csv (e.g. column labels 1-10, row labels
 *    a-t, or the other way round), run it. Then check every core: drag
 *    misplaced ones, add/remove rows/columns (right-click > TMA), and set empty
 *    positions to "missing".
 * 3. Set SLIDE_ID below to the slide_id in slides.csv and OUT_DIR to a folder,
 *    then Automate > Script editor > Run.
 * 4. Copy the CSV to <Scripts_tma>/config/core_grid/<slide_id>.csv and run 02.
 *
 * Columns written: core_name, x_um, y_um, diameter_um, missing
 * (centre and diameter in microns = Xenium coordinates).
 */

def SLIDE_ID = 'TMA_ORGAN'                          // <- slide_id in slides.csv
def OUT_DIR  = System.getProperty('user.home')      // <- folder for the CSV

def grid = getCurrentHierarchy().getTMAGrid()
if (grid == null) {
    print 'No TMA grid in this image -- run TMA > TMA dearrayer first.'
    return
}
def cal = getCurrentServer().getPixelCalibration()
double px = cal.getAveragedPixelSizeMicrons()
if (Double.isNaN(px)) {
    px = 0.2125                                     // Xenium morphology image: 0.2125 um / pixel
    print 'Image has no pixel size -- using the Xenium default 0.2125 um/pixel.'
}
def out = new File(OUT_DIR, SLIDE_ID + '.csv')
out.withPrintWriter { w ->
    w.println('core_name,x_um,y_um,diameter_um,missing')
    grid.getTMACoreList().each { core ->
        def roi = core.getROI()
        double d = 0.5 * (roi.getBoundsWidth() + roi.getBoundsHeight()) * px
        w.println(String.format(Locale.US, '%s,%.2f,%.2f,%.2f,%s', core.getName(),
                roi.getCentroidX() * px, roi.getCentroidY() * px, d, core.isMissing() ? 'TRUE' : 'FALSE'))
    }
}
print 'Wrote ' + grid.getTMACoreList().size() + ' cores (pixel size ' + px + ' um) -> ' + out.getAbsolutePath()
