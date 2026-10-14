/*
 * export_core_edits.groovy -- QuPath script: corrected core circles -> CSV for 02_tma_dearray.R
 *
 * For grids made by 02b_cellscape_grid_to_xenium.py (CellScape grid moved onto Xenium):
 * 1. Open the Xenium MORPHOLOGY image of the slide in QuPath
 *    (<xenium outs>/morphology_focus/ch0000_dapi.ome.tif, or morphology_focus.ome.tif).
 * 2. File > Import objects from file (or drag the file onto the image):
 *    <Scripts_tma>/config/core_grid/<slide_id>_grid_qupath.geojson
 *    Every core is a named circle; its class is its tissue status from 02b
 *    (both / xenium_only / cellscape_only / neither).
 * 3. Correct what is wrong:
 *      - move a circle: select it (double-click), then drag it with the Move tool
 *      - resize: draw a new ellipse (Ellipse tool) over the core, Set name = the core's
 *        name (e.g. 3-f, right-click > Annotations > Set properties), delete the old one
 *      - core without usable tissue on Xenium: set its class to "missing"
 *        (right-click > Set classification; create the class "missing" once if needed)
 *      - core that 02b called empty but has tissue: set its class to "both"
 *      - delete a circle: that core is left out completely
 *    Keep the names: they are how 02 knows which core (and patient) it is.
 * 4. Set SLIDE_ID and OUT_DIR below, then Automate > Script editor > Run.
 *    Writes <OUT_DIR>/<SLIDE_ID>_edited.csv. Put it in <Scripts_tma>/config/core_grid/ --
 *    02 then uses it instead of <slide_id>.csv (re-running 02b does not overwrite it).
 *
 * Columns written: core_name, x_um, y_um, diameter_um, missing
 */

def SLIDE_ID = 'TMA_ORGAN'                          // <- slide_id in slides.csv
def OUT_DIR  = System.getProperty('user.home')      // <- folder for the CSV

def cal = getCurrentServer().getPixelCalibration()
double px = cal.getAveragedPixelSizeMicrons()
if (Double.isNaN(px)) {
    px = 0.2125                                     // Xenium morphology image: 0.2125 um / pixel
    print 'Image has no pixel size -- using the Xenium default 0.2125 um/pixel.'
}
def NO_TISSUE = ['missing', 'cellscape_only', 'neither'] as Set

def cores = getAnnotationObjects().findAll { it.getName() != null && !it.getName().trim().isEmpty() }
def unnamed = getAnnotationObjects().size() - cores.size()
def names = cores.collect { it.getName().trim() }
def dups = names.findAll { n -> names.count(n) > 1 }.unique()
if (!dups.isEmpty())
    print 'WARNING: several annotations with the same name (fix before using): ' + dups.join(', ')

def out = new File(OUT_DIR, SLIDE_ID + '_edited.csv')
out.withPrintWriter { w ->
    w.println('core_name,x_um,y_um,diameter_um,missing')
    cores.each { c ->
        def roi = c.getROI()
        double d = 2 * Math.sqrt(roi.getArea() / Math.PI) * px          // equivalent circle diameter
        def cls = c.getPathClass() == null ? '' : c.getPathClass().getName().toLowerCase()
        w.println(String.format(Locale.US, '%s,%.2f,%.2f,%.2f,%s', c.getName().trim(),
                roi.getCentroidX() * px, roi.getCentroidY() * px, d, NO_TISSUE.contains(cls) ? 'TRUE' : 'FALSE'))
    }
}
print 'Wrote ' + cores.size() + ' cores (pixel size ' + px + ' um) -> ' + out.getAbsolutePath()
if (unnamed > 0) print 'Skipped ' + unnamed + ' annotations without a name.'
