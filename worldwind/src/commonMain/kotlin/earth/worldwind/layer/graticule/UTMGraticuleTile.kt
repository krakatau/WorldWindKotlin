package earth.worldwind.layer.graticule

import earth.worldwind.geom.Position
import earth.worldwind.geom.Sector
import earth.worldwind.geom.coords.Hemisphere
import earth.worldwind.geom.coords.UTMCoord
import earth.worldwind.layer.graticule.AbstractUTMGraticuleLayer.Companion.UTM_MAX_LATITUDE
import earth.worldwind.layer.graticule.AbstractUTMGraticuleLayer.Companion.UTM_MIN_LATITUDE
import earth.worldwind.layer.graticule.GridElement.Companion.TYPE_GRIDZONE_LABEL
import earth.worldwind.layer.graticule.GridElement.Companion.TYPE_LINE
import earth.worldwind.layer.graticule.UTMGraticuleLayer.Companion.UTM_ZONE_RESOLUTION
import earth.worldwind.render.RenderContext
import earth.worldwind.shape.PathType

internal class UTMGraticuleTile(layer: UTMGraticuleLayer, sector: Sector, private val zone: Int): AbstractGraticuleTile(layer, sector) {
    private val hemisphere = if (sector.centroidLatitude.degrees > 0) Hemisphere.N else Hemisphere.S
    private var squares: MutableList<UTMSquareZone>? = null
    override val layer get() = super.layer as UTMGraticuleLayer

    override fun selectRenderables(rc: RenderContext) {
        super.selectRenderables(rc)

        // Select tile grid elements
        val graticuleType = layer.getTypeFor(UTM_ZONE_RESOLUTION)
        for (ge in gridElements!!) if (ge.isInView(rc)) layer.addRenderable(ge.renderable, graticuleType)
        if (getSizeInPixels(rc) / 10 < MIN_CELL_SIZE_PIXELS * 2) return

        // Select child elements
        squares ?: createSquares()
        for (sz in squares!!) if (sz.isInView(rc)) sz.selectRenderables(rc) else sz.clearRenderables()
    }

    override fun clearRenderables() {
        super.clearRenderables()
        squares?.run {
            for (sz in this) sz.clearRenderables()
            clear()
        }.also { squares = null }
    }

    private fun createSquares() {
        // Find grid zone easting and northing boundaries
        var utm = UTMCoord.fromLatLon(sector.minLatitude, sector.centroidLongitude)
        val minNorthing = utm.northing
        utm = UTMCoord.fromLatLon(sector.maxLatitude, sector.centroidLongitude)
        var maxNorthing = utm.northing
        maxNorthing = if (maxNorthing == 0.0) 10e6 else maxNorthing
        utm = UTMCoord.fromLatLon(sector.minLatitude, sector.minLongitude)
        var minEasting = utm.easting
        utm = UTMCoord.fromLatLon(sector.maxLatitude, sector.minLongitude)
        minEasting = utm.easting.coerceAtMost(minEasting)
        val maxEasting = 1e6 - minEasting

        // Create squares
        squares = layer.createSquaresGrid(zone, hemisphere, sector, minEasting, maxEasting, minNorthing, maxNorthing)
    }

    override fun createRenderables() {
        super.createRenderables()
        val positions = mutableListOf(
            Position(sector.minLatitude, sector.minLongitude, 0.0),
            Position(sector.maxLatitude, sector.minLongitude, 0.0)
        )
        var polyline = layer.createLineRenderable(positions.toList(), PathType.LINEAR)
        var lineSector = Sector(sector.minLatitude, sector.maxLatitude, sector.minLongitude, sector.minLongitude)
        gridElements!!.add(GridElement(lineSector, polyline, TYPE_LINE, sector.minLongitude))

        // Generate south parallel at south pole and equator
        if (sector.minLatitude.degrees == UTM_MIN_LATITUDE || sector.minLatitude.degrees == 0.0) {
            positions.clear()
            positions.add(Position(sector.minLatitude, sector.minLongitude, 0.0))
            positions.add(Position(sector.minLatitude, sector.maxLongitude, 0.0))
            polyline = layer.createLineRenderable(ArrayList(positions), PathType.LINEAR)
            lineSector = Sector(sector.minLatitude, sector.minLatitude, sector.minLongitude, sector.maxLongitude)
            gridElements!!.add(GridElement(lineSector, polyline, TYPE_LINE, sector.minLatitude))
        }

        // Generate north parallel at north pole
        if (sector.maxLatitude.degrees == UTM_MAX_LATITUDE) {
            positions.clear()
            positions.add(Position(sector.maxLatitude, sector.minLongitude, 0.0))
            positions.add(Position(sector.maxLatitude, sector.maxLongitude, 0.0))
            polyline = layer.createLineRenderable(ArrayList(positions), PathType.LINEAR)
            lineSector = Sector(sector.maxLatitude, sector.maxLatitude, sector.minLongitude, sector.maxLongitude)
            gridElements!!.add(GridElement(lineSector, polyline, TYPE_LINE, sector.maxLatitude))
        }

        // Add label
        if (hasLabel()) {
            val text = layer.createTextRenderable(
                Position(sector.centroidLatitude, sector.centroidLongitude, 0.0),
                zone.toString() + hemisphere, 10e6
            )
            gridElements!!.add(GridElement(sector, text, TYPE_GRIDZONE_LABEL))
        }
    }

    private fun hasLabel(): Boolean {
        // Has label if it contains hemisphere mid latitude
        val southLat = UTM_MIN_LATITUDE / 2.0
        val southLabel = sector.minLatitude.degrees < southLat && southLat <= sector.maxLatitude.degrees
        val northLat = UTM_MAX_LATITUDE / 2.0
        val northLabel = (sector.minLatitude.degrees < northLat && northLat <= sector.maxLatitude.degrees)
        return southLabel || northLabel
    }

    companion object {
        private const val MIN_CELL_SIZE_PIXELS = 40 // TODO: make settable
    }
}