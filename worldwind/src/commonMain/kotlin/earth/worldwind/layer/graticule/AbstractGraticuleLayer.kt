package earth.worldwind.layer.graticule

import earth.worldwind.geom.*
import earth.worldwind.geom.Angle.Companion.fromDegrees
import earth.worldwind.geom.Angle.Companion.normalizeLatitude
import earth.worldwind.geom.Angle.Companion.normalizeLongitude
import earth.worldwind.geom.Angle.Companion.toDegrees
import earth.worldwind.layer.AbstractLayer
import earth.worldwind.render.Color
import earth.worldwind.render.Font
import earth.worldwind.render.RenderContext
import earth.worldwind.render.Renderable
import earth.worldwind.shape.Label
import earth.worldwind.shape.Path
import earth.worldwind.shape.PathType
import kotlin.math.abs
import kotlin.math.sign

/**
 * Displays a graticule.
 */
abstract class AbstractGraticuleLayer(name: String): AbstractLayer(name) {
    override var isPickEnabled = false

    // Helper variables to avoid memory leaks
    private val surfacePoint = Vec3()
    private val forwardRay = Line()
    private val lookAtPoint = Vec3()
    private val lookAtPos = Position()
    private val graticuleSupport = GraticuleSupport()

    // Update reference states
    private val lastCameraPoint = Vec3()
    private var lastCameraHeading = 0.0
    private var lastCameraTilt = 0.0
    private var lastFOV = 0.0
    private var lastVerticalExaggeration = 0.0

//    private var lastGlobe: Globe? = null
//    private var lastProjection: GeographicProjection? = null
//    private var frameTimeStamp = Instant.DISTANT_PAST // used only for 2D continuous globes to determine whether render is in same frame
//    private var terrainConformance = 50.0

    companion object {
//        /**
//         * Solid line rendering style. This style specifies that a line will be drawn without any breaks. <br></br>
//         * <pre>`_________`</pre>
//         * <br></br> is an example of a solid line.
//         */
//        val LINE_STYLE_SOLID = GraticuleRenderingParams.VALUE_LINE_STYLE_SOLID
//        /**
//         * Dashed line rendering style. This style specifies that a line will be drawn as a series of long strokes, with
//         * space in between. <br></br>
//         * <pre>`- - - - -`</pre>
//         * <br></br> is an example of a dashed line.
//         */
//        val LINE_STYLE_DASHED = GraticuleRenderingParams.VALUE_LINE_STYLE_DASHED
//        /**
//         * Dotted line rendering style. This style specifies that a line will be drawn as a series of evenly spaced "square"
//         * dots. <br></br>
//         * <pre>`. . . . .`</pre>
//         * is an example of a dotted line.
//         */
//        val LINE_STYLE_DOTTED = GraticuleRenderingParams.VALUE_LINE_STYLE_DOTTED
        private const val LOOK_AT_LATITUDE_PROPERTY = "look_at_latitude"
        private const val LOOK_AT_LONGITUDE_PROPERTY = "look_at_longitude"
        private const val GRATICULE_PIXEL_SIZE_PROPERTY = "graticule_pixel_size"
        private const val GRATICULE_LABEL_OFFSET_PROPERTY = "graticule_label_offset"
    }

    init {
        this.initRenderingParams()
    }

    protected abstract fun initRenderingParams()

    /**
     * Returns whether or not graticule lines will be rendered.
     *
     * @param key the rendering parameters key.
     *
     * @return true if graticule lines will be rendered; false otherwise.
     */
    fun isDrawGraticule(key: String) = getRenderingParams(key).isDrawLines

    /**
     * Sets whether or not graticule lines will be rendered.
     *
     * @param drawGraticule true to render graticule lines; false to disable rendering.
     * @param key           the rendering parameters key.
     */
    fun setDrawGraticule(drawGraticule: Boolean, key: String) { getRenderingParams(key).isDrawLines = drawGraticule }

    /**
     * Returns the graticule line Color.
     *
     * @param key the rendering parameters key.
     *
     * @return Color used to render graticule lines.
     */
    fun getGraticuleLineColor(key: String) = getRenderingParams(key).lineColor

    /**
     * Sets the graticule line Color.
     *
     * @param color Color that will be used to render graticule lines.
     * @param key   the rendering parameters key.
     */
    fun setGraticuleLineColor(color: Color, key: String) { getRenderingParams(key).lineColor = color }

    /**
     * Returns the graticule line width.
     *
     * @param key the rendering parameters key.
     *
     * @return width of the graticule lines.
     */
    fun getGraticuleLineWidth(key: String) = getRenderingParams(key).lineWidth

    /**
     * Sets the graticule line width.
     *
     * @param lineWidth width of the graticule lines.
     * @param key       the rendering parameters key.
     */
    fun setGraticuleLineWidth(lineWidth: Double, key: String) { getRenderingParams(key).lineWidth = lineWidth }

//    /**
//     * Returns the graticule line rendering style.
//     *
//     * @param key the rendering parameters key.
//     *
//     * @return rendering style of the graticule lines.
//     */
//    fun getGraticuleLineStyle(key: String) = getRenderingParams(key).lineStyle
//
//    /**
//     * Sets the graticule line rendering style.
//     *
//     * @param lineStyle rendering style of the graticule lines. One of LINE_STYLE_SOLID, LINE_STYLE_DASHED, or
//     * LINE_STYLE_DOTTED.
//     * @param key       the rendering parameters key.
//     */
//    fun setGraticuleLineStyle(lineStyle: String, key: String) { getRenderingParams(key).lineStyle = lineStyle }

    /**
     * Returns whether or not graticule labels will be rendered.
     *
     * @param key the rendering parameters key.
     *
     * @return true if graticule labels will be rendered; false otherwise.
     */
    fun isDrawLabels(key: String) = getRenderingParams(key).isDrawLabels

    /**
     * Sets whether or not graticule labels will be rendered.
     *
     * @param drawLabels true to render graticule labels; false to disable rendering.
     * @param key        the rendering parameters key.
     */
    fun setDrawLabels(drawLabels: Boolean, key: String) { getRenderingParams(key).isDrawLabels = drawLabels }

    /**
     * Returns the graticule label Color.
     *
     * @param key the rendering parameters key.
     *
     * @return Color used to render graticule labels.
     */
    fun getLabelColor(key: String) = getRenderingParams(key).labelColor

    /**
     * Sets the graticule label Color.
     *
     * @param color Color that will be used to render graticule labels.
     * @param key   the rendering parameters key.
     */
    fun setLabelColor(color: Color, key: String) { getRenderingParams(key).labelColor = color }

    /**
     * Returns the Font used for graticule labels.
     *
     * @param key the rendering parameters key.
     *
     * @return Font used to render graticule labels.
     */
    fun getLabelFont(key: String) = getRenderingParams(key).labelFont

    /**
     * Sets the Font used for graticule labels.
     *
     * @param font Font that will be used to render graticule labels.
     * @param key  the rendering parameters key.
     */
    fun setLabelFont(font: Font, key: String) { getRenderingParams(key).labelFont = font }

    fun getRenderingParams(key: String) = graticuleSupport.getRenderingParams(key)

    fun setRenderingParams(key: String, renderingParams: GraticuleRenderingParams) {
        graticuleSupport.setRenderingParams(key, renderingParams)
    }

    fun addRenderable(renderable: Renderable, paramsKey: String) { graticuleSupport.addRenderable(renderable, paramsKey) }

    private fun removeAllRenderables() { graticuleSupport.removeAllRenderables() }

    public override fun doRender(rc: RenderContext) {
//        if (rc.isContinuous2DGlobe) {
//            if (needsToUpdate(rc)) {
//                clear(rc)
//                selectRenderables(rc)
//            }
//
//            // If the frame time stamp is the same, then this is the second or third pass of the same frame. We continue
//            // selecting renderables in these passes.
//            if (rc.frameTimeStamp === frameTimeStamp) selectRenderables(rc)
//
//            frameTimeStamp = rc.frameTimeStamp
//        } else {
        if (needsToUpdate(rc)) {
            clear(rc)
            selectRenderables(rc)
        }
//        }

        // Render
        graticuleSupport.render(rc, opacity)
    }

    /**
     * Select the visible grid elements
     *
     * @param rc the current `RenderContext`.
     */
    protected abstract fun selectRenderables(rc: RenderContext)
    protected abstract val orderedTypes: List<String>
    abstract fun getTypeFor(resolution: Double): String

    /**
     * Determines whether the grid should be updated. It returns true if:   * the eye has moved more than 1% of its
     * altitude above ground  * the view FOV, heading or pitch have changed more than 1 degree  * vertical
     * exaggeration has changed  `RenderContext`.
     *
     * @return true if the graticule should be updated.
     */
    private fun needsToUpdate(rc: RenderContext): Boolean {
        if (lastVerticalExaggeration != rc.verticalExaggeration) return true
        if (abs(lastCameraHeading - rc.camera!!.heading.degrees) > 1) return true
        if (abs(lastCameraTilt - rc.camera!!.tilt.degrees) > 1) return true
        if (abs(lastFOV - rc.camera!!.fieldOfView.degrees) > 1) return true
        return rc.cameraPoint.distanceTo(lastCameraPoint) > computeAltitudeAboveGround(rc) / 100

        // We must test the globe and its projection to see if either changed. We can't simply use the globe state
        // key for this because we don't want a 2D globe offset change to cause an update. Offset changes don't
        // invalidate the current set of renderables.
//        if (rc.globe != lastGlobe) return true
//        if (rc.is2DGlobe) if ((rc.globe as Globe2D).projection != lastProjection) return true
    }

    protected open fun clear(rc: RenderContext) {
        removeAllRenderables()
        lastCameraPoint.copy(rc.cameraPoint)
        lastFOV = rc.camera!!.fieldOfView.degrees
        lastCameraHeading = rc.camera!!.heading.degrees
        lastCameraTilt = rc.camera!!.tilt.degrees
        lastVerticalExaggeration = rc.verticalExaggeration
//        lastGlobe = rc.globe
//        if (rc.is2DGlobe) lastProjection = (rc.globe as Globe2D).projection
//        terrainConformance = computeTerrainConformance(rc)
//        applyTerrainConformance()
    }

//    private fun computeTerrainConformance(rc: RenderContext): Double {
//        var value = 100
//        val alt = rc.camera!!.position.altitude
//        when {
//            alt < 10e3 -> value = 20
//            alt < 50e3 -> value = 30
//            alt < 100e3 -> value = 40
//            alt < 1000e3 -> value = 60
//        }
//        return value.toDouble()
//    }
//
//    private fun applyTerrainConformance() {
//        val graticuleType = getOrderedTypes()
//        for (type in graticuleType) {
//            getRenderingParams(type)[GraticuleRenderingParams.KEY_LINE_CONFORMANCE] = terrainConformance
//        }
//    }

    fun computeLabelOffset(rc: RenderContext): Location {
        return if (hasLookAtPos(rc)) {
            val labelOffsetDegrees = getLabelOffset(rc)
            val labelPos = Location(
                getLookAtLatitude(rc).minusDegrees(labelOffsetDegrees),
                getLookAtLongitude(rc).minusDegrees(labelOffsetDegrees)
            )
            labelPos.setDegrees(
                normalizeLatitude(labelPos.latitude.degrees).coerceIn(-70.0, 70.0),
                normalizeLongitude(labelPos.longitude.degrees)
            )
            labelPos
        } else rc.camera!!.position
    }

    fun createLineRenderable(positions: List<Position>, pathType: PathType) =
        Path(positions).apply {
            this.pathType = pathType
            isFollowTerrain = true
            // terrainConformance = 1.0 // TODO Why not terrainConformance?
            altitudeMode = AltitudeMode.CLAMP_TO_GROUND
        }

    fun createTextRenderable(position: Position, label: String, resolution: Double) =
        Label(position, label).apply {
            altitudeMode = AltitudeMode.CLAMP_TO_GROUND
            // priority = resolution * 1e6 // TODO Implement priority
        }

    fun hasLookAtPos(rc: RenderContext): Boolean {
        calculateLookAtProperties(rc)
        return rc.getUserProperty(LOOK_AT_LATITUDE_PROPERTY) != null && rc.getUserProperty(LOOK_AT_LONGITUDE_PROPERTY) != null
    }

    fun getLookAtLatitude(rc: RenderContext): Angle {
        calculateLookAtProperties(rc)
        return rc.getUserProperty(LOOK_AT_LATITUDE_PROPERTY) as Angle
    }

    fun getLookAtLongitude(rc: RenderContext): Angle {
        calculateLookAtProperties(rc)
        return rc.getUserProperty(LOOK_AT_LONGITUDE_PROPERTY) as Angle
    }

    fun getPixelSize(rc: RenderContext): Double {
        calculateLookAtProperties(rc)
        return rc.getUserProperty(GRATICULE_PIXEL_SIZE_PROPERTY) as Double
    }

    private fun getLabelOffset(rc: RenderContext): Double {
        calculateLookAtProperties(rc)
        return rc.getUserProperty(GRATICULE_LABEL_OFFSET_PROPERTY) as Double
    }

    fun getSurfacePoint(rc: RenderContext, latitude: Angle, longitude: Angle): Vec3 {
        if (!rc.terrain!!.surfacePoint(latitude, longitude, surfacePoint))
            rc.globe!!.geographicToCartesian(
                latitude, longitude, rc.globe!!.getElevationAtLocation(latitude, longitude)
                        * rc.verticalExaggeration, surfacePoint
            )
        return surfacePoint
    }

    fun computeAltitudeAboveGround(rc: RenderContext): Double {
        val surfacePoint = getSurfacePoint(rc, rc.camera!!.position.latitude, rc.camera!!.position.longitude)
        return rc.cameraPoint.distanceTo(surfacePoint)
    }

    fun computeTruncatedSegment(p1: Position, p2: Position, sector: Sector, positions: MutableList<Position>) {
        val p1In = sector.contains(p1.latitude, p1.longitude)
        val p2In = sector.contains(p2.latitude, p2.longitude)
        if (!p1In && !p2In) return  // whole segment is (likely) outside
        if (p1In && p2In) {
            // whole segment is (likely) inside
            positions.add(p1)
            positions.add(p2)
        } else {
            // segment does cross the boundary
            var outPoint = if (!p1In) p1 else p2
            val inPoint = if (p1In) p1 else p2
            for (i in 1..2) {
                // there may be two intersections
                var intersection: Location? = null
                if (outPoint.longitude.degrees > sector.maxLongitude.degrees
                    || sector.maxLongitude.degrees == 180.0 && outPoint.longitude.degrees < 0.0) {
                    // intersect with east meridian
                    intersection = greatCircleIntersectionAtLongitude(
                        inPoint, outPoint, sector.maxLongitude
                    )
                } else if (outPoint.longitude.degrees < sector.minLongitude.degrees
                    || sector.minLongitude.degrees == -180.0 && outPoint.longitude.degrees > 0.0) {
                    // intersect with west meridian
                    intersection = greatCircleIntersectionAtLongitude(
                        inPoint, outPoint, sector.minLongitude
                    )
                } else if (outPoint.latitude.degrees > sector.maxLatitude.degrees) {
                    // intersect with top parallel
                    intersection = greatCircleIntersectionAtLatitude(
                        inPoint, outPoint, sector.maxLatitude
                    )
                } else if (outPoint.latitude.degrees < sector.minLatitude.degrees) {
                    // intersect with bottom parallel
                    intersection = greatCircleIntersectionAtLatitude(
                        inPoint, outPoint, sector.minLatitude
                    )
                }
                outPoint = if (intersection != null) Position(
                    intersection.latitude,
                    intersection.longitude,
                    outPoint.altitude
                ) else break
            }
            positions.add(inPoint)
            positions.add(outPoint)
        }
    }

    /**
     * Computes the intersection point position between a great circle segment and a meridian.
     *
     * @param p1        the great circle segment start position.
     * @param p2        the great circle segment end position.
     * @param longitude the meridian longitude `Angle`
     *
     * @return the intersection `Position` or null if there was no intersection found.
     */
    private fun greatCircleIntersectionAtLongitude(p1: Location, p2: Location, longitude: Angle): Location? {
        if (p1.longitude == longitude) return p1
        if (p2.longitude == longitude) return p2
        var pos: Location? = null
        val deltaLon = getDeltaLongitude(p1, p2.longitude)
        if (getDeltaLongitude(p1, longitude) < deltaLon && getDeltaLongitude(p2, longitude) < deltaLon) {
            var count = 0
            val precision = 1.0 / 6378137.0 // 1m angle in radians
            var a = p1
            var b = p2
            var midPoint = greatCircleMidPoint(a, b)
            while (getDeltaLongitude(midPoint, longitude).radians > precision && count <= 20) {
                count++
                if (getDeltaLongitude(a, longitude) < getDeltaLongitude(b, longitude)) b = midPoint else a = midPoint
                midPoint = greatCircleMidPoint(a, b)
            }
            pos = midPoint
        }
        // Adjust final longitude for an exact match
        if (pos != null) pos = Location(pos.latitude, longitude)
        return pos
    }

    /**
     * Computes the intersection point position between a great circle segment and a parallel.
     *
     * @param p1       the great circle segment start position.
     * @param p2       the great circle segment end position.
     * @param latitude the parallel latitude `Angle`
     *
     * @return the intersection `Position` or null if there was no intersection found.
     */
    private fun greatCircleIntersectionAtLatitude(p1: Location, p2: Location, latitude: Angle): Location? {
        var pos: Location? = null
        if (sign(p1.latitude.degrees - latitude.degrees) != sign(p2.latitude.degrees - latitude.degrees)) {
            var count = 0
            val precision = 1.0 / 6378137.0 // 1m angle in radians
            var a = p1
            var b = p2
            var midPoint = greatCircleMidPoint(a, b)
            while (abs(midPoint.latitude.radians - latitude.radians) > precision && count <= 20) {
                count++
                if (sign(a.latitude.degrees - latitude.degrees) != sign(midPoint.latitude.degrees - latitude.degrees))
                    b = midPoint else a = midPoint
                midPoint = greatCircleMidPoint(a, b)
            }
            pos = midPoint
        }
        // Adjust final latitude for an exact match
        if (pos != null) pos = Location(latitude, pos.longitude)
        return pos
    }

    private fun greatCircleMidPoint(p1: Location, p2: Location): Location {
        val azimuth = p1.greatCircleAzimuth(p2)
        val distance = p1.greatCircleDistance(p2)
        return p1.greatCircleLocation(azimuth, distance / 2, Location())
    }

    private fun getDeltaLongitude(p1: Location, longitude: Angle): Angle {
        val deltaLon = abs(p1.longitude.degrees - longitude.degrees)
        return fromDegrees(if (deltaLon < 180) deltaLon else 360 - deltaLon)
    }

    private fun calculateLookAtProperties(rc: RenderContext) {
        if (!rc.hasUserProperty(LOOK_AT_LATITUDE_PROPERTY) || !rc.hasUserProperty(LOOK_AT_LONGITUDE_PROPERTY)) {
            //rc.modelview.extractEyePoint(forwardRay.origin)
            forwardRay.origin.copy(rc.cameraPoint)
            rc.modelview.extractForwardVector(forwardRay.direction)
            val range = if (rc.globe!!.intersect(forwardRay, lookAtPoint)) {
                rc.globe!!.cartesianToGeographic(
                    lookAtPoint.x, lookAtPoint.y, lookAtPoint.z, lookAtPos
                )
                rc.putUserProperty(LOOK_AT_LATITUDE_PROPERTY, lookAtPos.latitude)
                rc.putUserProperty(LOOK_AT_LONGITUDE_PROPERTY, lookAtPos.longitude)
                lookAtPoint.distanceTo(rc.cameraPoint)
            } else {
                rc.removeUserProperty(LOOK_AT_LATITUDE_PROPERTY)
                rc.removeUserProperty(LOOK_AT_LONGITUDE_PROPERTY)
                rc.horizonDistance
            }
            val pixelSizeMeters = rc.pixelSizeAtDistance(range)
            rc.putUserProperty(GRATICULE_PIXEL_SIZE_PROPERTY, pixelSizeMeters)
            val pixelSizeDegrees = toDegrees(pixelSizeMeters / rc.globe!!.equatorialRadius)
            rc.putUserProperty(
                GRATICULE_LABEL_OFFSET_PROPERTY, pixelSizeDegrees * rc.viewport.width / 4
            )
        }
    }
}