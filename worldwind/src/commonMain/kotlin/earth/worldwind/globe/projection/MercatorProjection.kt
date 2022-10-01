package earth.worldwind.globe.projection

import earth.worldwind.geom.*
import earth.worldwind.geom.Sector.Companion.fromDegrees
import earth.worldwind.globe.Globe
import earth.worldwind.util.Logger
import kotlin.math.*

open class MercatorProjection : GeographicProjection {
    override val displayName = "Mercator"
    override val is2D = true
    private val projectionLimits = fromDegrees(-78.0, -180.0, 156.0, 360.0)

    override fun geographicToCartesian(globe: Globe, latitude: Angle, longitude: Angle, altitude: Double, result: Vec3): Vec3 {
        val lat = when {
            latitude.degrees > projectionLimits.maxLatitude.degrees -> projectionLimits.maxLatitude
            latitude.degrees < projectionLimits.minLatitude.degrees -> projectionLimits.minLatitude
            else -> latitude
        }
        val lon = when {
            longitude.degrees > projectionLimits.maxLongitude.degrees -> projectionLimits.maxLongitude
            longitude.degrees < projectionLimits.minLongitude.degrees -> projectionLimits.minLongitude
            else -> longitude
        }

        // See "Map Projections: A Working Manual", page 44 for the source of the below formulas.
        val ecc = sqrt(globe.eccentricitySquared)
        val sinPhi = sin(lat.radians)
        val s = (1 + sinPhi) / (1 - sinPhi) * ((1 - ecc * sinPhi) / (1 + ecc * sinPhi)).pow(ecc)
        result.x = globe.equatorialRadius * lon.radians //+ origin?.x ?: 0.0
        result.y = 0.5 * globe.equatorialRadius * ln(s)
        result.z = altitude
        return result
    }

    override fun geographicToCartesianNormal(globe: Globe, latitude: Angle, longitude: Angle, result: Vec3): Vec3 {
        TODO("Not yet implemented")
    }

    override fun geographicToCartesianTransform(
        globe: Globe, latitude: Angle, longitude: Angle, altitude: Double, result: Matrix4
    ): Matrix4 {
        TODO("Not yet implemented")
    }

    override fun geographicToCartesianGrid(
        globe: Globe, sector: Sector, numLat: Int, numLon: Int, height: FloatArray?, verticalExaggeration: Float,
        origin: Vec3?, result: FloatArray, offset: Int, rowStride: Int
    ): FloatArray {
        require(numLat >= 1 && numLon >= 1) {
            Logger.logMessage(
                Logger.ERROR, "MercatorProjection", "geographicToCartesianGrid",
                "Number of latitude or longitude locations is less than one"
            )
        }
        require(height == null || height.size >= numLat * numLon) {
            Logger.logMessage(Logger.ERROR, "MercatorProjection", "geographicToCartesianGrid", "missingArray")
        }

        val eqr = globe.equatorialRadius
        val ecc = sqrt(globe.eccentricitySquared)
        val minLat = sector.minLatitude.radians
        val maxLat = sector.maxLatitude.radians
        val minLon = sector.minLongitude.radians
        val maxLon = sector.maxLongitude.radians
        val deltaLat = (maxLat - minLat) / if (numLat > 1) numLat - 1 else 1
        val deltaLon = (maxLon - minLon) / if (numLon > 1) numLon - 1 else 1
        val minLatLimit = projectionLimits.minLatitude.radians
        val maxLatLimit = projectionLimits.maxLatitude.radians
        val minLonLimit = projectionLimits.minLongitude.radians
        val maxLonLimit = projectionLimits.maxLongitude.radians
        var elevIndex = 0
        val xOffset = origin?.x ?: 0.0
        val yOffset = origin?.y ?: 0.0
        val zOffset = origin?.z ?: 0.0

        // Iterate over the latitude and longitude coordinates in the specified sector, computing the Cartesian point
        // corresponding to each latitude and longitude.
        var rowIndex = offset
        val stride = if (rowStride == 0) numLon * 3 else rowStride
        var lat = minLat
        for (latIndex in 0 until numLat) {
            if (latIndex == numLat - 1) lat = maxLat // explicitly set the last lat to the max latitude to ensure alignment
            lat = lat.coerceIn(minLatLimit, maxLatLimit) // limit lat to projection limits

            // Latitude is constant for each row. Values that are a function of latitude can be computed once per row.
            val sinLat = sin(lat)
            val s = (1 + sinLat) / (1 - sinLat) * ((1 - ecc * sinLat) / (1 + ecc * sinLat)).pow(ecc)
            val y = eqr * ln(s) * 0.5 - yOffset
            
            var lon = minLon
            var colIndex = rowIndex
            for (lonIndex in 0 until numLon) {
                if (lonIndex == numLon - 1) lon = maxLon // explicitly set the last lon to the max longitude to ensure alignment
                lon = lon.coerceIn(minLonLimit, maxLonLimit) // limit lon to projection limits
                result[colIndex++] = (eqr * lon - xOffset).toFloat()
                result[colIndex++] = y.toFloat()
                result[colIndex++] = if (height != null) (height[elevIndex++] * verticalExaggeration - zOffset).toFloat() else 0f
                lon += deltaLon
            }
            rowIndex += stride
            lat += deltaLat
        }

        return result
    }

    override fun geographicToCartesianBorder(
        globe: Globe, sector: Sector, numLat: Int, numLon: Int, height: Float, origin: Vec3?, result: FloatArray
    ): FloatArray {
        TODO("Not yet implemented")
    }

    override fun cartesianToGeographic(globe: Globe, x: Double, y: Double, z: Double, result: Position): Position {
        // See "Map Projections: A Working Manual", pages 45 and 19 for the source of the below formulas.
        val ecc2 = globe.eccentricitySquared
        val ecc4 = ecc2 * ecc2
        val ecc6 = ecc4 * ecc2
        val ecc8 = ecc6 * ecc2
        val t = E.pow(-y / globe.equatorialRadius)

        val A = PI / 2 - 2 * atan(t)
        val B = ecc2 / 2 + 5 * ecc4 / 24 + ecc6 / 12 + 13 * ecc8 / 360
        val C = 7 * ecc4 / 48 + 29 * ecc6 / 240 + 811 * ecc8 / 11520
        val D = 7 * ecc6 / 120 + 81 * ecc8 / 1120
        val E = 4279 * ecc8 / 161280

        val Ap = A - C + E
        val Bp = B - 3 * D
        val Cp = 2 * C - 8 * E
        val Dp = 4 * D
        val Ep = 8 * E

        val s2p = sin(2 * A)

        val lat = Ap + s2p * (Bp + s2p * (Cp + s2p * (Dp + Ep * s2p)))

        return result.setRadians(lat, x / globe.equatorialRadius, z)
    }

    override fun cartesianToLocalTransform(globe: Globe, x: Double, y: Double, z: Double, result: Matrix4): Matrix4 {
        TODO("Not yet implemented")
    }

    override fun intersect(globe: Globe, line: Line, result: Vec3): Boolean {
        TODO("Not yet implemented")
    }
}