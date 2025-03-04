package earth.worldwind.render.image

import android.graphics.Bitmap
import earth.worldwind.geom.Sector
import earth.worldwind.ogc.GpkgBitmapFactory
import earth.worldwind.util.DownloadPostprocessor
import earth.worldwind.util.Level
import earth.worldwind.util.Tile

actual open class ImageTile actual constructor(
    sector: Sector, level: Level, row: Int, column: Int
): Tile(sector, level, row, column), DownloadPostprocessor<Bitmap> {
    actual var imageSource: ImageSource? = null
    actual var cacheSource: ImageSource? = null

    override fun process(resource: Bitmap): Bitmap {
        val source = cacheSource?.asUnrecognized()
        if (source is GpkgBitmapFactory) source.saveBitmap(resource)
        return resource
    }
}