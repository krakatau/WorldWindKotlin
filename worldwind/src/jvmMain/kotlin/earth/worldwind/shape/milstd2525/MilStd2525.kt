package earth.worldwind.shape.milstd2525

import ArmyC2.C2SD.RendererPluginInterface.ISinglePointInfo
import ArmyC2.C2SD.Utilities.*
import sec.web.renderer.SinglePoint2525Renderer
import java.awt.Font

/**
 * This utility class generates MIL-STD-2525 symbols and tactical graphics using the MIL-STD-2525 Symbol Rendering Library
 * @see <a href="https://github.com/missioncommand/mil-sym-java">https://github.com/missioncommand/mil-sym-java</a>
 */
actual object MilStd2525 {
    const val GRAPHICS_OUTLINE_WIDTH = 3f
    const val DEFAULT_PIXEL_SIZE = 35

    actual var outlineWidth = GRAPHICS_OUTLINE_WIDTH
    /**
     * The actual rendering engine for the MIL-STD-2525 graphics.
     */
    private val renderer = SinglePoint2525Renderer()

    /**
    * Initializes the static MIL-STD-2525 symbol renderer.
    */
    init {
        // Establish the default rendering values.
        // See: https://github.com/missioncommand/mil-sym-android/blob/master/Renderer/src/main/java/armyc2/c2sd/renderer/utilities/RendererSettings.java
        val rs = RendererSettings.getInstance()
        rs.symbologyStandard = RendererSettings.Symbology_2525C

        // Depending on screen size and DPI you may want to change the font size.
        rs.setLabelFont("Arial", Font.BOLD, 8)
        rs.setMPLabelFont("Arial", Font.BOLD, 12)

        // Configure modifier text output
        rs.textBackgroundMethod = RendererSettings.TextBackgroundMethod_OUTLINE
        rs.textOutlineWidth = 4

        // Configure Single point symbol outline width
        rs.singlePointSymbolOutlineWidth = 0
    }

    /**
     * Creates an MIL-STD-2525 symbol from the specified symbol code, modifiers and attributes.
     *
     * @param symbolCode The MIL-STD-2525 symbol code.
     * @param modifiers  The MIL-STD-2525 modifiers. If null, a default (empty) modifier list will be used.
     * @param attributes The MIL-STD-2525 attributes. If null, a default (empty) attribute list will be used.
     *
     * @return An ImageInfo object containing the symbol's bitmap and metadata; may be null
     */
    @JvmStatic
    fun renderImage(symbolCode: String, modifiers: Map<String, String>?, attributes: Map<String, String>?): ISinglePointInfo? {
        val params = mutableMapOf<String, String>()
        if (modifiers != null) {
            params.putAll(modifiers)
        }
        if (attributes != null) {
            params.putAll(attributes)
        }
        return renderer.render(symbolCode, params)
    }

    /**
     * Get symbol text description and hierarchy reference.
     */
    @JvmStatic
    fun getSymbolDef(sidc: String): SymbolDef? = SymbolDefTable.getInstance()
        .getSymbolDef(SymbolUtilities.getBasicSymbolID(sidc), RendererSettings.getInstance().symbologyStandard)

    /**
     * Get symbol visual attributes like draw category, min points, max points, etc.
     */
    @JvmStatic
    fun getUnitDef(sidc: String): UnitDef? = UnitDefTable.getInstance()
        .getUnitDef(SymbolUtilities.getBasicSymbolID(sidc), RendererSettings.getInstance().symbologyStandard)

    @JvmStatic
    actual fun getSimplifiedSymbolID(sidc: String) =
        setAffiliation(SymbolUtilities.getBasicSymbolID(sidc), sidc.substring(1, 2))

    @JvmStatic
    actual fun setAffiliation(sidc: String, affiliation: String?) =
        // Weather symbols has no affiliation
        if (sidc.length >= 2 && !SymbolUtilities.isWeather(sidc) && affiliation != null && affiliation.length == 1) {
            val result = sidc.substring(0, 1) + affiliation.uppercase() + sidc.substring(2)
            if (SymbolUtilities.hasValidAffiliation(result)) result else sidc
        } else sidc

    @JvmStatic
    actual fun setStatus(sidc: String, status: String?) =
        // Weather symbols has no status
        if (sidc.length >= 4 && !SymbolUtilities.isWeather(sidc) && status != null && status.length == 1) {
            val result = sidc.substring(0, 3) + status.uppercase() + sidc.substring(4)
            if (SymbolUtilities.hasValidStatus(result)) result else sidc
        } else sidc

    @JvmStatic
    actual fun setEchelon(sidc: String, echelon: String?): String {
        val isTG = SymbolUtilities.isTacticalGraphic(sidc)
        return if (sidc.length >= 12 && (isTG && SymbolUtilities.canSymbolHaveModifier(sidc, ModifiersTG.B_ECHELON)
                    || !isTG && SymbolUtilities.canUnitHaveModifier(sidc, ModifiersUnits.B_ECHELON))
            && echelon != null && echelon.length == 1 && SymbolUtilities.getEchelonText(echelon).isNotEmpty()
        ) sidc.substring(0, 11) + echelon.uppercase() + sidc.substring(12) else sidc
    }

    @JvmStatic
    actual fun setSymbolModifier(
        sidc: String, hq: Boolean, taskForce: Boolean, feintDummy: Boolean, installation: Boolean, mobility: String?
    ): String {
        var result = sidc
        if (result.length >= 11 && !SymbolUtilities.isTacticalGraphic(result)) {
            // Check if mobility is applicable
            if (mobility != null && mobility.length == 2 && result.length >= 12
                && SymbolUtilities.canUnitHaveModifier(result, ModifiersUnits.R_MOBILITY_INDICATOR)
            ) {
                // Try applying mobility
                val sidcWithMobility = result.substring(0, 10) + mobility.uppercase() + result.substring(12)
                // Check if mobility is valid
                if (SymbolUtilities.isMobility(sidcWithMobility)) result = sidcWithMobility
            } else {
                // Check if HQ, TaskForce, Feint or Dummy symbol modifiers are applicable
                val modifier = if (hq && SymbolUtilities.canUnitHaveModifier(result, ModifiersUnits.S_HQ_STAFF_OR_OFFSET_INDICATOR)) "A" // Headquarters
                else if (installation && SymbolUtilities.canUnitHaveModifier(result, ModifiersUnits.AC_INSTALLATION)) "H" // Installation
                else if (taskForce && SymbolUtilities.canUnitHaveModifier(result, ModifiersUnits.D_TASK_FORCE_INDICATOR)
                    && feintDummy && SymbolUtilities.canUnitHaveModifier(result, ModifiersUnits.AB_FEINT_DUMMY_INDICATOR)) "D" // TaskForce/Feint/Dummy
                else if (taskForce && SymbolUtilities.canUnitHaveModifier(result, ModifiersUnits.D_TASK_FORCE_INDICATOR)) "B" // TaskForce
                else if (feintDummy && SymbolUtilities.canUnitHaveModifier(result, ModifiersUnits.AB_FEINT_DUMMY_INDICATOR)) "C" // Feint/Dummy
                else null
                // Apply symbol modifier
                if (modifier != null) result = result.substring(0, 10) + modifier + result.substring(11)
                // Fix Feint/Dummy modifier for installation
                if (feintDummy && result.length >= 12 && SymbolUtilities.hasInstallationModifier(result)
                    && SymbolUtilities.canUnitHaveModifier(result, ModifiersUnits.AB_FEINT_DUMMY_INDICATOR)
                ) result = result.substring(0, 11) + "B" + result.substring(12)
            }
        }
        return result
    }

    @JvmStatic
    actual fun setCountryCode(sidc: String, countryCode: String?) =
        if (sidc.length >= 14 && countryCode != null && countryCode.length == 2) {
            val result = sidc.substring(0, 12) + countryCode.uppercase() + sidc.substring(14)
            if (SymbolUtilities.hasValidCountryCode(result)) result else sidc
        } else sidc

    @JvmStatic
    actual fun getLineColor(sidc: String) = SymbolUtilities.getLineColorOfAffiliation(sidc)?.rgb
        ?: RendererSettings.getInstance().friendlyGraphicLineColor.rgb

    @JvmStatic
    actual fun getFillColor(sidc: String) = SymbolUtilities.getFillColorOfAffiliation(sidc)?.rgb
        ?: RendererSettings.getInstance().friendlyGraphicFillColor.rgb
}