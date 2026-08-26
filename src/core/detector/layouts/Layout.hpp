/**
 * @file Layout.hpp
 * @brief Interface every physical detector-board pixel layout must
 * implement (TILGAD, TREF, ...) - the mandatory contract that lets
 * GeometryManager and every consuming module (DUTAssociation,
 * AnalysisEfficiency, AnalysisTiming, PrealignmentCAEN, AlignmentCAEN)
 * work with any board shape without knowing what it looks like.
 */
#ifndef WARLOCK_LAYOUT_HPP
#define WARLOCK_LAYOUT_HPP

#include <string>
#include <vector>

namespace framework {

    struct DetectorGeo;
    class Histogram1D;

    /// Result of fitting one axis's flat-topped illuminated "plateau"
    /// (a board's own footprint in a track-intercept distribution) to a
    /// center position - see Layout::fitPlateauCenter().
    struct PlateauFit {
        double center;
        double derived_width;
        /// Non-empty whenever there's something worth logging at WARNING
        /// - either a severe "a channel looks missing" correction, or
        /// (independent of that) derived_width simply disagreeing with
        /// the configured width by more than the tolerance in either
        /// direction.
        std::string warning;
    };

    /**
     * @brief A detector board's physical pixel-grid layout: how column/row
     * position maps to physical spacing, whether/how rows split into
     * independent dies, and how to locate the board from track-intercept
     * data. Implementations are stateless - a detector's own pitch/
     * row-count fields live on DetectorGeo itself, passed in by reference,
     * so a Layout carries no per-detector data of its own.
     *
     * GeometryManager owns the one place that maps a `.geo` file's
     * `layout` string to a concrete Layout instance; every other module
     * only ever calls through DetectorGeo::dieOf()/hasDieSplit()/
     * effectivePitchX()/effectivePitchY()/getLayout()->prealign(...),
     * never touching a concrete Layout type directly unless it needs
     * something beyond this common interface, via DetectorGeo::getLayout()
     * + dynamic_cast.
     */
    class Layout {
    public:
        virtual ~Layout() = default;

        /// Column-to-column spacing actually used by DetectorGeo::
        /// getColumn()/getLocalX()/sizeX() for `det`.
        virtual double effectivePitchX(const DetectorGeo& det) const = 0;

        /// Row-to-row spacing actually used by DetectorGeo::getRow()/
        /// getLocalY()/sizeY() for `det`.
        virtual double effectivePitchY(const DetectorGeo& det) const = 0;

        /// "TOP"/"BOTTOM" (or whatever labels a future layout uses) for a
        /// (col, row) pixel/position on a split board, "" if this layout
        /// doesn't split `det` there. Takes both axes so a layout that
        /// splits by column (or by both) doesn't need the interface
        /// changed later. Every value this can return for a given `det`
        /// appears in dieLabels(det) - the two must never drift apart,
        /// which is exactly why hasDieSplit() below is derived from
        /// dieLabels() rather than declared separately.
        ///
        /// Not pure - defaults to "never splits" (returns ""), since die
        /// membership is normally config-driven (DetectorGeo::custom_dies,
        /// from a `.geo` file's `dies_file`) rather than computed from a
        /// hardcoded per-layout formula. A layout only needs to override
        /// this if it wants its own formula-based default independent of
        /// any dies_file - none of the current ones do.
        virtual std::string dieOf(const DetectorGeo& det, double col, double row) const;

        /// Every label dieOf() can return for `det` (e.g. {"TOP",
        /// "BOTTOM"}), in a stable, meaningful order (used to register/
        /// iterate per-die output) - empty if `det` doesn't split at all.
        /// Not pure - defaults to "never splits" (returns {}), same
        /// reasoning as dieOf() above. A future layout with more than two
        /// dies, or a column-split/irregular scheme, just returns however
        /// many labels it actually has - every consumer (DUTAssociation,
        /// AnalysisEfficiency, AnalysisTiming) iterates this list rather
        /// than assuming exactly two.
        virtual std::vector<std::string> dieLabels(const DetectorGeo& det) const;

        /// True iff dieOf() can ever return non-"" for `det`. Not
        /// virtual - always exactly !dieLabels(det).empty(), so it can
        /// never disagree with dieLabels() the way two independently
        /// hand-written implementations could.
        bool hasDieSplit(const DetectorGeo& det) const { return !dieLabels(det).empty(); }

        /// Locates `det` from its own track-intercept histograms (X/Y
        /// local-frame intercepts at det's Z, filled by PrealignmentCAEN)
        /// and updates det.position accordingly - each layout decides how
        /// to interpret its own physical footprint shape (e.g. "one
        /// combined rectangular plateau" vs. something else for a future
        /// layout).
        virtual void prealign(DetectorGeo& det, const Histogram1D& hist_x, const Histogram1D& hist_y) const = 0;

        /// Finds the center of a flat-topped "plateau" (a board's
        /// illuminated footprint in a track-intercept distribution) along
        /// one axis, by locating where the profile crosses half-max on
        /// either side of the peak rather than trusting a single maximum
        /// bin - the default centering strategy every current layout's
        /// prealign() uses. `width`/`pitch` are the detector's own
        /// geometry-derived expected footprint width and pixel pitch
        /// along this axis (DetectorGeo::sizeX()/sizeY(),
        /// pitch_x/pitch_y) - used both to size the search window and to
        /// sanity-check the result. `width_tolerance` is the relative
        /// disagreement (|derived_width - width| / width) above which
        /// PlateauFit::warning gets populated for the general mismatch
        /// warning, independent of the separate, more severe
        /// missing-channel check. Virtual (not pure) so a future layout
        /// whose footprint isn't a single flat plateau - several
        /// disconnected pads, an irregular shape - can override just this
        /// centering step without reimplementing prealign()'s own
        /// histogram/logging plumbing.
        virtual PlateauFit fitPlateauCenter(const Histogram1D& hist, double width, double pitch,
                                             double width_tolerance = 0.15) const;
    };
}
#endif
