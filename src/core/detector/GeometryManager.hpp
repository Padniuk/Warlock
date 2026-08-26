/**
 * @file GeometryManager.hpp
 * @brief Telescope/detector geometry: per-detector position, orientation,
 * pixel layout, masking, and the global registry that owns them all.
 */

#ifndef WARLOCK_GEOMETRYMANAGER_HPP
#define WARLOCK_GEOMETRYMANAGER_HPP

#include <string>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <toml++/toml.hpp>
#include "core/detector/layouts/Layout.hpp"
#include "core/detector/layouts/GenericLayout.hpp"

namespace framework {

    /**
     * @brief One detector's full geometric and physical description, plus
     * the local-frame <-> pixel-coordinate helpers every module needs.
     */
    struct DetectorGeo {
        std::string name, type, role;
        Eigen::Vector3d position{0,0,0}, orientation{0,0,0};
        double pitch_x{0}, pitch_y{0};
        /// Physical gap between column-split dies, mirroring #pitch_y_gap
        /// for the X axis. From the `.geo` file's `pitch_x_gap` key; 0
        /// (default) for every currently-known layout, since none splits
        /// by column - kept symmetric with #pitch_y_gap so a future
        /// column-split layout doesn't need a new field added to
        /// DetectorGeo itself, only its own Layout::effectivePitchX().
        double pitch_x_gap{0};
        /// Physical gap between the two dies of a multi-row layout (e.g.
        /// TILGAD), on top of pitch_y*n_pixels_y - not itself a pitch. The
        /// declared pitch_y is each die's own single-pad pitch; the two
        /// dies sit on the same carrier with an extra physical gap between
        /// them, so the board's real illuminated Y-extent is
        /// pitch_y*n_pixels_y + pitch_y_gap. From the `.geo` file's
        /// `pitch_y_gap` key; 0 (default) for a single-row board or any
        /// layout without a row gap. Only PrealignmentCAEN reads this
        /// directly (to derive the expected pre-split whole-board
        /// footprint); other per-die consumers use dieOf() instead.
        double pitch_y_gap{0};
        int n_pixels_x{0}, n_pixels_y{0};
        double spatial_res_x{0}, spatial_res_y{0};
        bool is_reference{false};

        /// If set, this detector's raw pixel hits are pulled from another
        /// detector's raw identity (`sensor_type + "_" + plane_id`, e.g.
        /// "CAEN_UZH_3") rather than its own name, filtered to raw hardware
        /// channels in #source_channels and remapped to row 0 in this
        /// detector's own local frame. This lets one physical multi-die raw
        /// sensor (a coarse DUT board with two independently-mounted pads
        /// sharing one DAQ channel group) be split into multiple Warlock
        /// detectors, each with its own accurate size and independently-fit
        /// position. Empty (default) means "not a split detector" -
        /// EventLoaderHDF5 uses this detector's own name as the raw
        /// identity instead, with no channel filter.
        std::string source_sensor;
        /// Which raw hardware channels (see minionRawConverter's own
        /// "/hits/channel" column) of #source_sensor's raw data belong to
        /// this detector; meaningless if #source_sensor is empty. Every
        /// listed channel is assumed to share the same raw row.
        std::vector<int> source_channels;

        /// Physical pad layout of this board, from the `.geo` file's
        /// `layout` key - lets a single generic CAEN pipeline
        /// (EventLoaderHDF5, WaveformProcessingCAEN, WaveformSelector,
        /// DUTAssociation, AnalysisEfficiency, AnalysisTiming) serve
        /// multiple board shapes without any of those modules knowing what
        /// a given layout looks like; they only ever consult #dieOf().
        /// Empty (default) means "not a recognized multi-die layout" -
        /// dieOf() returns "" for every row, same as a single-pad board.
        /// Known values: "OrthogonalLayout" (two independently-mounted
        /// pads stacked vertically, split by row - covers both TILGAD,
        /// which has a real gap between dies, and TREF, which doesn't).
        /// See src/core/detector/layouts/ for its actual implementation
        /// - GeometryManager::loadGeometry() is the only place that maps
        /// this string to #layout_impl.
        std::string layout;
        /// The concrete Layout answering #effectivePitchY()/#dieOf()/
        /// #hasDieSplit() below, set alongside #layout at parse time.
        /// Defaults to the no-split GenericLayout so a DetectorGeo is
        /// always safe to query even before parsing runs.
        const Layout* layout_impl = &GenericLayout::instance();

        Eigen::Matrix3d R{Eigen::Matrix3d::Identity()};      ///< Local-to-global rotation, from #orientation.
        Eigen::Matrix3d R_inv{Eigen::Matrix3d::Identity()};  ///< Global-to-local rotation (R's inverse).
        Eigen::Matrix2d V_inv{Eigen::Matrix2d::Identity()};  ///< Inverse spatial-resolution covariance, for chi2 weighting.

        /// Recomputes #R, #R_inv, and #V_inv from #orientation and
        /// #spatial_res_x/#spatial_res_y - call after changing either.
        void updateMatrices();

        /// The originally-parsed TOML table for this detector's own
        /// section, kept so saveGeometry() can round-trip every field -
        /// including ones DetectorGeo doesn't itself model (e.g.
        /// `coordinates`/`material_budget`, kept only for corry-format
        /// compatibility) - without bespoke bookkeeping. saveGeometry()
        /// only ever overwrites #position/#orientation in a copy of this
        /// before writing back out, since those are the only fields any
        /// alignment module actually mutates.
        toml::table raw_table;
        /// #raw_table's own keys, in original file order - toml::table
        /// doesn't preserve per-node source-position info through a copy
        /// (verified directly: a copied node's source position always
        /// comes back zeroed), so GeometryManager::loadGeometry() records
        /// this separately, from the freshly-parsed table, before it's
        /// ever copied into #raw_table. saveGeometry() writes keys in
        /// this order rather than #raw_table's own (unordered) iteration.
        std::vector<std::string> raw_key_order;

        /// @return Whether this detector has a pixel grid (as opposed to e.g. a TLU).
        bool isPixelDetector() const;

        /// Masked (col, row) pixels, loaded once from the `.geo`/`.toml`
        /// file's `mask_file` key at geometry-load time - shared by every
        /// module that needs them instead of each parsing its own copy.
        std::set<std::pair<int, int>> masked_pixels;
        /// Region-of-interest polygon vertices (col, row); empty means no restriction.
        std::vector<std::pair<int, int>> roi;

        /// One named group of (col, row) pixels - see #custom_dies.
        struct CustomDie {
            std::string name;
            std::set<std::pair<int, int>> pixels;
        };
        /// Explicit per-die pixel membership, from
        /// GeometryManager::loadDies() (a separate file, referenced by
        /// this geometry file's own `[General].dies_file` key) - if
        /// non-empty, authoritative and used INSTEAD of #layout_impl's
        /// own dieOf()/dieLabels() (which, for every current layout, just
        /// means "no split"). Empty (default, and every detector before/
        /// without a dies_file) means "use the layout's own answer."
        std::vector<CustomDie> custom_dies;

        /// Column-to-column spacing actually used by #getColumn()/
        /// #getLocalX() - see #layout_impl's own effectivePitchX() for the
        /// formula.
        double effectivePitchX() const { return layout_impl->effectivePitchX(*this); }
        /// Row-to-row spacing actually used by #getRow()/#getLocalY() -
        /// see #layout_impl's own effectivePitchY() for the formula.
        double effectivePitchY() const { return layout_impl->effectivePitchY(*this); }

        // (column, row) <-> local-plane-position (mm) conversions, matching
        // corry's PixelDetector::getColumn/getRow/getLocalPosition (col/row
        // centered so that pixel 0 sits at -((n-1)/2)*pitch). Both axes use
        // #effectivePitchX()/#effectivePitchY() so a multi-die board's real
        // spacing is reflected on whichever axis its layout splits.
        double getColumn(double local_x) const { return local_x / effectivePitchX() + (n_pixels_x - 1) / 2.0; }
        double getRow(double local_y) const { return local_y / effectivePitchY() + (n_pixels_y - 1) / 2.0; }
        double getLocalX(double column) const { return (column - (n_pixels_x - 1) / 2.0) * effectivePitchX(); }
        double getLocalY(double row) const { return (row - (n_pixels_y - 1) / 2.0) * effectivePitchY(); }

        /// Real total physical extent of the pixel grid (mm), using
        /// #effectivePitchX()/#effectivePitchY() so a multi-die board's
        /// real footprint (pitch plus any gap between dies) is reflected,
        /// not just n_pixels_*pitch_*. Use these instead of
        /// n_pixels_*pitch_* directly anywhere the result means "how big
        /// is this sensor really" (plot ranges, bounds sanity checks, ...).
        double sizeX() const { return (n_pixels_x - 1) * effectivePitchX() + pitch_x; }
        double sizeY() const { return (n_pixels_y - 1) * effectivePitchY() + pitch_y; }

        /// @return Whether (col, row) is in #masked_pixels.
        bool isMasked(int col, int row) const { return masked_pixels.count({col, row}) > 0; }
        /// @return Whether a local-plane position lies on the physical
        /// sensor, with an edge margin of `pixel_tolerance` pixels (corry's
        /// `hasIntercept()`).
        bool hasIntercept(double local_x, double local_y, double pixel_tolerance) const;
        /// @return Whether a local-plane position is on/near (within
        /// `tolerance` pixels, Chebyshev distance) a masked pixel (corry's
        /// `hitMasked()`).
        bool hitMasked(double local_x, double local_y, int tolerance) const;
        /// @return Whether a local-plane position falls inside #roi
        /// (always true if #roi is empty) - corry's `isWithinROI()`.
        bool isWithinROI(double local_x, double local_y) const;

        /// @return "" if this (col, row) position doesn't belong to a
        /// split die, else a die label ("TOP"/"BOTTOM" for a typical
        /// two-die board, matching the convention established by
        /// AnalysisEfficiency's `z_rotation_row` - positive = top =
        /// highest row index, negative = bottom = row 0). Checks
        /// #custom_dies first (rounding col/row to the nearest pixel,
        /// same convention as #hitMasked()) - if that's non-empty, it's
        /// authoritative; otherwise falls back to #layout_impl's own
        /// dieOf() (always "" for every current layout, see
        /// Layout::dieOf()'s own docs).
        std::string dieOf(double col, double row) const {
            if (!custom_dies.empty()) {
                std::pair<int, int> px{static_cast<int>(std::floor(col + 0.5)), static_cast<int>(std::floor(row + 0.5))};
                for (auto const& die : custom_dies) {
                    if (die.pixels.count(px)) return die.name;
                }
                return "";
            }
            return layout_impl->dieOf(*this, col, row);
        }

        /// Every label #dieOf() can return for this detector (e.g.
        /// {"TOP", "BOTTOM"}), empty if it doesn't split - the list a
        /// caller should register/iterate per-die output over
        /// (DUTAssociation, AnalysisEfficiency, AnalysisTiming), rather
        /// than assuming any fixed number of dies. Mirrors #dieOf()'s own
        /// #custom_dies-first, #layout_impl-fallback order.
        std::vector<std::string> dieLabels() const {
            if (!custom_dies.empty()) {
                std::vector<std::string> names;
                for (auto const& die : custom_dies) names.push_back(die.name);
                return names;
            }
            return layout_impl->dieLabels(*this);
        }

        /// True iff #dieOf() can ever return non-"" for this detector.
        /// Equivalent to !#dieLabels().empty() (so it reflects a
        /// #custom_dies override too, unlike delegating straight to
        /// #layout_impl->hasDieSplit()) - callers that need to know
        /// "should I register/fill per-die plot variants for this
        /// detector at all" should check this instead of deriving the
        /// condition from #n_pixels_y alone, since that doesn't account
        /// for #layout: a plain 2D-grid sensor (e.g. MIMOSA26, RD53B) can
        /// satisfy n_pixels_y >= 2 despite being a single die.
        bool hasDieSplit() const { return !dieLabels().empty(); }

        /// True only for a genuinely unrecognized, non-empty #layout
        /// string - lets a module warn about a likely typo/unsupported
        /// value while still falling back to non-split behavior via
        /// #dieOf()/#hasDieSplit(). A single shared check so
        /// DUTAssociation/AnalysisEfficiency/AnalysisTiming/
        /// PrealignmentCAEN don't each re-derive it (they still word their
        /// own WARNING).
        bool hasUnrecognizedLayout() const { return !layout.empty() && layout_impl == &GenericLayout::instance(); }

        /// Escape hatch for a module that needs something beyond the
        /// common Layout interface: dynamic_cast<const OrthogonalLayout*>
        /// (det.getLayout()) to reach a layout-specific extra.
        const Layout* getLayout() const { return layout_impl; }
    };

    /**
     * @brief Singleton registry owning every DetectorGeo for the current run.
     */
    class GeometryManager {
    public:
        static GeometryManager& getInstance();
        /// Parses a geometry TOML file into the detector registry. If the
        /// file has a `[General]` table with a `dies_file` key, also
        /// calls loadDies() on it automatically - a geometry file
        /// self-describes where its own dies live, so every pipeline
        /// stage that loads it gets dies for free with no per-config
        /// opt-in.
        void loadGeometry(const std::string& path);
        /// Parses `path` (a small standalone TOML file, keyed by detector
        /// name, each holding a `[[Detector.dies]]` array-of-tables) and
        /// populates DetectorGeo::custom_dies for whichever detectors
        /// appear in it and already exist in this registry. Normally
        /// called automatically by loadGeometry(); exposed separately in
        /// case a caller wants to load/replace dies independently.
        void loadDies(const std::string& path);
        /// Writes the current detector set back out as TOML, including a
        /// `[General]` block if a dies_file was loaded.
        void saveGeometry(const std::string& path);
        /// @return Whether a detector named `name` exists.
        bool hasDetector(const std::string& name) const;
        /// @return The named detector. Undefined if it doesn't exist - check hasDetector() first.
        DetectorGeo& getDetector(const std::string& name);
        /// @return Every detector name, in file order.
        std::vector<std::string> getDetectorNames() const;
        /// @return The name of the detector marked as reference (`is_reference`).
        std::string getReferenceName() const;

    private:
        GeometryManager() = default;
        std::map<std::string, DetectorGeo> detectors_;
        std::vector<std::string> detector_names_;
        std::string reference_name_;
        /// The `[General].dies_file` path this registry's geometry was
        /// loaded with, if any - not per-detector (unlike
        /// DetectorGeo::raw_table), so saveGeometry() re-emits it from
        /// here rather than from any one detector's own table.
        std::string dies_file_;
    };
}
#endif
