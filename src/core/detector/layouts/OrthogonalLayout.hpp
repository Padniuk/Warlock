/**
 * @file OrthogonalLayout.hpp
 * @brief Two independently-mounted dies stacked vertically along an axis
 * (row-split), possibly with a real physical gap between them
 * (`pitch_y_gap`) - covers both TILGAD (real gap) and TREF (no gap)
 * since the two turned out to need identical pitch/alignment handling in
 * practice. Die membership itself is never hardcoded here - see
 * DetectorGeo::custom_dies/`dies_file` (Layout::dieOf()/dieLabels()'s own
 * default, inherited unmodified, is "never splits" until config says
 * otherwise). A future layout with genuinely different pitch/alignment
 * math (not just a different die split) gets its own Layout subclass.
 */
#ifndef WARLOCK_ORTHOGONALLAYOUT_HPP
#define WARLOCK_ORTHOGONALLAYOUT_HPP

#include "core/detector/layouts/Layout.hpp"

namespace framework {

    class OrthogonalLayout : public Layout {
    public:
        double effectivePitchX(const DetectorGeo& det) const override;
        double effectivePitchY(const DetectorGeo& det) const override;
        void prealign(DetectorGeo& det, const Histogram1D& hist_x, const Histogram1D& hist_y) const override;

        /// Process-lifetime singleton - OrthogonalLayout carries no
        /// per-detector state, so one instance serves every board using
        /// it.
        static const OrthogonalLayout& instance();
    };
}
#endif
