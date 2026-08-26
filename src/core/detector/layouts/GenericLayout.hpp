/**
 * @file GenericLayout.hpp
 * @brief Fallback Layout for an empty or unrecognized `.geo` `layout`
 * string - a genuine 2D pixel grid (MIMOSA26, RD53B, ...). Never splits
 * into dies on its own (inherits Layout::dieOf()/dieLabels()'s "never
 * splits" default unmodified) - only DetectorGeo::custom_dies (from a
 * `dies_file`) can split a detector, regardless of which Layout it uses.
 */
#ifndef WARLOCK_GENERICLAYOUT_HPP
#define WARLOCK_GENERICLAYOUT_HPP

#include "core/detector/layouts/Layout.hpp"

namespace framework {

    class GenericLayout : public Layout {
    public:
        double effectivePitchX(const DetectorGeo& det) const override;
        double effectivePitchY(const DetectorGeo& det) const override;
        void prealign(DetectorGeo& det, const Histogram1D& hist_x, const Histogram1D& hist_y) const override;

        /// Process-lifetime singleton - GenericLayout carries no
        /// per-detector state, so one instance serves every detector
        /// that falls back to it.
        static const GenericLayout& instance();
    };
}
#endif
