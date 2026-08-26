/**
 * @file PrealignmentCAEN.hpp
 * @brief X/Y prealignment for CAEN digitizer DUTs from track-intercept edge profiles.
 */

#ifndef WARLOCK_PREALIGNMENTCAEN_HPP
#define WARLOCK_PREALIGNMENTCAEN_HPP

#include "core/Module.hpp"
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <map>

namespace framework {
    /**
     * @brief Aligns CAEN digitizer DUTs (coarse-pixel boards with no
     * corryvreckan Prealignment counterpart) by profiling track intercepts
     * at each DUT's Z position and locating the half-max edges of the
     * resulting plateau - the DUT's illuminated physical footprint - rather
     * than correlating against a reference plane's own hit distribution.
     * A CAEN DUT's footprint is a hard-edged rectangle set by its physical
     * extent (for the layouts that fit that model), not a track-correlated
     * peak, so standard Prealignment (which looks for a correlation peak)
     * has nothing to lock onto here. This module only builds the per-DUT
     * intercept histograms (run()) and hands them to each DUT's own
     * DetectorGeo::getLayout()->prealign() (finalize()) - the actual
     * locating strategy, and whether it's edge-detection at all, is a
     * property of that DUT's Layout (see src/core/detector/layouts/ -
     * OrthogonalLayout and GenericLayout both use their own independent
     * half-max plateau finder).
     *
     * Example configuration:
     * @code
     * [PrealignmentCAEN]
     * dut_names = ["CAEN_IJS_0", "CAEN_IJS_4"]
     * @endcode
     */
    class PrealignmentCAEN : public Module {
    public:
        PrealignmentCAEN(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

    private:
        std::vector<std::string> dut_names_; ///< CAEN DUTs to align, in config order.
    };
}
#endif
