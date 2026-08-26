/**
 * @file Correlations.hpp
 * @brief Hit/cluster correlation diagnostics between every detector and the reference plane.
 */

#ifndef WARLOCK_CORRELATIONS_HPP
#define WARLOCK_CORRELATIONS_HPP

#include "core/Module.hpp"

namespace framework {
    /**
     * @brief Fills raw hitmaps and pixel/cluster correlation histograms
     * (position and time, local and global) between every non-reference
     * detector and the reference plane - the standard first diagnostic run
     * after loading data, used to check gross alignment/timing before any
     * actual alignment module runs. Mirrors corryvreckan's Correlations
     * module.
     */
    class Correlations : public Module {
    public:
        Correlations(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;

    private:
        double time_cut_abs_; ///< Cluster-pair time cut, in ns (only applied when #do_time_cut_ is set).
        bool do_time_cut_;    ///< Whether #time_cut_abs_ actually gates the position-correlation fills.
        double time_binning_; ///< Bin width for the correlationTime histogram, in ns.
    };
}
#endif
