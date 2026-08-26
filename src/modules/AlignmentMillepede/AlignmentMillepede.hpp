/**
 * @file AlignmentMillepede.hpp
 * @brief Global track-based alignment using the Millepede-II algorithm.
 */

#ifndef WARLOCK_ALIGNMENTMILLEPEDE_HPP
#define WARLOCK_ALIGNMENTMILLEPEDE_HPP

#include "core/Module.hpp"
#include <mutex>
#include <vector>
#include <string>
#include <unordered_map>

namespace framework {
    struct DetectorGeo;

    /// One linear measurement equation (a single X or Y cluster residual)
    /// as consumed by Millepede's local/global least-squares fit: local
    /// (track) derivatives, global (alignment) derivatives, and the
    /// nonlinear-correction terms needed to re-linearize them each outer
    /// iteration.
    struct Equation {
        double rmeas;              ///< Measured residual (cluster global X or Y).
        double weight;             ///< 1/sigma^2 for this measurement.
        std::vector<int> indL;     ///< Indices into the local (per-track) parameter vector.
        std::vector<double> derL;  ///< Local-parameter derivatives, parallel to #indL.
        std::vector<int> indG;     ///< Indices into the global (alignment) parameter vector.
        std::vector<double> derG;  ///< Global-parameter derivatives, parallel to #indG.
        std::vector<int> indNL;    ///< Which local-parameter slot (0 or 2) each nonlinear term corrects, parallel to #indG.
        std::vector<double> derNL; ///< Nonlinear derivative terms, parallel to #indG.
        std::vector<double> slopes; ///< Track slope (tx or ty) used by the nonlinear correction, parallel to #indG.
    };

    /// One linear constraint row (e.g. "sum of all X shifts = 0") added to
    /// the global normal-equations matrix to fix the otherwise-degenerate
    /// overall alignment (translation/rotation/shear of the whole setup).
    struct Constraint {
        double rhs;                      ///< Right-hand side of the constraint equation.
        std::vector<double> coefficients; ///< One coefficient per global parameter.
    };

    /// Per-thread accumulators for the parallel track loops in finalize()
    /// and fitGlobal() - each worker thread gets its own ThreadWorkspace so
    /// tracks can be processed without a lock, merging into the shared
    /// m_bgvec/m_cgmat only once per chunk.
    struct ThreadWorkspace {
        std::vector<double> bgvec;                 ///< Global parameter vector contribution (size m_nagb).
        std::vector<std::vector<double>> cgmat;     ///< Global normal-equations matrix contribution (m_nagb x m_nagb).
        std::vector<double> corrv;                  ///< Scratch: local-to-global correction vector.
        std::vector<std::vector<double>> corrm;     ///< Scratch: local-to-global correction matrix.
        std::vector<std::vector<double>> clcmat;    ///< Scratch: local/global cross-correlation matrix, reused per track.

        std::vector<double> dergb;   ///< Scratch: per-global-parameter derivative accumulator, cleared by addEquation().
        std::vector<double> dernl;   ///< Scratch: nonlinear derivative term, parallel to #dergb.
        std::vector<int> dernli;     ///< Scratch: nonlinear correction target slot, parallel to #dergb.
        std::vector<double> dernls;  ///< Scratch: track slope for the nonlinear term, parallel to #dergb.

        ThreadWorkspace(size_t nagb, size_t nalc) {
            bgvec.assign(nagb, 0.);
            cgmat.assign(nagb, std::vector<double>(nagb, 0.));
            corrv.assign(nagb, 0.);
            corrm.assign(nagb, std::vector<double>(nagb, 0.));
            clcmat.assign(nagb, std::vector<double>(nalc, 0.));
            dergb.assign(nagb, 0.);
            dernl.assign(nagb, 0.);
            dernli.assign(nagb, 0);
            dernls.assign(nagb, 0.);
        }
    };

    /**
     * @brief Global track-based alignment via the classic Millepede-II
     * algorithm (as used by corryvreckan's own AlignmentMillepede): every
     * track contributes local (track-parameter) and global (alignment
     * -parameter) measurement equations, local parameters are eliminated
     * per track, and the resulting reduced global normal equations are
     * solved once per outer iteration to get a simultaneous correction for
     * every plane at once - unlike AlignmentTrackChi2, which optimizes one
     * detector at a time. Linear-algebra routines (invertMatrix,
     * multiplyAX, ...) are direct ports of Millepede-II's own Gauss-Jordan
     * elimination, kept structurally close to the original rather than
     * replaced with a library solver so behavior matches exactly.
     */
    class AlignmentMillepede : public Module {
    public:
        AlignmentMillepede(Configuration cfg, Configuration g_cfg, ThreadPool* pool);
        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

        double getProgress() const override {
            if (max_tracks_ > 0) return static_cast<double>(total_evaluated_tracks_) / max_tracks_;
            return -1.0;
        }

    private:
        /// Builds the constraint set fixing the overall (degenerate)
        /// alignment of the whole setup - see the note in the .cpp on why
        /// Z-translation/scale constraints are always skipped.
        void setConstraints(size_t nPlanes);
        void addConstraint(const std::vector<double>& dercs, double rhs);
        /// Resets per-iteration global-fit state (matrices, vectors,
        /// which DOFs are free) ahead of a new outer iteration.
        bool reset(size_t nPlanes, double startfact);
        /// Builds this track's local+global measurement equations (two per
        /// plane hit: X and Y) and performs an initial local fit to reject
        /// it as an outlier before it contributes to the global matrices.
        bool putTrack(Track* track, size_t nPlanes, std::vector<Equation>& equations, ThreadWorkspace& ws);
        /// Converts one row of local/global derivatives into an Equation,
        /// draining ws.dergb (see ThreadWorkspace::dergb) as a side effect.
        void addEquation(std::vector<Equation>& equations, const std::vector<double>& derlc, ThreadWorkspace& ws, double rmeas, double sigma);
        /// Fits one track's local parameters by eliminating them from its
        /// equations, applying the residual cut, and (unless singlefit)
        /// folding its contribution into ws's global matrices.
        bool fitTrackLocal(const std::vector<Equation>& equations, std::vector<double>& trackParams, bool singlefit, unsigned int iteration, ThreadWorkspace& ws);
        /// Solves the accumulated global normal equations for one outer
        /// iteration, re-fitting every track's local parameters against
        /// the updated global parameters up to nMaxIterations times.
        bool fitGlobal();
        /// Applies the just-solved global parameter deltas (m_dparm) to
        /// every aligned detector's position/orientation.
        void updateGeometry();
        /// Millepede-II's original symmetric-matrix Gauss-Jordan inversion
        /// (row/column scaling + pivoting), operating on the full
        /// global+constraints matrix.
        int invertMatrix(std::vector<std::vector<double>>& v, std::vector<double>& b, size_t n);
        /// Same algorithm as invertMatrix() but for the small per-track
        /// local-parameter system (no row/column scaling needed there).
        int invertMatrixLocal(std::vector<std::vector<double>>& v, std::vector<double>& b, size_t n);
        /// Chi2 cut for `number_of_stddev_`-sigma outlier rejection,
        /// looked up from a fixed table for ndf <= 30, asymptotic formula
        /// otherwise - matches Millepede-II's own table exactly.
        double chi2Limit(int n, int nd) const;
        bool multiplyAX(const std::vector<std::vector<double>>& a, const std::vector<double>& x, std::vector<double>& y, unsigned int n, unsigned int m);
        bool multiplyAVAt(const std::vector<std::vector<double>>& v, const std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& w, unsigned int n, unsigned int m);
        void printResults();

        /// True if `det` should be left out of alignment - matched against
        /// #exclude_detectors_ by exact detector name OR by detector type,
        /// mirroring Tracking4D::isExcluded() so a plane can be included
        /// or excluded from the fit purely via config, without touching
        /// its role in the geometry file.
        bool isExcluded(const DetectorGeo& det) const;

        int iterations_;                ///< Number of outer (global) Millepede iterations.
        std::vector<std::string> exclude_detectors_; ///< Detector names and/or types left out of alignment; see #isExcluded().
        std::vector<bool> dofs_;        ///< 6 flags (Tx,Ty,Tz,Rx,Ry,Rz) selecting which DOFs are free.
        std::vector<double> sigmas_;    ///< 6 prior sigmas, parallel to #dofs_, regularizing each free DOF.
        size_t max_tracks_;             ///< Framework-wide track target (corry's `number_of_tracks`), used only for #getProgress().
        double residual_cut_;           ///< Per-equation residual cut used from the 2nd inner iteration onward.
        double residual_cut_init_;      ///< Looser residual cut used for the very first inner iteration.
        int number_of_stddev_;          ///< If nonzero, tracks whose local chi2 exceeds this many sigma are rejected.
        double convergence_cut_;        ///< Outer-iteration loop stops once mean |m_dparm| drops below this.

        std::vector<std::shared_ptr<Track>> global_tracks_; ///< Tracks accumulated across all batches, aligned against in finalize().
        size_t total_evaluated_tracks_{0}; ///< Count of every track seen in run().
        std::mutex track_mtx_; ///< Guards #global_tracks_ and #total_evaluated_tracks_ across concurrent run() calls.

        std::unordered_map<std::string, unsigned int> m_millePlanes; ///< Aligned detector name -> plane index (0..nPlanes-1).
        std::vector<Constraint> m_constraints; ///< Constraint rows built by setConstraints(), refreshed every outer iteration.
        std::vector<std::vector<Equation>> m_equations; ///< Per-track equation lists, indexed like #global_tracks_.

        unsigned int m_nagb{0};       ///< Number of global (alignment) parameters, 6 * nPlanes.
        const unsigned int m_nalc{4}; ///< Number of local (per-track) parameters: intercept + slope, for X and Y.
        std::vector<double> m_bgvec;  ///< Global parameter vector / right-hand side of the normal equations.
        std::vector<std::vector<double>> m_cgmat; ///< Global normal-equations matrix (accumulated, then inverted in place).
        std::vector<double> m_dparm;  ///< Current global parameter corrections (position/orientation deltas), one entry per #m_nagb.
        std::vector<bool> m_fixed;    ///< Per-global-parameter flag: true if held fixed (DOF disabled) this iteration.
        std::vector<double> m_psigm;  ///< Per-global-parameter prior sigma, parallel to #m_fixed.
        std::vector<double> m_diag;   ///< Diagonal of #m_cgmat captured before inversion, used for pull/correlation reporting.

        double m_cfactref{1.}; ///< Reference chi2 scale factor for outlier cuts.
        double m_cfactr{1.};   ///< Current (possibly annealed) chi2 scale factor for outlier cuts.
    };
}
#endif
