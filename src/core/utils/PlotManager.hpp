/**
 * @file PlotManager.hpp
 * @brief Thread-safe histogram/graph storage, and the singleton registry
 * every module fills into.
 */

#ifndef WARLOCK_PLOTMANAGER_HPP
#define WARLOCK_PLOTMANAGER_HPP

#include <H5Cpp.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <map>
#include <mutex>
#include <memory>

namespace framework {

    /// @brief A fixed-range, fixed-bin-count 1D histogram, safe to fill from multiple threads.
    class Histogram1D {
    public:
        Histogram1D(std::string name, size_t bins, double min, double max);
        /// Fills one entry.
        void fill(double value, double weight = 1.0);
        /// Fills many entries at once, taking the internal lock only once.
        void fillBatch(const std::vector<double>& values);
        /// Replaces this histogram's bin contents outright (e.g. when reloading saved data).
        void setData(const std::vector<double>& data);
        /// Optional descriptive title, distinct from #name_ (the
        /// registration key/HDF5 dataset name) - e.g. a per-channel plot
        /// whose real-world identity (which raw hardware channel it is)
        /// is only known once the first waveform routed to it has been
        /// seen, not at registration time. Empty (default) means
        /// H5ToRootConverter falls back to using #name_ as the ROOT
        /// object's title too, same as if this were never called. Safe to
        /// call repeatedly with the same value from multiple threads
        /// (e.g. once per waveform until #hasTitle() starts returning
        /// true) - redundant identical sets are harmless, just wasted
        /// work, so callers doing that don't need their own synchronization.
        void setTitle(const std::string& title);
        bool hasTitle() const { return !title_.empty(); }
        /// Writes this histogram to an HDF5 group as a 1D dataset (see PlotManager::saveToFile).
        void write(H5::Group& group);

        size_t getEntries() const { return entries_; }
        /// Overrides the entries counter directly - for callers that
        /// rescale bin content via setData() after real fill() calls
        /// already happened (setData()'s own recompute-from-sum is only
        /// correct for a pure derived/ratio plot with no prior fills).
        void setEntries(size_t entries) { std::lock_guard<std::mutex> lock(mtx_); entries_ = entries; }
        size_t getBins() const { return data_.size(); }
        double getMaxBinCenter() const;
        double getPrecisePeak() const;
        double getMean() const;
        double getMeanInRange(double lo, double hi) const;
        double getStdDev() const;
        double getBinContent(size_t i) const;
        double getBinCenter(size_t i) const;
        double getBinWidth() const { return bin_width_; }

    private:
        std::string name_;
        std::string title_; ///< See #setTitle() - empty means "no override, use name_".
        double min_, max_, bin_width_;
        size_t entries_{0};
        std::vector<double> data_;
        /// True (unbinned) first/second moments of every filled value,
        /// tracked alongside the bin content in fill()/fillBatch(). Unlike
        /// deriving mean/stddev from bin centers alone, this preserves the
        /// exact fractional value of a fill (e.g. a charge-weighted cluster
        /// centroid), so getMean()/getStdDev() - and the sum_w/sum_wx/
        /// sum_wx2 attributes write() persists for H5ToRootConverter's
        /// TH1::PutStats() - match ROOT's native TH1::Fill() statistics
        /// exactly rather than only to within one bin width.
        double sum_w_{0.0}, sum_wx_{0.0}, sum_wx2_{0.0};
        mutable std::mutex mtx_;
    };

    /// @brief A fixed-range, fixed-bin-count 2D histogram, safe to fill from multiple threads.
    class Histogram2D {
    public:
        Histogram2D(std::string name, size_t bx, double minx, double maxx,
                                     size_t by, double miny, double maxy);
        /// Fills one entry.
        void fill(double x, double y, double weight = 1.0);
        /// Adds a whole matrix of bin counts at once, taking the internal
        /// lock only once. True moments aren't available from a pre-binned
        /// matrix, so this overload approximates them from bin centers -
        /// exact for a plot whose fill values are always bin-integer (e.g.
        /// per-pixel occupancy), but only approximate for fractional fill
        /// values (e.g. a charge-weighted cluster centroid); see the other
        /// overload for that case.
        void fillBatch(const Eigen::MatrixXd& block);
        /// Same, but for a caller that already tracked the block's true
        /// (unbinned) moments and cross-moment locally - stores those
        /// directly instead of approximating from bin centers, so the
        /// sum_w/sum_x/.../sum_xy attributes write() persists for
        /// H5ToRootConverter's TH2::PutStats() match ROOT's native
        /// TH2::Fill() exactly.
        void fillBatch(const Eigen::MatrixXd& block, double sum_w, double sum_x, double sum_y,
                       double sum_x2, double sum_y2, double sum_xy);
        /// Replaces this histogram's bin contents outright (e.g. when reloading saved data).
        void setData(const Eigen::MatrixXd& matrix);
        /// Writes this histogram to an HDF5 group as a 2D dataset (see PlotManager::saveToFile).
        void write(H5::Group& group);

        size_t getEntries() const { return entries_; }
        size_t getBinsX() const { return bx_; }
        size_t getBinsY() const { return by_; }

    private:
        std::string name_;
        size_t bx_, by_;
        double minx_, maxx_, miny_, maxy_, w_x_, w_y_;
        size_t entries_{0};
        Eigen::MatrixXd data_;
        /// True (unbinned) moments - see the two-overload fillBatch() docs
        /// above. Left at 0 (matching zero entries) unless the exact-moment
        /// fillBatch() overload is used.
        double sum_w_{0.0}, sum_x_{0.0}, sum_y_{0.0}, sum_x2_{0.0}, sum_y2_{0.0}, sum_xy_{0.0};
        mutable std::mutex mtx_;
    };

    /// @brief An (x, y) point series (e.g. an alignment correction vs.
    /// iteration), safe to append to from multiple threads.
    class Graph1D {
    public:
        Graph1D(std::string name) : name_(std::move(name)) {}
        /// Appends one (x, y) point.
        void addPoint(double x, double y) {
            std::lock_guard<std::mutex> lock(mtx_);
            x_.push_back(x);
            y_.push_back(y);
        }
        /// Writes this graph to an HDF5 group as a 2D (N x 2) dataset (see PlotManager::saveToFile).
        void write(H5::Group& group);
    private:
        std::string name_;
        std::vector<double> x_, y_;
        mutable std::mutex mtx_;
    };

    /**
     * @brief Singleton registry of every histogram/graph in the run,
     * keyed by `"<module>/<name>"`.
     *
     * Modules register their plots in initialize() and fill them from
     * run() (safe to do concurrently - see Histogram1D/Histogram2D/Graph1D);
     * FrameworkManager calls saveToFile() once, after every module's
     * finalize(), to write the whole registry out to a single HDF5 file.
     */
    class PlotManager {
    public:
        static PlotManager& getInstance();

        void registerPlot1D(const std::string& mod, const std::string& name, size_t bins, double min, double max);
        void registerPlot2D(const std::string& mod, const std::string& name, size_t bx, double minx, double maxx, size_t by, double miny, double maxy);
        void registerGraph(const std::string& mod, const std::string& name);

        /// Discards a previously registered 2D histogram/graph so
        /// saveToFile() never writes it - for a caller that has to
        /// pre-register a safe upper bound of slots before the real count
        /// is known (e.g. a bounded iteration loop that may stop early),
        /// so the unused tail doesn't show up as empty output. No-op if
        /// not registered.
        void unregisterPlot2D(const std::string& mod, const std::string& name);
        void unregisterGraph(const std::string& mod, const std::string& name);

        /// Fills `"<mod>/<name>"`'s 1D histogram. A no-op if it was never registered.
        void fill1D(const std::string& mod, const std::string& name, double value);
        /// Fills `"<mod>/<name>"`'s 2D histogram. A no-op if it was never registered.
        void fill2D(const std::string& mod, const std::string& name, double x, double y);
        /// Appends a point to `"<mod>/<name>"`'s graph. A no-op if it was never registered.
        void addPoint(const std::string& mod, const std::string& name, double x, double y);

        /// @return Whether `"<mod>/<name>"` is a registered 1D histogram.
        bool hasPlot1D(const std::string& mod, const std::string& name) const;
        /// @return The registered 1D histogram. Undefined if it doesn't exist - check hasPlot1D() first.
        Histogram1D& getPlot1D(const std::string& mod, const std::string& name);
        /// @return The registered 2D histogram. Undefined if it doesn't exist.
        Histogram2D& getPlot2D(const std::string& mod, const std::string& name);
        /// @return The registered graph. Undefined if it doesn't exist.
        Graph1D& getGraph(const std::string& mod, const std::string& name);

        /// Writes every registered histogram/graph to one HDF5 file,
        /// mirroring the `"<module>/<name>"` key structure as nested groups.
        void saveToFile(const std::string& filename);

    private:
        PlotManager() = default;
        std::map<std::string, std::unique_ptr<Histogram1D>> plots1d_;
        std::map<std::string, std::unique_ptr<Histogram2D>> plots2d_;
        std::map<std::string, std::unique_ptr<Graph1D>> graphs_;
        mutable std::mutex global_mtx_;
    };
}
#endif
