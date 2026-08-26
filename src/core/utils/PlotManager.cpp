#include "core/utils/PlotManager.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace framework {

    // --- Histogram1D ---
    Histogram1D::Histogram1D(std::string name, size_t bins, double min, double max)
        : name_(std::move(name)), min_(min), max_(max), data_(bins, 0.0) {
        bin_width_ = (max - min) / static_cast<double>(bins);
    }
    void Histogram1D::fill(double value, double weight) {
        std::lock_guard<std::mutex> lock(mtx_);
        entries_++;
        // Matches ROOT TH1::Fill()'s default (StatOverflows() off): entries_
        // counts every fill unconditionally, but the moment sums used for
        // GetMean()/GetStdDev() only include IN-RANGE values - an out-of-
        // range fill still lands in the (conceptual) under/overflow bin but
        // is excluded from the statistics, same as corry's native histograms.
        if (value < min_ || value >= max_) return;
        sum_w_ += weight; sum_wx_ += weight * value; sum_wx2_ += weight * value * value;
        size_t bin = static_cast<size_t>((value - min_) / bin_width_);
        data_[bin] += weight;
    }
    void Histogram1D::fillBatch(const std::vector<double>& values) {
        if (values.empty()) return;
        std::lock_guard<std::mutex> lock(mtx_);
        entries_ += values.size();
        for (double val : values) {
            // See fill()'s comment: only in-range values feed the moments.
            if (val < min_ || val >= max_) continue;
            sum_w_ += 1.0; sum_wx_ += val; sum_wx2_ += val * val;
            size_t bin = static_cast<size_t>((val - min_) / bin_width_);
            data_[bin] += 1.0;
        }
    }
    void Histogram1D::setTitle(const std::string& title) {
        std::lock_guard<std::mutex> lock(mtx_);
        title_ = title;
    }
    void Histogram1D::setData(const std::vector<double>& data) {
        std::lock_guard<std::mutex> lock(mtx_);
        data_ = data;
        double sum = 0.0;
        for (double v : data_) sum += v;
        entries_ = static_cast<size_t>(sum);
        // Not a raw-value fill (a derived/ratio histogram being written
        // outright) - true moments aren't meaningful here, only the bin
        // content is. Leave sum_w_/sum_wx_/sum_wx2_ at whatever they were
        // (0 for a freshly-constructed histogram); write() falls back to
        // marking them unavailable when sum_w_ doesn't match entries_.
    }
    double Histogram1D::getBinContent(size_t i) const {
        std::lock_guard<std::mutex> lock(mtx_);
        return (i < data_.size()) ? data_[i] : 0.0;
    }
    double Histogram1D::getBinCenter(size_t i) const {
        return min_ + (static_cast<double>(i) + 0.5) * bin_width_;
    }
    double Histogram1D::getMaxBinCenter() const {
        std::lock_guard<std::mutex> lock(mtx_);
        if (data_.empty()) return 0.0;
        auto it = std::max_element(data_.begin(), data_.end());
        size_t bin = std::distance(data_.begin(), it);
        return min_ + (static_cast<double>(bin) + 0.5) * bin_width_;
    }
    double Histogram1D::getMean() const {
        std::lock_guard<std::mutex> lock(mtx_);
        // True (exact-value) mean when this histogram was populated via
        // fill()/fillBatch() - matches corry's native ROOT TH1::GetMean()
        // instead of only approximating it to within one bin width via bin
        // centers. Falls back to the bin-center approximation for a
        // setData()-populated (derived/ratio) histogram, where sum_w_
        // never got updated and so can't be trusted.
        if (sum_w_ > 0) return sum_wx_ / sum_w_;
        double sum = 0, count = 0;
        for (size_t i = 0; i < data_.size(); ++i) {
            double center = min_ + (i + 0.5) * bin_width_;
            sum += data_[i] * center;
            count += data_[i];
        }
        return (count > 0) ? (sum / count) : 0.0;
    }
    double Histogram1D::getMeanInRange(double lo, double hi) const {
        std::lock_guard<std::mutex> lock(mtx_);
        double sum = 0, count = 0;
        for (size_t i = 0; i < data_.size(); ++i) {
            double center = min_ + (i + 0.5) * bin_width_;
            if (center < lo || center > hi) continue;
            sum += data_[i] * center;
            count += data_[i];
        }
        return (count > 0) ? (sum / count) : 0.0;
    }
    double Histogram1D::getStdDev() const {
        std::lock_guard<std::mutex> lock(mtx_);
        if (entries_ == 0) return 0.0;
        // Same true-vs-bin-approximated distinction as getMean() above.
        if (sum_w_ > 0) {
            double mean = sum_wx_ / sum_w_;
            double variance = (sum_wx2_ / sum_w_) - (mean * mean);
            return std::sqrt(std::max(0.0, variance));
        }
        double sum_w = 0, sum_wx = 0, sum_wx2 = 0;
        for (size_t i = 0; i < data_.size(); ++i) {
            double x = min_ + (i + 0.5) * bin_width_;
            sum_w   += data_[i];  sum_wx  += data_[i] * x;  sum_wx2 += data_[i] * x * x;
        }
        if (sum_w <= 0) return 0.0;
        double mean = sum_wx / sum_w;
        double variance = (sum_wx2 / sum_w) - (mean * mean);
        return std::sqrt(std::max(0.0, variance));
    }
    double Histogram1D::getPrecisePeak() const {
        std::lock_guard<std::mutex> lock(mtx_);
        if (data_.empty()) return 0.0;
        auto it = std::max_element(data_.begin(), data_.end());
        size_t bin = std::distance(data_.begin(), it);
        if (bin > 0 && bin < data_.size() - 1) {
            double y1 = *(it - 1), y2 = *it, y3 = *(it + 1);
            double denom = (y1 - 2.0 * y2 + y3);
            if (std::abs(denom) > 1e-9) {
                double offset = 0.5 * (y1 - y3) / denom;
                return min_ + (static_cast<double>(bin) + 0.5 + offset) * bin_width_;
            }
        }
        return min_ + (static_cast<double>(bin) + 0.5) * bin_width_;
    }
    void Histogram1D::write(H5::Group& group) {
        if (data_.empty()) return;
        hsize_t dims[1] = { data_.size() };
        H5::DataSpace dataspace(1, dims);
        H5::DataSet dataset = group.createDataSet(name_, H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(data_.data(), H5::PredType::NATIVE_DOUBLE);
        dataset.createAttribute("min", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &min_);
        dataset.createAttribute("max", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &max_);
        double ent = static_cast<double>(entries_);
        dataset.createAttribute("entries", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &ent);
        // Only meaningful (and only ever nonzero) for a fill()/fillBatch()-
        // populated histogram - see getMean()'s docs. H5ToRootConverter
        // only trusts these (calls TH1::PutStats()) when sum_w > 0.
        dataset.createAttribute("sum_w", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &sum_w_);
        dataset.createAttribute("sum_wx", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &sum_wx_);
        dataset.createAttribute("sum_wx2", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &sum_wx2_);
        if (!title_.empty()) {
            H5::StrType st(H5::PredType::C_S1, title_.size() + 1);
            dataset.createAttribute("title", st, H5S_SCALAR).write(st, title_);
        }
    }

    // --- Histogram2D ---
    Histogram2D::Histogram2D(std::string name, size_t bx, double minx, double maxx, size_t by, double miny, double maxy)
        : name_(std::move(name)), bx_(bx), by_(by), minx_(minx), maxx_(maxx), miny_(miny), maxy_(maxy) {
        data_ = Eigen::MatrixXd::Zero(bx, by);
        w_x_ = (maxx - minx) / static_cast<double>(bx);
        w_y_ = (maxy - miny) / static_cast<double>(by);
    }
    void Histogram2D::fill(double x, double y, double weight) {
        std::lock_guard<std::mutex> lock(mtx_);
        entries_++;
        // Same in-range-only moment convention as Histogram1D::fill() -
        // matches ROOT TH2::Fill()'s default (StatOverflows() off).
        if (x < minx_ || x >= maxx_ || y < miny_ || y >= maxy_) return;
        sum_w_ += weight; sum_x_ += weight * x; sum_y_ += weight * y;
        sum_x2_ += weight * x * x; sum_y2_ += weight * y * y; sum_xy_ += weight * x * y;
        size_t ix = static_cast<size_t>((x - minx_) / w_x_);
        size_t iy = static_cast<size_t>((y - miny_) / w_y_);
        data_(ix, iy) += weight;
    }
    void Histogram2D::fillBatch(const Eigen::MatrixXd& block) {
        std::lock_guard<std::mutex> lock(mtx_);
        data_ += block;
        entries_ += static_cast<size_t>(block.sum());
        // No true per-value moments available from a pre-binned block -
        // approximate them from bin centers (exact for an always-integer
        // fill value like MaskCreator's occupancy, since bin centers land
        // exactly on integers there; only an approximation otherwise - see
        // the other overload for that case).
        for (size_t ix = 0; ix < bx_; ++ix) {
            double x = minx_ + (ix + 0.5) * w_x_;
            for (size_t iy = 0; iy < by_; ++iy) {
                double w = block(ix, iy);
                if (w == 0.0) continue;
                double y = miny_ + (iy + 0.5) * w_y_;
                sum_w_ += w; sum_x_ += w * x; sum_y_ += w * y;
                sum_x2_ += w * x * x; sum_y2_ += w * y * y; sum_xy_ += w * x * y;
            }
        }
    }
    void Histogram2D::fillBatch(const Eigen::MatrixXd& block, double sum_w, double sum_x, double sum_y,
                                double sum_x2, double sum_y2, double sum_xy) {
        std::lock_guard<std::mutex> lock(mtx_);
        data_ += block;
        entries_ += static_cast<size_t>(block.sum());
        sum_w_ += sum_w; sum_x_ += sum_x; sum_y_ += sum_y;
        sum_x2_ += sum_x2; sum_y2_ += sum_y2; sum_xy_ += sum_xy;
    }
    void Histogram2D::setData(const Eigen::MatrixXd& matrix) {
        std::lock_guard<std::mutex> lock(mtx_);
        data_ = matrix;
        entries_ = static_cast<size_t>(matrix.sum());
    }
    void Histogram2D::write(H5::Group& group) {
        if (data_.size() == 0 || bx_ == 0 || by_ == 0) return;
        hsize_t dims[2] = { static_cast<hsize_t>(bx_), static_cast<hsize_t>(by_) };
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = group.createDataSet(name_, H5::PredType::NATIVE_DOUBLE, dataspace);
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> row_major = data_;
        dataset.write(row_major.data(), H5::PredType::NATIVE_DOUBLE);
        double ent = static_cast<double>(entries_);
        dataset.createAttribute("entries", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &ent);
        dataset.createAttribute("min_x", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &minx_);
        dataset.createAttribute("max_x", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &maxx_);
        dataset.createAttribute("min_y", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &miny_);
        dataset.createAttribute("max_y", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &maxy_);
        // See Histogram1D::write()'s equivalent attributes - H5ToRootConverter
        // only trusts these (calls TH2::PutStats()) when sum_w > 0.
        dataset.createAttribute("sum_w", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &sum_w_);
        dataset.createAttribute("sum_x", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &sum_x_);
        dataset.createAttribute("sum_y", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &sum_y_);
        dataset.createAttribute("sum_x2", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &sum_x2_);
        dataset.createAttribute("sum_y2", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &sum_y2_);
        dataset.createAttribute("sum_xy", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR).write(H5::PredType::NATIVE_DOUBLE, &sum_xy_);
    }

    // --- Graph1D ---
    void Graph1D::write(H5::Group& group) {
        if (x_.empty()) return;
        hsize_t dims[2] = { x_.size(), 2 };
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = group.createDataSet(name_, H5::PredType::NATIVE_DOUBLE, dataspace);
        std::vector<double> flat(x_.size() * 2);
        for(size_t i = 0; i < x_.size(); ++i) { flat[2*i] = x_[i]; flat[2*i+1] = y_[i]; }
        dataset.write(flat.data(), H5::PredType::NATIVE_DOUBLE);
        int is_graph = 1;
        dataset.createAttribute("is_graph", H5::PredType::NATIVE_INT, H5S_SCALAR).write(H5::PredType::NATIVE_INT, &is_graph);
    }

    // --- PlotManager Implementation ---
    PlotManager& PlotManager::getInstance() { static PlotManager instance; return instance; }

    bool PlotManager::hasPlot1D(const std::string& mod, const std::string& name) const {
        std::lock_guard<std::mutex> lock(global_mtx_);
        return plots1d_.find(mod + "/" + name) != plots1d_.end();
    }
    void PlotManager::registerPlot1D(const std::string& mod, const std::string& name, size_t bins, double min, double max) {
        std::lock_guard<std::mutex> lock(global_mtx_);
        plots1d_[mod + "/" + name] = std::make_unique<Histogram1D>(name, bins, min, max);
    }
    void PlotManager::registerPlot2D(const std::string& mod, const std::string& name, size_t bx, double minx, double maxx, size_t by, double miny, double maxy) {
        std::lock_guard<std::mutex> lock(global_mtx_);
        plots2d_[mod + "/" + name] = std::make_unique<Histogram2D>(name, bx, minx, maxx, by, miny, maxy);
    }
    void PlotManager::registerGraph(const std::string& mod, const std::string& name) {
        std::lock_guard<std::mutex> lock(global_mtx_);
        graphs_[mod + "/" + name] = std::make_unique<Graph1D>(name);
    }
    void PlotManager::fill1D(const std::string& mod, const std::string& name, double value) {
        auto it = plots1d_.find(mod + "/" + name);
        if (it != plots1d_.end()) it->second->fill(value);
    }
    void PlotManager::fill2D(const std::string& mod, const std::string& name, double x, double y) {
        auto it = plots2d_.find(mod + "/" + name);
        if (it != plots2d_.end()) it->second->fill(x, y);
    }
    void PlotManager::addPoint(const std::string& mod, const std::string& name, double x, double y) {
        auto it = graphs_.find(mod + "/" + name);
        if (it != graphs_.end()) it->second->addPoint(x, y);
    }
    void PlotManager::unregisterPlot2D(const std::string& mod, const std::string& name) {
        std::lock_guard<std::mutex> lock(global_mtx_);
        plots2d_.erase(mod + "/" + name);
    }
    void PlotManager::unregisterGraph(const std::string& mod, const std::string& name) {
        std::lock_guard<std::mutex> lock(global_mtx_);
        graphs_.erase(mod + "/" + name);
    }
    Histogram1D& PlotManager::getPlot1D(const std::string& mod, const std::string& name) { return *plots1d_.at(mod + "/" + name); }
    Histogram2D& PlotManager::getPlot2D(const std::string& mod, const std::string& name) { return *plots2d_.at(mod + "/" + name); }
    Graph1D& PlotManager::getGraph(const std::string& mod, const std::string& name) { return *graphs_.at(mod + "/" + name); }

    void PlotManager::saveToFile(const std::string& filename) {
        H5::Exception::dontPrint();

        std::unique_ptr<H5::H5File> file;
        try {
            file = std::make_unique<H5::H5File>(filename, H5F_ACC_TRUNC);
        } catch (const H5::Exception& e) {
            std::cerr << "[Warlock PlotManager] Could not create " << filename << ": " << e.getDetailMsg() << std::endl;
            return;
        }

        std::map<std::string, std::shared_ptr<H5::Group>> groups;

        // Each plot's write() is wrapped in its own try/catch (see
        // write_plot below) rather than one try/catch around the whole
        // loop, so a single plot's HDF5 failure is reported and skipped
        // without aborting every plot after it in map order.
        std::function<H5::Group*(const std::string&)> get_or_create = [&](const std::string& path) -> H5::Group* {
            if (path.empty()) return nullptr;
            auto it = groups.find(path);
            if (it != groups.end()) return it->second.get();

            size_t pos = path.find_last_of('/');
            if (pos != std::string::npos) {
                std::string parent_path = path.substr(0, pos);
                std::string child_name = path.substr(pos + 1);
                H5::Group* parent_grp = get_or_create(parent_path);
                if (!parent_grp) return nullptr;
                auto grp = std::make_shared<H5::Group>(parent_grp->createGroup(child_name));
                groups[path] = grp;
                return grp.get();
            }
            auto grp = std::make_shared<H5::Group>(file->createGroup(path));
            groups[path] = grp;
            return grp.get();
        };

        size_t n_failed = 0;
        auto write_plot = [&](const std::string& full_path, auto& plot) {
            size_t pos = full_path.find_last_of('/');
            if (pos == std::string::npos) return;
            try {
                std::string group_path = full_path.substr(0, pos);
                H5::Group* grp = get_or_create(group_path);
                if (grp) plot->write(*grp);
            } catch (const H5::Exception& e) {
                std::cerr << "[Warlock PlotManager] Failed to write '" << full_path << "': "
                           << e.getDetailMsg() << " - skipping this plot, continuing with the rest." << std::endl;
                ++n_failed;
            }
        };

        for (auto const& [path, plot] : plots1d_) write_plot(path, plot);
        for (auto const& [path, plot] : plots2d_) write_plot(path, plot);
        for (auto const& [path, plot] : graphs_) write_plot(path, plot);

        if (n_failed > 0) {
            std::cerr << "[Warlock PlotManager] " << n_failed << " plot(s) failed to write - see above." << std::endl;
        }
    }
}
