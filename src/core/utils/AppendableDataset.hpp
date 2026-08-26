/**
 * @file AppendableDataset.hpp
 * @brief A chunked HDF5 dataset that can be grown incrementally, batch by batch.
 */

#ifndef WARLOCK_CORE_APPENDABLE_DATASET_HPP
#define WARLOCK_CORE_APPENDABLE_DATASET_HPP

#include <H5Cpp.h>
#include <vector>
#include <string>

namespace framework {

    /**
     * @brief One chunked, growable HDF5 column (or fixed-width row) that
     * modules append to batch by batch, instead of buffering an entire
     * run's output in memory before a single write.
     *
     * Avoids two performance traps: an explicit `reserve()` call before
     * every append would force a full reallocation-and-copy of everything
     * appended so far (`push_back`'s own amortized growth doesn't have this
     * problem, so this class never calls `reserve()` itself), and an
     * undersized per-dataset chunk cache (HDF5's 1 MiB default) forces disk
     * I/O on every write to a chunk that doesn't fit in it - sized here via
     * `setChunkCache()` off the actual chunk footprint instead of a fixed
     * guess.
     */
    class AppendableDataset {
    public:
        /**
         * @param group Parent HDF5 group this dataset is created in.
         * @param name Dataset name within `group`.
         * @param type Element datatype.
         * @param rank 1 for a plain scalar column; 2 for a fixed-width row
         * (`width` columns per row).
         * @param width Row width, only meaningful when `rank == 2`.
         * @param chunk_rows HDF5 chunk size, in rows.
         */
        AppendableDataset(H5::Group& group, const std::string& name, H5::PredType type,
                           int rank = 1, hsize_t width = 1, hsize_t chunk_rows = 10000);

        /// Appends `data` as new rows of a rank-1 (scalar-column) dataset.
        template <typename T>
        void append(const std::vector<T>& data) {
            if (data.empty()) return;
            hsize_t new_dims[1] = { current_size_ + data.size() };
            dataset_.extend(new_dims);
            H5::DataSpace filespace = dataset_.getSpace();
            hsize_t offset[1] = { current_size_ };
            hsize_t dims[1] = { data.size() };
            filespace.selectHyperslab(H5S_SELECT_SET, dims, offset);
            H5::DataSpace memspace(1, dims);
            dataset_.write(data.data(), datatype_, memspace, filespace);
            current_size_ += data.size();
        }

        /// Appends `data` (flattened, row-major) as new rows of a rank-2
        /// (fixed-`width`-per-row) dataset.
        void append2D(const std::vector<double>& data, hsize_t width);

        /// @return Number of rows written so far.
        size_t size() const { return current_size_; }

    private:
        H5::DataSet dataset_;
        size_t current_size_;
        H5::PredType datatype_;
    };
}
#endif
