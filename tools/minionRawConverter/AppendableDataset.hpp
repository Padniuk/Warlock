/**
 * @file AppendableDataset.hpp
 * @brief A chunked, growable HDF5 dataset that minionRawConverter streams
 * decoded hits to batch by batch.
 */

#ifndef APPENDABLE_DATASET_HPP
#define APPENDABLE_DATASET_HPP

#include <H5Cpp.h>
#include <vector>
#include <string>

namespace converter {
    /**
     * @brief Standalone counterpart to framework::AppendableDataset (this
     * tool doesn't link against WarlockCore) - one extendable HDF5 column
     * appended to as the converter decodes each batch of raw events,
     * instead of buffering the whole run in memory before a single write.
     */
    class AppendableDataset {
    public:
        /**
         * @param group Parent HDF5 group this dataset is created in.
         * @param name Dataset name within `group`.
         * @param type Element datatype.
         *
         * `rank` is documented on the .cpp's constructor definition (1 for
         * hits metadata, 2 for raw waveforms; see also #append2D).
         */
        AppendableDataset(H5::Group& group, const std::string& name, H5::PredType type, int rank = 1);

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

        /// Appends `data` (flattened, row-major) as new rows of a rank-2 (fixed-`width`-per-row) dataset.
        void append2D(const std::vector<double>& data, hsize_t width);

    private:
        H5::DataSet dataset_;
        size_t current_size_;
        H5::PredType datatype_;
        int rank_;
    };
}
#endif