#include "core/utils/AppendableDataset.hpp"
#include <algorithm>

namespace framework {
    AppendableDataset::AppendableDataset(H5::Group& group, const std::string& name, H5::PredType type,
                                          int rank, hsize_t width, hsize_t chunk_rows)
        : current_size_(0), datatype_(type) {

        H5::DataSpace dataspace;
        H5::DSetCreatPropList prop;
        H5::DSetAccPropList dapl;

        if (rank == 1) {
            hsize_t dims[1] = {0};
            hsize_t max_dims[1] = {H5S_UNLIMITED};
            hsize_t chunk_dims[1] = {chunk_rows};
            dataspace = H5::DataSpace(1, dims, max_dims);
            prop.setChunk(1, chunk_dims);
            hsize_t chunk_bytes = chunk_rows * type.getSize();
            dapl.setChunkCache(10007, std::max<hsize_t>(chunk_bytes * 4, 1024 * 1024), 0.75);
        } else {
            hsize_t dims[2] = {0, width};
            hsize_t max_dims[2] = {H5S_UNLIMITED, width};
            hsize_t chunk_dims[2] = {chunk_rows, width};
            dataspace = H5::DataSpace(2, dims, max_dims);
            prop.setChunk(2, chunk_dims);
            hsize_t chunk_bytes = chunk_rows * width * type.getSize();
            dapl.setChunkCache(10007, std::max<hsize_t>(chunk_bytes * 4, 1024 * 1024), 0.75);
        }

        dataset_ = group.createDataSet(name, datatype_, dataspace, prop, dapl);
    }

    void AppendableDataset::append2D(const std::vector<double>& data, hsize_t width) {
        if (data.empty()) return;
        hsize_t num_rows = data.size() / width;
        hsize_t new_dims[2] = { current_size_ + num_rows, width };
        dataset_.extend(new_dims);
        H5::DataSpace filespace = dataset_.getSpace();
        hsize_t offset[2] = { current_size_, 0 };
        hsize_t dims[2] = { num_rows, width };
        filespace.selectHyperslab(H5S_SELECT_SET, dims, offset);
        H5::DataSpace memspace(2, dims);
        dataset_.write(data.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
        current_size_ += num_rows;
    }
}
