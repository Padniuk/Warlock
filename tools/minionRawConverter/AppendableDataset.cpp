
/**
 * @file AppendableDataset.cpp
 * @brief Implementation of the HDF5 dataset management for Warlock data.
 */

#include "AppendableDataset.hpp"

namespace converter {
    /**
     * @brief Constructs a new dataset with chunking and compression.
     * @param rank Dimensionality (1 for hits metadata, 2 for raw waveforms).
     */
    AppendableDataset::AppendableDataset(H5::Group& group, const std::string& name, H5::PredType type, int rank)
        : current_size_(0), datatype_(type), rank_(rank) {
        
        H5::DataSpace dataspace;
        H5::DSetCreatPropList prop;

        // Per-dataset chunk cache, sized to comfortably hold several whole
        // chunks. HDF5's default is only 1 MiB, too small to hold even a
        // single chunk of the 2D samples dataset below (1000 rows x 1024
        // doubles = ~8 MiB/chunk); an undersized cache forces every write
        // to fetch/flush its chunk from/to disk regardless of file size.
        H5::DSetAccPropList dapl;

        if (rank == 1) {
            hsize_t dims[1] = {0};
            hsize_t max_dims[1] = {H5S_UNLIMITED};
            hsize_t chunk_dims[1] = {10000};
            dataspace = H5::DataSpace(1, dims, max_dims);
            prop.setChunk(1, chunk_dims);
            // 8 bytes/element x 10000 = ~80 KiB/chunk; a few MiB comfortably
            // holds many chunks resident at once.
            dapl.setChunkCache(10007, 8 * 1024 * 1024, 0.75);
        } else {
            // Rank 2 optimized for waveform data: [Hits] x [1024 Samples]
            hsize_t dims[2] = {0, 1024};
            hsize_t max_dims[2] = {H5S_UNLIMITED, 1024};
            hsize_t chunk_dims[2] = {1000, 1024};
            dataspace = H5::DataSpace(2, dims, max_dims);
            prop.setChunk(2, chunk_dims);
            // 1000 rows x 1024 doubles = ~8 MiB/chunk - cache several of
            // them so an in-progress (not-yet-full) chunk stays resident
            // across flushes instead of round-tripping to disk every time.
            dapl.setChunkCache(10007, 64 * 1024 * 1024, 0.75);
        }

        // No compression: for a chunked dataset appended to incrementally, a
        // write that doesn't exactly fill a chunk forces HDF5 to
        // decompress the existing chunk, merge in the new rows, and
        // recompress the whole chunk on every flush until it fills up.
        // That cost scales with flush count, not data volume, and outweighs
        // compression's benefit for a write-heavy, incrementally-appended
        // file like this.

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