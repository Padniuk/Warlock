/**
 * @file H5ToRootConverter.hpp
 * @brief Converts Warlock's own HDF5 plot output into ROOT TH1/TH2 histograms.
 */

#ifndef WARLOCK_H5TOROOTCONVERTER_HPP
#define WARLOCK_H5TOROOTCONVERTER_HPP

#include <string>
#include <H5Cpp.h>
#include <TFile.h>
#include <TDirectory.h>

namespace framework {
    /**
     * @brief Walks an HDF5 file written by PlotManager (Histogram1D::write()
     * / Histogram2D::write()) and rebuilds the same directory tree of ROOT
     * TH1D/TH2D histograms, so Warlock's output can be inspected with
     * standard ROOT tools without Warlock itself ever linking against ROOT.
     *
     * Each converted histogram's Mean/StdDev come from PutStats() using the
     * source's true (unbinned) moment attributes (sum_w/sum_wx/sum_wx2) when
     * available, not recomputed from the binned array - see PlotManager's
     * own docs on the true-vs-bin-center-approximated distinction. Bin
     * content also gets an optional "title" attribute (Histogram1D::
     * setTitle()) carried over via #readStringAttribute.
     */
    class H5ToRootConverter {
    public:
        H5ToRootConverter() = default;
        ~H5ToRootConverter() = default;

        /// Converts every dataset/group in `input_h5` into `output_root`, preserving the directory structure.
        bool convert(const std::string& input_h5, const std::string& output_root);

    private:
        void processGroup(H5::Group& h5_group, TDirectory* root_dir);
        void convertDataset(H5::Group& h5_group, const std::string& name, TDirectory* root_dir);
        
        double readAttribute(const H5::H5Object& obj, const std::string& name);
        /// @return The named string attribute's value, or "" if absent -
        /// used for Histogram1D::write()'s optional "title" attribute
        /// (see PlotManager.hpp's docs on Histogram1D::setTitle()).
        std::string readStringAttribute(const H5::H5Object& obj, const std::string& name);
    };
}

#endif