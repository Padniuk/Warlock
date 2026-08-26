/**
 * @file H5ToRootConverter.cpp
 * @brief H5ToRootConverter implementation - see H5ToRootConverter.hpp for the class docs.
 */

#include "H5ToRootConverter.hpp"
#include <algorithm>
#include <iostream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TAxis.h>

namespace framework {

    bool H5ToRootConverter::convert(const std::string& input_h5, const std::string& output_root) {
        try {
            H5::Exception::dontPrint();
            H5::H5File file(input_h5, H5F_ACC_RDONLY);
            TFile* rFile = new TFile(output_root.c_str(), "RECREATE");

            if (!rFile || rFile->IsZombie()) {
                std::cerr << "[minion] Error: Could not create ROOT file " << output_root << std::endl;
                return false;
            }

            std::cout << "[minion] Starting conversion: " << input_h5 << " -> " << output_root << std::endl;

            H5::Group root_group = file.openGroup("/");
            processGroup(root_group, rFile);

            // No top-level rFile->Write() here - every histogram/graph
            // already calls its own .Write() inside convertDataset(). ROOT
            // histograms auto-register into the current TDirectory at
            // construction time (via gDirectory), so a further
            // TFile::Write() would re-write everything still attached to
            // the file, producing a duplicate ";2" cycle for every
            // histogram in the output.
            rFile->Close();
            std::cout << "[minion] Conversion successful!" << std::endl;
            return true;

        } catch (const H5::Exception& e) {
            std::cerr << "[minion] HDF5 Error: " << e.getDetailMsg() << std::endl;
            return false;
        }
    }

    void H5ToRootConverter::processGroup(H5::Group& h5_group, TDirectory* root_dir) {
        for (hsize_t i = 0; i < h5_group.getNumObjs(); ++i) {
            std::string obj_name = h5_group.getObjnameByIdx(i);
            H5G_obj_t type = h5_group.getObjTypeByIdx(i);

            if (type == H5G_GROUP) {
                H5::Group child_h5 = h5_group.openGroup(obj_name);
                TDirectory* child_root = root_dir->mkdir(obj_name.c_str());
                processGroup(child_h5, child_root);
            } 
            else if (type == H5G_DATASET) {
                convertDataset(h5_group, obj_name, root_dir);
            }
        }
    }

    void H5ToRootConverter::convertDataset(H5::Group& h5_group, const std::string& name, TDirectory* root_dir) {
        H5::DataSet ds = h5_group.openDataSet(name);
        H5::DataSpace space = ds.getSpace();
        int rank = space.getSimpleExtentNdims();
        hsize_t dims[2];
        space.getSimpleExtentDims(dims);

        root_dir->cd();

        if (rank == 1) {
            double min = readAttribute(ds, "min");
            double max = readAttribute(ds, "max");
            
            double entries = 0;
            try { entries = readAttribute(ds, "entries"); } catch(...) { entries = dims[0]; }

            size_t nBins = dims[0];
            std::vector<double> buffer(nBins);
            ds.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);

            std::string title = readStringAttribute(ds, "title");
            TH1D* h1 = new TH1D(name.c_str(), (title.empty() ? name : title).c_str(),
                                 static_cast<int>(nBins), min, max);
            for (size_t i = 0; i < nBins; ++i) h1->SetBinContent(static_cast<int>(i + 1), buffer[i]);
            h1->SetEntries(entries);

            // Building a histogram via SetBinContent() (as above) rather
            // than Fill() per raw value leaves ROOT to derive GetMean()/
            // GetStdDev() from bin contents x bin centers, only an
            // approximation for non-bin-aligned fill values (e.g. a
            // fractional cluster centroid). Warlock's Histogram1D separately
            // tracks the true (exact-value) first/second moments (see
            // PlotManager.cpp) and persists them as these attributes
            // whenever populated via fill()/fillBatch() (sum_w > 0) rather
            // than setData() (a derived/ratio histogram, where they aren't
            // meaningful); PutStats() injects them so GetMean()/GetStdDev()
            // match corryvreckan's native TH1::Fill()-tracked statistics
            // exactly instead of only to within one bin width.
            try {
                double sum_w = readAttribute(ds, "sum_w");
                if (sum_w > 0) {
                    double sum_wx = readAttribute(ds, "sum_wx");
                    double sum_wx2 = readAttribute(ds, "sum_wx2");
                    Double_t stats[4] = {sum_w, sum_w, sum_wx, sum_wx2};
                    h1->PutStats(stats);
                }
            } catch (...) {}

            h1->Write();
        } 
        else if (rank == 2) {
            int is_graph = 0;
            try { 
                H5::Attribute attr = ds.openAttribute("is_graph");
                attr.read(H5::PredType::NATIVE_INT, &is_graph);
            } catch(...) {}

            if (is_graph == 1) {
                size_t nPoints = dims[0];
                std::vector<double> buffer(nPoints * 2);
                ds.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);

                std::vector<double> vx(nPoints), vy(nPoints);
                for (size_t i = 0; i < nPoints; ++i) {
                    vx[i] = buffer[i * 2];
                    vy[i] = buffer[i * 2 + 1];
                }

                TGraph* gr = new TGraph(static_cast<int>(nPoints), vx.data(), vy.data());
                gr->SetName(name.c_str());
                gr->SetTitle(name.c_str());

                // Axis titles depend on the graph's actual kind, since not
                // every graph is a correction-vs-iteration plot (e.g.
                // AnalysisEfficiency's "efficiency_vs_event" has a real
                // event number on X and a unitless fraction on Y).
                if (name.find("alignment_correction") != std::string::npos) {
                    gr->GetXaxis()->SetTitle("# iteration");
                    if (name.find("rotation") != std::string::npos) {
                        gr->GetYaxis()->SetTitle("correction [deg]");
                    } else {
                        gr->GetYaxis()->SetTitle("correction [#mum]");
                    }
                } else if (name.find("efficiency_vs_event") != std::string::npos) {
                    gr->GetXaxis()->SetTitle("event number");
                    gr->GetYaxis()->SetTitle("efficiency");
                } else if (name.find("z_rotation_convergence") != std::string::npos) {
                    gr->GetXaxis()->SetTitle("# iteration");
                    gr->GetYaxis()->SetTitle("cumulative correction [deg]");
                } else if (name.find("nominal_pixels") != std::string::npos ||
                           name.find("measured_pixels") != std::string::npos) {
                    // AlignmentCAEN: one point per trusted pixel, local
                    // frame - overlay the nominal and measured graphs to
                    // see the tilt the rotation fit solved for directly. A
                    // discrete scatter of pixel positions has no meaningful
                    // "curve" between points, so use markers only (ROOT's
                    // default line-connects consecutive points, which draws
                    // a meaningless zigzag for insertion-ordered positions).
                    // Matches both the top-level z_rotation_nominal_pixels/
                    // z_rotation_measured_pixels (final state) and the
                    // per-iteration z_rotation_steps/iterNN/nominal_pixels/
                    // measured_pixels (bare names, no z_rotation_ prefix,
                    // since that context already lives in the group path).
                    gr->GetXaxis()->SetTitle("local x [mm]");
                    gr->GetYaxis()->SetTitle("local y [mm]");
                    gr->SetMarkerStyle(20);
                    gr->SetDrawOption("AP");

                    // A single-row or single-column board has every point
                    // sharing the exact same Y or X by construction. ROOT's
                    // default auto-range collapses a zero-variance axis to
                    // zero width, rendering the graph invisible, so pad any
                    // degenerate axis to a fixed minimum span.
                    auto padAxis = [](TAxis* axis, double lo, double hi) {
                        constexpr double kMinSpan = 0.1; // mm - comfortably larger than any of these boards' pitch
                        double span = hi - lo;
                        if (span < kMinSpan) {
                            double center = 0.5 * (lo + hi);
                            axis->SetLimits(center - kMinSpan / 2.0, center + kMinSpan / 2.0);
                        }
                    };
                    if (nPoints > 0) {
                        double xlo = *std::min_element(vx.begin(), vx.end());
                        double xhi = *std::max_element(vx.begin(), vx.end());
                        double ylo = *std::min_element(vy.begin(), vy.end());
                        double yhi = *std::max_element(vy.begin(), vy.end());
                        padAxis(gr->GetXaxis(), xlo, xhi);
                        padAxis(gr->GetYaxis(), ylo, yhi);
                    }
                } else {
                    std::cerr << "[minion] Warning: graph '" << name
                              << "' doesn't match any known kind - leaving its axes untitled "
                                 "rather than guessing wrong. Add a branch for it here." << std::endl;
                }

                gr->Write(name.c_str());
            } 
            else {
                double minX = readAttribute(ds, "min_x");
                double maxX = readAttribute(ds, "max_x");
                double minY = readAttribute(ds, "min_y");
                double maxY = readAttribute(ds, "max_y");
                
                double entries = 0;
                try { entries = readAttribute(ds, "entries"); } catch(...) { entries = dims[0]*dims[1]; }

                size_t binsX = dims[0];
                size_t binsY = dims[1];
                std::vector<double> buffer(binsX * binsY);
                ds.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);

                TH2D* h2 = new TH2D(name.c_str(), name.c_str(), 
                                    static_cast<int>(binsX), minX, maxX, 
                                    static_cast<int>(binsY), minY, maxY);
                
                for (size_t i = 0; i < binsX; ++i) {
                    for (size_t j = 0; j < binsY; ++j) {
                        h2->SetBinContent(static_cast<int>(i + 1), static_cast<int>(j + 1), buffer[i * binsY + j]);
                    }
                }
                h2->SetEntries(entries);

                // Same reasoning as the TH1D case above - inject the true
                // (exact-value) moments Warlock's Histogram2D tracked
                // alongside its binned content, when available.
                try {
                    double sum_w = readAttribute(ds, "sum_w");
                    if (sum_w > 0) {
                        double sum_x = readAttribute(ds, "sum_x");
                        double sum_y = readAttribute(ds, "sum_y");
                        double sum_x2 = readAttribute(ds, "sum_x2");
                        double sum_y2 = readAttribute(ds, "sum_y2");
                        double sum_xy = readAttribute(ds, "sum_xy");
                        Double_t stats[7] = {sum_w, sum_w, sum_x, sum_x2, sum_y, sum_y2, sum_xy};
                        h2->PutStats(stats);
                    }
                } catch (...) {}

                h2->Write();
            }
        }
    }

    double H5ToRootConverter::readAttribute(const H5::H5Object& obj, const std::string& name) {
        double val = 0;
        H5::Attribute attr = obj.openAttribute(name);
        attr.read(H5::PredType::NATIVE_DOUBLE, &val);
        return val;
    }

    std::string H5ToRootConverter::readStringAttribute(const H5::H5Object& obj, const std::string& name) {
        try {
            H5::Attribute attr = obj.openAttribute(name);
            std::string val;
            attr.read(attr.getStrType(), val);
            return val;
        } catch (...) {
            return "";
        }
    }
}