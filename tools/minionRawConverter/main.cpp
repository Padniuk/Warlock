/**
 * @file main.cpp
 * @brief CLI entry point for minionRawConverter - see RawConverter.hpp and
 * the usage text below for what each flag controls.
 */

#include "RawConverter.hpp"
#include <iostream>
#include <unistd.h>

int main(int argc, char** argv) {
    framework::ConverterConfig cfg;
    cfg.output_file = "output.h5";
    cfg.max_events = 0;

    int opt;
    while ((opt = getopt(argc, argv, "i:o:n:x:s:w:t:m:")) != -1) {
        switch (opt) {
            case 'i': cfg.input_file = optarg; break;
            case 'o': cfg.output_file = optarg; break;
            case 'n': cfg.max_events = std::stoul(optarg); break;
            case 'x': cfg.skip_descriptions.push_back(optarg); break;
            case 's': cfg.serial_descriptions.push_back(optarg); break;
            case 'w': cfg.warmup_descriptions.push_back(optarg); break;
            case 't': cfg.threads = std::stoul(optarg); break;
            case 'm': cfg.master_description = optarg; break;
            default:
                std::cerr << "Usage: " << argv[0]
                          << " -i <in.raw> [-o <out.h5>] [-n <limit>]\n"
                          << "  [-x <description>]...  skip decoding this producer entirely\n"
                          << "  [-s <description>]...  force this producer to decode sequentially,\n"
                          << "                         PERMANENTLY, even under -t. Use for a converter\n"
                          << "                         whose thread-safety problem is ongoing/per-event\n"
                          << "                         (e.g. CMSIT: compares each event's trigger ID to\n"
                          << "                         the previous one - needs strict order forever).\n"
                          << "  [-w <description>]...  force this producer to decode sequentially for the\n"
                          << "                         file's FIRST decode batch only, then allow it into\n"
                          << "                         the thread pool for every batch after that. Use for\n"
                          << "                         a converter whose only shared-state problem is\n"
                          << "                         one-time lazy setup (e.g. CAEN DT5742: its\n"
                          << "                         correction-table load is unsynchronized, but stable\n"
                          << "                         and read-only once loaded, and the actual per-\n"
                          << "                         sample correction math is reentrant).\n"
                          << "  [-t <threads>]         decode threads for everything not forced serial;\n"
                          << "                         defaults to 1 (fully sequential) so a plain run\n"
                          << "                         with no flags is always safe for any producer mix\n"
                          << "  [-m <description>]     anchor window formation on this raw producer\n"
                          << "                         description, matching corry's sync_by_event\n"
                          << "                         (whichever detector loader is declared first in\n"
                          << "                         its .conf defines each event's boundary). Every\n"
                          << "                         other producer's (sub-)event is matched against\n"
                          << "                         that boundary instead of all producers being\n"
                          << "                         treated symmetrically. Omit to keep the old\n"
                          << "                         symmetric windowing (all streams' numbers must\n"
                          << "                         agree directly).\n"
                          << "  The converter logs every distinct raw event description it encounters once,\n"
                          << "  at startup - run without -s/-w/-x first to see what strings this file\n"
                          << "  actually uses before deciding what to pass.\n";
                return 1;
        }
    }

    if (cfg.input_file.empty()) {
        std::cerr << "[minion] Error: Input file (-i) is required." << std::endl;
        return 1;
    }

    try {
        framework::RawConverter converter(cfg);
        converter.run();
    } catch (const std::exception& e) {
        std::cerr << "[minion] Fatal error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}