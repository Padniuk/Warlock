#include "core/FrameworkManager.hpp"
#include "core/utils/Logger.hpp"
#include <iostream>
#include <vector>
#include <string>

#ifdef _WIN32
#include <windows.h>
namespace {
    // The dynamic progress bar (FrameworkManager) writes raw ANSI/VT100
    // escape codes to move the cursor and clear lines. Modern Windows
    // Terminal understands these out of the box, but the legacy conhost.exe
    // console only does if a process explicitly opts in via this flag - a
    // straight `\033[...` write there otherwise just prints garbled escape
    // sequences instead of animating the bar. No-op on Linux/macOS.
    void enableAnsiOnWindows() {
        HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
        if (handle == INVALID_HANDLE_VALUE) return;
        DWORD mode = 0;
        if (!GetConsoleMode(handle, &mode)) return;
        SetConsoleMode(handle, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
    }
}
#endif

int main(int argc, char** argv) {
#ifdef _WIN32
    enableAnsiOnWindows();
#endif
    std::string config_file;
    std::vector<std::string> overrides;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-c" && i + 1 < argc) {
            config_file = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            overrides.push_back(argv[++i]);
        } else if (config_file.empty() && arg[0] != '-') {
            config_file = arg;
        }
    }

    if (config_file.empty()) {
        std::cerr << "Usage: ./warlock -c <config.toml> [-o Module.key=value]\n";
        return 1;
    }

    try {
        framework::FrameworkManager manager;
        manager.loadConfiguration(config_file, overrides);
        
        std::cout << "\n"
                  << "   _      __         __         __ \n"
                  << "  | | /| / /__ _____/ /__  ____/ /__\n"
                  << "  | |/ |/ / _ `/ __/ / _ \\/ __/  '_/\n"
                  << "  |__/|__/\\_,_/_/ /_/\\___/\\__/_/\\_\\ \n"
                  << "  Dark magic for the reconstruction\n" << std::endl;

        manager.run();
        std::cout << "\n[Warlock] Reconstruction completed successfully." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\n[FATAL ERROR] " << e.what() << std::endl;
        return 1;
    }

    return 0;
}