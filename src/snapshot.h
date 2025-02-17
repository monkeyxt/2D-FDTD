#ifndef EMSOLVER_SNAPSHOT_H
#define EMSOLVER_SNAPSHOT_H

#include <vector>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <iterator>

#include "constants.h"

///============================================================================
/// Snapshot class definition
///============================================================================
template <FloatingPoint T>
class Snapshot {
public:
    Snapshot(const std::string& directoryName_ = ".");
    ~Snapshot() = default;

    Snapshot(const Snapshot&) = delete;
    Snapshot(Snapshot&&) = delete;
    Snapshot& operator=(const Snapshot&) = delete;
    Snapshot& operator=(Snapshot&&) = delete;

    template <typename... Vectors>
    void write(const std::vector<std::string>& names,
               const Vectors&... vectors);

private:
    std::string directoryName;
};

///============================================================================
/// Snapshot constructor
///============================================================================
template <FloatingPoint T>
Snapshot<T>::Snapshot(const std::string& directoryName_) 
    : directoryName(directoryName_) {
    std::filesystem::path dir(directoryName);
    if(!std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
    }
}

///============================================================================
/// Snapshot write function
///============================================================================
template <FloatingPoint T>
template <typename... Vectors>
void Snapshot<T>::write(const std::vector<std::string>& names, 
                        const Vectors&... vectors) {
    static int frameNumber = 0;
    const int currentFrame = frameNumber; // Create local copy
    std::cout << "Writing snapshot " << currentFrame << "..." << std::endl;
    auto writeArray = [this, currentFrame](const auto& vector, 
                                           const std::string& name) {
        std::string filename = directoryName + "/" + name + "_" 
            + std::to_string(currentFrame) + ".dat";
        std::ofstream file(filename);
        
        if (!file.is_open()) {
            throw std::runtime_error("Could not open output file " + filename + 
                                   " for frame " + std::to_string(currentFrame));
        }

        // Handle 2D vectors by writing each row on a new line
        if constexpr (std::is_same_v<typename std::decay_t<decltype(vector)>::value_type, 
                                    std::vector<T>>) {
            for (const auto& row : vector) {
                std::copy(row.begin(), row.end(), 
                         std::ostream_iterator<T>(file, " "));
                file << "\n";
            }
        } else {
            // Original 1D vector handling
            std::copy(vector.begin(), vector.end(), 
                     std::ostream_iterator<T>(file, " "));
            file << "\n";
        }
        file.close();
    };

    // Verify we have the right number of names
    if (names.size() != sizeof...(vectors)) {
        throw std::runtime_error("Number of names must match number of vectors");
    }

    // Unpack and write each array with its corresponding name
    int arrayIndex = 0;
    (writeArray(vectors, names[arrayIndex++]), ...);

    frameNumber++;
}

#endif