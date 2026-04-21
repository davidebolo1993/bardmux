#pragma once
#include <string>

struct FastqRecord {
    std::string name;
    std::string sequence;
    std::string quality;
};

class FastqReader {
public:
    explicit FastqReader(const std::string& filename);
    ~FastqReader();

    bool next(FastqRecord& record);

private:
    struct Impl;
    Impl* impl_;
};
