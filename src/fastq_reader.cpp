#include "fastq_reader.h"
#include <zlib.h>
#include <stdexcept>
#include <cstring>

// ---------------------------------------------------------------------------
// Buffered line reader on top of gzFile
// ---------------------------------------------------------------------------
struct FastqReader::Impl {
    gzFile  gz      = nullptr;
    char*   buf     = nullptr;
    int     buf_sz  = 1 << 20;   // 1 MiB read buffer for gzread
    char*   line    = nullptr;
    int     line_sz = 1 << 17;   // 128 KiB per-line buffer (longest ONT read)

    bool getline(std::string& out) {
        out.clear();
        for (;;) {
            // Read one char at a time via gzgetc (gzFile is already buffered)
            int c = gzgetc(gz);
            if (c == -1) return !out.empty();  // EOF
            if (c == '\r') continue;           // skip Windows CR
            if (c == '\n') return true;
            out += static_cast<char>(c);
        }
    }
};

FastqReader::FastqReader(const std::string& filename) : impl_(new Impl()) {
    impl_->gz = gzopen(filename.c_str(), "rb");
    if (!impl_->gz)
        throw std::runtime_error("Cannot open FASTQ: " + filename);
    gzbuffer(impl_->gz, impl_->buf_sz);
}

FastqReader::~FastqReader() {
    if (impl_->gz) gzclose(impl_->gz);
    delete impl_;
}

bool FastqReader::next(FastqRecord& rec) {
    std::string line;

    // --- Find '@' header line (skip blank lines / stray content) ---
    for (;;) {
        if (!impl_->getline(line)) return false;   // EOF
        if (!line.empty() && line[0] == '@') break;
    }

    // Parse name: everything after '@' up to first space is the read ID;
    // the rest (if any) is the comment, appended after a space.
    {
        const char* p = line.c_str() + 1;  // skip '@'
        const char* sp = std::strchr(p, ' ');
        if (!sp) sp = std::strchr(p, '\t');
        if (sp) {
            rec.name.assign(p, sp - p);
            // skip the space, append comment
            rec.name += ' ';
            rec.name += (sp + 1);
        } else {
            rec.name.assign(p);
        }
    }

    // --- Sequence (may span multiple lines, until '+' line) ---
    rec.sequence.clear();
    for (;;) {
        if (!impl_->getline(line)) return false;
        if (!line.empty() && line[0] == '+') break;
        rec.sequence += line;
    }

    // --- Quality (read exactly as many chars as sequence length) ---
    rec.quality.clear();
    while (rec.quality.size() < rec.sequence.size()) {
        if (!impl_->getline(line)) break;
        rec.quality += line;
    }

    return true;
}
