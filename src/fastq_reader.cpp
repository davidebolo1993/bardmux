#include "fastq_reader.h"
#include <zlib.h>
#include <stdexcept>
#include <cstring>
#include <vector>
#include <unistd.h>

// ---------------------------------------------------------------------------
// Buffered line reader on top of gzFile
// ---------------------------------------------------------------------------
struct FastqReader::Impl {
    gzFile              gz      = nullptr;
    int                 buf_sz  = 1 << 20;   // 1 MiB gzip buffer
    std::vector<char>   in      = std::vector<char>(1 << 20);
    std::size_t         pos     = 0;
    std::size_t         end     = 0;

    bool refill() {
        int n = gzread(gz, in.data(), static_cast<unsigned>(in.size()));
        if (n < 0) {
            int err = Z_OK;
            const char* msg = gzerror(gz, &err);
            throw std::runtime_error(std::string("FASTQ read error: ") +
                                     (msg ? msg : "unknown zlib error"));
        }
        pos = 0;
        end = static_cast<std::size_t>(n);
        return end > 0;
    }

    bool getline(std::string& out) {
        out.clear();
        for (;;) {
            if (pos == end) {
                if (!refill()) return !out.empty();
            }

            char* start = in.data() + pos;
            std::size_t n = end - pos;
            char* nl = static_cast<char*>(std::memchr(start, '\n', n));
            if (nl) {
                std::size_t len = static_cast<std::size_t>(nl - start);
                out.append(start, len);
                pos += len + 1;
                if (!out.empty() && out.back() == '\r') out.pop_back();
                return true;
            }

            out.append(start, n);
            pos = end;
        }
    }
};

FastqReader::FastqReader(const std::string& filename) : impl_(new Impl()) {
    if (filename == "-") {
        int fd = ::dup(STDIN_FILENO);
        if (fd < 0)
            throw std::runtime_error("Cannot duplicate stdin for FASTQ input");
        impl_->gz = gzdopen(fd, "rb");
        if (!impl_->gz) {
            ::close(fd);
            throw std::runtime_error("Cannot open FASTQ from stdin");
        }
    } else {
        impl_->gz = gzopen(filename.c_str(), "rb");
    }
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
        if (!impl_->getline(line))
            throw std::runtime_error("Malformed FASTQ: truncated quality for read " + rec.name);
        rec.quality += line;
    }
    if (rec.quality.size() != rec.sequence.size())
        throw std::runtime_error("Malformed FASTQ: quality length mismatch for read " + rec.name);

    return true;
}
