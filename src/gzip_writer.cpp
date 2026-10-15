#include "gzip_writer.h"
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <memory>
#include <cstring>

namespace cosmic {

// GzipStreambuf implementation
GzipOutputStream::GzipStreambuf::GzipStreambuf(const std::string& filename, int compression_level, bool append)
    : filename_(filename) {
    std::string mode = std::string(append ? "ab" : "wb") + std::to_string(compression_level);
    gz_file_ = gzopen(filename.c_str(), mode.c_str());
    if (!gz_file_) {
        throw std::runtime_error("Failed to open gzip file for writing: " + filename);
    }
    setp(buffer_, buffer_ + sizeof(buffer_));
}

GzipOutputStream::GzipStreambuf::~GzipStreambuf() {
    close();
}

bool GzipOutputStream::GzipStreambuf::close() {
    if (gz_file_) {
        const bool flushed = sync() == 0;
        const int status = gzclose(gz_file_);
        gz_file_ = nullptr;
        return flushed && status == Z_OK;
    }
    return true;
}

GzipOutputStream::GzipStreambuf::int_type GzipOutputStream::GzipStreambuf::overflow(int_type c) {
    if (c == traits_type::eof()) {
        return sync() == 0 ? traits_type::eof() : traits_type::eof();
    }
    // Flush the full buffer before writing the new character
    if (sync() != 0) return traits_type::eof();
    *pptr() = traits_type::to_char_type(c);
    pbump(1);
    return c;
}

int GzipOutputStream::GzipStreambuf::sync() {
    if (!gz_file_) return -1;
    std::streamsize n = pptr() - pbase();
    if (n > 0) {
        int written = gzwrite(gz_file_, pbase(), static_cast<unsigned int>(n));
        if (written != n) return -1;
        setp(buffer_, buffer_ + sizeof(buffer_));
    }
    return 0;
}

GzipOutputStream::GzipOutputStream(const std::string& filename, int compression_level, bool append)
    : std::ostream(nullptr), streambuf_(std::make_unique<GzipStreambuf>(filename, compression_level, append)) {
    rdbuf(streambuf_.get());
}

GzipOutputStream::~GzipOutputStream() {
    // The streambuf will be closed automatically via unique_ptr destructor
}

GzipOutputStream::GzipOutputStream(GzipOutputStream&& other) noexcept 
    : std::ostream(nullptr), streambuf_(std::move(other.streambuf_)) {
    rdbuf(streambuf_.get());
}

GzipOutputStream& GzipOutputStream::operator=(GzipOutputStream&& other) noexcept {
    if (this != &other) {
        streambuf_ = std::move(other.streambuf_);
        rdbuf(streambuf_.get());
    }
    return *this;
}

void GzipOutputStream::write(const char* data, size_t length) {
    if (!rdbuf() || rdbuf()->sputn(data, length) != static_cast<std::streamsize>(length))
        setstate(std::ios::badbit);
}

void GzipOutputStream::write(const std::string& data) {
    write(data.data(), data.size());
}

void GzipOutputStream::flush() {
    if (!rdbuf() || rdbuf()->pubsync() != 0) setstate(std::ios::badbit);
}

void GzipOutputStream::close() {
    if (!streambuf_ || !streambuf_->close()) setstate(std::ios::badbit);
    if (!*this) throw std::runtime_error("Failed while closing gzip output");
}

bool GzipOutputStream::isOpen() const {
    return streambuf_ != nullptr && streambuf_->isOpen();
}

void GzipWriter::writeFile(const std::string& filename, const std::string& content, int compression_level) {
    GzipOutputStream out(filename, compression_level);
    out << content;
    out.close();
}

void GzipWriter::writeFileLines(const std::string& filename, const std::vector<std::string>& lines, int compression_level) {
    GzipOutputStream out(filename, compression_level);
    for (const auto& line : lines) {
        out << line << "\n";
    }
    out.close();
}

void GzipWriter::compressFile(const std::string& input_file, const std::string& output_file, int compression_level) {
    GzipOutputStream out(output_file, compression_level);
    std::ifstream in(input_file, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open input file: " + input_file);
    
    constexpr size_t buffer_size = 8192;
    char buffer[buffer_size];
    
    while (in.read(buffer, buffer_size) || in.gcount() > 0) {
        out.write(buffer, in.gcount());
    }
    out.close();
}

} // namespace cosmic
