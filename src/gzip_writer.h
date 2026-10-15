#pragma once
#include <string>
#include <fstream>
#include <memory>
#include <zlib.h>
#include <vector>
#include <sstream>
#include <cstring>

namespace cosmic {

class GzipOutputStream : public std::ostream {
public:
    explicit GzipOutputStream(const std::string& filename, int compression_level = Z_DEFAULT_COMPRESSION,
                              bool append = false);
    ~GzipOutputStream();
    
    GzipOutputStream(const GzipOutputStream&) = delete;
    GzipOutputStream& operator=(const GzipOutputStream&) = delete;
    
    GzipOutputStream(GzipOutputStream&& other) noexcept;
    GzipOutputStream& operator=(GzipOutputStream&& other) noexcept;
    
    void write(const char* data, size_t length);
    void write(const std::string& data);
    void flush();
    void close();
    
    bool isOpen() const;
    
    // Stream-like interface
    template<typename T>
    GzipOutputStream& operator<<(const T& value) {
        std::ostringstream oss;
        oss << value;
        write(oss.str().c_str(), oss.str().size());
        return *this;
    }
    
    GzipOutputStream& operator<<(std::ostream& (*manip)(std::ostream&)) {
        std::ostringstream oss;
        oss << manip;
        write(oss.str().c_str(), oss.str().size());
        return *this;
    }
    
    GzipOutputStream& operator<<(const char* str) {
        write(str, std::strlen(str));
        return *this;
    }
    
    GzipOutputStream& operator<<(const std::string& str) {
        write(str.c_str(), str.size());
        return *this;
    }
    
    GzipOutputStream& operator<<(char c) {
        write(&c, 1);
        return *this;
    }
    
    private:
    class GzipStreambuf : public std::streambuf {
    public:
        GzipStreambuf(const std::string& filename, int compression_level, bool append);
        ~GzipStreambuf();
        
        bool close();
        
    protected:
        virtual int_type overflow(int_type c) override;
        virtual int sync() override;
        
    public:
        bool isOpen() const { return gz_file_ != nullptr; }
        
    private:
        gzFile gz_file_ = nullptr;
        std::string filename_;
        char buffer_[8192];
    };
    
    std::unique_ptr<GzipStreambuf> streambuf_;
};

class GzipWriter {
public:
    static void writeFile(const std::string& filename, const std::string& content, int compression_level = Z_DEFAULT_COMPRESSION);
    static void writeFileLines(const std::string& filename, const std::vector<std::string>& lines, int compression_level = Z_DEFAULT_COMPRESSION);
    static void compressFile(const std::string& input_file, const std::string& output_file, int compression_level = Z_DEFAULT_COMPRESSION);
};

} // namespace cosmic
