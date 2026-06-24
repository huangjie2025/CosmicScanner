#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "input/input_common.h"
#include "input/input_source_reader.h"
#include "input/time_normalizer.h"
#include "input/wide_to_long.h"

using namespace cosmic;

static void test_common() {
    assert(trim("  abc \t") == "abc");
    assert(isMissing("NA"));
    assert(isMissing("."));
    assert(!isMissing("0"));
    double x = 0.0;
    assert(parseDouble("12.5", x) && std::abs(x - 12.5) < 1e-12);
    assert(!parseDouble("12x", x));
    int y = 0, m = 0, d = 0;
    assert(parseDate("2024/01/03", y, m, d));
    assert(y == 2024 && m == 1 && d == 3);
}

static void test_source_reader() {
    const std::string path = "/tmp/opencode/input_units.csv";
    {
        std::ofstream out(path);
        out << "id,date,value\n";
        out << "I1,2024-01-01,10\n";
        out << "I2,2024-01-02,11\n";
    }
    Table t = readTable(path, "csv");
    assert(t.header.size() == 3);
    assert(t.rows.size() == 2);
    assert(getCell(t, t.rows[0], "id") == "I1");
    requireColumn(t, "date", "--time-col");
    bool threw = false;
    try { requireColumn(t, "missing", "--trait-cols"); } catch (const std::runtime_error& e) { threw = std::string(e.what()).find("MISSING_REQUIRED_COLUMN") != std::string::npos; }
    assert(threw);
}

static void test_time_normalizer() {
    auto dates = normalizeBaseTimes({"2024-01-01", "2024-01-03"}, "date");
    assert(dates.size() == 2 && dates[0] == 0.0 && dates[1] == 2.0);
    auto ints = normalizeBaseTimes({"3", "4"}, "integer");
    assert(ints[0] == 3.0 && ints[1] == 4.0);
    auto batches = normalizeBaseTimes({"B", "A", "B"}, "batch");
    assert(batches[0] == 0.0 && batches[1] == 1.0 && batches[2] == 0.0);
    assert(resolveTimeMultiplier("auto", 6) == 10);
    assert(resolveTimeMultiplier("100", 6) == 100);
    bool threw = false;
    try { (void)resolveTimeMultiplier("2", 6); } catch (const std::runtime_error& e) { threw = std::string(e.what()).find("INVALID_TIME_MULTIPLIER") != std::string::npos; }
    assert(threw);
}

static void test_wide_to_long() {
    auto a = parseWideSpec("volume_d12");
    assert(a.trait == "volume" && a.raw_time == "12");
    auto b = parseWideSpec("concentration_3");
    assert(b.trait == "concentration" && b.raw_time == "3");
    auto c = parseWideSpec("day5_volume");
    assert(c.trait == "volume" && c.raw_time == "5");
    auto d = parseWideSpec("motility_T2");
    assert(d.trait == "motility" && d.raw_time == "2");
    bool threw = false;
    try { (void)parseWideSpec("badcolumn"); } catch (const std::runtime_error& e) { threw = std::string(e.what()).find("WIDE_COLUMN_PARSE_FAILED") != std::string::npos; }
    assert(threw);
}

int main() {
    test_common();
    test_source_reader();
    test_time_normalizer();
    test_wide_to_long();
    std::cout << "Input adaptation unit tests passed.\n";
    return 0;
}
