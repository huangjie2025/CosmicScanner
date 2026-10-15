#include "pheno_reader.h"
#include <cmath>
#include <fstream>
#include <stdexcept>
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::string path = argv[1];
    { std::ofstream f(path); f << "FID IID Y T\nF I 2.5 7\nF J 3.5 9\n"; }
    cosmic::PhenoReader reader;
    auto plain = reader.loadAllRecords(path, 1);
    if (plain.size()!=2 || plain[0].value!=2.5 || !std::isnan(plain[0].time)
        || plain[0].line_num!=2 || plain[1].line_num!=3)
        throw std::runtime_error("Non-longitudinal record fields are incorrect");
    auto longitudinal = reader.loadLongitudinalRecords(path, 2, 1);
    if (longitudinal.size()!=2 || longitudinal[0].value!=2.5 || longitudinal[0].time!=7
        || longitudinal[0].line_num!=2 || longitudinal[1].time!=9)
        throw std::runtime_error("Longitudinal record fields are incorrect");
}
