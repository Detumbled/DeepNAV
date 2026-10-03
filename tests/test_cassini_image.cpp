#include "opnav/image/CassiniImage.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace img = fd::opnav::image;
namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F operation) {
    try { operation(); } catch (const std::exception&) { return; }
    throw std::runtime_error("Invalid image was accepted.");
}
struct Fixture {
    std::filesystem::path directory = std::filesystem::temp_directory_path()/
        ("deepnav-cassini-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::path image = directory/"test.IMG", label = directory/"test.LBL";
    Fixture() { std::filesystem::create_directory(directory); }
    ~Fixture() { std::error_code ignored; std::filesystem::remove_all(directory,ignored); }
    void write(std::string format = "HIGH", bool truncated = false, int labelSamples = 3) {
        const std::string header = "LBLSIZE=256 FORMAT='HALF' ORG='BSQ' INTFMT='"+format+
            "' NB=1 NL=2 NS=3 NBB=122 NLB=1 EOL=0 RECSIZE=128 ";
        std::vector<unsigned char> bytes(5*128,0xab);
        std::fill(bytes.begin(),bytes.begin()+256,' ');
        std::copy(header.begin(),header.end(),bytes.begin());
        const int values[] = {-2,1000,4095,32767,-32768,42};
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 3; ++x) {
                const unsigned v = unsigned(values[y*3+x]) & 0xffff;
                const auto offset = 3*128+y*128+122+2*x;
                bytes[offset] = v >> 8; bytes[offset+1] = v & 255;
            }
        std::ofstream binary(image,std::ios::binary);
        binary.write(reinterpret_cast<const char*>(bytes.data()),bytes.size()-(truncated ? 1 : 0));
        std::ofstream text(label);
        text << "PDS_VERSION_ID = PDS3\nRECORD_TYPE = FIXED_LENGTH\nRECORD_BYTES = 128\nFILE_RECORDS = 5\n"
             << "^IMAGE = (\"test.IMG\",4)\nMISSING_LINES = 0\nMISSING_PACKET_FLAG = \"NO\"\n"
             << "DATA_CONVERSION_TYPE = \"12BIT\"\nIMAGE_MID_TIME = 2004-284T18:13:06.046\n"
             << "TARGET_DESC = \"ENCELADUS\"\nGAIN_MODE_ID = \"29 ELECTRONS PER DN\"\nEXPOSURE_DURATION = 60\n"
             << "OBJECT = IMAGE\nLINES = 2\nLINE_SAMPLES = " << labelSamples
             << "\nSAMPLE_BITS = 16\nSAMPLE_TYPE = SUN_INTEGER\nLINE_PREFIX_BYTES = 122\nEND_OBJECT = IMAGE\nEND\n";
    }
};
}
int main() {
    try {
        Fixture fixture;
        fixture.write();
        const auto loaded = img::loadCassiniImage(fixture.image);
        require(loaded.dn.rows()==2 && loaded.dn.cols()==3,"Wrong image shape.");
        require(loaded.dn(0,0)==-2 && loaded.dn(0,1)==1000 && loaded.dn(0,2)==4095 &&
                loaded.dn(1,0)==32767 && loaded.dn(1,1)==-32768 && loaded.dn(1,2)==42,
                "Endian decoding, signed samples, prefix skipping, or row order failed.");
        require(loaded.electronsPerDn==29 && loaded.exposureMilliseconds==60 &&
                loaded.targetDescription=="ENCELADUS" && loaded.imageMidTimeUtc=="2004-284T18:13:06.046",
                "Metadata decoding failed.");
        std::cout << "PASS: big-endian signed pixels, header/telemetry/prefix skipping, metadata\n";
        fixture.write("LOW"); rejects([&]{ (void)img::loadCassiniImage(fixture.image); });
        fixture.write("HIGH",true); rejects([&]{ (void)img::loadCassiniImage(fixture.image); });
        fixture.write("HIGH",false,4); rejects([&]{ (void)img::loadCassiniImage(fixture.image); });
        fixture.write();
        { std::fstream bytes(fixture.image,std::ios::in|std::ios::out|std::ios::binary); bytes << "JFIF JPEG preview"; }
        rejects([&]{ (void)img::loadCassiniImage(fixture.image); });
        rejects([&]{ (void)img::loadCassiniImage(fixture.image,fixture.directory/"missing.LBL"); });
        std::cout << "PASS: unsupported endian layout, truncation, inconsistent label, renamed JPEG, missing label\n";
        const auto real = std::filesystem::path(DEEPNAV_SOURCE_DIR)/"Input Optical/N1476124698_2.IMG";
        if (std::filesystem::exists(real)) {
            const auto image = img::loadCassiniImage(real);
            require(image.dn.rows()==1024 && image.dn.cols()==1024 && image.dn.minCoeff()==4 &&
                    image.dn.maxCoeff()==1670 && image.dn(569,687)==1670,"Real Cassini pixels do not match inspection.");
            std::cout << "PASS: original Cassini fixture, 1024x1024 DN and known peak\n";
        }
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n'; return 1;
    }
}
