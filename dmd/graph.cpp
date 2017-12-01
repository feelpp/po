#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

using floating_precision_t = double;

const auto pi = std::acos(-1.0);

struct dmd_mode_t {
    using eigenvalue_t = std::complex<floating_precision_t>;
    eigenvalue_t eigenvalue;
    floating_precision_t frequency;
    floating_precision_t growth_rate;
    floating_precision_t energy;
    floating_precision_t power;

    dmd_mode_t() = delete;
    // constructor from eigenvalue and time step
    // no member initializer list for clarity
    dmd_mode_t(eigenvalue_t val, floating_precision_t delta_time) {
        eigenvalue = val;
        frequency = (std::log(val)).imag() / (2. * pi * delta_time);
        growth_rate = (std::log(val)).real() / delta_time;
        energy = std::norm(val);
        power = (val * std::conj(val)).real();
    }
};

using modes_t = std::vector<dmd_mode_t>;
using vector_t = std::vector<floating_precision_t>;

void usage(char* arg0) { std::cout << "usage: " << arg0 << " dt" << std::endl; }

auto parse_line(const std::string& line, floating_precision_t time_step)
    -> std::optional<dmd_mode_t> {
    std::string current;
    std::istringstream stream{line};

    // "Mode"
    stream >> current;
    if (current.compare("Mode") != 0)
        return std::nullopt;

    // mode number
    stream >> current;

    // eigenvalue string
    stream >> current;
    current = current.substr(std::strlen("eigenValues=("), current.size());

    std::string buffer_real;
    std::string buffer_imaginary;

    // extract real part substring
    for (const char c : current) {
        if (c == ',')
            break;
        buffer_real += c;
    }

    // extract imaginary part substring
    current = current.substr(buffer_real.size() + 1, current.size());
    for (const char c : current) {
        if (c == ')')
            break;
        buffer_imaginary += c;
    }

    // deduce timestep from frequency
    stream >> current;
    current = current.substr(std::strlen("frequency="), current.size());

    std::string buffer_frequency;
    for (const char c : current) {
        if (c == 'H') // Hz;
            break;
        buffer_frequency += c;
    }

    floating_precision_t real, imaginary;

    try {
        real = std::stod(buffer_real);
        imaginary = std::stod(buffer_imaginary);
    } catch (...) {
        std::cerr << "error converting string to number" << std::endl;
        return std::nullopt;
    }

    return dmd_mode_t{std::complex<floating_precision_t>{real, imaginary},
                      time_step};
}

void generate_output(const modes_t& modes_to_plot) {
    std::cout << "import numpy as np\n"
              << "import matplotlib.pyplot as plt\n";

    // array of frequencies
    std::cout << "frequencies = [";
    for (size_t i = 0; i < modes_to_plot.size(); i++) {
        std::cout << modes_to_plot[i].frequency;
        if (i != modes_to_plot.size() - 1)
            std::cout << ", ";
    }
    std::cout << "]\n";

    // array of powers
    std::cout << "powers = [";
    for (size_t i = 0; i < modes_to_plot.size(); i++) {
        std::cout << modes_to_plot[i].power;
        if (i != modes_to_plot.size() - 1)
            std::cout << ", ";
    }
    std::cout << "]\n";

    // array of growth rates
    std::cout << "growth_rates = [";
    for (size_t i = 0; i < modes_to_plot.size(); i++) {
        std::cout << modes_to_plot[i].growth_rate;
        if (i != modes_to_plot.size() - 1)
            std::cout << ", ";
    }
    std::cout << "]\n"
                 "fig, power_axis = plt.subplots()\n"
                 "growth_axis = power_axis.twinx()\n"
                 "power_axis.plot(frequencies, powers, color='r')\n"
                 "power_axis.set_xlabel('frequency')\n"
                 "power_axis.set_ylabel('power')\n"
                 "growth_axis.plot(frequencies, growth_rates, color='b')\n"
                 "growth_axis.set_ylabel('growth')\n"
                 "def color_y(ax, color):\n"
                 "\tfor t in ax.get_yticklabels():\n"
                 "\t\tt.set_color(color)\n"
                 "\treturn None\n"
                 "color_y(power_axis, 'r')\n"
                 "color_y(growth_axis, 'b')\n"
                 "plt.show()\n"
              << std::endl;
}

int main(int argc, char** argv) {
    modes_t modes;
    std::string buffer;
    double time_step;

    if (argc != 2)
        usage(argv[0]);

    try {
        time_step = std::stod(argv[1]);
    } catch (...) {
        std::cout << "error while reading the time step, make sure you "
                     "entered a valid double precision floating point number"
                  << std::endl;
        return -1;
    }

    // skip standard input until the part we want
    for (std::getline(std::cin, buffer);
         buffer.compare("Modes projection and normalization");
         std::getline(std::cin, buffer))
        ;

    for (std::string s; std::getline(std::cin, s);) {
        if (auto mode = parse_line(s, time_step); mode.has_value()) {
            modes.push_back(*mode);
            std::cerr << "adding mode (" << mode->eigenvalue.real() << ", "
                      << mode->eigenvalue.imag()
                      << ") with f = " << mode->frequency
                      << ", sigma = " << mode->growth_rate
                      << ", alpha = " << mode->energy
                      << " and P = " << mode->power << '\n';
        }
    }

    modes_t sorted_modes;

    std::copy_if(std::begin(modes), std::end(modes),
                 std::back_inserter(sorted_modes),
                 [](const auto& mode) { return mode.frequency > -1e-6; });
    // sort modes by frequency
    const auto compare_by_frequency = [](const auto& mode1, const auto& mode2) {
        return mode1.frequency < mode2.frequency;
    };
    std::sort(std::begin(sorted_modes), std::end(sorted_modes),
              compare_by_frequency);

    generate_output(sorted_modes);

    return 0;
}
