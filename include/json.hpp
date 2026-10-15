#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>
#include <type_traits>
#include <iomanip>


struct Json {
    std::stringstream s;
    bool first = true;

    Json() {
        s << "{\n";
    }

    template<typename T>
    Json& add(const std::string& name, const T& x) {
        if (first) {
            first = false;
        }
        else {
            s << ",\n";
        }

        indent(s, 1);
        write_string(s, name);
        s << ": ";
        write_val(s, x, 1);
        return *this;
    }

    std::string str() const {
        return s.str() + "\n}";
    }

    private:
        static void indent(std::stringstream& out, int depth) {
            for (int i = 0; i < depth; ++i) {
                out << "  ";
            }
        }

        static void write_string(std::stringstream& out, const std::string& x) {
            out << "\"";
            for (char c : x) {
                switch (c) {
                    case '"' : out << "\\\""; break;
                    case '\\': out << "\\\\"; break;
                    case '\n': out << "\\n"; break;
                    case '\t': out << "\\t"; break;
                    case '\r': out << "\\r"; break;
                    default  : out << c;      break;
                }
            }
            out << "\"";
        }

        static void write_val(std::stringstream& out, const std::string& x, int) {
            write_string(out, x);
        }

        static void write_val(std::stringstream& out, const char *x, int) {
            write_string(out, x);
        }

        static void write_val(std::stringstream& out, const Json& x, int depth) {
            std::string str = x.str();
            for (size_t i = 0; i < str.size(); ++i) {
                out << str[i];
                if ((str[i] == '\n') && (i != str.size() - 1)) {
                    indent(out, depth+1);
                }
            }
        }

        static void write_val(std::stringstream& out, bool x, int) {
            out << (x ? "true" : "false");
        }

        template<typename T>
        static std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool>, void>
        write_val(std::stringstream& out, T x, int) {
            out << x;
        }

        template<typename T>
        static std::enable_if_t<std::is_floating_point_v<T>, void>
        write_val(std::stringstream& out, const T& x, int) {
            if (!std::isfinite(x)) { throw std::runtime_error("NaN and infinity are not valid JSON"); }
            out << std::setprecision(std::numeric_limits<T>::max_digits10) << x;
        }

        template<typename T>
        static void write_val(std::stringstream& out, const std::vector<T>& x, int depth) {
            out << '[';
            for (size_t i = 0; i < x.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                write_val(out, x[i], depth+1);
            }
            out << ']';
        }

        template<typename V>
        static void write_val(std::stringstream& out, const std::unordered_map<std::string, V>& xs, int depth) {
            out << '{';

            bool first = true;
            for (const auto& [k, v] : xs) {
                if (first) {
                    first = false;
                }
                else {
                    out << ", ";
                }
                write_string(out, k);
                out << ": ";
                write_val(out, v, depth+1);
            }
            out << '}';
        }
};

