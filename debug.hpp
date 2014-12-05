#include "Point.hpp"
#include<iostream>
#include<string>

/* Debug color
 * c: 30 - 37
 * black, red, green, yellow, blue, magenta, cyan, white
 */
template <typename S1>
void dbc(const S1 s1, int c) {
	std::cout << "\x1b[" << c << "m" << s1 << "\x1b[0m" << std::endl;
}
template <typename S1> // debug red
void dbr(const S1 s1) {
	std::cout << "\x1b[31m" << s1 << "\x1b[0m" << std::endl;
}
template <typename S1> // debug green
void dbg(const S1 s1) {
	std::cout << "\x1b[32m" << s1 << "\x1b[0m" << std::endl;
}
template <typename S1, typename S2 = std::string>
void db(const S1 s1,const S2 s2 = "") {
	std::cout << s1 << " " << s2 << std::endl;
}
template <typename S1>
void db(const S1 s1, const Point p) {
	std::cout << s1 << " " << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;
}
