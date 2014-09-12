#include "CS207/Util.hpp"
#include <forward_list>
#define FIRST_PRIME 2

/** Return true iff @a n is prime.
 * @pre @a n >= 0
 */
bool is_prime(int n)
{
  assert(n >= 0);
  static std::vector<int> known_primes = {FIRST_PRIME};
  auto it = known_primes.begin();
  int test_devisor = *it;
  while (1) {
	if (test_devisor * test_devisor > n) {
		// Place the element at the end of the vector
		if (test_devisor > FIRST_PRIME)
			known_primes.push_back(test_devisor);
		return true;
	}
  	else if (n % test_devisor == 0) {
		return false;
	}
	else {
		++it;
		if (it != known_primes.end()) {
			test_devisor = *it;
		} 
		else {
			test_devisor = test_devisor + 1;
		}
	}
  }
}

int main()
{
  while (!std::cin.eof()) {
    // How many primes to test? And should we print them?
    std::cerr << "Input Number: ";
    int n = 0;
    CS207::getline_parsed(std::cin, n);
    if (n <= 0)
      break;

    std::cerr << "Print Primes (y/n): ";
    char confirm = 'n';
    CS207::getline_parsed(std::cin, confirm);
    bool print_primes = (confirm == 'y' || confirm == 'Y');

    CS207::Clock timer;

    // Loop and count primes from 2 up to n
    int num_primes = 0;
    for (int i = 2; i <= n; ++i) {
      if (is_prime(i)) {
        ++num_primes;
        if (print_primes)
          std::cout << i << std::endl;
      }
    }

    double elapsed_time = timer.seconds();

    std::cout << "There are " << num_primes
              << " primes less than or equal to " << n << ".\n"
              << "Found in " << (1000 * elapsed_time) << " milliseconds.\n\n";
  }

  return 0;
}
