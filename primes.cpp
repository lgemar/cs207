#include "CS207/Util.hpp"
#include <forward_list>

/** Return true iff @a n is prime.
 * @pre @a n >= 0
 */
bool is_prime(int n)
{
  assert(n >= 0);
  static std::forward_list<int> known_primes = { 2 };
  int latest_prime;
  auto it = known_primes.begin();
  latest_prime = *it;
  while (1) {
	if (latest_prime * latest_prime > n) {
		known_primes.insert_after(it, n);
		return true;
	}
  	else if (n % latest_prime == 0) {
		return false;
	}
	else {
		if (it != known_primes.end()) {
			++it;
			latest_prime = *it;
		} 
		else {
			latest_prime = latest_prime + 1;
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
