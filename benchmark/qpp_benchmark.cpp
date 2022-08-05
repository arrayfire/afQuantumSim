///////////////////////////////////////////////////////////
// Code for the algorithms taken from Quantum++ examples
// See License in docs/licenses/Quantum++
//////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>
#include <qpp/qpp.h>

#include <chrono>

void grover(uint32_t qubits)
{
    qpp::idx n = qubits; // number of qubits

    std::vector<qpp::idx> dims(n, 2); // local dimensions
    std::vector<qpp::idx> subsys(n);  // ordered subsystems
    std::iota(std::begin(subsys), std::end(subsys), 0);

    // number of elements in the database
    auto N = static_cast<qpp::idx>(std::llround(std::pow(2, n)));

    // mark an element randomly
    qpp::idx marked = qpp::randidx(0, N - 1);

    qpp::ket psi = qpp::mket(std::vector<qpp::idx>(n, 0)); // computational |0>^\otimes n

    // apply H^\otimes n, no aliasing
    psi = (qpp::kronpow(qpp::gt.H, n) * psi).eval();

    qpp::cmat G = 2 * qpp::prj(psi) - qpp::gt.Id(N); // Diffusion operator

    // number of queries
    auto nqueries = static_cast<qpp::idx>(std::ceil(qpp::pi / 4 * std::sqrt(N)));

    for (qpp::idx i = 0; i < nqueries; ++i) {
        psi(marked) = -psi(marked); // apply the oracle first, no aliasing
        psi = (G * psi).eval();     // then the diffusion operator, no aliasing
    }

    // we now measure the state in the computational basis, destructively
    auto measured = qpp::measure_seq(psi, subsys, dims);

    // sample
    auto result = std::get<qpp::RES>(measured);
}

void qft(uint32_t qubits)
{
    using namespace qpp;

    std::vector<idx> init_qubits(qubits); // initial state

    ket psi = mket(init_qubits);
    ket result = psi;

    idx n = init_qubits.size();                                   // number of qubits
    auto D = static_cast<idx>(std::llround(std::pow(2, n))); // dimension 2^n

    for (idx i = 0; i < n; ++i) {
        result = apply(result, gt.H, {i}); // apply Hadamard on qubit 'i'
        // apply controlled rotations
        for (idx j = 2; j <= n - i; ++j) {
            cmat Rj(2, 2);
            auto pow_j = static_cast<idx>(std::llround(std::pow(2, j)));
            Rj << 1, 0, 0, omega(pow_j);
            result = applyCTRL(result, Rj, {i + j - 1}, {i});
        }
    }

    // we have the qubits in reversed order, we must swap them
    for (idx i = 0; i < n / 2; ++i)
        result = apply(result, gt.SWAP, {i, n - i - 1});
}

void shor(uint32_t qubits)
{
    using namespace qpp;

    bigint N = 15;                                  // number to factor
    bigint a = rand(static_cast<bigint>(3), N - 1); // random co-prime with N
    while (gcd(a, N) != 1) {
        a = rand(static_cast<bigint>(3), N - 1);
    }
    // qubits required for half of the circuit, in total we need 2n qubits
    // if you know the order 'r' of 'a', then you can take the smallest 'n' s.t.
    // 2^n >= 2 * r^2, i.e., n = ceil(log2(2 * r^2))
    auto n = static_cast<idx>(std::ceil(2 * std::log2(N)));
    auto D = static_cast<idx>(std::llround(std::pow(2, n))); // dimension 2^n

    // vector with labels of the first half of the qubits
    std::vector<idx> first_subsys(n);
    std::iota(std::begin(first_subsys), std::end(first_subsys), 0);

    // vector with labels of the second half of the qubits
    std::vector<idx> second_subsys(n);
    std::iota(std::begin(second_subsys), std::end(second_subsys), n);

    // QUANTUM STAGE
    // prepare the initial state |0>^\otimes n \otimes |0...01>
    ket psi = kron(st.zero(2 * n - 1), 1_ket);

    // apply Hadamards H^\otimes n on first half of the qubits
    for (idx i = 0; i < n; ++i) {
        psi = apply(psi, gt.H, {i});
    }

    // perform the modular exponentiation as a sequence of
    // modular multiplications
    for (idx i = 0; i < n; ++i) {
        // compute 2^(n-i-1) mod N
        idx j = static_cast<idx>(std::llround(std::pow(2, n - i - 1)));
        // compute the a^(2^(n-i-1)) mod N
        idx aj = modpow(a, j, N);
        // apply the controlled modular multiplication
        psi = applyCTRL(psi, gt.MODMUL(aj, N, n), {i}, second_subsys);
    }

    // apply inverse QFT on first half of the qubits
    psi = applyTFQ(psi, first_subsys);
    // END QUANTUM STAGE

    // FIRST MEASUREMENT STAGE
    auto measured1 = measure_seq(psi, first_subsys); // measure first n qubits
    std::vector<idx> vect_results1 = std::get<RES>(measured1); // results
    double prob1 = std::get<PROB>(measured1); // probability of the result
    idx n1 = multiidx2n(vect_results1, std::vector<idx>(n, 2)); // binary to int
    auto x1 =
        static_cast<double>(n1) / static_cast<double>(D); // multiple of 1/r

    bool failed = true;
    idx r1 = 0, c1 = 0;
    for (auto&& elem : convergents(x1, 10)) {
        std::tie(c1, r1) = elem;
        auto c1r1 = static_cast<double>(c1) / r1;
        if (abs(x1 - c1r1) < 1. / std::pow(2, (n - 1.) / 2.)) {
            failed = false;
            break;
        }
    }
    if (failed) {
        return;
    }
    // END FIRST MEASUREMENT STAGE

    // SECOND MEASUREMENT STAGE
    auto measured2 = measure_seq(psi, first_subsys); // measure first n qubits
    std::vector<idx> vect_results2 = std::get<RES>(measured2); // results
    double prob2 = std::get<PROB>(measured2); // probability of the result
    idx n2 = multiidx2n(vect_results2, std::vector<idx>(n, 2)); // binary to int
    auto x2 = static_cast<double>(n2) / D; // multiple of 1/r

    failed = true;
    idx r2 = 0, c2 = 0;
    for (auto&& elem : convergents(x2, 10)) {
        std::tie(c2, r2) = elem;
        auto c2r2 = static_cast<double>(c2) / r2;
        if (abs(x2 - c2r2) < 1. / std::pow(2, (n - 1.) / 2.)) {
            failed = false;
            break;
        }
    }
    if (failed) {
        return;
    }
    // END SECOND MEASUREMENT STAGE

    // THIRD POST-PROCESSING STAGE
    idx r = lcm(r1, r2); // candidate order of a mod N
    if (r % 2 == 0 && modpow(a, r / 2, N) != static_cast<bigint>(N - 1)) {
        // at least one of those is a non-trivial factor
        bigint p = gcd(modpow(a, r / 2, N) - 1, N);
        bigint q = gcd(modpow(a, r / 2, N) + 1, N);
        if (p == 1)
            p = N / q;
        if (q == 1)
            q = N / p;
    } else {
        return;
    }
}

void benchmark(uint32_t qubits, uint32_t count, void (*func)(uint32_t), const char* name)
{
    int64_t total = 0;
    int64_t total_squared = 0;

    for (uint32_t i = 0; i < count; ++i)
    {
        auto start = std::chrono::high_resolution_clock::now();
        func(qubits);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        total += duration;
        total_squared += duration * duration;
    }

    float average = total / (float)count;
    float squared_average = total_squared / (float)count;

    std::cout << name << "\n";
    std::cout << "Total time: " << total / 1000.f << " ms; Average Time: " << average / 1000.f
              << " ms; Std. Dev: " << count / (float)(count - 1) * std::sqrt(squared_average - average * average) / 1000.f << " ms; Count: " << count << '\n' << std::endl;
}

int main()
{
    const uint32_t qubits = 10;
    const uint32_t count = 100;

    benchmark(qubits, count, grover, "Grover");
    benchmark(qubits, count, qft, "QFT");
    benchmark(qubits, count, shor, "Shor");

    return 0;    
}