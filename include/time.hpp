#ifndef TIME_HPP
#define TIME_HPP

#include <chrono>
#include "string.hpp"

// Return the current time [seconds]
double now() {
    return std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now().time_since_epoch()
    ).count()*1e-6;
}

// Converts an ammount of time [seconds] to a formatted string hh:mm:ss
template<typename T>
std::string date_str(T t) {
    std::string date;
    
    if (t < 1.0) {
        int mag = floor(log10(t));
        if (mag >= -3) {
            return strn(round(t*1e3), 3)+"ms";
        } else if (mag >= -6) {
            return strn(round(t*1e6), 3)+"us";
        } else {
            return strn(round(t*1e9), 3)+"ns";
        }
    } else {
        std::size_t day  = floor(t/(24*60*60));
        std::size_t hour = floor(t/(60*60)) - day*24;
        std::size_t min  = floor(t/60) - day*24*60 - hour*60;
        std::size_t sec  = floor(t) - day*24*60*60 - hour*60*60 - min*60;
        
        if (day  != 0) date += strn(day)+'d';
        if (hour != 0) date += strn(hour,2)+'h';
        if (min  != 0) date += strn(min,2)+'m';
        date += strn(sec,2)+'s';
        
        if (date[0] == '0' && date.size() != 2) {
            date.erase(0,1);
        }
    }
    
    return date;
}

// Execute the provided code and return the time it took to execute [seconds]
template<typename F>
double profile(F func) {
    auto start = now();
    func();
    return now() - start;
}

template<typename T>
void progress_(std::size_t i, std::size_t n, T start) {
    double total = now() - start;
    double remaining = total*double(n)/(i+1) - total;

    const std::size_t ndash = 50;
    std::cout << "\r[" << std::string(floor(ndash*(i+1)/double(n)),'-')
        << std::string(ndash - floor(ndash*(i+1)/double(n)),' ') << "]";
    std::cout << " " << strn((i+1), floor(log10(double(n))) + 1);
    std::cout << " " << strn(std::size_t(floor(100.0*(i+1)/double(n))), 3)
        << "%, " << date_str(total) << " elapsed, " << date_str(remaining) << " left" << std::flush;
}

// Display/updates a progress bar, for 'i' iterations among 'n', where the computation
// started at time 'start' [seconds].
template<typename I, typename T>
void progress(I& i, std::size_t n, T start) {
    progress_(i, n, start);

    ++i;

    if (std::size_t(i) >= n) {
        std::cout << std::endl;
    }
}

// Display/updates a progress bar, for 'i' iterations among 'n', where the computation
// started at time 'start' [seconds].
template<typename I, typename T>
void print_progress(const I& ti, std::size_t n, T start) {
    std::size_t i = ti;
    if (i >= n) i = n-1;

    progress_(i, n, start);

    if (i == n-1) {
        std::cout << std::endl;
    }
}

#endif
