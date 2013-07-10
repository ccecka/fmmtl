#pragma once

#include <iostream>
#include <iomanip>
#include <cassert>

#include <map>
#include <string>

#include <sys/time.h>

/** Clock class, useful when timing code.
 */
struct Clock {
  /** Construct a Clock and start timing. */
  Clock() {
    start();
  }
  /** Start the clock. */
  inline void start() {
    time_ = now();
  }
  /** Return the seconds elapsed since the last start. */
  inline double elapsed() const {
    return sec_diff(time_, now());
  }
  /** Return the seconds difference between two Clocks */
  inline double operator-(const Clock& clock) const {
    return sec_diff(time_, clock.time_);
  }
 private:
  timeval time_;
  inline static timeval now() {
    timeval tv;
    gettimeofday(&tv, nullptr);
    return tv;
  }
  // Return the time difference (t2 - t1) in seconds
  inline static double sec_diff(const timeval& t1, const timeval& t2) {
    timeval dt;
    timersub(&t1, &t2, &dt);
    return seconds(dt);
  }
  inline static double seconds(const timeval& tv) {
    return tv.tv_sec + 1e-6 * tv.tv_usec;
  }
};


/** Timer class using RAII to easily time sections of code
 * E.g.
 * { Timer("Critial code") t;
 *   // Code to time
 * }
 */
class Timer {
  Clock clock;
  std::string msg;
 public:
  Timer(const std::string& s) : msg(s) {
    clock.start();
  }
  ~Timer() {
    double secs = clock.elapsed();
    std::cerr << "Timer [" << msg << "]: " << secs << "s" << std::endl;
  }
  std::string& string() { return msg; }
};



//! Timer and Trace logger
class Logger {

  struct EventData {
    Clock start_time;
    double total_time;
    int hit;
    EventData()
        : start_time(), total_time(0), hit(0) {
    }
    void start() {
      start_time.start();
    }
    double log(const Clock& end_time) {
      double elapsed = (end_time - start_time);
      total_time += elapsed;
      ++hit;
      return elapsed;
    }
    friend std::ostream& operator<<(std::ostream& s, const EventData& e) {
      return s << e.hit << " (calls) * "
               << e.total_time/e.hit << " (sec/call) = "
               << e.total_time << " (secs)";
    }
  };

  typedef std::map<std::string, EventData> EventMap;
  EventMap data_;

  // Output the event name and EventData
  friend std::ostream& operator<<(std::ostream& s,
                                  const typename EventMap::value_type& event) {
    return s << std::setw(20) << std::left << event.first
             << " : " << event.second;
  }

 public:
  //! Start a clock for an event
  inline void start(const std::string& event) {
    data_[event].start();
  }

  //! Return the elasped time for given event
  double stop(const std::string& event, bool print_event = false) {
    Clock end_time;      // Stop the clock
    typename EventMap::iterator it = data_.find(event);
    assert(it != data_.end());
    double elapsed = it->second.log(end_time);
    if (print_event) std::cout << *it << std::endl;
    return elapsed;
  }

  //! Erase entry in timer
  inline void clear(const std::string& event) {
    data_.erase(event);
  }

  //! Erase all events in timer
  inline void clear() {
    data_.clear();
  }

  // Get an event's data? RAII with start(event)?
  //PublicEventData operator[](const std::string& event) { ... }

  //! Print all events and timing to an ostream
  friend std::ostream& operator<<(std::ostream& s, const Logger& log) {
    for (const typename EventMap::value_type& event : log.data_)
      s << event << std::endl;
    return s;
  }
};

//! Global static logger rather than a singleton for efficiency/consistency
static Logger fmm_logger;


#if 0
int main() {
  Logger log;
  log.start("Outer");
  log.start("Inner");
  log.stop("Inner");
  log.stop("Outer");
  std::cout << log << std::endl;

  global_logger.start("Outer");
  global_logger.start("Inner");
  global_logger.stop("Inner");
  global_logger.stop("Outer");
  std::cout << global_logger << std::endl;
  return 0;
}
#endif
